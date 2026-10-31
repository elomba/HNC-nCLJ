#!/usr/bin/env python3
"""Evaluate the prespecified rehabilitation protocol on existing velocity records."""
from pathlib import Path
import argparse,hashlib,json,os
import numpy as np
from scipy.fft import rfft,irfft,next_fast_len
from netCDF4 import Dataset
from two_pt import spectrum,evaluate

HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
CACHE=HERE/'cache';RESULTS=HERE/'replicas'
DT=.002
PROTOCOL=json.loads((HERE/'protocol.json').read_text())

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def key(s):return (s['q'],s['eps'],s['sigma'],s['mass'])

def biased_vacf(v):
    n,na,_=v.shape;nfft=next_fast_len(2*n-1)
    c=np.zeros(n)
    for j in range(0,na,32):
        f=rfft(v[:,j:j+32,:],n=nfft,axis=0)
        c+=irfft(f*f.conjugate(),n=nfft,axis=0)[:n].sum(axis=(1,2))
    return c/(3*na*n)

def cache_replica(folder,force=False):
    m=json.loads((folder/'manifest.json').read_text());name=m['composition']['name'];rep=m['replica'];stem=f'{name}_r{rep}'
    zfile=CACHE/(stem+'.npz');jfile=CACHE/(stem+'.json')
    if zfile.exists() and jfile.exists() and not force:
        meta=json.loads(jfile.read_text())
        assert meta['cache_version']==1
        return meta,dict(np.load(zfile))
    traj=folder/'transport/velocities.nc';cf=folder/'transport/global_vcom.dat'
    with Dataset(traj) as f:
        f.set_auto_mask(False)
        time=np.asarray(f.variables['time'][:],float)
        ids=np.asarray(f.variables['id'][0],int);typ=np.asarray(f.variables['type'][0],int)
        assert np.all(np.diff(ids)>0)
        assert np.array_equal(ids,f.variables['id'][-1]) and np.array_equal(typ,f.variables['type'][-1])
        vel=np.asarray(f.variables['velocities'][:],float)
    com=np.loadtxt(cf);steps=np.rint(time/DT).astype(int)
    common,iv,ic=np.intersect1d(steps,com[:,0].astype(int),return_indices=True)
    assert len(common)>=len(time)-2
    vel=vel[iv]-com[ic,1:4,None].transpose(0,2,1)
    interval=int(np.median(np.diff(common)));assert np.all(np.diff(common)==interval)
    dt=interval*DT
    species=m['composition']['species'];merged={}
    for i,s in enumerate(species):merged.setdefault(key(s),[]).append(i)
    groups=[]
    for ii in merged.values():
        sp=species[ii[0]];indices=np.flatnonzero(np.isin(typ,np.array(ii)+1))
        if sp['q']==0:indices=indices[:200]
        groups.append(dict(label=sp['label'] if len(ii)==1 else 'merged_identical_solvent',
                           x=sum(species[j]['x'] for j in ii),key=list(key(sp)),indices=indices,
                           source_labels=[species[j]['label'] for j in ii],source_fractions=[species[j]['x'] for j in ii],
                           primary=True,secondary=sp['q']==0 and sp['eps']!=1,probe=abs(sp['q'])==2))
    if name=='dispersion_split_d000':
        for i,s in enumerate(species):
            if s['label'] in ['weak_A','weak_B']:
                groups.append(dict(label='label_diagnostic_'+s['label'],x=s['x'],key=list(key(s)),
                    indices=np.flatnonzero(typ==i+1),source_labels=[s['label']],source_fractions=[s['x']],primary=False,secondary=True,probe=False))
    arrays={};info=[];n=len(vel)
    segments={'full':slice(None),'first_half':slice(0,n//2),'last_half':slice(n//2,None)}
    for gi,g in enumerate(groups):
        indices=g.pop('indices');v=vel[:,indices,:]
        for seg,sl in segments.items():arrays[f'g{gi}_{seg}']=biased_vacf(v[sl])
        info.append(dict(**g,sampled_ids=ids[indices].tolist(),sample_count=len(indices),
                    corrected_velocity_sha256=hashlib.sha256(v.astype('<f8').tobytes()).hexdigest()))
    meta=dict(cache_version=1,name=name,replica=rep,dt=dt,nframes=n,groups=info,
              source=str(traj.relative_to(STUDY)),source_bytes=traj.stat().st_size,source_mtime_ns=traj.stat().st_mtime_ns,
              manifest_sha256=digest(folder/'manifest.json'),common_com_sha256=digest(cf),segments={k:len(vel[v]) for k,v in segments.items()})
    CACHE.mkdir(exist_ok=True);np.savez_compressed(zfile,**arrays);jfile.write_text(json.dumps(meta,indent=2)+'\n')
    return meta,arrays


def run(folder,force=False):
    meta,arrays=cache_replica(folder,force);name=meta['name'];rep=meta['replica'];groups=meta['groups']
    tr=json.loads((folder/'replica_result.json').read_text());out=dict(name=name,replica=rep,TI_s_ex=tr['s_ex'],nframes=meta['nframes'],dt=meta['dt'],
       method='classical R2PT; one-fluid equal-partial-volume approximation; no entropy fit',variants={},
       source_cache=str((CACHE/f'{name}_r{rep}.npz').relative_to(HERE)),protocol_sha256=digest(HERE/'protocol.json'),
       code_sha256={x:digest(HERE/x) for x in ['two_pt.py','run_rehabilitation.py']})
    main=[i for i,g in enumerate(groups) if g['primary']];assert abs(sum(groups[i]['x'] for i in main)-1)<1e-12
    for seg in ['full','first_half','last_half']:
        out['variants'][seg]={}
        for tau in PROTOCOL['window_sensitivity_tau']:
            spectra={i:spectrum(arrays[f'g{i}_{seg}'],meta['dt'],tau,2.) for i in range(len(groups))}
            for exponent in [1.,1.5]:
                label=f"tau{tau if tau is not None else 'none'}_d{exponent:g}"
                per={};sx=0
                for i,g in enumerate(groups):
                    nu,G,diag=spectra[i];a=evaluate(nu,G,exponent=exponent)
                    mD=sum(tr['transport']['D_species'][l]*x for l,x in zip(g['source_labels'],g['source_fractions']))/g['x']
                    a.update(**diag,label=g['label'],x=g['x'],primary=g['primary'],D_MSD=mD,
                        spectral_D_relative_error=diag['D_spectral_equipartition']/mD-1,
                        reference_h_invariance_span=float(np.ptp([evaluate(nu,G,exponent=exponent,h=h,diagnostics=False)['s_ex'] for h in [.1,1.,10.]])))
                    per[g['label']]=a
                    if g['primary']:sx+=g['x']*a['s_ex']
                mixG=sum(groups[i]['x']*spectra[i][1] for i in main);mixnu=spectra[main[0]][0]
                pooled=evaluate(mixnu,mixG,exponent=exponent)
                out['variants'][seg][label]=dict(s_ex=sx,species=per,pooled_spectrum_s_ex=pooled['s_ex'],
                    probe_spectral_D=sum(groups[i]['x']*spectra[i][2]['D_spectral_equipartition'] for i in main if groups[i]['probe'])/.06,
                    D_probe_MSD=tr['transport']['D_probe'])
                if seg=='full' and tau==10 and exponent==1:
                    # Mechanical merger removes arbitrary-label nonlinearities; unmerged finite-sample estimate is diagnostic only.
                    naive=sx
                    if name=='dispersion_split_d000':
                        gi=next(i for i in main if len(groups[i]['source_labels'])==2)
                        naive-=groups[gi]['x']*per[groups[gi]['label']]['s_ex']
                        naive+=sum(g['x']*per[g['label']]['s_ex'] for g in groups if not g['primary'])
                    out['label_unmerged_finite_sample_s_ex']=naive
                    out['volume_sensitivity']={}
                    for which in ['probe','secondary']:
                        sel=[i for i in main if groups[i][which]]
                        if not sel:continue
                        rest=[i for i in main if i not in sel]
                        xs=sum(groups[i]['x'] for i in sel);xr=1-xs
                        for shift in [-.1,.1]:
                            value=0.;weighted_volume=0
                            for i in main:
                                vf=1+shift if i in sel else 1-xs*shift/xr
                                weighted_volume+=groups[i]['x']*vf/.7
                                nu,G,_=spectra[i]
                                value+=groups[i]['x']*evaluate(nu,G,rho_eff=.7/vf,diagnostics=False)['s_ex']
                            assert abs(weighted_volume-1/.7)<1e-12
                            out['volume_sensitivity'][f'{which}_{shift:+.1f}']=value
                    out['historical_lnZ_diagnostic_s_ex']=sum(groups[i]['x']*evaluate(spectra[i][0],spectra[i][1],hs_lnZ=True)['s_ex'] for i in main)
                    out['frequency_resolution_sensitivity']={}
                    for mode in ['padding2','padding4','decimate2','decimate4']:
                        stride=int(mode[-1]) if mode.startswith('decimate') else 1
                        padding=int(mode[-1]) if mode.startswith('padding') else 1
                        sex=0.
                        for i in main:
                            nu2,G2,_=spectrum(arrays[f'g{i}_full'][::stride],meta['dt']*stride,10,2.,padding_factor=padding)
                            sex+=groups[i]['x']*evaluate(nu2,G2,diagnostics=False)['s_ex']
                        out['frequency_resolution_sensitivity'][mode]=sex
    RESULTS.mkdir(exist_ok=True);file=RESULTS/f'{name}_r{rep}.json'
    tmp=file.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2)+'\n');tmp.replace(file)
    sx=out['variants']['full']['tau10_d1']['s_ex'];print(f'{name} r{rep}: TI={tr["s_ex"]:.6f}, R2PT={sx:.6f}, error={sx-tr["s_ex"]:+.6f}',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--only');p.add_argument('--replica',type=int);p.add_argument('--force-cache',action='store_true');a=p.parse_args()
    folders=sorted((STUDY/'md').glob('*/replica_*'))
    for folder in folders:
        if a.only and folder.parent.name!=a.only:continue
        if a.replica and folder.name!=f'replica_{a.replica}':continue
        run(folder,a.force_cache)
