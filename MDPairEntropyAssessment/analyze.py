#!/usr/bin/env python3
"""Analyze one immutable frame-histogram cache; TI is never used to tune RDFs."""
import argparse,json
from pathlib import Path
import numpy as np
from pair_entropy import ROOT,STUDY,ORDER,DR,merge_map,merge_hist,evaluate,scalar_metrics,digest,pairs

def jsonable(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:jsonable(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [jsonable(x) for x in v]
    return v

def tau_initial_positive(v):
    v=np.asarray(v,float);v-=v.mean()
    if np.var(v)<1e-20:return .5
    n=len(v);f=np.fft.rfft(v,n=2*n)
    ac=np.fft.irfft(f*f.conj(),n=2*n)[:n]/np.arange(n,0,-1)
    ac/=ac[0];total=.5
    for j in range(1,n-1,2):
        pair=ac[j]+ac[j+1]
        if pair<=0:break
        total+=pair
        if j>10*total:break
    return float(total)

def analyze(index,ideal=False):
    code,name=ORDER[index//3];replica=index%3+1
    tag=f'{code}_r{replica}'+('_ideal' if ideal else '')
    cache=ROOT/'cache'/(tag+'.npz');metafile=cache.with_suffix('.json')
    meta=json.loads(metafile.read_text())
    assert meta['nframes']==2001 and meta['limit'] is None
    assert digest(cache)==meta['histogram_sha256']
    with np.load(cache) as d:
        labels=d['counts'];frames=d['frames'];times=d['times'];box=d['boxes'][0];nlabels=d['particle_counts']
    assert np.array_equal(frames,np.arange(2001))
    group,names=merge_map(meta['species'])
    h,counts=merge_hist(labels,nlabels,group)
    physical=[next(s for j,s in enumerate(meta['species']) if group[j]==i) for i in range(len(counts))]
    charges=[s['q'] for s in physical]
    volume=float(np.prod(box));late=np.arange(1001,2001);early=np.arange(1,1001);allframes=np.arange(1,2001)
    def measure(ix,stride=2,rmax=8.88,normalization='ideal'):
        result=evaluate(h[ix].sum(axis=0,dtype=np.uint64),counts,len(ix),volume,
                        stride=stride,rmax=rmax,normalization=normalization)
        return scalar_metrics(result,counts,charges)
    full=measure(late)
    primary=evaluate(h[late].sum(axis=0,dtype=np.uint64),counts,len(late),volume,stride=2)
    result=dict(code=code,name=name,replica=replica,ideal=ideal,
        source_metadata=meta,cache_sha256=meta['histogram_sha256'],
        analysis_sha256=digest(Path(__file__)),estimator_sha256=digest(ROOT/'pair_entropy.py'),
        species=physical,counts=counts,volume=volume,rho=sum(counts)/volume,
        primary_frame_indices=[int(late[0]),int(late[-1])],primary_frames=len(late),
        primary_time_interval=[float(times[late[0]]),float(times[late[-1]])],
        primary=full,pair_contributions=primary['pair_contributions'],
        pair_species_indices=pairs(len(counts)),species_centered=primary['species_centered'])
    result['radial_sensitivity']={str(r):measure(late,rmax=r) for r in [4.,5.,6.,7.,8.,8.88]}
    result['bin_sensitivity']={str(DR*k):measure(late,stride=k) for k in [1,2,4,8]}
    result['sample_sensitivity']={str(n):measure(late[np.rint(np.linspace(0,999,n)).astype(int)]) for n in [50,100,200,500,1000]}
    result['time_windows']={'early':measure(early),'late':full,'all':measure(allframes)}
    result['normalization_sensitivity']={'ideal':full,'bulk':measure(late,normalization='bulk')}
    # Unequal blocks are avoided: exactly 1000 late frames and B=10,20,40.
    total=h[late].sum(axis=0,dtype=np.uint64)
    result['jackknife']={};result['jackknife_bin_sensitivity']={}
    for blocks in [10,20,40]:
        size=len(late)//blocks; leave=[]
        for ix in np.array_split(late,blocks):
            removed=h[ix].sum(axis=0,dtype=np.uint64)
            ev=evaluate(total-removed,counts,len(late)-size,volume,stride=2)
            leave.append(scalar_metrics(ev,counts,charges))
        corrected={k:float(blocks*full[k]-(blocks-1)*np.mean([v[k] for v in leave])) for k in full}
        stderr={k:float(np.sqrt((blocks-1)/blocks*np.sum((np.array([v[k] for v in leave])-np.mean([v[k] for v in leave]))**2))) for k in full}
        result['jackknife'][str(blocks)]={'corrected':corrected,'raw_minus_corrected':{k:full[k]-corrected[k] for k in full},
            'within_replica_se':stderr,'block_frames':size,'block_duration_tau':size*.5}
        if blocks==20:
            for stride in [1,2,4,8]:
                base=measure(late,stride=stride);variants=[]
                for ix in np.array_split(late,blocks):
                    ev=evaluate(total-h[ix].sum(axis=0,dtype=np.uint64),counts,
                                len(late)-size,volume,stride=stride)
                    variants.append(scalar_metrics(ev,counts,charges))
                result['jackknife_bin_sensitivity'][str(DR*stride)]={
                    k:float(blocks*base[k]-(blocks-1)*np.mean([v[k] for v in variants])) for k in base}
    if len(names)!=len(nlabels):
        unmerged=evaluate(labels[late].sum(axis=0,dtype=np.uint64),nlabels,len(late),volume,stride=2)
        result['unmerged_label_sensitivity']={'s2':unmerged['s2'],'minus_merged':unmerged['s2']-full['s2']}
    # Coordination-like observables diagnose temporal correlation without using D.
    shell=h[late,:,:150].sum(axis=-1,dtype=np.uint64)
    taus=np.array([tau_initial_positive(shell[:,k]) for k in range(shell.shape[1])])
    result['correlation_diagnostics']={'observable':'pair counts within 1.5 sigma per frame',
        'integrated_times_frames':taus,'max_tau':float(taus.max()*.5),
        'effective_frames_min':float(len(late)/(2*taus.max()))}
    if not ideal:
        ti_file=STUDY/'md'/name/f'replica_{replica}'/'replica_result.json'
        ti=json.loads(ti_file.read_text());result['TI_s_ex']=ti['s_ex']
        result['remainder_raw']=ti['s_ex']-full['s2']
        result['remainder_jackknife20']=ti['s_ex']-result['jackknife']['20']['corrected']['s2']
        result['pair_fraction_raw']=full['s2']/ti['s_ex']
        result['reference_replica_sha256']=digest(ti_file)
        result['transport']=ti['transport']
    # Full RDFs and frame counts remain reproducible local artifacts, not bloated JSON.
    rdf=ROOT/'cache'/(tag+'_rdf.npz')
    np.savez_compressed(rdf,g=primary['g'],edges=primary['edges'],counts=counts,
                        pair_contributions=primary['pair_contributions'])
    dest=ROOT/('ideal_results' if ideal else 'replicas');dest.mkdir(exist_ok=True)
    (dest/(tag+'.json')).write_text(json.dumps(jsonable(result),indent=2)+'\n')
    print(tag,'s2',full['s2'],'jackknife20',result['jackknife']['20']['corrected']['s2'],flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',required=True,type=int);p.add_argument('--ideal',action='store_true')
    a=p.parse_args();analyze(a.index,a.ideal)
