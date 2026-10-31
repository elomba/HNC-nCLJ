#!/usr/bin/env python3
"""Automatic per-replica entropy/MSD/VACF analysis and locked study aggregation."""
from pathlib import Path
import argparse,fcntl,json,os,sys
import numpy as np
from scipy.interpolate import CubicSpline,PchipInterpolator
from scipy.integrate import cumulative_trapezoid
from scipy.stats import t as student_t
from netCDF4 import Dataset
ROOT=Path(__file__).resolve().parent
DT=.002


def entropy(rho,z,u,q2,pchip=False):
    x=np.sqrt(rho);C=(4*np.pi*.5*q2)**1.5/(12*np.pi)
    cls=PchipInterpolator if pchip else CubicSpline
    integ=cls(np.r_[0,x],np.r_[-C,2*(z-1)/x]).antiderivative()
    return u/2-(integ(x)-integ(0))


def configurational_z(mean_rows,rho,N):
    # Default LAMMPS all-atom temperature has 3N-3 degrees of freedom.
    # P_kin = rho*(1-1/N)*T_inst. Remove its sampled contribution and add
    # the analytic target-temperature ideal term for the excess EOS.
    return 1+(mean_rows[:,2]-rho*(1-1/N)*mean_rows[:,1])/(rho*2.)


def transport(folder,c):
    with Dataset(folder/'positions.nc') as f:
        f.set_auto_mask(False)
        time=np.asarray(f.variables['time'][:],float)
        pos=np.asarray(f.variables['unwrapped_coordinates'][:],float)
        typ=np.asarray(f.variables['type'][0],int)
        L=np.asarray(f.variables['cell_lengths'][0],float)
    assert len(time)>100 and np.all(np.diff(time)>0)
    pos-=pos.mean(axis=1,keepdims=True)  # all equal masses; ONE global frame
    n=len(time);maxlag=n//4
    lags=np.unique(np.rint(np.geomspace(1,maxlag,80)).astype(int))
    lagtime=np.array([np.mean(time[lag:]-time[:-lag]) for lag in lags])
    msd=[]
    for lag in lags:
        origins=np.arange(0,n-lag,max(1,(n-lag)//200))
        dr=pos[origins+lag]-pos[origins]
        per_atom=np.mean(np.sum(dr*dr,axis=2),axis=0)
        msd.append([per_atom.mean()]+[per_atom[typ==i+1].mean() for i in range(len(c['species']))])
    msd=np.array(msd)
    x=np.array([s['x'] for s in c['species']])
    assert np.allclose(msd[:,0],msd[:,1:]@x,rtol=1e-10,atol=1e-10)
    windows={}
    for lo,hi in [(10,50),(20,100),(50,200)]:
        ok=(lagtime>=lo)&(lagtime<=hi)
        assert ok.sum()>=5
        windows[f'{lo}_{hi}']=[float(np.polyfit(lagtime[ok],column[ok],1)[0]/6) for column in msd.T]
    np.savetxt(folder/'msd_multi_origin.csv',np.column_stack([lagtime,msd]),delimiter=',',
               header='lag_tau,net,'+','.join(s['label'] for s in c['species']))
    # Low-k concentration fluctuations are a diagnostic, not a phase-stability proof.
    fluctuations={}
    take=np.linspace(0,n-1,min(n,200)).astype(int)
    phases=np.exp(2j*np.pi*pos[take]/L)
    total=phases.sum(axis=1)
    for i,s in enumerate(c['species']):
        signal=phases[:,typ==i+1].sum(axis=1)-s['x']*total
        fluctuations[s['label']]=float(np.mean(abs(signal)**2)/(len(typ)*s['x']*(1-s['x'])))
    with Dataset(folder/'velocities.nc') as f:
        f.set_auto_mask(False)
        vt=np.asarray(f.variables['time'][:],float)
        vel=np.asarray(f.variables['velocities'][:],float)
        vtyp=np.asarray(f.variables['type'][0],int)
    com=np.loadtxt(folder/'global_vcom.dat')
    steps=np.rint(vt/DT).astype(int)
    common,iv,ic=np.intersect1d(steps,com[:,0].astype(int),return_indices=True)
    assert len(common)>=len(vt)-2
    vel=vel[iv]-com[ic,1:4,None].transpose(0,2,1)
    vt=vt[iv];dt=float(np.median(np.diff(vt)));nv=len(vt)
    maxcv=min(int(20/dt),nv//4);nfft=1<<int(2*nv-1).bit_length()
    running=[];gk={};cv_all=[]
    for i,s in enumerate(c['species']):
        indices=np.flatnonzero(vtyp==i+1);assert len(indices)>0
        cv=np.zeros(maxcv+1)
        for begin in range(0,len(indices),32):
            v=vel[:,indices[begin:begin+32],:]
            spec=np.fft.rfft(v,n=nfft,axis=0)
            cv+=np.fft.irfft(spec*spec.conjugate(),n=nfft,axis=0)[:maxcv+1].sum(axis=(1,2))
        cv/=3*len(indices)*(nv-np.arange(maxcv+1))
        rt=np.arange(maxcv+1)*dt
        integral=cumulative_trapezoid(cv,dx=dt,initial=0)
        gk[s['label']]={str(cut):float(np.interp(cut,rt,integral)) for cut in [2,5,10,20] if cut<=rt[-1]}
        cv_all.append(cv);running.append(integral)
    np.savetxt(folder/'vacf_and_integrals.csv',np.column_stack([rt]+cv_all+running),delimiter=',')
    probe=[i for i,s in enumerate(c['species']) if s['label'].startswith('probe_')]
    primary=np.array(windows['20_100'])
    return dict(D_net=float(primary[0]),D_probe=float(primary[1+np.array(probe)].mean()),
                D_species={s['label']:float(primary[i+1]) for i,s in enumerate(c['species'])},
                fitting_windows=windows,Green_Kubo=gk,concentration_fluctuations=fluctuations,
                n_position_frames=n,n_velocity_frames=nv,
                note='MSD confidence intervals come from independent replica contrasts, not overlapping origins')


def analyze_replica(folder):
    m=json.loads((folder/'manifest.json').read_text());c=m['composition']
    assert not m.get('smoke'),folder
    rhos=np.array(m['densities']);r=[];means=[];windows=[]
    for rho in rhos:
        d=folder/('transport' if np.isclose(rho,.7) else f'rho_{rho:.4f}')
        assert (d/'verified.json').exists(),d
        a=np.loadtxt(d/'thermo_avg.dat')
        retained=a[len(a)//2:]
        means.append(retained.mean(axis=0))
        windows.append([a[int(f*len(a)):].mean(axis=0) for f in [0.,.25,.5]])
    a=np.array(means)
    # Columns: step, instantaneous temperature, pressure, target-T Z, potential energy.
    z=configurational_z(a,rhos,m['N'])
    sex=entropy(rhos,z,a[:,4],c['q2'])
    assert np.isfinite(sex).all()
    transport_result=transport(folder/'transport',c)
    sensitivity=[]
    windows=np.array(windows)
    for k in range(3):
        vals=windows[:,k,:]
        sensitivity.append(float(entropy(rhos,configurational_z(vals,rhos,m['N']),vals[:,4],c['q2'])[-1]))
    np.savetxt(folder/'entropy_density.csv',np.column_stack([rhos,z,a[:,4]/2,sex]),delimiter=',',
               header='rho,Z,u_ex_over_kBT,s_ex_over_kB')
    result=dict(name=c['name'],replica=m['replica'],seed=m['seed'],s_ex=float(sex[-1]),
                s_ex_pchip=float(entropy(rhos,z,a[:,4],c['q2'],True)[-1]),
                averaging_window_entropy=sensitivity,transport=transport_result,
                s_ex_from_total_pressure=float(entropy(rhos,a[:,2]/(rhos*2),a[:,4],c['q2'])[-1]),
                pressure_reference='configurational virial plus target-T ideal term; total kinetic DOF=3N-3',
                HNC_s_ex=json.loads((ROOT/'hnc/fine'/c['name']/'result.json').read_text())['s_ex'])
    target=folder/'replica_result.json';temp=folder/f'.result.{os.getpid()}.tmp'
    temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(target)
    print('Analyzed',folder,flush=True)


def estimate(values):
    v=np.array(values,float)
    if len(v)<2:return dict(n=len(v),mean=float(v.mean()),CI=None)
    error=student_t.ppf(.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))
    return dict(n=len(v),mean=float(v.mean()),CI=[float(v.mean()-error),float(v.mean()+error)])


def aggregate():
    (ROOT/'results').mkdir(exist_ok=True)
    with (ROOT/'results/.aggregate.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        rows=[json.loads(p.read_text()) for p in sorted((ROOT/'md').glob('*/replica_*/replica_result.json'))]
        by={}
        for row in rows:by.setdefault(row['name'],{})[row['replica']]=row
        summary={name:{key:estimate([r['s_ex'] if key=='s_ex' else r['transport'][key] for r in reps.values()])
                           for key in ['s_ex','D_net','D_probe']} for name,reps in by.items()}
        contrasts={}
        comparisons=[('label_null','dispersion_split_d000','weak_f030_e025'),
                     ('moderate_split','dispersion_split_d025','weak_f030_e025'),
                     ('high_split','dispersion_split_d040','weak_f030_e025'),
                     ('high_vs_baseline','dispersion_split_d040','baseline_3C'),
                     ('mono_high_split','mono_weak_split7C','mono_weak_uniform6C'),
                     ('mono_high_vs_baseline','mono_weak_split7C','mono_f010'),
                     ('salt_addition_baseline','mono_f010','baseline_3C')]
        for label,a,b in comparisons:
            common=sorted(set(by.get(a,{}))&set(by.get(b,{})))
            if not common:continue
            diff=[by[a][r]['s_ex']-by[b][r]['s_ex'] for r in common]
            net=[np.log(by[a][r]['transport']['D_net']/by[b][r]['transport']['D_net']) for r in common]
            probe=[np.log(by[a][r]['transport']['D_probe']/by[b][r]['transport']['D_probe']) for r in common]
            contrasts[label]=dict(delta_s_ex=estimate(diff),delta_log_D_net=estimate(net),delta_log_D_probe=estimate(probe))
        result=dict(completed_replicas=len(rows),expected_replicas=24,complete=len(rows)==24,
                    statistics='paired replica differences and Student t intervals; 3 replicas per mixture',
                    systems=summary,contrasts=contrasts)
        out=ROOT/'results/summary.json';tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(out)
        text=f'# Component entropy study: {len(rows)}/24 replicas analyzed\n\n'
        text+='| Mixture | Replicas | s_ex/kB | D_net | D_probe |\n|---|---:|---:|---:|---:|\n'
        for name,r in summary.items():text+=f"| {name} | {r['s_ex']['n']} | {r['s_ex']['mean']:.5f} | {r['D_net']['mean']:.5f} | {r['D_probe']['mean']:.5f} |\n"
        text+='\nSee summary.json for paired confidence intervals and individual replica diagnostics. Incomplete replica sets cannot establish a transport trend.\n'
        (ROOT/'results/progress.md').write_text(text)
        print('Aggregated',len(rows),'/24 replicas',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--replica',type=Path);a=p.parse_args()
    if a.replica:analyze_replica(a.replica)
    aggregate()
