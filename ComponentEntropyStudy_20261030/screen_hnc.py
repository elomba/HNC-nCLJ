#!/usr/bin/env python3
"""HNC screening for composition/entropy controls; all probe fractions fixed."""
from pathlib import Path
import argparse,csv,json,math,os,subprocess,time,shutil
import numpy as np
ROOT=Path(__file__).resolve().parent
BINARY=ROOT.parent/'bin/hncr_nclj'


def species(label,q,x,eps=1.):return dict(label=label,q=q,x=x,eps=eps,mass=1.,sigma=1.)


def candidate(name,background,description,family):
    neutral=[p for p in background if p['q']==0]
    charged=[p for p in background if p['q']!=0]
    entries=neutral+charged+[species('probe_plus2',2,.03),species('probe_minus2',-2,.03)]
    entries=[e for e in entries if e['x']>1e-10]
    x=np.array([e['x'] for e in entries]);q=np.array([e['q'] for e in entries])
    assert abs(x.sum()-1)<1e-12 and abs(x@q)<1e-12
    assert len(entries)<=10
    return dict(name=name,species=entries,description=description,family=family,
                n_components=len(entries),I=float(.5*x@(q*q)),ion_fraction=float(x[q!=0].sum()),
                mixing_entropy=float(-x@np.log(x)),q2=float(x@(q*q)),q4=float(x@(q**4)))


def candidates():
    c=[candidate('baseline_3C',[species('solvent',0,.94)],'Original 3C probe electrolyte','baseline'),
       candidate('labels_4C',[species('solvent',0,.64),species('solvent_label',0,.30)],
                 'Mechanically identical neutral-label split; exact null control','labels')]
    for f in [.05,.10,.20]:
        c.append(candidate(f'mono_f{round(f*100):03d}',[species('solvent',0,.94-f),
                 species('background_plus1',1,f/2),species('background_minus1',-1,f/2)],
                 'Symmetric monovalent salt addition; not an ionic-strength control','monovalent'))
    for f in [.04,.10,.20]:
        for t in [0.,.25,.5,1.]:
            c.append(candidate(f'asym13_f{round(f*100):03d}_t{round(t*100):03d}',[
                species('solvent',0,.94-f),species('background_plus1',1,.75*f*(1-t)),
                species('background_minus3',-3,.25*f*(1-t)),
                species('background_minus1',-1,.75*f*t),species('background_plus3',3,.25*f*t)],
                'Mix conjugate 1:3 salts at fixed solvent, probes, ion fraction and ALL even charge moments','conjugate_salts'))
    for f in [.10,.30]:
        for eps in [.25,.50,.75]:
            c.append(candidate(f'weak_f{round(f*100):03d}_e{round(eps*100):03d}',[
                species('solvent',0,.94-f),species('weak_solvent',0,f,eps)],
                'Weak neutral solvent: fixed probes/I, changed dispersion interactions','weak_solvent'))
    return c


def write_input(path,c,fine=False,densities=None,restart=False):
    entries=c['species'];n=len(entries)
    nr,dr,num=(4096,.005,100) if fine else (2048,.01,40)
    densities=densities or [.7]
    betas=[.5] if restart else [.1,.2,.3,.4,.5]
    lines=[c['name'],f'{nr} {dr} 1.0d-7 0.5 300.0 0.05',f'{num} 1.0d-7',str(n)]
    lines+=[' '.join(['1.0']*n) for _ in entries]
    lines+=[' '.join(f'{math.sqrt(a["eps"]*b["eps"]):.12f}' for b in entries) for a in entries]
    lines+=[' '.join(str(e['q']) for e in entries),f'{len(densities)} {len(betas)}',
             ' '.join(str(r) for r in densities),' '.join(str(b) for b in betas),'.f.' if restart else '.t.',
             ' '.join(f'{e["x"]:.12f}' for e in entries[:-2])]
    lines+=[' '.join(['2.5']*n) for _ in entries]+['1.0']
    path.write_text('\n'.join(lines)+'\n')


def solve(c,fine=False):
    folder=ROOT/'hnc'/('fine' if fine else 'screen')/c['name']
    folder.mkdir(parents=True,exist_ok=True)
    output=folder/'result.json'
    if output.exists():return json.loads(output.read_text())
    coarse=ROOT/'hnc'/'screen'/c['name']/'lsol.hncr'
    restart=fine and coarse.exists()
    if restart:shutil.copy2(coarse,folder/'solin.hncr')
    write_input(folder/'input.hncr',c,fine,restart=restart)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    start=time.monotonic()
    try:
        with (folder/'stdout.log').open('w') as f:
            result=subprocess.run([str(BINARY),'input.hncr'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,
                                  env=env,timeout=600,check=True)
        log=(folder/'stdout.log').read_text()
        arr=np.loadtxt(folder/'Thermohncr.dat',comments='#',ndmin=2)
        assert len(arr)==(1 if restart else 5) and abs(arr[-1,2]-.5)<1e-8
        assert 'Aborting' not in log and np.isfinite(arr).all()
        vals=arr[-1]
        assert vals[11]>0 and vals[3]>0
        r=dict(name=c['name'],status='converged',fine=fine,s_ex=float(vals[7]),s2=float(vals[9]),
               a_ex=float(vals[6]),u_ex_over_kBT=float(vals[4]),Z=float(vals[3]),chi=float(vals[11]),
               n_components=c['n_components'],I=c['I'],mixing_entropy=c['mixing_entropy'],
               family=c['family'],elapsed_seconds=time.monotonic()-start)
    except Exception as exc:
        r=dict(name=c['name'],status='failed',error=str(exc),fine=fine,elapsed_seconds=time.monotonic()-start)
    output.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r),flush=True)
    return r


def main():
    p=argparse.ArgumentParser();p.add_argument('--fine',action='store_true');p.add_argument('--names',nargs='+');p.add_argument('--limit',type=int)
    args=p.parse_args();cs=candidates()
    extra=ROOT/'additional_candidates.json'
    if extra.exists():cs+=json.loads(extra.read_text())
    (ROOT/'candidates.json').write_text(json.dumps(cs,indent=2)+'\n')
    if args.names:cs=[c for c in cs if c['name'] in args.names]
    if args.limit:cs=cs[:args.limit]
    results=[solve(c,args.fine) for c in cs]
    (ROOT/('hnc_fine_summary.json' if args.fine else 'hnc_screen_summary.json')).write_text(json.dumps(results,indent=2)+'\n')
    good=[r for r in results if r['status']=='converged']
    if good:
        with (ROOT/('hnc_fine_summary.csv' if args.fine else 'hnc_screen_summary.csv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(good[0]));w.writeheader();w.writerows(good)
    print('Completed:',len(results),'converged:',len(good),flush=True)


if __name__=='__main__':main()
