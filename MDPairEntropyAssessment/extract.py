#!/usr/bin/env python3
"""Extract reusable frame-level RDF histograms without changing MD files."""
import argparse,json,os,time,socket
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from pair_entropy import ROOT,STUDY,ORDER,DR,NBINS,Kernel,compositions,digest

def run(index,limit=None,ideal=False):
    code,name=ORDER[index//3]; replica=index%3+1
    c=compositions()[name]; ns=len(c['species'])
    source=STUDY/'md'/name/f'replica_{replica}'/'transport'/'positions.nc'
    tag=f'{code}_r{replica}'+('_ideal' if ideal else '')+(f'_pilot{limit}' if limit else '')
    cache=ROOT/'cache';cache.mkdir(exist_ok=True)
    out=cache/(tag+'.npz');meta=cache/(tag+'.json')
    if out.exists() and meta.exists():
        print('Already extracted',tag,flush=True);return
    started=time.monotonic();kernel=Kernel(); frame_seconds=[]
    with Dataset(source) as f:
        f.set_auto_mask(False)
        n=len(f.dimensions['atom']);nf=len(f.dimensions['frame'])
        expected=np.array([round(n*s['x']) for s in c['species']])
        types0=np.asarray(f['type'][0],np.int32)-1
        ids0=np.asarray(f['id'][0],np.int32)
        assert np.array_equal(np.sort(ids0),np.arange(1,n+1))
        assert np.array_equal(np.bincount(types0,minlength=ns),expected)
        box0=np.asarray(f['cell_lengths'][0],float)
        assert NBINS*DR<min(box0)/2
        frames=np.unique(np.rint(np.linspace(0,nf-1,limit)).astype(int)) if limit else np.arange(nf)
        counts=np.empty((len(frames),ns*(ns+1)//2,NBINS),np.uint32)
        times=np.empty(len(frames)); boxes=np.empty((len(frames),3))
        rng=np.random.default_rng(731000+index)
        for start in range(0,len(frames),32):
            ix=frames[start:start+32]
            if ideal:
                # No MD coordinates enter the independent random-point test.
                xyz=rng.random((len(ix),n,3))*box0
                ids=np.broadcast_to(ids0,(len(ix),n))
                typ=np.broadcast_to(types0,(len(ix),n))
                cells=np.broadcast_to(box0,(len(ix),3))
                tt=ix*.5
            else:
                # One chunk bounds memory and reduces repeated Python/NetCDF calls.
                xyz=np.asarray(f['unwrapped_coordinates'][ix],float)
                ids=np.asarray(f['id'][ix],np.int32)
                typ=np.asarray(f['type'][ix],np.int32)-1
                cells=np.asarray(f['cell_lengths'][ix],float)
                angles=np.asarray(f['cell_angles'][ix],float)
                tt=np.asarray(f['time'][ix],float)
                assert np.max(abs(angles-90))<1e-5
            assert np.max(abs(cells-box0))<1e-5
            for j in range(len(ix)):
                assert np.array_equal(np.sort(ids[j]),np.arange(1,n+1))
                order=np.argsort(ids[j]); ref=np.argsort(ids0)
                assert np.array_equal(typ[j,order],types0[ref])
                assert np.isfinite(xyz[j]).all()
                tick=time.monotonic()
                counts[start+j]=kernel.histogram(xyz[j],typ[j],cells[j],ns=ns)
                frame_seconds.append(time.monotonic()-tick)
            times[start:start+len(ix)]=tt; boxes[start:start+len(ix)]=cells
            if start%256==0: print(tag,'frames',min(start+32,len(frames)),'/',len(frames),flush=True)
        assert np.all(np.diff(times)>0)
    temp=cache/(tag+'.tmp.npz')
    np.savez_compressed(temp,counts=counts,frames=frames,times=times,boxes=boxes,particle_counts=expected,types=types0,ids=ids0)
    temp.replace(out)
    info=dict(code=code,name=name,replica=replica,ideal=ideal,limit=limit,nframes=len(frames),n=n,
        nspecies=ns,dr=DR,nbins=NBINS,rmax=DR*NBINS,box=box0.tolist(),
        seed=731000+index if ideal else None,source=str(source),source_bytes=source.stat().st_size,
        source_mtime_ns=source.stat().st_mtime_ns,source_sha256=digest(source) if not ideal else None,
        kernel_sha256=digest(ROOT/'histogram.c'),extract_sha256=digest(Path(__file__)),
        histogram_sha256=digest(out),host=socket.gethostname(),threads=int(os.environ.get('OMP_NUM_THREADS','1')),
        frame_histogram_seconds_mean=float(np.mean(frame_seconds)),elapsed_seconds=time.monotonic()-started,
        species=c['species'])
    meta.write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps({k:info[k] for k in ['code','replica','ideal','nframes','elapsed_seconds','frame_histogram_seconds_mean']}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);p.add_argument('--limit',type=int);p.add_argument('--ideal',action='store_true')
    a=p.parse_args();run(a.index,a.limit,a.ideal)
