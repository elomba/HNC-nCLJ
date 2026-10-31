#!/usr/bin/env python3
"""Resolve structural binning error below 3.2 sigma; retain measured outer RDFs."""
import argparse,json,time
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from pair_entropy import ROOT,ORDER,Kernel,digest,merge_map,merge_hist,evaluate,scalar_metrics

def run(index,ideal=False):
    code,name=ORDER[index//3];rep=index%3+1;tag=f'{code}_r{rep}'+('_ideal' if ideal else '')
    meta=json.loads((ROOT/'cache'/(tag+'.json')).read_text())
    with np.load(ROOT/'cache'/(tag+'.npz')) as d:
        old_labels=d['counts'];nlabels=d['particle_counts'];types=d['types'];box=d['boxes'][0]
    cache=ROOT/'cache'/(tag+'_core0025.npz');dr=.0025;nbins=1280;late=np.arange(1001,2001)
    start=time.monotonic();kernel=Kernel();ns=len(nlabels);n=sum(nlabels)
    if cache.exists():
        with np.load(cache) as d:hlabels=d['counts']
    else:
        hlabels=np.empty((1000,ns*(ns+1)//2,nbins),np.uint32)
        if ideal:
            rng=np.random.default_rng(meta['seed'])
            for f in range(2001):
                pos=rng.random((n,3))*box
                if f>=1001:hlabels[f-1001]=kernel.histogram(pos,types,box,ns=ns,dr=dr,nbins=nbins)
        else:
            with Dataset(meta['source']) as d:
                d.set_auto_mask(False)
                for first in range(0,1000,32):
                    ix=late[first:first+32]
                    pos=np.asarray(d['unwrapped_coordinates'][ix],float)
                    tt=np.asarray(d['type'][ix],np.int32)-1
                    for j in range(len(ix)):
                        hlabels[first+j]=kernel.histogram(pos[j],tt[j],box,ns=ns,dr=dr,nbins=nbins)
        tmp=cache.with_name(cache.stem+'.tmp.npz');np.savez_compressed(tmp,counts=hlabels);tmp.replace(cache)
    # Reproducible matching of trajectories/random ensembles and exact histogram edges.
    recoarse=hlabels.reshape(*hlabels.shape[:-1],320,4).sum(axis=-1)
    assert np.array_equal(recoarse,old_labels[late,:,:320]),tag+' mismatch with original core counts'
    del recoarse
    group,names=merge_map(meta['species']);old,counts=merge_hist(old_labels[late],nlabels,group)
    fine,_=merge_hist(hlabels,nlabels,group);del old_labels,hlabels
    physical=[next(s for j,s in enumerate(meta['species']) if group[j]==i) for i in range(len(counts))]
    q=[s['q'] for s in physical];volume=float(np.prod(box));old_total=old.sum(axis=0,dtype=np.uint64);new_total=fine.sum(axis=0,dtype=np.uint64)
    def measurements(oh,fh,F,width):
        outer=scalar_metrics(evaluate(oh,counts,F,volume,stride=1),counts,q)
        core=scalar_metrics(evaluate(oh,counts,F,volume,stride=1,rmax=3.2),counts,q)
        new=scalar_metrics(evaluate(fh,counts,F,volume,dr=dr,stride=round(width/dr)),counts,q)
        return {k:outer[k]-core[k]+new[k] for k in outer}
    widths=[.0025,.005,.01,.02,.04,.08]
    raw={str(w):measurements(old_total,new_total,1000,w) for w in widths}
    jk={};block_checks={}
    for blocks in [10,20,40]:
        old_blocks=old.reshape(blocks,1000//blocks,*old.shape[1:]).sum(axis=1,dtype=np.uint64)
        new_blocks=fine.reshape(blocks,1000//blocks,*fine.shape[1:]).sum(axis=1,dtype=np.uint64)
        for width in (widths if blocks==20 else [.0025,.005]):
            vals=[measurements(old_total-old_blocks[b],new_total-new_blocks[b],1000-1000//blocks,width) for b in range(blocks)]
            corrected={k:float(blocks*raw[str(width)][k]-(blocks-1)*np.mean([v[k] for v in vals])) for k in raw[str(width)]}
            if blocks==20:jk[str(width)]=corrected
            block_checks[f'{blocks}_{width}']=corrected
    # The 0.01-core hybrid must reproduce the original full 0.01 estimator.
    original=json.loads((ROOT/('ideal_results' if ideal else 'replicas')/(tag+'.json')).read_text())
    assert abs(raw['0.01']['s2']-original['bin_sensitivity']['0.01']['s2'])<1e-11
    assert abs(jk['0.01']['s2']-original['jackknife_bin_sensitivity']['0.01']['s2'])<1e-10
    result={'code':code,'replica':rep,'ideal':ideal,'raw_hybrid':raw,'jackknife20_hybrid':jk,'jackknife_block_sensitivity':block_checks,
            'core_dr':dr,'core_radius':3.2,'tail_dr':.01,'total_radius':8.88,'frames':[1001,2000],
            'original_core_histograms_match_exactly':True,'cache_sha256':digest(cache),
            'script_sha256':digest(Path(__file__)),'elapsed_seconds':time.monotonic()-start,
            'note':'Raw original dr=.02 remains primary; refinement was triggered by bin-width sensitivity, not TI agreement. Hybrid uses fine core plus sampled .01 tail; no fitted tail.'}
    out=ROOT/'refined';out.mkdir(exist_ok=True)
    (out/(tag+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(tag,'refined raw',raw['0.005']['s2'],'corrected',jk['0.005']['s2'],'finest',jk['0.0025']['s2'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);p.add_argument('--ideal',action='store_true');a=p.parse_args();run(a.index,a.ideal)
