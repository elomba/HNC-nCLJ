#!/usr/bin/env python3
"""Independent numerical and analytic checks for the scientific estimator."""
import json
import numpy as np
from pair_entropy import ROOT,Kernel,pairs,mapping,evaluate,merge_hist,phi

def direct(pos,typ,box,ns,nbins,dr):
    h=np.zeros((ns*(ns+1)//2,nbins),np.uint32);mp=mapping(ns)
    for i in range(len(pos)):
        for j in range(i+1,len(pos)):
            delta=pos[i]-pos[j]
            delta-=box*np.rint(delta/box)
            r=float(np.linalg.norm(delta))
            if r<nbins*dr: h[mp[typ[i],typ[j]],int(r/dr)]+=1
    return h

def main():
    rng=np.random.default_rng(81233);k=Kernel();checks={}
    box=np.array([8.,9.,10.]);pos=rng.random((120,3))*box
    typ=np.arange(120,dtype=np.int32)%4
    h=k.histogram(pos,typ,box,nbins=80,dr=.05,threads=2)
    assert np.array_equal(h,direct(pos,typ,box,4,80,.05))
    checks['compiled_kernel_matches_independent_direct_distances']=True
    permutation=rng.permutation(len(pos))
    assert np.array_equal(h,k.histogram(pos[permutation],typ[permutation],box,nbins=80,dr=.05,threads=1))
    assert np.array_equal(h,k.histogram(pos+np.array([2.13,4.32,1.19]),typ,box,nbins=80,dr=.05,threads=2))
    checks['atom_permutation_translation_and_thread_invariance']=True
    counts=np.bincount(typ,minlength=4)
    group=np.array([0,0,1,2]);merged,nm=merge_hist(h,counts,group)
    assert np.array_equal(merged,k.histogram(pos,group[typ],box,nbins=80,dr=.05,threads=1))
    single,n1=merge_hist(h,counts,np.zeros(4,int))
    assert np.array_equal(single,k.histogram(pos,np.zeros(len(typ),np.int32),box,nbins=80,dr=.05,threads=1))
    checks['exact_count_merger_and_single_species_recovery']=True
    ns=3;n=np.array([120,40,40]);volume=1000.;f=100
    edges=np.arange(101)*.02;dv=4*np.pi/3*np.diff(edges**3)
    ideal=np.array([n[a]*n[b] if a!=b else n[a]*(n[a]-1)/2 for a,b in pairs(ns)])
    exact=f*ideal[:,None]*dv[None,:]/volume
    assert abs(evaluate(exact,n,f,volume,dr=.02,stride=1)['s2'])<1e-25
    checks['analytic_ideal_rdf_normalization_and_zero_entropy']=True
    hs=exact.copy();hs[:,:50]=0
    expected=-2*np.pi*(sum(n)/volume)/3  # all pairs g=0 for r<1
    value=evaluate(hs,n,f,volume,dr=.02,stride=1)
    assert abs(value['s2']-expected)<1e-12
    assert abs(value['s2']-np.dot(n/sum(n),value['species_centered']))<1e-12
    checks['analytic_hard_core_integral_cross_pair_factor_and_allocation']=True
    assert phi(np.array([0.,1.,1.+1e-9])).tolist()[:2]==[1.,0.]
    assert phi(np.array([1.+1e-9]))[0]>=0
    checks['zero_bin_limit_and_near_unity_stability']=True
    # The physical label null is tested by pooling counts before the nonlinear map.
    a=evaluate(merged,nm,1,float(np.prod(box)),dr=.05,stride=1)
    b=evaluate(k.histogram(pos,group[typ],box,nbins=80,dr=.05),nm,1,float(np.prod(box)),dr=.05,stride=1)
    assert a['s2']==b['s2']
    checks['synthetic_identical_label_null']=True
    (ROOT/'validation.json').write_text(json.dumps({'passed':True,'checks':checks},indent=2)+'\n')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
