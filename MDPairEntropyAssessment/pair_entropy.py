"""Pair bookkeeping, finite-count RDF normalization, and entropy allocations."""
from pathlib import Path
import ctypes
import hashlib
import json
import os
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parent
STUDY = ROOT.parent / 'ComponentEntropyStudy_20261030'
ORDER = [('B3','baseline_3C'),('B5','mono_f010'),('U4','weak_f030_e025'),
         ('L5','dispersion_split_d000'),('M5','dispersion_split_d025'),
         ('H5','dispersion_split_d040'),('U6','mono_weak_uniform6C'),('H7','mono_weak_split7C')]
DR = .01
NBINS = 888

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def pairs(ns):
    return [(a,b) for a in range(ns) for b in range(a,ns)]

def mapping(ns):
    m=np.empty((ns,ns),np.int32)
    for k,(a,b) in enumerate(pairs(ns)): m[a,b]=m[b,a]=k
    return m

def build_kernel():
    out=ROOT/'build'/'histogram.so'
    out.parent.mkdir(exist_ok=True)
    src=ROOT/'histogram.c'
    if not out.exists() or out.stat().st_mtime<src.stat().st_mtime:
        subprocess.run(['gcc','-std=c11','-O3','-fopenmp','-fPIC','-shared',str(src),'-lm','-o',str(out)],check=True)
    return out

class Kernel:
    def __init__(self):
        self.lib=ctypes.CDLL(str(build_kernel()))
        self.call=self.lib.pair_histogram
        self.call.argtypes=[np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS'),ctypes.c_int,
            np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),
            np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS'),ctypes.c_int,
            ctypes.c_int,ctypes.c_double,ctypes.c_int,
            np.ctypeslib.ndpointer(np.uint32,flags='C_CONTIGUOUS')]
        self.call.restype=ctypes.c_int
    def histogram(self,pos,types,box,ns=None,nbins=NBINS,dr=DR,threads=None):
        pos=np.ascontiguousarray(pos,dtype=np.float64)
        box=np.ascontiguousarray(box,dtype=np.float64)
        pos=np.ascontiguousarray(np.mod(pos,box))
        types=np.ascontiguousarray(types,dtype=np.int32)
        ns=int(types.max()+1) if ns is None else ns
        assert np.isfinite(pos).all() and np.all(box>0)
        assert types.shape==(len(pos),) and np.all((types>=0)&(types<ns))
        assert nbins*dr<=box.min()*.5
        h=np.empty((ns*(ns+1)//2,nbins),np.uint32)
        err=self.call(pos,types,len(pos),box,mapping(ns),ns,nbins,dr,
                      threads or int(os.environ.get('OMP_NUM_THREADS','1')),h)
        if err: raise RuntimeError(f'histogram kernel returned {err}')
        return h

def compositions():
    return {c['name']:c for c in json.loads((STUDY/'selected.json').read_text())}

def merge_map(species):
    keys=[]; labels=[]; group=[]
    for s in species:
        key=tuple(s[k] for k in ['q','eps','mass','sigma'])
        if key not in keys: keys.append(key); labels.append(s['label'])
        group.append(keys.index(key))
    return np.asarray(group,int),labels

def merge_hist(hist,counts,group):
    ns=len(counts); ng=int(max(group)+1)
    out=np.zeros((*hist.shape[:-2],ng*(ng+1)//2,hist.shape[-1]),hist.dtype)
    mp=mapping(ng)
    for k,(a,b) in enumerate(pairs(ns)):
        out[...,mp[group[a],group[b]],:]+=hist[...,k,:]
    return out,np.bincount(group,weights=counts,minlength=ng).astype(int)

def phi(g):
    # log1p is stable near g=1; define 0*log(0)=0 exactly.
    g=np.asarray(g,float)
    out=np.ones_like(g)
    good=g>0
    out[good]=g[good]*np.log(g[good])-g[good]+1
    near=(g>.95)&(g<1.05)
    h=g[near]-1
    out[near]=(1+h)*np.log1p(h)-h
    assert out.min()>-1e-14
    return np.maximum(out,0)

def evaluate(hist,counts,nframes,volume,dr=DR,stride=2,rmax=None,normalization='ideal'):
    counts=np.asarray(counts,int); ns=len(counts); n=counts.sum();rho=n/volume
    h=np.asarray(hist,float)
    assert h.ndim==2 and h.shape[0]==ns*(ns+1)//2
    assert h.shape[-1]%stride==0
    h=h.reshape(len(h),-1,stride).sum(-1)
    edges=np.arange(h.shape[-1]+1)*dr*stride
    if rmax is not None:
        keep=np.flatnonzero(edges[1:]<=rmax+1e-10)
        h=h[:,keep];edges=edges[:len(keep)+1]
    dv=4*np.pi/3*np.diff(edges**3)
    ideal=np.array([counts[a]*counts[b] if a!=b else counts[a]*(counts[a]-1)/2 for a,b in pairs(ns)])
    if normalization=='bulk':
        ideal=np.array([counts[a]*counts[b] if a!=b else counts[a]**2/2 for a,b in pairs(ns)])
    g=h/(nframes*ideal[:,None]*dv[None,:]/volume)
    integrals=(phi(g)*np.diff(edges**3)[None,:]/3).sum(-1)
    x=counts/n
    weighted=np.array([-2*np.pi*rho*x[a]*x[b]*(1 if a==b else 2)*integrals[k]
                       for k,(a,b) in enumerate(pairs(ns))])
    centered=np.zeros(ns)
    for k,(a,b) in enumerate(pairs(ns)):
        centered[a]+=-2*np.pi*rho*x[b]*integrals[k]
        if a!=b: centered[b]+=-2*np.pi*rho*x[a]*integrals[k]
    assert np.isclose(weighted.sum(),x@centered,atol=1e-12)
    assert np.all(weighted<=1e-12)
    return dict(s2=float(weighted.sum()),pair_contributions=weighted,
                species_centered=centered,g=g,edges=edges,integrals=integrals)

def scalar_metrics(result,counts,charges):
    x=np.asarray(counts)/sum(counts); q=np.asarray(charges)
    p=result['pair_contributions']; cen=result['species_centered']
    div=np.abs(q)==2
    groups={'neutral_neutral':0.,'neutral_ion':0.,'ion_ion':0.,'divalent_union':0.,'monovalent_union':0.}
    for k,(a,b) in enumerate(pairs(len(q))):
        category='neutral_neutral' if q[a]==q[b]==0 else ('ion_ion' if q[a]!=0 and q[b]!=0 else 'neutral_ion')
        groups[category]+=float(p[k])
        if div[a] or div[b]: groups['divalent_union']+=float(p[k])
        if abs(q[a])==1 or abs(q[b])==1: groups['monovalent_union']+=float(p[k])
    assert np.isclose(sum(groups[k] for k in ['neutral_neutral','neutral_ion','ion_ion']),result['s2'])
    return dict(s2=result['s2'],probe_centered_mean=float(cen[div].mean()),
                probe_weighted_centered=float(x[div]@cen[div]),
                plus2_centered=float(cen[q==2][0]),minus2_centered=float(cen[q==-2][0]),**groups)
