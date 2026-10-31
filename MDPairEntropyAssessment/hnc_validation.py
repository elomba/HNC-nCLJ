#!/usr/bin/env python3
"""Check HNC s2 exports and compare selected HNC/MD probe-solvent RDFs."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pair_entropy import ROOT, STUDY, merge_map, merge_hist, evaluate

SELECT = ['B3', 'U4', 'H5', 'U6', 'H7']

def hnc_integral(name, species):
    path=STUDY/'hnc'/'fine'/name/'RDFhncr.dat'
    data=np.loadtxt(path)
    r=data[:,0]; g=np.maximum(data[:,1:],0.0)
    pairs=[(a,b) for a in range(len(species)) for b in range(a+1)]
    assert g.shape[1] == len(pairs)
    total=0.0
    for k,(a,b) in enumerate(pairs):
        y=g[:,k]
        f=np.ones_like(y); good=y>0
        f[good]=y[good]*np.log(y[good])-y[good]+1
        # RDF starts at r=0.55; the omitted repulsive core has g~0.
        integral=trapezoid(f*r*r,r)+r[0]**3/3
        total-=2*np.pi*.7*species[a]['x']*species[b]['x']*(1 if a==b else 2)*integral
    return data,total

def main():
    compositions={x['name']:x for x in json.loads((STUDY/'selected.json').read_text())}
    hnc={x['name']:x for x in json.loads((STUDY/'hnc_fine_summary.json').read_text())}
    names={c:n for c,n in __import__('pair_entropy').ORDER}
    validation={}; curves={}
    for code,name in names.items():
        c=compositions[name]; rdf,s2=hnc_integral(name,c['species'])
        stored=float(hnc[name]['s2'])
        validation[code]={'hnc_name':name,'stored_s2':stored,'rdf_integrated_s2':s2,
                          'difference':s2-stored,'rdf_file':str(STUDY/'hnc'/'fine'/name/'RDFhncr.dat')}
        if code not in SELECT: continue
        # Pair (probe +2, solvent eps=1) is lower-triangle column (1,0).
        species=c['species']; probe=next(i for i,x in enumerate(species) if x['label']=='probe_plus2')
        solvent=next(i for i,x in enumerate(species) if x['label']=='solvent')
        pair_index=probe*(probe+1)//2+solvent
        rh=rdf[:,0];gh=np.maximum(rdf[:,1+pair_index],0)
        pooled=None; counts=None; nframes=0; volume=None
        for rep in [1,2,3]:
            with np.load(ROOT/'cache'/f'{code}_r{rep}.npz') as z:
                meta=json.loads((ROOT/'cache'/f'{code}_r{rep}.json').read_text())
                group,_=merge_map(meta['species'])
                h,nn=merge_hist(z['counts'][1001:2001].sum(axis=0),z['particle_counts'],group)
                pooled=h if pooled is None else pooled+h
                if counts is None: counts=nn
                else: assert np.array_equal(counts,nn)
                nframes+=1000
                box=z['boxes'][-1]
                volume=float(np.prod(box))
        merged_species=[]
        seen={}
        for s in species:
            key=tuple(s[k] for k in ['q','eps','mass','sigma'])
            if key not in seen: seen[key]=len(merged_species);merged_species.append(s)
        group=[seen[tuple(s[k] for k in ['q','eps','mass','sigma'])] for s in species]
        result=evaluate(pooled,counts,nframes,volume,dr=.01,stride=1,rmax=8.88)
        labels=[x['label'] for x in merged_species]
        pm={tuple(sorted((labels[a],labels[b]))):k for k,(a,b) in enumerate([(a,b) for a in range(len(labels)) for b in range(a,len(labels))])}
        md_idx=pm[tuple(sorted((merged_species[group[probe]]['label'],merged_species[group[solvent]]['label'])))]
        mdg=result['g'][md_idx]
        mdr=(result['edges'][:-1]+result['edges'][1:])/2
        curves[code]={'hnc_r':rh.tolist(),'hnc_g':gh.tolist(),'md_r':mdr.tolist(),'md_g':mdg.tolist(),
                      'pair':'probe +2 / solvent (epsilon=1)'}
    out=ROOT/'results';out.mkdir(exist_ok=True)
    (out/'hnc_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    (out/'hnc_md_probe_rdf.json').write_text(json.dumps(curves,indent=2)+'\n')
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,code in zip(axs,['B3','H5']):
        d=curves[code]
        ax.plot(d['hnc_r'],d['hnc_g'],label='HNC',lw=1.5)
        ax.plot(d['md_r'],d['md_g'],label='MD (3-replica pooled)',lw=1.2)
        ax.set_xlim(.6,4);ax.set_ylim(bottom=0);ax.set_xlabel(r'$r/\sigma$');ax.set_ylabel(r'$g_{+2,solvent}(r)$')
        ax.set_title(code);ax.legend()
    fig.savefig(out/'hnc_md_probe_rdf.pdf');fig.savefig(out/'hnc_md_probe_rdf.png',dpi=180);plt.close(fig)
    print(json.dumps({'max_abs_hnc_s2_export_error':max(abs(x['difference']) for x in validation.values()),
                      'validations':validation,'rdf_plot':str(out/'hnc_md_probe_rdf.pdf')},indent=2))

if __name__=='__main__': main()
