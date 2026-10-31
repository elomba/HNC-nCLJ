#!/usr/bin/env python3
"""Joint replica-level TI, HNC, MD pair-entropy, and transport assessment."""
import csv,json,itertools
from pathlib import Path
import numpy as np
from scipy.stats import t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pair_entropy import ROOT,STUDY,ORDER,digest

CONTRASTS=[('label_null','L5','U4'),('moderate_split','M5','U4'),
    ('high_split','H5','U4'),('mono_high_split','H7','U6'),
    ('high_vs_baseline','H5','B3'),('mono_high_vs_baseline','H7','B5'),('salt_addition','B5','B3')]

def estimate(v):
    v=np.asarray(v,float);n=len(v);mean=v.mean()
    hw=t.ppf(.975,n-1)*v.std(ddof=1)/np.sqrt(n)
    return {'mean':float(mean),'CI':[float(mean-hw),float(mean+hw)],'n':n}

def window_probe(row,w):
    tr=row['transport'];ids=[j+1 for j,k in enumerate(tr['D_species']) if k.startswith('probe_')]
    return float(np.mean(np.asarray(tr['fitting_windows'][w])[ids]))

def diffusion_value(row,target):
    if target=='D_plus': return float(row['transport']['D_species']['probe_plus2'])
    if target=='D_minus': return float(row['transport']['D_species']['probe_minus2'])
    return float(row['transport'][target])

def fit(x,y):
    x=np.array(x,float);y=np.array(y,float);slope,intercept=np.polyfit(x,y,1)
    residual=y-intercept-slope*x
    predictions=[]
    for i in range(len(x)):
        take=np.arange(len(x))!=i;b,a=np.polyfit(x[take],y[take],1)
        predictions.append(y[i]-a-b*x[i])
    return {'slope':float(slope),'intercept':float(intercept),
        'R2':float(1-np.sum(residual**2)/np.sum((y-y.mean())**2)),
        'leave_one_state_out_RMSE_lnD':float(np.sqrt(np.mean(np.array(predictions)**2))),
        'n_states':len(x),'residuals_lnD':residual.tolist()}

def main():
    rows={};ideals={}
    for code,name in ORDER:
        rows[code]=[];ideals[code]=[]
        for rep in [1,2,3]:
            rows[code].append(json.loads((ROOT/'replicas'/f'{code}_r{rep}.json').read_text()))
            ideals[code].append(json.loads((ROOT/'ideal_results'/f'{code}_r{rep}_ideal.json').read_text()))
            for rr,tag in [(rows[code][-1],f'{code}_r{rep}'),(ideals[code][-1],f'{code}_r{rep}_ideal')]:
                fine=json.loads((ROOT/'refined'/(tag+'.json')).read_text())
                assert fine['original_core_histograms_match_exactly']
                rr['refinement']=fine
                rr['resolved_metrics']=fine['jackknife20_hybrid']['0.005']
    # Reject stale results after estimator changes, rather than silently mixing versions.
    for d in [rows,ideals]:
        assert all(r['analysis_sha256']==digest(ROOT/'analyze.py') for rs in d.values() for r in rs)
        assert all(r['estimator_sha256']==digest(ROOT/'pair_entropy.py') for rs in d.values() for r in rs)
    hnc={x['name']:x for x in json.loads((STUDY/'hnc_fine_summary.json').read_text())}
    metrics=list(rows['B3'][0]['primary'])
    result={'scope':'Complete all-pair MD RDF extraction; independent TI and ideal controls; no entropy fit to TI or D.',
        'complete_MD':24,'complete_ideal':24,'protocol_sha256':digest(ROOT/'protocol.json'),
        'statistical_note':'Nominal paired 95% Student-t intervals across 3 independent replicas (df=2); no multiplicity correction or finite-size bias included.',
        'primary':'Raw ideal-normalized pooled RDF, dr=.02, R=8.88, last 1000 frames.',
        'jackknife':'20 equal time blocks, 25 tau per block; sensitivity, not selected to match TI.',
        'systems':{},'contrasts':{},'fits':{},'diagnostics':{}}
    for code,name in ORDER:
        rs=rows[code];ii=ideals[code]
        d={'s2_raw':estimate([r['primary']['s2'] for r in rs]),
            's2_jackknife20':estimate([r['jackknife']['20']['corrected']['s2'] for r in rs]),
            's2_resolved':estimate([r['resolved_metrics']['s2'] for r in rs]),
            's2_resolved_finest':estimate([r['refinement']['jackknife20_hybrid']['0.0025']['s2'] for r in rs]),
            'remainder_resolved':estimate([r['TI_s_ex']-r['resolved_metrics']['s2'] for r in rs]),
            'pair_fraction_resolved':estimate([r['resolved_metrics']['s2']/r['TI_s_ex'] for r in rs]),
            'resolved_metrics':{k:estimate([r['resolved_metrics'][k] for r in rs]) for k in metrics},
            'TI_s_ex':estimate([r['TI_s_ex'] for r in rs]),
            'remainder_raw':estimate([r['remainder_raw'] for r in rs]),
            'remainder_jackknife20':estimate([r['remainder_jackknife20'] for r in rs]),
            'pair_fraction_raw':estimate([r['pair_fraction_raw'] for r in rs]),
            'pair_fraction_jackknife20':estimate([r['jackknife']['20']['corrected']['s2']/r['TI_s_ex'] for r in rs]),
            'D_net':estimate([r['transport']['D_net'] for r in rs]),
            'D_probe':estimate([r['transport']['D_probe'] for r in rs]),
            'D_plus':estimate([diffusion_value(r,'D_plus') for r in rs]),
            'D_minus':estimate([diffusion_value(r,'D_minus') for r in rs]),
            'HNC_s2':hnc[name]['s2'],'HNC_s_ex':hnc[name]['s_ex'],
            'MD_s2_minus_HNC':estimate([r['primary']['s2']-hnc[name]['s2'] for r in rs]),
            'metrics':{k:estimate([r['primary'][k] for r in rs]) for k in metrics},
            'corrected_metrics':{k:estimate([r['jackknife']['20']['corrected'][k] for r in rs]) for k in metrics},
            'ideal_raw':estimate([r['primary']['s2'] for r in ii]),
            'ideal_jackknife20':estimate([r['jackknife']['20']['corrected']['s2'] for r in ii]),
            'ideal_resolved':estimate([r['resolved_metrics']['s2'] for r in ii]),
            'ideal_resolved_finest':estimate([r['refinement']['jackknife20_hybrid']['0.0025']['s2'] for r in ii]),
            'refined_core_bins':{w:estimate([r['refinement']['jackknife20_hybrid'][w]['s2'] for r in rs]) for w in rs[0]['refinement']['jackknife20_hybrid']},
            'radial':{radius:estimate([r['radial_sensitivity'][radius]['s2'] for r in rs]) for radius in rs[0]['radial_sensitivity']},
            'bins_raw':{width:estimate([r['bin_sensitivity'][width]['s2'] for r in rs]) for width in rs[0]['bin_sensitivity']},
            'bins_corrected':{width:estimate([r['jackknife_bin_sensitivity'][width]['s2'] for r in rs]) for width in rs[0]['jackknife_bin_sensitivity']},
            'samples':{n:estimate([r['sample_sensitivity'][n]['s2'] for r in rs]) for n in rs[0]['sample_sensitivity']},
            'time_windows':{w:estimate([r['time_windows'][w]['s2'] for r in rs]) for w in rs[0]['time_windows']},
            'normalization_bulk_minus_ideal':estimate([r['normalization_sensitivity']['bulk']['s2']-r['primary']['s2'] for r in rs]),
            'jackknife_block_sensitivity':{b:estimate([r['jackknife'][b]['corrected']['s2'] for r in rs]) for b in ['10','20','40']},
            'max_coordination_integrated_tau':max(r['correlation_diagnostics']['max_tau'] for r in rs)}
        result['systems'][code]=d
    for key,a,b in CONTRASTS:
        aa=rows[a];bb=rows[b]
        d={'a':a,'b':b,'delta_TI':estimate([x['TI_s_ex']-y['TI_s_ex'] for x,y in zip(aa,bb)]),
           'delta_s2_raw':estimate([x['primary']['s2']-y['primary']['s2'] for x,y in zip(aa,bb)]),
           'delta_s2_corrected':estimate([x['jackknife']['20']['corrected']['s2']-y['jackknife']['20']['corrected']['s2'] for x,y in zip(aa,bb)]),
           'delta_s2_resolved':estimate([x['resolved_metrics']['s2']-y['resolved_metrics']['s2'] for x,y in zip(aa,bb)]),
           'delta_remainder_resolved':estimate([(x['TI_s_ex']-x['resolved_metrics']['s2'])-(y['TI_s_ex']-y['resolved_metrics']['s2']) for x,y in zip(aa,bb)]),
           'resolved_metric_contrasts':{m:estimate([x['resolved_metrics'][m]-y['resolved_metrics'][m] for x,y in zip(aa,bb)]) for m in metrics},
           'resolved_bin_sensitivity':{w:estimate([x['refinement']['jackknife20_hybrid'][w]['s2']-y['refinement']['jackknife20_hybrid'][w]['s2'] for x,y in zip(aa,bb)]) for w in aa[0]['refinement']['jackknife20_hybrid']},
           'delta_remainder_corrected':estimate([x['remainder_jackknife20']-y['remainder_jackknife20'] for x,y in zip(aa,bb)]),
           'delta_log_D_probe':estimate([np.log(x['transport']['D_probe']/y['transport']['D_probe']) for x,y in zip(aa,bb)]),
           'delta_log_D_net':estimate([np.log(x['transport']['D_net']/y['transport']['D_net']) for x,y in zip(aa,bb)]),
           'delta_log_D_plus':estimate([np.log(diffusion_value(x,'D_plus')/diffusion_value(y,'D_plus')) for x,y in zip(aa,bb)]),
           'delta_log_D_minus':estimate([np.log(diffusion_value(x,'D_minus')/diffusion_value(y,'D_minus')) for x,y in zip(aa,bb)]),
           'metric_contrasts':{m:estimate([x['primary'][m]-y['primary'][m] for x,y in zip(aa,bb)]) for m in metrics},
           'bin_sensitivity_corrected':{w:estimate([x['jackknife_bin_sensitivity'][w]['s2']-y['jackknife_bin_sensitivity'][w]['s2'] for x,y in zip(aa,bb)]) for w in aa[0]['jackknife_bin_sensitivity']},
           'radial_sensitivity':{r:estimate([x['radial_sensitivity'][r]['s2']-y['radial_sensitivity'][r]['s2'] for x,y in zip(aa,bb)]) for r in aa[0]['radial_sensitivity']},
           'sample_sensitivity':{n:estimate([x['sample_sensitivity'][n]['s2']-y['sample_sensitivity'][n]['s2'] for x,y in zip(aa,bb)]) for n in aa[0]['sample_sensitivity']},
           'time_window_sensitivity':{w:estimate([x['time_windows'][w]['s2']-y['time_windows'][w]['s2'] for x,y in zip(aa,bb)]) for w in aa[0]['time_windows']},
           'normalization_bulk_minus_ideal':estimate([(x['normalization_sensitivity']['bulk']['s2']-y['normalization_sensitivity']['bulk']['s2'])-(x['primary']['s2']-y['primary']['s2']) for x,y in zip(aa,bb)])}
        d['pair_fraction_of_TI_increment']=d['delta_s2_resolved']['mean']/d['delta_TI']['mean'] if key!='label_null' else None
        result['contrasts'][key]=d
    groups={'no_monovalent':['B3','U4','M5','H5'],'monovalent':['B5','U6','H7']}
    def coordinate(r,k):
        if k=='TI':return r['TI_s_ex']
        if k=='s2_corrected':return r['jackknife']['20']['corrected']['s2']
        if k=='s2_resolved':return r['resolved_metrics']['s2']
        if k=='probe_centered_resolved':return r['resolved_metrics']['probe_centered_mean']
        if k=='divalent_union_resolved':return r['resolved_metrics']['divalent_union']
        return r['primary'][k]
    for family,codes in groups.items():
        result['fits'][family]={}
        for target in ['D_net','D_probe','D_plus','D_minus']:
            ff={}
            for coord in ['TI','s2','s2_corrected','s2_resolved','probe_centered_mean','probe_centered_resolved','divalent_union','divalent_union_resolved']:
                x=[np.mean([coordinate(r,coord) for r in rows[c]]) for c in codes]
                y=[np.mean([np.log(diffusion_value(r,target)*.7**(1/3)/np.sqrt(2)) for r in rows[c]]) for c in codes]
                fit_result=fit(x,y)
                slopes=[]
                for excluded in range(3):
                    take=[i for i in range(3) if i!=excluded]
                    xx=[np.mean([coordinate(rows[c][i],coord) for i in take]) for c in codes]
                    yy=[np.mean([np.log(diffusion_value(rows[c][i],target)*.7**(1/3)/np.sqrt(2)) for i in take]) for c in codes]
                    slopes.append(fit(xx,yy)['slope'])
                fit_result['leave_one_replica_out_slopes']=slopes
                ff[coord]=fit_result
            result['fits'][family][target]=ff
    # Direct cache-level mechanical label test, beyond matching rounded entropy.
    label_checks=[]
    from pair_entropy import merge_map,merge_hist
    for rep in [1,2,3]:
        with np.load(ROOT/'cache'/f'U4_r{rep}.npz') as a,np.load(ROOT/'cache'/f'L5_r{rep}.npz') as b:
            ma=json.loads((ROOT/'cache'/f'L5_r{rep}.json').read_text())
            g,_=merge_map(ma['species']);hh,nn=merge_hist(b['counts'],b['particle_counts'],g)
            label_checks.append(bool(np.array_equal(a['counts'],hh) and np.array_equal(a['particle_counts'],nn)))
    result['diagnostics']['U4_L5_exact_frame_histograms']=label_checks
    result['diagnostics']['unmerged_label_minus_merged']=estimate([r['unmerged_label_sensitivity']['minus_merged'] for r in rows['L5']])
    out=ROOT/'results';out.mkdir(exist_ok=True)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    with (out/'states.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['state','TI_s_ex','MD_s2_raw','MD_s2_jackknife20','MD_s2_resolved','HNC_s2','residual_resolved','D_probe'])
        for c,_ in ORDER:
            d=result['systems'][c];w.writerow([c,d['TI_s_ex']['mean'],d['s2_raw']['mean'],d['s2_jackknife20']['mean'],d['s2_resolved']['mean'],d['HNC_s2'],d['remainder_resolved']['mean'],d['D_probe']['mean']])
    with (out/'contrasts.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['contrast','delta_TI','delta_s2_raw','delta_s2_corrected','delta_s2_resolved','resolved_CI_low','resolved_CI_high','delta_residual_resolved','fraction_of_TI_increment'])
        for k,_,_ in CONTRASTS:
            d=result['contrasts'][k];w.writerow([k,d['delta_TI']['mean'],d['delta_s2_raw']['mean'],d['delta_s2_corrected']['mean'],d['delta_s2_resolved']['mean'],*d['delta_s2_resolved']['CI'],d['delta_remainder_resolved']['mean'],d['pair_fraction_of_TI_increment']])
    plot(result,out,groups)
    print(json.dumps({'systems':{c:{k:v[k]['mean'] for k in ['TI_s_ex','s2_raw','s2_jackknife20','remainder_jackknife20']} for c,v in result['systems'].items()},'contrasts':{k:v['delta_s2_corrected'] for k,v in result['contrasts'].items()}},indent=2))

def plot(s,out,groups):
    plt.rcParams.update({'font.size':10})
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    codes=[c for c,_ in ORDER if c!='L5']
    for ax,kind in zip(axs,['s2_raw','s2_resolved']):
        for c in codes:
            d=s['systems'][c];x=d['TI_s_ex'];y=d[kind]
            ax.errorbar(x['mean'],y['mean'],xerr=[[x['mean']-x['CI'][0]],[x['CI'][1]-x['mean']]],yerr=[[y['mean']-y['CI'][0]],[y['CI'][1]-y['mean']]],fmt='o',color='tab:blue')
            ax.annotate(c,(x['mean'],y['mean']),xytext=(4,4),textcoords='offset points')
        ax.plot([-2.1,-1.5],[-2.1,-1.5],'k:',label='equality');ax.set_xlabel('TI excess entropy / $k_B$');ax.set_ylabel('MD pair entropy / $k_B$');ax.set_title(kind.replace('_',' '))
    fig.savefig(out/'entropy_benchmark.pdf');fig.savefig(out/'entropy_benchmark.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,target in zip(axs,['D_net','D_probe']):
        for family,codes in groups.items():
            color='tab:blue' if family=='no_monovalent' else 'tab:red'
            for c in codes:
                d=s['systems'][c];x=d['s2_resolved'];y=d[target];factor=.7**(1/3)/np.sqrt(2)
                ax.errorbar(x['mean'],y['mean']*factor,xerr=[[x['mean']-x['CI'][0]],[x['CI'][1]-x['mean']]],yerr=[[(y['mean']-y['CI'][0])*factor],[(y['CI'][1]-y['mean'])*factor]],fmt='o',color=color)
                ax.annotate(c,(x['mean'],y['mean']*factor),xytext=(3,3),textcoords='offset points')
        ax.set_yscale('log');ax.set_xlabel('MD pair entropy / $k_B$ (fine-core resolved)');ax.set_ylabel('Reduced '+target);ax.set_title(target)
    fig.savefig(out/'pair_entropy_diffusion.pdf');fig.savefig(out/'pair_entropy_diffusion.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,target in zip(axs,['D_plus','D_minus']):
        for family,codes in groups.items():
            color='tab:blue' if family=='no_monovalent' else 'tab:red'
            for c in codes:
                d=s['systems'][c];x=d['s2_resolved'];y=d[target];factor=.7**(1/3)/np.sqrt(2)
                ax.errorbar(x['mean'],y['mean']*factor,xerr=[[x['mean']-x['CI'][0]],[x['CI'][1]-x['mean']]],yerr=[[(y['mean']-y['CI'][0])*factor],[(y['CI'][1]-y['mean'])*factor]],fmt='o',color=color)
                offsets={'B3':(5,5 if target=='D_plus' else -12),'B5':(5,-13 if target=='D_plus' else 5),
                         'U4':(5,7),'U6':(5,-14),'H5':(5,7),'H7':(5,-14)}
                ax.annotate(c,(x['mean'],y['mean']*factor),xytext=offsets.get(c,(3,3)),textcoords='offset points')
        ax.set_yscale('log');ax.set_xlabel('MD pair entropy / $k_B$ (fine-core resolved)');ax.set_ylabel('Reduced '+target);ax.set_title(target.replace('_',' '))
    fig.savefig(out/'pair_entropy_species_diffusion.pdf');fig.savefig(out/'pair_entropy_species_diffusion.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for code in ['B3','U4','H5','U6','H7']:
        d=s['systems'][code]
        axs[0].plot([float(r) for r in d['radial']],[v['mean'] for v in d['radial'].values()],'o-',label=code)
        axs[1].plot([int(n) for n in d['samples']],[v['mean'] for v in d['samples'].values()],'o-',label=code)
    axs[0].set_xlabel(r'RDF integration radius / $\sigma$');axs[1].set_xlabel(r'Frames in final 500 $\tau$');axs[1].set_xscale('log')
    for ax in axs:ax.set_ylabel('Raw MD $s_2/k_B$');ax.legend()
    fig.savefig(out/'convergence.pdf');fig.savefig(out/'convergence.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
