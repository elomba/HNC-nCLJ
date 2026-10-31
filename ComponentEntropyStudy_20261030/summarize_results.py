#!/usr/bin/env python3
"""Summarize completed replicas without rereading or changing MD trajectories."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
ORDER = [('B3','baseline_3C'),('B5','mono_f010'),('U4','weak_f030_e025'),
         ('L5','dispersion_split_d000'),('M5','dispersion_split_d025'),
         ('H5','dispersion_split_d040'),('U6','mono_weak_uniform6C'),('H7','mono_weak_split7C')]

def estimate(v):
    v=np.asarray(v,float)
    m=v.mean(); h=t.ppf(.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))
    return dict(mean=float(m),CI=[float(m-h),float(m+h)])

def probe(row, window):
    tr=row['transport']; ids=[i+1 for i,k in enumerate(tr['D_species']) if k.startswith('probe_')]
    return float(np.mean(np.asarray(tr['fitting_windows'][window])[ids]))

def percent(e):
    return f"{100*np.expm1(e['mean']):+.2f}% [{100*np.expm1(e['CI'][0]):+.2f}, {100*np.expm1(e['CI'][1]):+.2f}]"

rows={}
for f in sorted((ROOT/'md').glob('*/replica_*/replica_result.json')):
    a=json.loads(f.read_text()); rows.setdefault(a['name'],{})[a['replica']]=a
assert len(rows)==8 and all(set(rs)=={1,2,3} for rs in rows.values())
summary=json.loads((OUT/'summary.json').read_text())
assert summary['complete']
assert len(list((ROOT/'md').glob('*/replica_*/*/verified.json')))==456

diag={}
for code,name in ORDER:
    rs=list(rows[name].values())
    diag[code]={'probe_MSD_windows':{w:estimate([probe(r,w) for r in rs]) for w in ['10_50','20_100','50_200']},
       'probe_GK':{str(c):estimate([np.mean([v[str(c)] for k,v in r['transport']['Green_Kubo'].items() if k.startswith('probe_')]) for r in rs]) for c in [5,10]},
       'sign_difference':estimate([r['transport']['D_species']['probe_plus2']-r['transport']['D_species']['probe_minus2'] for r in rs]),
       'max_abs_pchip_shift':max(abs(r['s_ex_pchip']-r['s_ex']) for r in rs),
       'max_abs_averaging_window_shift':max(abs(s-r['s_ex']) for r in rs for s in r['averaging_window_entropy']),
       'concentration_fluctuations':{k:estimate([r['transport']['concentration_fluctuations'][k] for r in rs]) for k in rs[0]['transport']['concentration_fluctuations']},
       'position_frames':[r['transport']['n_position_frames'] for r in rs],
       'velocity_frames':[r['transport']['n_velocity_frames'] for r in rs]}
comparisons=[('H5/U4','dispersion_split_d040','weak_f030_e025'),('H7/U6','mono_weak_split7C','mono_weak_uniform6C')]
window_contrasts={}
for label,a,b in comparisons:
    window_contrasts[label]={w:estimate([np.log(probe(rows[a][r],w)/probe(rows[b][r],w)) for r in [1,2,3]]) for w in ['10_50','20_100','50_200']}
fits={}
for family,names in [('no_monovalent',['baseline_3C','weak_f030_e025','dispersion_split_d025','dispersion_split_d040']),('monovalent',['mono_f010','mono_weak_uniform6C','mono_weak_split7C'])]:
    fits[family]={}
    for key in ['D_net','D_probe']:
        x=np.array([summary['systems'][n]['s_ex']['mean'] for n in names]); y=np.log([summary['systems'][n][key]['mean']*.7**(1/3)/np.sqrt(2) for n in names])
        B,A=np.polyfit(x,y,1); r2=1-np.sum((y-A-B*x)**2)/np.sum((y-y.mean())**2)
        fits[family][key]=dict(A=float(A),B=float(B),R2=float(r2),n_states=len(names),note='descriptive means fit; no universal slope or independent validation inferred')
(OUT/'diagnostics.json').write_text(json.dumps(dict(systems=diag,window_contrasts=window_contrasts,descriptive_Rosenfeld_fits=fits),indent=2)+'\n')

lines=['# Completed HNC-guided MD study: conclusions','',
'All 24 Slurm allocations (100669–100692) completed with exit code 0; all 456 stages have verification records. Eight mixtures, three independent replicas each, N=4000, T*=2 and rho*=0.7. The LJ cutoff is fixed at 2.5 throughout the entropy integration. Reported entropy is excess entropy per particle relative to the ideal mixture of the same composition; ideal mixing entropy is excluded. Diffusion coefficients use LJ units, with a common whole-system COM frame.','',
'**Conclusion:** the selected interaction-diverse solvents raise excess entropy and divalent-ion mobility relative to the original baselines. The controlled diversification effect on the ions is smaller (about 5%) and remains provisional: its significance is sensitive to fitting window, and it is unresolved with monovalent background. Increasing the number of labels alone has exactly zero transport effect. The results support an entropy–mobility association within these model families, not a causal law based on component count.','',
'## State averages','',
'Composition definitions are in ../README.md. D_probe averages the +2 and −2 self-diffusion coefficients; D_net is the composition-weighted mean self-diffusion, not mutual diffusion or conductivity.','',
'| State | MD s_ex/kB | HNC s_ex/kB | D_net | D_probe |','|---|---:|---:|---:|---:|']
for code,n in ORDER:
    s=summary['systems'][n]; hn=rows[n][1]['HNC_s_ex']; lines.append(f"| {code} | {s['s_ex']['mean']:.5f} | {hn:.5f} | {s['D_net']['mean']:.5f} | {s['D_probe']['mean']:.5f} |")
lines+=['','## Paired comparisons','',
'Percent changes are exp(mean paired log ratio)−1; brackets are nominal 95% Student-t intervals across three independent replicas (2 degrees of freedom), without multiple-comparison correction. They do not include finite-size or estimator bias. Entropy differences are in kB per particle.','',
'| Comparison | Delta s_ex/kB | D_net change [95% CI] | D_probe change [95% CI] |','|---|---:|---:|---:|']
labels={'label_null':'L5 − U4: identical labels','moderate_split':'M5 − U4: moderate diversification','high_split':'H5 − U4: strong diversification','high_vs_baseline':'H5 − B3: full design change','mono_high_split':'H7 − U6: diversification with ±1','mono_high_vs_baseline':'H7 − B5: full design with ±1','salt_addition_baseline':'B5 − B3: add ±1 salt'}
for key,label in labels.items():
    a=summary['contrasts'][key]; s=a['delta_s_ex']; lines.append(f"| {label} | {s['mean']:+.5f} [{s['CI'][0]:+.5f}, {s['CI'][1]:+.5f}] | {percent(a['delta_log_D_net'])} | {percent(a['delta_log_D_probe'])} |")
lines+=['','## Interpretation in order of evidential strength','',
'1. **Component labels alone do not accelerate ions.** U4 and L5 have identical mechanical interactions, and their transport trajectories give exactly identical diffusion in every replica. Their small entropy difference (+0.00042, interval including zero) is consistent with thermodynamic sampling. Formal ideal mixing rises by 0.30 ln(2)=0.20794 kB if labels are counted as distinguishable; that contribution is absent from s_ex and produces no force change.',
'2. **HNC successfully selected higher-excess-entropy liquids.** MD reproduces the selected-state entropy ordering. HNC predicts H5−U4=+0.1107 versus MD +0.1032, and H7−U6=+0.1108 versus MD +0.1033. Absolute HNC entropies are more negative by roughly 0.064–0.093 kB; agreement of increments is stronger than agreement of absolute values.',
'3. **Large design changes reproducibly accelerate both ion signs.** H5 versus B3 increases pooled divalent diffusion by about 20%, and H7 versus B5 by about 18%. At the primary window, the separate +2 and −2 changes are +16.8% and +23.4% for H5/B3 and +18.0% and +17.6% for H7/B5; all four paired nominal 95% intervals exclude zero. These interventions also change mean interaction strength and local solvent friction, so they cannot isolate entropy as the cause.',
'4. **Controlled physical diversification has a smaller ion response than the mean response.** U4→H5 preserves composition of charged particles, masses, nominal diameters, ionic strength and the mean bare LJ pair amplitude. It raises s_ex by 0.1032 kB, mean diffusion by 9.47%, and ion diffusion by 4.81%. The primary-window ion interval barely excludes zero. With ±1 background the analogous means are +0.1033 kB, +8.37%, and +5.20%, but the ion interval includes zero. The moderate split raises entropy clearly but its ~1.5% ion response is unresolved.',
'5. **The original monovalent-addition effect is not established.** B3→B5 slightly lowers excess entropy (−0.00320 kB) and gives no resolved diffusion change. More charged components alone therefore do not ensure higher excess entropy or faster divalent ions.',
'6. **Rosenfeld scaling is descriptive, not universal here.** At common T, rho and mass the reduced diffusion D^R=D rho^(1/3)/sqrt(T) differs from D by a constant. The broad state ordering is consistent with a positive ln(D^R)–s_ex trend for both mean and probe diffusion. The controlled H5/U4 local response is B≈0.88 for the mean and B≈0.46 for the probes, substantially below the prior provisional B≈1.4 forecast; H7/U6 gives about 0.78 and 0.49. A single bulk entropy does not determine species mobility. Descriptive fits of state means give B=0.80 (mean) and 0.62 (probes) without monovalent background, and 0.78 and 0.55 with it; R² exceeds 0.99 in all four fits. These three- or four-state fits are dominated by the broad interaction changes and do not validate a universal relation. See diagnostics.json for the fit parameters.','',
'## Convergence and limitations','',
'### MSD fitting-window sensitivity','',
'| Comparison | 10–50 tau | 20–100 tau (primary) | 50–200 tau |','|---|---:|---:|---:|']
for label,d in window_contrasts.items():lines.append('| '+label+' | '+' | '.join(percent(d[w]) for w in ['10_50','20_100','50_200'])+' |')
lines+=['',
'The positive ~4–6% ion response survives in point estimates, but H5/U4 loses nominal significance in the later window. This is evidence for a trend, not a converged precise 5% effect. VACF integrals at 5 and 10 tau support the broad mobility ordering, but the final-100-tau velocity records are too short for a precise independent small-effect determination. No reliable 20-tau plateau is claimed.','',
'### Charge-sign symmetry','',
'The sign-symmetric B3 Hamiltonian should have equal long-time +2 and −2 self-diffusion. The observed difference D(+2)−D(−2)=0.01039, nominal replica interval [0.00291,0.01786], flags remaining sampling/uncertainty limitations. The other states have sign-difference intervals containing zero. Inspection of generated inputs found equal masses/LJ parameters for both probes, random probe assignments and one all-particle thermostat/common COM frame; no sign-dependent input parameter was found. This does not establish that the residual asymmetry is harmless. Pooling both signs is physically appropriate but does not make the convergence flag disappear. Both signs independently speed up strongly in the full baseline contrasts.','',
'### Entropy and homogeneity diagnostics','']
lines.append(f"Across all replicas, the largest change from cubic-spline to PCHIP integration is {max(x['max_abs_pchip_shift'] for x in diag.values()):.5f} kB; the largest change when retaining all/75% versus 50% of production is {max(x['max_abs_averaging_window_shift'] for x in diag.values()):.5f} kB. These are small compared with the ~0.10 kB diversification increments. They do not test the unsampled dilute-density endpoint or finite-size bias.")
lines+=['',
'Low-k concentration fluctuation diagnostics remain of order unity (individual species values approximately 0.48–1.80), without a large fluctuation growth across these runs. One box size and short density-stage trajectories cannot establish single-phase stability or exclude slower segregation. Full species diagnostics are saved in diagnostics.json.','',
'## Most useful next checks','',
'1. Add independent, longer matched U4/H5 and U6/H7 replicas to resolve the ~5% effect; extend B3 to test the sign-symmetry flag. Keep the same estimators and prespecified windows. More overlapping time origins are not extra independent replicas.',
'2. Repeat at a larger N at fixed density/composition (e.g. 8000 or 16000), matching all cutoffs and PPPM accuracy. Compare both signs and solvent mobility; the expected finite-size transport correction need not cancel when viscosity differs.',
'3. From existing trajectories, measure probe–solvent/ion RDFs, association populations and residence times, and species-resolved pair excess entropy. These may distinguish weaker solvation/caging from a bulk entropy proxy. A pair-entropy approximation must not be relabeled as total thermodynamic excess entropy.',
'4. After convergence checks, raise divalent concentration with matched uniform/diverse solvents at each concentration. Re-screen with HNC; maintain neutrality and compare within the same ionic strength at each matched pair. Higher concentration can change pairing and introduce species-specific decoupling rather than strengthen the present trend.',
'5. Test several densities/temperatures within each fixed composition to separate state-point Rosenfeld scaling from composition-induced trends. Varying nominal particle size is a distinct perturbation of packing and friction and should be a separate controlled series.','',
'No new MD jobs were submitted during this analysis. All raw trajectories and original automatic results are preserved. Reproduce this report and plots with `python3 ComponentEntropyStudy_20261030/summarize_results.py`.']
(OUT/'conclusions.md').write_text('\n'.join(lines)+'\n')

fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
colors={'B3':'#2166ac','U4':'#2166ac','M5':'#2166ac','H5':'#2166ac','B5':'#b2182b','U6':'#b2182b','H7':'#b2182b'}
for ax,key,title in zip(axs,['D_net','D_probe'],['Mean self-diffusion','Divalent-ion self-diffusion']):
    for code,n in ORDER:
        if code=='L5':continue
        s=summary['systems'][n];x=s['s_ex']['mean'];y=s[key]['mean'];xe=(s['s_ex']['CI'][1]-s['s_ex']['CI'][0])/2;ye=(s[key]['CI'][1]-s[key]['CI'][0])/2
        ax.errorbar(x,y,xerr=xe,yerr=ye,fmt='o',color=colors[code],capsize=3)
        ax.annotate(code,(x,y),xytext=((-27,-13) if code=='B5' else (5,5 if code not in ['U6','H7'] else -13)),textcoords='offset points',fontsize=9)
    ax.set(xlabel=r'$s_{\rm ex}/k_B$',ylabel=r'$D$ (LJ units)',title=title)
    ax.set_yscale('log');ax.margins(x=.10,y=.12);ax.grid(alpha=.2)
fig.suptitle('HNC-guided MD: state averages and nominal 95% replica intervals',fontsize=11)
for ext in ['png','pdf']:fig.savefig(OUT/f'entropy_diffusion.{ext}',dpi=200)
plt.close(fig)
print('Wrote conclusions.md, diagnostics.json and entropy_diffusion.png/pdf')
