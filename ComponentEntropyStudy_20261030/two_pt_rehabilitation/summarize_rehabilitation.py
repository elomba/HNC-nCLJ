#!/usr/bin/env python3
"""Aggregate independent-replica rehabilitation checks; never refit the model."""
from pathlib import Path
import csv,hashlib,json,platform
import numpy as np
import scipy
from scipy.stats import t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from two_pt import spectrum,evaluate

HERE=Path(__file__).resolve().parent
ORDER=[('B3','baseline_3C'),('B5','mono_f010'),('U4','weak_f030_e025'),('L5','dispersion_split_d000'),('M5','dispersion_split_d025'),('H5','dispersion_split_d040'),('U6','mono_weak_uniform6C'),('H7','mono_weak_split7C')]
PAIRS=[('label_null','L5/U4','dispersion_split_d000','weak_f030_e025'),('moderate','M5/U4','dispersion_split_d025','weak_f030_e025'),('controlled','H5/U4','dispersion_split_d040','weak_f030_e025'),('baseline','H5/B3','dispersion_split_d040','baseline_3C'),('mono_controlled','H7/U6','mono_weak_split7C','mono_weak_uniform6C'),('mono_baseline','H7/B5','mono_weak_split7C','mono_f010'),('salt','B5/B3','mono_f010','baseline_3C')]
PRIMARY='tau10_d1';ALTERNATE='tau10_d1.5'

def est(v):
    v=np.array(v,float);m=v.mean();h=t.ppf(.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))
    return dict(n=len(v),mean=float(m),CI=[float(m-h),float(m+h)])

def fmt(e,n=5):return f"{e['mean']:+.{n}f} [{e['CI'][0]:+.{n}f}, {e['CI'][1]:+.{n}f}]"

def sx(row,method=PRIMARY,segment='full'):return row['variants'][segment][method]['s_ex']

def contains_zero(e):return e['CI'][0]<=0<=e['CI'][1]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

rows={}
for f in (HERE/'replicas').glob('*.json'):
    a=json.loads(f.read_text());rows.setdefault(a['name'],{})[a['replica']]=a
assert set(rows)=={n for c,n in ORDER} and all(set(v)=={1,2,3} for v in rows.values())
for rs in rows.values():
    for r in rs.values():
        assert 'frequency_resolution_sensitivity' in r
        assert r['protocol_sha256']==digest(HERE/'protocol.json')
        assert r['code_sha256']['two_pt.py']==digest(HERE/'two_pt.py')
        assert r['code_sha256']['run_rehabilitation.py']==digest(HERE/'run_rehabilitation.py')

states={};global_h=[];global_mode=[];global_reconstruct=[];raw_modes=[];negative=[];species_D=[]
for code,name in ORDER:
    rs=[rows[name][r] for r in [1,2,3]]
    variants={k:est([sx(r,k) for r in rs]) for k in rs[0]['variants']['full']}
    bias=est([sx(r)-r['TI_s_ex'] for r in rs])
    pd=[];nd=[];fprobe=[]
    dynamics=[]
    for r in rs:
        v=r['variants']['full'][PRIMARY];sp=[g for g in v['species'].values() if g['primary']]
        pd.append(v['probe_spectral_D']/v['D_probe_MSD']-1)
        dnet=sum(g['x']*g['D_MSD'] for g in sp)
        nd.append(sum(g['x']*g['D_spectral_equipartition'] for g in sp)/dnet-1)
        fprobe.append(np.mean([g['f'] for g in sp if g['label'].startswith('probe_')]))
        meta=json.loads((HERE/'cache'/f'{name}_r{r["replica"]}.json').read_text());z=np.load(HERE/'cache'/f'{name}_r{r["replica"]}.npz')
        ps=[i for i,g in enumerate(meta['groups']) if g['primary'] and g['probe']]
        cb=np.mean([z[f'g{i}_full'] for i in ps],axis=0);tm=np.arange(len(cb))*meta['dt'];c=cb/(1-np.arange(len(cb))/len(cb))/cb[0];short=tm<=1.
        dynamics.append(dict(probe_VACF_min_0_1=float(c[short].min()),probe_negative_VACF_area_0_1=float(np.trapezoid(np.minimum(c[short],0),tm[short]))))
        for seg,d in r['variants'].items():
            for key,value in d.items():
                for g in value['species'].values():
                    global_h.append(g['reference_h_invariance_span']);global_mode.append(abs(g['normalized_modes']-3));global_reconstruct.append(g['spectrum_reconstruction_error'])
        for g in sp:
            raw_modes.append(g['raw_targetT_modes']);negative.append(g['negative_residual_fraction']);species_D.append(g['spectral_D_relative_error'])
    windowvals=[variants[w]['mean'] for w in ['tau5_d1','tau10_d1','tau20_d1']]
    states[code]=dict(name=name,TI=est([r['TI_s_ex'] for r in rs]),primary=variants[PRIMARY],alternate=variants[ALTERNATE],
        primary_bias_vs_TI=bias,alternate_bias_vs_TI=est([sx(r,ALTERNATE)-r['TI_s_ex'] for r in rs]),
        variants=variants,probe_D_relative_error=est(pd),net_D_relative_error=est(nd),
        probe_fluidicity=est(fprobe),probe_VACF_min_0_1=est([d['probe_VACF_min_0_1'] for d in dynamics]),
        probe_negative_VACF_area_0_1=est([d['probe_negative_VACF_area_0_1'] for d in dynamics]),
        window_span_tau5_10_20=float(np.ptp(windowvals)),
        half_record_sensitivity={seg:est([sx(r,segment=seg)-sx(r) for r in rs]) for seg in ['first_half','last_half']},
        partial_volume_sensitivity={k:est([r['volume_sensitivity'][k]-sx(r) for r in rs]) for k in rs[0]['volume_sensitivity']},
        frequency_resolution_sensitivity={k:est([r['frequency_resolution_sensitivity'][k]-sx(r) for r in rs]) for k in rs[0]['frequency_resolution_sensitivity']},
        historical_lnZ_bias_vs_TI=est([r['historical_lnZ_diagnostic_s_ex']-r['TI_s_ex'] for r in rs]),
        absolute_threshold_pass=abs(bias['mean'])<=.05,
        window_threshold_pass=bool(np.ptp(windowvals)<=.02),
        probe_D_threshold_pass=bool(abs(np.mean(pd))<=.05),net_D_threshold_pass=bool(abs(np.mean(nd))<=.05))

contrasts={}
for key,label,a,b in PAIRS:
    diffTI=np.array([rows[a][r]['TI_s_ex']-rows[b][r]['TI_s_ex'] for r in [1,2,3]])
    d={}
    for method in [PRIMARY,ALTERNATE]:
        delta=np.array([sx(rows[a][r],method)-sx(rows[b][r],method) for r in [1,2,3]])
        err=est(delta-diffTI)
        d[method]=dict(delta=est(delta),bias_vs_TI=err,practical_threshold_pass=abs(err['mean'])<=.02,zero_bias_in_interval=contains_zero(err))
    contrasts[key]=dict(label=label,TI=est(diffTI),methods=d,window_sensitivity={w:dict(delta=est([sx(rows[a][r],w)-sx(rows[b][r],w) for r in [1,2,3]]),bias_vs_TI=est([sx(rows[a][r],w)-sx(rows[b][r],w)-diffTI[r-1] for r in [1,2,3]])) for w in ['tau5_d1','tau10_d1','tau20_d1']})

null=[]
for rep in [1,2,3]:
    u='weak_f030_e025';l='dispersion_split_d000'
    um=json.loads((HERE/'cache'/f'{u}_r{rep}.json').read_text());lm=json.loads((HERE/'cache'/f'{l}_r{rep}.json').read_text())
    uz=np.load(HERE/'cache'/f'{u}_r{rep}.npz');lz=np.load(HERE/'cache'/f'{l}_r{rep}.npz')
    cg=[g for g in lm['groups'] if g['primary']];same=True;err=0.
    for i,(a,b) in enumerate(zip(um['groups'],cg)):
        same &= a['key']==b['key'] and a['sampled_ids']==b['sampled_ids'] and a['corrected_velocity_sha256']==b['corrected_velocity_sha256']
        for seg in ['full','first_half','last_half']:err=max(err,float(np.max(abs(uz[f'g{i}_{seg}']-lz[f'g{i}_{seg}']))))
    ent=max(abs(sx(rows[u][rep],k,seg)-sx(rows[l][rep],k,seg)) for seg in ['full','first_half','last_half'] for k in rows[u][rep]['variants'][seg])
    null.append(dict(replica=rep,identical_particle_ids_and_corrected_velocity_hashes=same,max_VACF_difference=err,max_entropy_difference=ent))
null_unmerged=est([rows['dispersion_split_d000'][r]['label_unmerged_finite_sample_s_ex']-sx(rows['weak_f030_e025'][r]) for r in [1,2,3]])

small_salt={}
for label,a,b in [('H7/H5','mono_weak_split7C','dispersion_split_d040'),('U6/U4','mono_weak_uniform6C','weak_f030_e025')]:
    ti=np.array([rows[a][r]['TI_s_ex']-rows[b][r]['TI_s_ex'] for r in [1,2,3]])
    vv=np.array([sx(rows[a][r])-sx(rows[b][r]) for r in [1,2,3]])
    small_salt[label]=dict(TI=est(ti),primary=est(vv),bias_vs_TI=est(vv-ti),scope='exploratory small-salt diagnostic, not a prespecified primary contrast')

unique=[c for c,n in ORDER if c!='L5']
checks=dict(analytic_checks=json.loads((HERE/'analytic_tests.json').read_text()),
    all_mode_reference_reconstruction_checks_pass=max(global_h+global_mode+global_reconstruct)<1e-10,
    max_mode_error=max(global_mode),max_h_reference_span=max(global_h),max_reconstruction_error=max(global_reconstruct),
    raw_targetT_mode_range=[min(raw_modes),max(raw_modes)],
    primary_negative_residual_fraction_range=[min(negative),max(negative)],
    primary_clipping_h01_to_h10_spurious_entropy_span_range=[min(negative)*3*np.log(100),max(negative)*3*np.log(100)],
    exact_mechanical_label_null_pass=all(r['identical_particle_ids_and_corrected_velocity_hashes'] and r['max_entropy_difference']<1e-10 for r in null),
    exact_label_null_details=null,unmerged_label_finite_sample_difference=null_unmerged,
    absolute_entropy_pass_codes=[c for c in unique if states[c]['absolute_threshold_pass']],
    absolute_entropy_fail_codes=[c for c in unique if not states[c]['absolute_threshold_pass']],
    absolute_bias_interval_contains_zero_codes=[c for c in unique if contains_zero(states[c]['primary_bias_vs_TI'])],
    window_pass_codes=[c for c in unique if states[c]['window_threshold_pass']],
    window_fail_codes=[c for c in unique if not states[c]['window_threshold_pass']],
    probe_D_pass_codes=[c for c in unique if states[c]['probe_D_threshold_pass']],
    probe_D_fail_codes=[c for c in unique if not states[c]['probe_D_threshold_pass']])
checks['both_controlled_increments_validated']=all(contrasts[k]['methods'][PRIMARY]['practical_threshold_pass'] and contrasts[k]['methods'][PRIMARY]['zero_bias_in_interval'] for k in ['controlled','mono_controlled'])
checks['quantitative_rehabilitation_pass']=not checks['absolute_entropy_fail_codes'] and checks['both_controlled_increments_validated'] and not checks['window_fail_codes'] and checks['exact_mechanical_label_null_pass']
summary=dict(complete=True,completed_replicas=24,primary_model=PRIMARY,alternate_model=ALTERNATE,states=states,contrasts=contrasts,exploratory_small_salt_effects=small_salt,validation=checks,
    provenance=dict(protocol_sha256=digest(HERE/'protocol.json'),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,code_sha256={p.name:digest(p) for p in HERE.glob('*.py')}),
    interpretation='internally consistent approximate classical estimator; independent quantitative TI validation required; no parameters calibrated to TI')
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (HERE/'state_comparison.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['code','TI','R2PT_delta1','R2PT_delta1p5','primary_bias','bias_CI_low','bias_CI_high','probe_spectral_D_error_percent','window_span_tau5_20'])
    for c,n in ORDER:
        a=states[c];w.writerow([c,a['TI']['mean'],a['primary']['mean'],a['alternate']['mean'],a['primary_bias_vs_TI']['mean'],*a['primary_bias_vs_TI']['CI'],100*a['probe_D_relative_error']['mean'],a['window_span_tau5_10_20']])
with (HERE/'contrast_comparison.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['comparison','TI_delta','primary_delta','primary_increment_bias','bias_CI_low','bias_CI_high','alternate_delta','alternate_increment_bias'])
    for key,label,a,b in PAIRS:
        d=contrasts[key];v=d['methods'][PRIMARY];w.writerow([label,d['TI']['mean'],v['delta']['mean'],v['bias_vs_TI']['mean'],*v['bias_vs_TI']['CI'],d['methods'][ALTERNATE]['delta']['mean'],d['methods'][ALTERNATE]['bias_vs_TI']['mean']])

lines=['# 2PT rehabilitation: completed tests','',
'**Assessment: the implementation can be rehabilitated as an internally consistent classical spectral estimator, but it does not pass the prespecified quantitative excess-entropy validation for these mixtures.** It captures the broad entropy ordering and one controlled increment reasonably well. Its absolute bias is composition dependent, and it overestimates the controlled increment with monovalent background. TI remains the primary entropy estimate.','',
'Figures: [absolute-entropy validation](entropy_validation.pdf), [entropy-increment validation](increment_validation.pdf), and [probe dynamics](probe_dynamics.pdf).','',
'All 24 existing velocity records (eight mixtures × three independent replicas, N=4000) were analyzed. No new MD jobs were required. The primary model and practical thresholds were recorded in protocol.json before computing the baseline result. No parameter, entropy offset, partial volume, or partition exponent was fitted to TI.','',
'State codes follow the [study compositions](../README.md#selected-compositions): B3/B5 are baselines without/with monovalent background; U4 has 30% uniform secondary solvent (epsilon=0.25); L5 is its label-only split; M5/H5 have physically split secondary solvents; U6/H7 repeat uniform/strong diversification with monovalent background. All contain 3% of each divalent sign.','',
'## Formulation and scope','',
'The new calculation follows the classical signed-residual R2PT construction of [Sun et al., JCP 147, 194505 (2017)](https://arxiv.org/abs/1711.02041), with partition exponent δ=1 primary and the literature δ=1.5 as an unfitted sensitivity. Species use the equal-partial-volume approximation, Vbar_i=1/rho; this is the one-fluid volume approximation discussed by [Lai et al., PCCP 14, 15206 (2012)](https://doi.org/10.1039/C2CP42011B). Equal nominal diameters do not prove equal actual partial volumes in charged, dispersion-diverse liquids. These are allocated spectral species contributions, not independently measured chemical-potential derivatives.','',
'The density-scaled diffusivity is Δ=(8/3)(36/π²)^(1/3) D sqrt(πm/T) rho_eff^(1/3). The archived implementation instead used rho^(-1/3) and a different prefactor. Gas entropy is weighted by its fraction f. Classical oscillator weights, ideal-mixture subtraction, and the analytical Lorentzian tail are treated consistently. Quantum weights are excluded when comparing against classical TI.','',
'The entire biased finite-record VACF is used, without subtracting each particle’s time-averaged velocity. One recorded system COM velocity is removed. Gaussian tapering preserves nonnegative total spectra; the full Nyquist band is retained. The primary taper is 10 tau, with 3, 5, 20 tau and untapered sensitivities. Spectral normalization uses measured C(0); raw target-temperature mode weights are reported, so the enforced sum rule is not mistaken for independent accuracy evidence.','',
'The gas Lorentzian integrates analytically to 3f modes. Its entropy tail is integrated to infinity, and the residual G−Ggas remains signed as R2PT bookkeeping. It is not clipped to create a nonnegative “solid population.” Ideal mixing is added to the approximate absolute mixture entropy and subtracted in the ideal reference, cancelling exactly; arbitrary classical phase-space normalization h also cancels.','',
'## Absolute entropy benchmark','',
'Values are kB per particle. Parenthesized model columns are means, not independent thermodynamic determinations. Bias intervals are nominal paired 95% Student-t intervals across three replicas, df=2; no multiplicity correction or systematic error is included.','',
'| State | TI | Primary R2PT, δ=1 | R2PT, δ=1.5 | Primary − TI [95% CI] |','|---|---:|---:|---:|---:|']
for c,n in ORDER:
    a=states[c];lines.append(f"| {c} | {a['TI']['mean']:.5f} | {a['primary']['mean']:.5f} | {a['alternate']['mean']:.5f} | {fmt(a['primary_bias_vs_TI'])} |")
lines+=['',
'The old baseline 2PT result was about −5.376, compared with corrected TI near −1.97. The corrected primary baseline is now about −2.074: the discrepancy is reduced from roughly 3.4 to 0.10 kB. This is substantial rehabilitation of the implementation, but a 0.10 kB error is as large as the controlled entropy signal. It is not an ensemble effect.','',
'Against the prespecified absolute-bias tolerance 0.05 kB, the independent physical states that fail are '+', '.join(checks['absolute_entropy_fail_codes'])+'. A constant entropy-offset correction would not solve the composition dependence. Changing the partition exponent moves absolute entropies substantially and does not make every mixture agree with TI. All primary absolute-bias intervals lie below zero at nominal 95% confidence, including states whose point bias meets the practical tolerance. Meeting that tolerance is not a statistical demonstration of unbiasedness.','',
'## Entropy increments','',
'| Comparison | TI increment | Primary increment | Primary increment bias [95% CI] | δ=1.5 increment |','|---|---:|---:|---:|---:|']
for key,label,a,b in PAIRS:
    d=contrasts[key];v=d['methods'][PRIMARY];lines.append(f"| {label} | {d['TI']['mean']:+.5f} | {v['delta']['mean']:+.5f} | {fmt(v['bias_vs_TI'])} | {d['methods'][ALTERNATE]['delta']['mean']:+.5f} |")
lines+=['',
'For H5/U4, the primary increment is 0.11818 versus TI 0.10322. Its +0.01496 point bias meets the practical 0.02 kB threshold and its interval contains zero, but the wide interval does not establish equivalence within 0.02 kB. The δ=1.5 increment is 0.10397, but this good agreement in one comparison does not validate that partition generally. For H7/U6, both partitions overestimate the TI increment: the primary bias is +0.03596 kB, interval [0.01501,0.05691]; δ=1.5 gives +0.02755, interval [0.00705,0.04804]. Both fail the controlled-increment test with monovalent background. The broad baseline-to-design entropy changes are also overestimated.','',
'An exploratory small-salt comparison further illustrates the precision limit: H5→H7 has TI entropy change −0.00385 kB (interval [−0.00669,−0.00100]), whereas primary 2PT gives +0.01081 (interval [−0.03736,+0.05898]). The 2PT estimate cannot resolve this small effect; its opposite point-estimate sign is not a significant contradiction given that wide interval. These comparisons are marked exploratory in summary.json.','',
'The controlled primary increments also remain above TI across the 5/10/20-tau smoothing variants. See the window_sensitivity fields in summary.json for paired increments and bias intervals; changing smoothing is not used to select the closest agreement.','',
'## Numerical and label checks','',
'All seven analytic tests pass: the standard fluidicity polynomial; invariance of Δ to length/time units; an Ornstein–Uhlenbeck mode-count/Green–Kubo limit; the analytical Lorentzian entropy tail; the dilute ideal-gas limit; classical h cancellation; and exact synthetic label/mixing cancellation. See analytic_tests.json and test_two_pt.py.','']
lines.append(f"Across all model/window/segment calculations, maximum mode-count error is {max(global_mode):.2e}, reconstruction error {max(global_reconstruct):.2e}, and entropy span under h=0.1,1,10 is {max(global_h):.2e} kB. Raw target-T mode weights for primary species range from {min(raw_modes):.5f} to {max(raw_modes):.5f}; they are not all exactly three before equipartition normalization.")
lines+=['',
'**Actual label-only null:** all three U4/L5 pairs use exactly the same sampled particle IDs and have identical hashes of the corrected velocity data for every mechanical group. Their cached VACFs are identical and entropy differences are at floating-point roundoff for every tested model, taper, and time segment. The exactly identical neutral interactions are pooled before entropy estimation; this avoids nonlinear estimator dependence on arbitrary labels. This merger is explicit, not a hidden fitting parameter.','']
lines.append(f"An unmerged L5 calculation using its separately retained label samples gives a finite-sample difference of {fmt(null_unmerged)} kB relative to U4. Its samples and nonlinear spectral transformation differ from the matched pooled null; it is a diagnostic, not an additional independent physical state.")
lines.append(f"The signed residual has negative spectral weight of {100*min(negative):.2f}–{100*max(negative):.2f}% of the total three modes across primary species. Clipping it while retaining the full gas contribution would add modes and reintroduce an arbitrary h dependence; over h=0.1 to 10 the spurious entropy span would be {min(negative)*3*np.log(100):.3f}–{max(negative)*3*np.log(100):.3f} kB. This is why the positive-part clipping in the archived code is not retained.")
lines+=['','## Sampling and model sensitivity','',
'Full 100-tau records, first/last 50-tau halves, and several taper widths were analyzed. Halves are sensitivity checks, not extra independent replicas. The table gives the span of state-mean primary entropy over taper widths 5, 10 and 20 tau, and the primary mean probe spectral-D discrepancy relative to the full-production MSD coefficient.','',
'| State | Entropy span, tau5/10/20 | Probe spectral D error (%) | Probe fluidicity |','|---|---:|---:|---:|']
for c,n in ORDER:
    a=states[c];lines.append(f"| {c} | {a['window_span_tau5_10_20']:.5f} | {100*a['probe_D_relative_error']['mean']:+.2f} | {a['probe_fluidicity']['mean']:.5f} |")
padding=max(abs(a['frequency_resolution_sensitivity'][k]['mean']) for a in states.values() for k in ['padding2','padding4'])
decim=max(abs(a['frequency_resolution_sensitivity'][k]['mean']) for a in states.values() for k in ['decimate2','decimate4'])
vol=max(abs(v['mean']) for a in states.values() for v in a['partial_volume_sensitivity'].values())
half=max(abs(v['mean']) for a in states.values() for v in a['half_record_sensitivity'].values())
lines+=['',f"Doubling/quadrupling frequency resolution changes state means by at most {padding:.6f} kB; doubling/quadrupling the sampling interval (a Nyquist/aliasing stress test) changes them by at most {decim:.6f} kB. The largest state-mean shift for either half record is {half:.5f} kB. These numerical checks do not remove the partition-model bias.",
 f"Illustrative, volume-conserving ±10% changes in probe or secondary-solvent partial volume shift primary state means by at most {vol:.5f} kB. These perturbations are sensitivity tests, not measured partial-volume errors or a proof that equal partial volumes are adequate. Converged partial volumes would be needed for a stronger mixture-2PT claim.",
 '',
'The 5% spectral-D practical check fails for '+(', '.join(checks['probe_D_fail_codes']) if checks['probe_D_fail_codes'] else 'none of the state-mean probe estimates')+'. This comparison uses a smoothed 100-tau velocity estimator versus a 1000-tau MSD estimator; agreement does not independently establish the entropy formula. Full replica intervals and all species diagnostics are in summary.json and replicas/*.json.','',
'## What dynamical insight survives?','',
'Probe fluidicity rises from approximately 0.530 in B3 to 0.569 in H5; U4→H5 increases it from about 0.554 to 0.569, and U6→H7 from about 0.552 to 0.568. This is consistent with increased low-frequency motion, but diffusivity already enters the construction, so it is not independent proof of an entropy-driven speedup. Fluidicity is not the fraction of unpaired ions.','',
'The mean probe VACFs remain positive over the first tau in these records. Thus a pronounced early negative backscattering minimum, as previously described for the archived figures, is not supported by this new calculation. A local association or cage-escape mechanism still requires coordination and residence-time measurements.','',
'## Recommendation','',
'Keep TI as the entropy coordinate in the Rosenfeld analysis. Retain the corrected 2PT implementation as a transparent approximation and the spectra as supporting dynamics. Do not restore the old numerical entropy tables, call the species allocations partial molar entropies, or use the improved 2PT–diffusion correlation as an independent confirmation. For quantitative screening, validation must hold across compositions; it currently fails for the monovalent-background increment and for some absolute entropies. Extending velocity records and obtaining partial volumes can test sampling/reference assumptions, but neither should be tuned to force agreement with TI.','',
'Reproduction: run test_two_pt.py, run_rehabilitation.py, then summarize_rehabilitation.py in this directory. The cache contains COM-corrected VACFs, particle IDs, velocity hashes and input provenance. Original MD data, TI results, and archived 2PT outputs were not overwritten.']
(HERE/'report.md').write_text('\n'.join(lines)+'\n')

plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for ax,method,title in zip(axs,['primary','alternate'],['Primary partition: δ=1','Literature sensitivity: δ=1.5']):
    lim=[]
    for c,n in ORDER:
        if c=='L5':continue
        a=states[c];x=a['TI']['mean'];y=a[method]['mean'];color='#b2182b' if c in ['B5','U6','H7'] else '#2166ac'
        ax.errorbar(x,y,xerr=np.diff(a['TI']['CI'])/2,yerr=np.diff(a[method]['CI'])/2,fmt='o',color=color,capsize=3)
        ax.annotate(c,(x,y),xytext=((-24,-10) if c=='B5' else (4,-13 if c in ['U6','H7'] else 5)),textcoords='offset points',fontsize=9);lim.extend([*a['TI']['CI'],*a[method]['CI']])
    lo=min(lim)-.025;hi=max(lim)+.025;ax.plot([lo,hi],[lo,hi],'k--',lw=.8);ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel=r'TI $s_{ex}/k_B$',ylabel=r'R2PT $s_{ex}/k_B$',title=title);ax.grid(alpha=.2)
fig.suptitle('Nominal 95% replica intervals; dashed line: equality',fontsize=10)
for ext in ['png','pdf']:fig.savefig(HERE/f'entropy_validation.{ext}',dpi=210)
plt.close(fig)

fig,ax=plt.subplots(figsize=(7.5,4),layout='constrained');labs=['H5/U4','H7/U6','H5/B3','H7/B5'];keys=['controlled','mono_controlled','baseline','mono_baseline']
for j,(method,label,color) in enumerate([('TI','TI','#222222'),(PRIMARY,'R2PT δ=1','#2166ac'),(ALTERNATE,'R2PT δ=1.5','#b2182b')]):
    vals=[contrasts[k]['TI'] if method=='TI' else contrasts[k]['methods'][method]['delta'] for k in keys]
    y=np.array([a['mean'] for a in vals]);lo=np.array([a['CI'][0] for a in vals]);hi=np.array([a['CI'][1] for a in vals]);ax.errorbar(np.arange(4)+(j-1)*.18,y,yerr=[y-lo,hi-y],fmt='o',capsize=3,color=color,label=label)
ax.set_xticks(range(4),labs);ax.set(ylabel=r'Entropy increment $\Delta s_{ex}/k_B$');ax.grid(axis='y',alpha=.2);ax.legend(frameon=False)
fig.suptitle('Paired entropy increments with nominal 95% replica intervals',fontsize=10)
for ext in ['png','pdf']:fig.savefig(HERE/f'increment_validation.{ext}',dpi=210)
plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for code,name,color in [('B3','baseline_3C','#555555'),('U4','weak_f030_e025','#2166ac'),('H5','dispersion_split_d040','#1b7837'),('U6','mono_weak_uniform6C','#e08214'),('H7','mono_weak_split7C','#b2182b')]:
    cv=[];Gs=[]
    for rep in [1,2,3]:
        m=json.loads((HERE/'cache'/f'{name}_r{rep}.json').read_text());z=np.load(HERE/'cache'/f'{name}_r{rep}.npz');idx=[i for i,g in enumerate(m['groups']) if g['primary'] and g['probe']]
        cb=np.mean([z[f'g{i}_full'] for i in idx],axis=0);tm=np.arange(len(cb))*m['dt'];cv.append(cb/(1-np.arange(len(cb))/len(cb))/cb[0]);Gs.append(np.mean([spectrum(z[f'g{i}_full'],m['dt'],10,2.)[1] for i in idx],axis=0))
    cv=np.array(cv);Gs=np.array(Gs);take=tm<=1.5;axs[0].plot(tm[take],cv.mean(axis=0)[take],color=color,label=code);nu=spectrum(cb,m['dt'],10,2.)[0];take=nu<=5;axs[1].plot(nu[take],Gs.mean(axis=0)[take],color=color,label=code)
axs[0].axhline(0,color='k',lw=.6);axs[0].set(xlabel=r'Time ($\tau$)',ylabel=r'Probe VACF $C(t)/C(0)$',title='Short-time dynamics');axs[1].set(xlabel=r'Frequency ($\tau^{-1}$)',ylabel=r'Probe spectrum $G(\nu)$',title='Gaussian-smoothed spectrum (10τ)')
for ax in axs:ax.legend(frameon=False);ax.grid(alpha=.2)
fig.suptitle('Probe dynamics: means of three independent replicas',fontsize=10)
for ext in ['png','pdf']:fig.savefig(HERE/f'probe_dynamics.{ext}',dpi=210)
plt.close(fig)

print('Completed 24-replica report, CSV/JSON summaries and three figures.')
print('Validation:',json.dumps({k:checks[k] for k in ['quantitative_rehabilitation_pass','absolute_entropy_fail_codes','window_fail_codes','probe_D_fail_codes','exact_mechanical_label_null_pass']},indent=2))
