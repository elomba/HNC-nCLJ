#!/usr/bin/env python3
"""Reproduce the fixed-LJ cutoff analysis, paired sampling errors, and figures.

Uses the second half of each trajectory, matching the previous analysis.
Bootstrap intervals quantify finite-time sampling conditional on independent
density simulations; they exclude finite-size, PPPM, and extrapolation errors.
Original simulation results and analysis CSVs are never overwritten.
"""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import os
import re
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
os.environ.setdefault('MPLCONFIGDIR','/tmp/entropy-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

DOC=Path(__file__).resolve().parent
BASE=DOC.parent/'lammps_study'
RUN=BASE/'fixed_lj_20261030'
OUT=RUN/'results'
SYSTEMS=['system_3C','system_4C','system_4C_eps075','system_5C']
LABELS=['3C',r'4C ($\epsilon_2=1.25$)',r'4C ($\epsilon_2=0.75$)','5C']
COLORS=['#1f77b4','#2ca02c','#ff7f0e','#d62728']
RHO=np.array([.01,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7])
TEMP=2.
NBOOT=5000


def legacy(r,z):
    y=(z-1)/r
    intercept=np.polyfit(r[:4],y[:4],2)[-1]
    spline=CubicSpline(np.r_[0,r],np.r_[intercept,y])
    return spline.antiderivative()(r)-spline.antiderivative()(0)


def dh_coefficient(q2):
    # units lj, dielectric 1, Coulomb energy prefactor 1, T=2: l_B=1/T.
    return (4*np.pi*q2/TEMP)**1.5/(12*np.pi)


def screened(r,z,q2,method='cubic'):
    # x=sqrt(rho) removes the Debye-Huckel integrable singularity.
    x=np.sqrt(r)
    f=2*(z-1)/x
    cls=CubicSpline if method=='cubic' else PchipInterpolator
    spline=cls(np.r_[0,x],np.r_[-dh_coefficient(q2),f])
    integral=spline.antiderivative()
    return integral(x)-integral(0)


def anchored(r,z,q2):
    # Z-1=-(C/2)sqrt(rho)+B*rho below the lowest sampled density.
    c=dh_coefficient(q2)
    a0=(z[0]-1)-.5*c*np.sqrt(r[0])
    spline=CubicSpline(r,(z-1)/r).antiderivative()
    return a0+spline(r)-spline(r[0])


def linear_weights(func):
    z=np.ones(len(RHO))
    constant=func(RHO,z)
    weights=np.column_stack([func(RHO,z+np.eye(len(RHO))[i])-constant for i in range(len(RHO))])
    offset=constant-weights@z
    assert np.allclose(func(RHO,z+.3),offset+weights@(z+.3))
    return offset,weights


def ci(x):
    return [float(a) for a in np.percentile(x,[2.5,97.5])]


def fit(x,d):
    y=np.log(d)
    slope,intercept=np.polyfit(x,y,1)
    return dict(B=float(slope),intercept=float(intercept),
                R2=float(1-np.sum((y-(slope*x+intercept))**2)/np.sum((y-y.mean())**2)))


def bootstrap_means(a,block,rng):
    # Circular moving blocks; columns from old and new trajectories stay paired.
    n=len(a)
    starts=rng.integers(0,n,size=(NBOOT,int(np.ceil(n/block))))
    indices=((starts[:,:,None]+np.arange(block))%n).reshape(NBOOT,-1)[:,:n]
    return a[indices].mean(axis=1)


def write_csv(path,rows):
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def write_report(report):
    results=report['systems']
    text="""# Fixed-LJ sweep: completed analysis

All **60 runs** (four mixtures × 15 densities, jobs 100493–100552) completed normally and reached their expected production length. Inputs and configurations were checked against the submission manifest. Original outputs remain untouched.

**The cutoff error was material for absolute excess entropy, but fixing it does not establish enhanced divalent-ion diffusion through Rosenfeld scaling.** It lowers dense-state excess entropy by about 0.085–0.091 kB per particle, removes approximately half of the former HNC–MD entropy discrepancy, and leaves the main qualitative interpretation intact.

## Isolating the cutoff correction

The following comparison uses the *same original postprocessing* on both datasets: final half of production, stored mean instantaneous Z, and a cubic spline with a finite fitted low-density intercept. This isolates the changed LJ potential from subsequent estimator improvements.

| System | Original s_ex/kB | Fixed LJ, same estimator | Cutoff-only shift |
|---|---:|---:|---:|
"""
    for key,label in zip(SYSTEMS,['3C','4C, ε2=1.25','4C, ε2=0.75','5C']):
        r=results[key]['dense']
        text+=f"| {label} | {r['s_ex_old_legacy']:.5f} | {r['s_ex_fixed_legacy']:.5f} | {r['cutoff_shift_same_method']:+.5f} |\n"
    text+="""
For rho≥0.20, the original LJ cutoff was already 2.5 sigma. Those reruns reproduce the original thermodynamic output, and the dense MSD arrays are numerically identical. They are **not independent transport replicas**. At rho=0.01–0.15, the corrected pressure is higher because the spurious longer-range LJ attraction is removed. Its integral raises excess free energy and consequently lowers excess entropy at all higher densities. The effect exceeds the conditional sampling errors by a large margin.

## Refined thermodynamic estimates

The primary calculation additionally uses **Z = mean(P)/(rho T_target)** with T_target=2, rather than mean[P/(rho T_instantaneous)]. The latter ratio is a different finite-sample estimator; converting to target-temperature pressure changes dense entropy by up to 0.0023 kB. This matters for small composition differences.

The integration variable is x=sqrt(rho), with F(x)=2[Z(x²)−1]/x. Debye–Hückel asymptotics give F(0)=−C, C=(4π lB Σxi qi²)^(3/2)/(12π). A cubic spline in x then integrates a finite endpoint. In the actual LAMMPS and HNC inputs the Coulomb **energy** prefactor is 1, so lB=1/T=0.5; calling the energy prefactor itself the Bjerrum length would be incorrect. The manuscript definition was clarified. See the [LAMMPS reduced-unit convention](https://docs.lammps.org/units.html).

| System | Primary s_ex/kB | 95% sampling interval | HNC s_ex/kB |
|---|---:|---:|---:|
"""
    for key,label in zip(SYSTEMS,['3C','4C, ε2=1.25','4C, ε2=0.75','5C']):
        r=results[key]['dense']
        text+=f"| {label} | {r['s_ex_fixed_screened']:.5f} | [{r['s_ex_CI_lo']:.5f}, {r['s_ex_CI_hi']:.5f}] | {r['s_ex_HNC']:.5f} |\n"
    text+="""
Intervals use 5000 circular moving-block bootstrap resamples, retaining pressure/energy pairing, with 10 tau blocks. Checks at 4 and 20 tau give similar error scales. These are **within-trajectory sampling intervals**, conditional on independent simulations at different densities; they are not total uncertainty bounds or substitutes for independent seeds.

Shape-preserving interpolation and a dilute anchor with a fixed Debye–Hückel term plus a fitted B*rho correction change the primary estimates by less than 0.0008 kB. They do not prove convergence below the lowest sampled density (0.01). Applied to the HNC equation of state, the same screened integration recovers its tabulated dense free energy to within 0.0006–0.0008 kBT, providing a useful numerical crosscheck. The analytic Debye–Hückel limiting curve is integrated exactly by the implementation.

Using all production samples, the final 75%, or the final half changes entropy by up to 0.0045 kB. This window sensitivity is comparable to some composition differences and should accompany the bootstrap intervals. The cutoff shift of roughly 0.09 kB is much larger. A denser grid below 0.01 and independent longer runs would be the useful next thermodynamic checks.

## Which scientific conclusions survive?

- **Solvent attraction:** weak-minus-strong 4C gives Δs_ex/kB=+0.01063, with conditional 95% interval [+0.00504,+0.01639]. Mean diffusion increases by 2.084%. This is consistent with the expected direction of a solvent entropy–mobility association. It does not establish causality or a universal exponential relation.
- **Additional components:** 5C minus 3C gives Δs_ex/kB=−0.00895, with conditional 95% interval [−0.01513,−0.00282]. Adding 1:1 salt does not raise excess entropy here. Ideal mixing entropy increases by +0.38787 kB, producing a total entropy increase of about +0.37892 kB; this mixing contribution is excluded from Rosenfeld scaling.
- **HNC comparison:** the original entropy gap of approximately 0.18 kB becomes approximately 0.092–0.098 kB. The correction explains about half, not all, of the discrepancy. Dense HNC pressure discrepancies are unaffected.
- **Divalent transport:** the original/rerun MSD diffusivities remain 0.21429, 0.19686, 0.16660, and 0.22887. The separate trajectory's GK estimates remain about 0.2038, 0.2017, 0.2057, and 0.2037. Correcting the entropy coordinate cannot reconcile these datasets. Group-specific COM subtraction, finite minority-ion sampling, and matched estimator convergence still require investigation.

## Rosenfeld fits

Unweighted ln(D_reduced) versus s_ex/kB fits, with diffusion fixed to the dense MSD values:

| Entropy analysis | Solvent-series B | Solvent-series R² | All-mixture B | All-mixture R² |
|---|---:|---:|---:|---:|
"""
    for method,label in [('s_ex_old_legacy','Original'),('s_ex_fixed_legacy','Fixed LJ, original estimator'),('s_ex_screened_archived_Z','Fixed LJ + screened integration, stored Z'),('s_ex_fixed_screened','Fixed LJ + screened integration, target-T pressure')]:
        r=report['fit_method_comparison'][method];a=r['solvent_series'];b=r['all_systems']
        text+=f"| {label} | {a['B']:.3f} | {a['R2']:.3f} | {b['B']:.3f} | {b['R2']:.3f} |\n"
    text+="""
The solvent-only fit improves, but it has just three points and is sensitive to entropy-estimator choices. The four-mixture fit remains poor. Fitting errors in both coordinates and obtaining independent transport uncertainties would be necessary before treating any slope as predictive. The negative descriptive slope for the ion MSD values is not evidence of a validated inverse entropy law.

## Next analysis priorities

1. Recompute species MSD and VACF in one common global reference frame, with multiple time origins and independent trajectory seeds. Verify cation/anion symmetry and species/pool additivity. The cutoff reruns do not supply new independent evidence for ion transport.
2. Extend the dilute grid and check independent thermodynamic replicas, especially for 5C, whose averaging-window sensitivity is largest.
3. Only then test whether salt addition changes association fractions and lifetimes, followed by the concentration, particle-diameter and box-size controls described in the manuscript.

The 2PT normalization/reference discrepancy is independent of this LJ cutoff correction and remains unresolved.

## Files and reproduction

Run `MPLCONFIGDIR=/tmp/entropy-mpl OPENBLAS_NUM_THREADS=1 python3 doc/analyze_fixed_lj_sweep.py` from the repository. The script validates all 60 runs and reproduces the calculations and figures. Full density results, paired-bootstrap samples, and machine-readable summaries are in `lammps_study/fixed_lj_20261030/results/`. Compact copies are in `doc/fixed_lj_dense_comparison.csv` and `doc/fixed_lj_analysis.json`.

The manuscript source, abstract, thermodynamic table, figures and conclusions have been updated. Old results and figures remain available. Build the PDFs with `bash doc/build_pdfs.sh`, which loads `texlive/20260310-GCC-13.2.0` and uses the system zlib with direct pdfLaTeX passes.
"""
    (DOC/'fixed_lj_results.md').write_text(text)
    (OUT/'analysis_report.md').write_text(text)


def main():
    OUT.mkdir(exist_ok=True)
    manifest=json.loads((RUN/'manifest.json').read_text())
    with (BASE/'results/rosenfeld_diffusion_coefficients.csv').open() as f:
        old_diff={r['System']:r for r in csv.DictReader(f)}
    dense_rows=[];sweep_rows=[];boot_rows=[];all_results={};raw={}
    rng=np.random.default_rng(20261030)
    for system,label in zip(SYSTEMS,LABELS):
        matrix=[];old_matrix=[];joined=[];checks=[]
        for rho in RHO:
            directory=RUN/system/f'rho_{rho:.2f}'
            record=next(r for r in manifest['runs'] if r['system']==system and np.isclose(r['rho'],rho))
            assert hashlib.sha256((directory/'in.run').read_bytes()).hexdigest()==record['input_sha256']
            assert hashlib.sha256((directory/'data.mixture').read_bytes()).hexdigest()==record['data_sha256']
            text=(directory/'lammps_run.log').read_text()
            assert 'Total wall time:' in text and 'ERROR' not in text,directory
            new=np.loadtxt(directory/'thermo_avg.dat',skiprows=2)
            old=np.loadtxt(BASE/system/directory.name/'thermo_avg.dat',skiprows=2)
            expected=150000 if rho==.7 else 50000
            assert new[-1,0]==expected and len(new)==expected//1000,directory
            assert np.array_equal(new[:,0],old[:,0]),directory
            new,old=new[len(new)//2:],old[len(old)//2:]
            matrix.append(new.mean(axis=0));old_matrix.append(old.mean(axis=0))
            joined.append(np.column_stack([new[:,2]/(rho*TEMP),new[:,4],old[:,3],old[:,4],new[:,3]]))
            checks.append(dict(rho=float(rho),production_steps=int(expected),retained_samples=len(new),
                               old_new_thermo_identical=bool(np.array_equal(new,old))))
        new=np.array(matrix);old=np.array(old_matrix)
        q2=.34 if system=='system_5C' else .24
        methods={'legacy':legacy,'screened_cubic':lambda r,z:screened(r,z,q2),
                 'screened_pchip':lambda r,z:screened(r,z,q2,'pchip'),
                 'legacy_fixed_T_pressure':legacy,
                 'DH_B_anchor':lambda r,z:anchored(r,z,q2)}
        new_z=new[:,2]/(RHO*TEMP)
        af={name:func(RHO,new[:,3] if name=='legacy' else new_z) for name,func in methods.items()}
        a_old=legacy(RHO,old[:,3]);s_old=old[:,4]/TEMP-a_old
        ss={name:new[:,4]/TEMP-a for name,a in af.items()}
        hnc_file=DOC.parent/'EntropyComparison/results'/('Thermo_'+system.removeprefix('system_')+'_density_sweep.dat')
        hnc=np.loadtxt(hnc_file);hnc=hnc[np.argsort(hnc[:,1])]
        assert np.allclose(hnc[:,1],RHO)
        dfile=RUN/system/'rho_0.70/msd.dat'
        old_msd=BASE/system/'rho_0.70/msd.dat'
        msd=np.loadtxt(dfile,skiprows=2); old_msd_data=np.loadtxt(old_msd,skiprows=2)
        assert np.array_equal(msd,old_msd_data),dfile
        diff={}
        for column,key in [(1,'net'),(2,'z2'),(3,'cat'),(4,'an')]:
            start=int(.3*len(msd));diff[key]=float(np.polyfit(msd[start:,0]*.002,msd[start:,column],1)[0]/6)
            assert abs(diff[key]-float(old_diff[system]['D_'+key]))<1e-6
        sensitivities={};boot=None
        for block in [2,5,10]:
            values=np.stack([bootstrap_means(a,block,rng) for a in joined],axis=1)
            c,w=linear_weights(methods['screened_cubic'])
            a_boot=values[:,:,0]@w.T+c
            s_boot=values[:,:,1]/TEMP-a_boot
            c_legacy,w_legacy=linear_weights(legacy)
            old_a_boot=values[:,:,2]@w_legacy.T+c_legacy
            old_s_boot=values[:,:,3]/TEMP-old_a_boot
            new_legacy_boot=values[:,:,1]/TEMP-(values[:,:,4]@w_legacy.T+c_legacy)
            sensitivities[str(block*2)]=dict(s_ex_ci=ci(s_boot[:,-1]),
                s_ex_SE=float(s_boot[:,-1].std(ddof=1)),
                cutoff_shift_same_method_ci=ci((new_legacy_boot-old_s_boot)[:,-1]))
            if block==5:
                boot=s_boot[:,-1]
                boot_rows.extend(dict(system=system,replicate=i,s_ex=float(b)) for i,b in enumerate(boot))
        raw[system]=dict(new=new,old=old,hnc=hnc,a=af,s=ss,s_old=s_old,boot=boot,diff=diff)
        row=dict(system=system,s_ex_old_legacy=float(s_old[-1]),s_ex_fixed_legacy=float(ss['legacy'][-1]),
                 cutoff_shift_same_method=float(ss['legacy'][-1]-s_old[-1]),
                 s_ex_fixed_screened=float(ss['screened_cubic'][-1]),
                 s_ex_screened_archived_Z=float(new[-1,4]/TEMP-screened(RHO,new[:,3],q2)[-1]),
                 s_ex_CI_lo=ci(boot)[0],s_ex_CI_hi=ci(boot)[1],
                 a_ex_fixed_screened=float(af['screened_cubic'][-1]),
                 u_ex_over_kBT=float(new[-1,4]/TEMP),Z=float(new_z[-1]),Z_archived=float(new[-1,3]),
                 s_ex_HNC=float(hnc[-1,7]),D_net=diff['net'],D_z2=diff['z2'])
        dense_rows.append(row)
        for i,rho in enumerate(RHO):
            sweep_rows.append(dict(system=system,rho=float(rho),T=float(new[i,1]),
                 Z_old=float(old[i,3]),Z_fixed=float(new_z[i]),
                 u_ex_over_kBT=float(new[i,4]/TEMP),a_ex_fixed=float(af['screened_cubic'][i]),
                 s_ex_fixed=float(ss['screened_cubic'][i]),s_ex_old=float(s_old[i]),
                 s_ex_fixed_legacy=float(ss['legacy'][i]),s_ex_HNC=float(hnc[i,7])))
        all_results[system]=dict(dense=row,completion=checks,
             integration_sensitivity={name:float(s[-1]) for name,s in ss.items()},
             block_bootstrap_sensitivity=sensitivities,
             HNC_free_energy_crosscheck=dict(functional=float(hnc[-1,6]),screened_integral=float(screened(RHO,hnc[:,3],q2)[-1])),
             fitting_window_entropy={})
        for fraction in [0.,.25,.5]:
            window_z=[];window_u=[]
            for rho in RHO:
                a=np.loadtxt(RUN/system/f'rho_{rho:.2f}'/'thermo_avg.dat',skiprows=2)
                a=a[int(fraction*len(a)):]
                window_z.append(a[:,2].mean()/(rho*TEMP));window_u.append(a[:,4].mean()/TEMP)
            all_results[system]['fitting_window_entropy'][str(fraction)]=float(window_u[-1]-screened(RHO,np.array(window_z),q2)[-1])
    ent=np.array([raw[k]['s']['screened_cubic'][-1] for k in SYSTEMS])
    d=np.array([raw[k]['diff']['net'] for k in SYSTEMS]);ion=np.array([raw[k]['diff']['z2'] for k in SYSTEMS])
    slopes={}
    for n,label in [(3,'solvent_series'),(4,'all_systems')]:
        for key,ds in [('net',d),('ion',ion)]:
            slopes[f'{label}_{key}']=fit(ent[:n],ds[:n])
    comparisons={}
    for method in ['s_ex_old_legacy','s_ex_fixed_legacy','s_ex_screened_archived_Z','s_ex_fixed_screened']:
        x=np.array([r[method] for r in dense_rows])
        comparisons[method]=dict(solvent_series=fit(x[:3],d[:3]),all_systems=fit(x,d))
    b=np.column_stack([raw[k]['boot'] for k in SYSTEMS])
    contrasts={}
    for i,j,label in [(2,1,'weak_minus_strong_4C'),(3,0,'5C_minus_3C'),(0,1,'3C_minus_strong_4C')]:
        contrasts[label]=dict(delta_s=float(ent[i]-ent[j]),CI=ci(b[:,i]-b[:,j]),
                              diffusion_ratio_net=float(d[i]/d[j]),diffusion_ratio_ion=float(ion[i]/ion[j]))
    write_csv(OUT/'dense_comparison.csv',dense_rows)
    write_csv(OUT/'density_sweep.csv',sweep_rows)
    write_csv(OUT/'entropy_bootstrap.csv',boot_rows)
    report=dict(n_runs=60,n_bootstrap=NBOOT,primary_method='screened cubic spline in sqrt(rho)',
                bootstrap='circular blocks, primary length 10 tau, final half of production; sampling only',
                bjerrum_length=1/TEMP,systems=all_results,fits=slopes,contrasts=contrasts,fit_method_comparison=comparisons)
    # Primary pressure is <P>/(rho*T_target), not <P/(rho*T_inst)>.
    report['pressure_definition']='mean pressure / (density * target temperature 2.0)'
    report['screened_integrand']='x=sqrt(rho), f=2*(Z-1)/x, f(0)=-(4*pi*lB*sum(xi*qi^2))^(3/2)/(12*pi)'
    (OUT/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    (DOC/'fixed_lj_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    write_csv(DOC/'fixed_lj_dense_comparison.csv',dense_rows)
    write_report(report)
    # Density-dependent thermodynamics and entropy comparison.
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for key,label,color in zip(SYSTEMS,LABELS,COLORS):
        r=raw[key]
        axes[0].plot(RHO,r['new'][:,2]/(RHO*TEMP),'o-',c=color,label=label)
        axes[0].plot(RHO,r['old'][:,2]/(RHO*TEMP),':',c=color,alpha=.6)
        axes[0].plot(RHO,r['hnc'][:,3],'--',c=color,alpha=.8)
        axes[1].plot(RHO,r['a']['screened_cubic'],'o-',c=color)
        axes[1].plot(RHO,r['hnc'][:,6],'--',c=color)
    axes[0].set_ylabel('Z');axes[1].set_ylabel(r'$A_{ex}/Nk_BT$')
    axes[0].legend(fontsize=8)
    axes[0].set_title('Fixed LJ: solid; original: dotted; HNC: dashed')
    axes[1].set_title('Corrected free energy (HNC: dashed)')
    for ax in axes:ax.set_xlabel(r'$\rho^*$');ax.grid(alpha=.2)
    fig.savefig(DOC/'figures/compressibility_and_ti_fixed_lj.png',dpi=220);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for key,label,color in zip(SYSTEMS,LABELS,COLORS):
        r=raw[key]
        axes[0].plot(RHO,r['s']['screened_cubic'],'o-',c=color,label=label)
        axes[0].plot(RHO,r['s_old'],':',c=color)
        axes[0].plot(RHO,r['hnc'][:,7],'--',c=color)
        axes[1].plot(RHO,r['hnc'][:,9],'-',c=color)
    axes[0].legend(fontsize=8)
    axes[0].set_ylabel(r'$s_{ex}/k_B$');axes[1].set_ylabel(r'$s^{(2)}/k_B$ (HNC)')
    axes[0].set_title('Fixed LJ: solid; original: dotted; HNC: dashed')
    axes[1].set_title('Pair structural entropy')
    for ax in axes:ax.set_xlabel(r'$\rho^*$');ax.grid(alpha=.2)
    fig.savefig(DOC/'figures/entropy_comparison_fixed_lj.png',dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
    baseline=raw['system_3C']
    for key,label,color in zip(SYSTEMS[1:],LABELS[1:],COLORS[1:]):
        ax.plot(RHO,raw[key]['s']['screened_cubic']-baseline['s']['screened_cubic'],'o-',c=color,label=label)
        ax.plot(RHO,raw[key]['hnc'][:,7]-baseline['hnc'][:,7],'--',c=color)
    ax.axhline(0,c='grey',lw=.8);ax.set_xlabel(r'$\rho^*$');ax.set_ylabel(r'$\Delta s_{ex}/k_B$ relative to 3C')
    ax.set_title('Fixed-LJ MD: solid; HNC: dashed');ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(DOC/'figures/delta_entropy_fixed_lj.png',dpi=220);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    scale=.7**(1/3)/np.sqrt(TEMP);gk=np.array([.2038,.2017,.2057,.2037])
    for i,(label,color) in enumerate(zip(LABELS,COLORS)):
        errors=np.array([[ent[i]-dense_rows[i]['s_ex_CI_lo']],[dense_rows[i]['s_ex_CI_hi']-ent[i]]])
        axes[0].errorbar(ent[i],d[i]*scale,xerr=errors,fmt='o',c=color,label=label,capsize=2)
        axes[1].errorbar(ent[i],ion[i]*scale,xerr=errors,fmt='o',c=color,capsize=2)
        axes[1].scatter(ent[i],gk[i]*scale,edgecolors=color,marker='s',facecolors='none')
    xx=np.linspace(ent.min()-.002,ent.max()+.002,100);f=slopes['solvent_series_net']
    axes[0].plot(xx,np.exp(f['intercept']+f['B']*xx)*scale,'k--',label=f"3-point fit: B={f['B']:.2f}, $R^2$={f['R2']:.2f}")
    axes[1].scatter([],[],c='black',label='Original / rerun MSD')
    axes[1].scatter([],[],edgecolors='black',facecolors='none',marker='s',label='Separate-trajectory GK')
    for ax in axes:
        ax.set_yscale('log');ax.set_xlabel(r'Corrected MD $s_{ex}/k_B$');ax.set_ylabel(r'$D\rho^{1/3}/\sqrt{k_BT/m}$')
        ax.xaxis.set_major_locator(MaxNLocator(4));ax.legend(fontsize=8);ax.grid(alpha=.2)
    axes[0].set_title('(a) Mean self-diffusion')
    axes[1].set_title(r'(b) Divalent-ion self-diffusion')
    fig.savefig(DOC/'figures/rosenfeld_scaling_fixed_lj.png',dpi=220);plt.close(fig)
    print(json.dumps(dict(dense=dense_rows,fits=slopes,contrasts=contrasts),indent=2))


if __name__=='__main__': main()
