# Review of entropy and divalent-ion diffusion

> **Update after the fixed-LJ reruns:** All 60 corrected runs completed. The
> [new results and uncertainty analysis](fixed_lj_results.md) supersede the
> numerical entropy estimates below. The original review is retained as a
> record of the issues that motivated the rerun; the cutoff-path issue is
> now corrected. Transport and 2PT issues remain unresolved.

The current results do **not** establish that increasing the number of components accelerates the +2/−2 ions through Rosenfeld scaling. There is a plausible solvent entropy–mobility association, but even the mean-diffusion scaling is preliminary. Ion-specific changes need reproducible transport estimates before a physical mechanism can be assigned.

This review concerns `paper_pre.tex` and the archived analysis inputs/results. `entropy_theory.tex` contains a historical transport comparison with different compositions (I=0.40 and 0.19); a note now distinguishes it from the current manuscript. Numerical values in the tables were preserved, with limitations made explicit. No new MD simulations were run.

## Entropy: what increases when components are added?

For these equal-mass particles at fixed temperature and total density,

\[
\Delta s/k_B=\Delta[-\sum_i x_i\ln x_i]+\Delta s_{ex}/k_B.
\]

Rosenfeld uses excess entropy relative to an ideal gas **of the same composition**, so the ideal mixing term is removed.

| Mixture | Ideal mixing entropy / kB | Reported MD excess entropy / kB | Mean D | Divalent D, first MSD dataset | Divalent D, second GK dataset |
|---|---:|---:|---:|---:|---:|
| 3C | 0.26856 | −1.87985 | 0.214795 | 0.214291 | 0.2038 |
| 4C, ε2=1.25 | 0.46390 | −1.89091 | 0.214924 | 0.196861 | 0.2017 |
| 4C, ε2=0.75 | 0.46390 | −1.87730 | 0.219403 | 0.166599 | 0.2057 |
| 5C | 0.65642 | −1.88535 | 0.217894 | 0.228874 | 0.2037 |

The 3C→5C total entropy increase is about 0.38237 kB per particle, largely the ideal mixing increase of 0.38787 kB. Excess entropy actually decreases by 0.00550 kB. A positive-slope Rosenfeld law with fixed ion prefactor would therefore predict a slight ion slowdown, not the 6.8% increase in the first MSD dataset. The second trajectory set gives essentially no increase: +0.4% by MSD and −0.05% by GK.

An exact control is to split the solvent into two labels with identical masses and interactions. Formal mixing entropy increases if the labels are treated as distinguishable components, while forces, excess entropy and diffusion remain unchanged. The ε2=1 limit of the proposed 4C series supplies this control.

The physical 3C→5C comparison also increases ion fraction from 6% to 16% and mole-fraction ionic strength from 0.12 to 0.17. It does not isolate component count. Both 4C fluids have exactly the same ideal mixing entropy, making their dispersion comparison more interpretable. Geometric energy mixing changes solvent–ion interactions as well as solvent–solvent ones.

## How convincing is the scaling?

Unweighted fits to the archived CSV give:

| Fit of ln(D reduced) | Points | B | R² |
|---|---:|---:|---:|
| Mean D vs MD excess entropy, solvent series | 3 | 1.045 | 0.391 |
| Mean D vs HNC excess entropy, solvent series | 3 | 1.402 | 0.780 |
| Mean D vs MD excess entropy, all mixtures | 4 | 0.867 | 0.251 |
| Solvent D vs partial pair entropy, original figure constants | 6 | 1.031 | 0.874 |

These are descriptive fits without entropy uncertainties. The manuscript's MSD linear-fit R²≈0.999 measures a different relationship. It does not validate an entropy law or provide independent-sample uncertainty. The solvent fit previously claimed R²=0.91, which was hard-coded in the plot legend; the six constants give 0.874. Their extraction needs documentation, and the 3C value uses mean diffusion as a solvent proxy.

Mean self-diffusion is a solvent-dominated weighted average, not collective diffusion, mutual diffusion or conductivity. An approximate law for the average need not hold for each species. Species-dependent kinetic prefactors and local environments can matter. Generalized tracer scaling is worth testing, but its dilute-limit construction cannot simply be transferred to unscreened bare Coulomb interactions without addressing electroneutrality and screening. [Krekelberg et al., PRE 80, 061205 (2009)](https://arxiv.org/abs/0910.0280) both develops mixture tracer scaling and documents limitations when mobilities decouple.

## Issues to resolve before interpreting ion trends

1. **Reproducibility and reference frame.** The original MSD ion range, 0.1666–0.2289, contrasts with GK values of 0.2017–0.2057. The second MSD dataset reverses the apparent weak-solvent slowdown. Each sign has only 60 ions. Charge inversion requires equal ensemble-mean cation/anion diffusivities; large differences are a convergence diagnostic. Pooling signs improves sampling but is not the diffusion of neutral bound pairs.

   Archived `in.run` files use `msd com yes` separately for every group. This removes each group's center-of-mass displacement, as specified in the [LAMMPS documentation](https://docs.lammps.org/compute_msd.html), rather than a common global drift. Recompute all species in the same frame and verify weighted sums. For example, Table IV's 3C ionic values average to 0.2012, whereas its pooled value is 0.2027; for 5C the corresponding values are 0.19955 and 0.2036. These differences exceed rounding. Group-specific corrections can contribute and should be measured rather than assumed to explain everything.

2. **Fixed Hamiltonian for TI.** At low density the archived input uses `pair_style lj/cut/coul/long 8.0` and `pair_coeff * * 1.0 1.0 8.0`, versus 2.5 at high density. The LJ Hamiltonian therefore changes along the integration path. Either repeat the density sweep with a fixed LJ cutoff or add explicit free-energy corrections. Changing only the real-space split of a converged Ewald electrostatic sum would be a different matter.

3. **Electrolyte dilute limit.** The TI scripts extrapolate `(Z−1)/rho` to a finite polynomial intercept. Screened-electrolyte Debye–Hückel asymptotics instead give an integrable `rho^(−1/2)` behavior. Test an appropriate low-density treatment and its effect on composition differences, not only absolute entropies. The current fit range in entropy is only about 0.014 kB.

4. **2PT and thermodynamic consistency.** A 2PT species decomposition is not automatically a partial derivative of a chemical potential. Entropy at the same state does not change merely because one uses a canonical calculation. The −5.38 versus −1.88 kB discrepancy needs an ideal-reference, mode-weight and normalization audit. As written, the manuscript's VACF definition and Fourier prefactor yield unit integrated spectral weight, although three is claimed. Small response of a computed entropy to a varied diffusion input does not establish transport convergence.

   Check partial-property sum rules using consistent ideal references and derivative constraints. The reported 3C partial volumes give Σxi Vi=1.44344 rather than 1/rho=1.42857. The pair entropies in its partial table give Σxi s2i≈−1.61405 rather than the bulk −1.5706. These discrepancies exceed rounding and suggest mixed sources or definitions. Temperature derivatives must hold physical interactions fixed and use dimensional chemical potentials correctly. Single-ion chemical potentials also require a specified electrostatic convention; neutral salt combinations are the safer thermodynamic benchmark.

5. **Mechanism.** Neutral LJ solvent has no explicit dielectric response. Weakening its LJ energy does not directly change the Coulomb dielectric prefactor. Structural reorganization can alter screening and association, but a contact peak is not a pair population or a lifetime. Avoid asserting that poorer dielectric shielding or persistent Bjerrum pairs caused the first MSD trend.

## Prioritized further checks

| Priority | Comparison | What would resolve the question |
|---|---|---|
| 1 | Independent seeds, longer trajectories, multiple time origins; common-frame MSD and GK | Replica/block confidence intervals; stable fitting windows and VACF integrals; charge symmetry; species/pool sum rules |
| 1 | Fixed-LJ-cutoff density sweep with screened dilute-limit treatment | Reliable composition-dependent excess entropy and its uncertainty |
| 2 | Identical-solvent label split, ε2=1 | Direct test that adding a label does not accelerate ions |
| 2 | Add symmetric 1:1 salt at fixed x+2=x−2=0.03, T and rho | Dose-response for the proposed screening/association mechanism; separately compare constant-I replacement states |
| 2 | Coordination, cluster distributions, continuous/intermittent association lifetimes, solvent residence times | Evidence that any reproducible speedup accompanies weakened or shorter-lived association; sensitivity to shell definition |
| 3 | Increase each divalent fraction from 0.03 to 0.06 and 0.10, replacing solvent | Greater ionic weight and sampling; quantify competing screening, crowding and association; check phase stability |
| 3 | Equal-sign ion diameters / solvent diameter = 0.8, 1.0, 1.2 at fixed masses | Separate packing and local coupling; report packing fraction and add matched contact-coupling controls |
| 3 | N=2000, 4000, 8000 at fixed physical state | Extrapolate D vs 1/L; distinguish box effects from changes in particle size |
| 4 | Several densities and temperatures; matched-excess-entropy composition pairs | Test predictiveness on held-out states, allowing distinct species prefactors; distinguish density trends from composition trends |

For cubic boxes the leading hydrodynamic correction is `D∞ ≈ DL + 2.837297 kBT/(6πηL)`, using measured mixture viscosity and checking the observed size trend. See [Yeh and Hummer (2004)](https://doi.org/10.1021/jp0477147). Thermostat sensitivity and run lengths relative to measured association lifetimes should be checked alongside finite size.

The revised conclusion order is: mixing versus excess entropy; the solvent dispersion trend and limited mean scaling evidence; unresolved divalent transport; supporting thermodynamic results and remaining audits. This places the scientific question before methodological achievements.

## Reproduction and files

Run `MPLCONFIGDIR=/tmp/entropy-mpl python3 doc/audit_entropy_diffusion.py` to reproduce the fit and mixing calculations, `entropy_diffusion_audit.json`, and `figures/rosenfeld_scaling_review.png`. The revised manuscript uses this new figure; the original figure is preserved. The audit uses reported values, not a new trajectory analysis.

The cluster provides `texlive/20260310-GCC-13.2.0`. Build the revised PDFs with `bash doc/build_pdfs.sh`. The historical theoretical report is otherwise preserved.
