# Fixed-LJ sweep: completed analysis

All **60 runs** (four mixtures × 15 densities, jobs 100493–100552) completed normally and reached their expected production length. Inputs and configurations were checked against the submission manifest. Original outputs remain untouched.

**The cutoff error was material for absolute excess entropy, but fixing it does not establish enhanced divalent-ion diffusion through Rosenfeld scaling.** It lowers dense-state excess entropy by about 0.085–0.091 kB per particle, removes approximately half of the former HNC–MD entropy discrepancy, and leaves the main qualitative interpretation intact.

## Isolating the cutoff correction

The following comparison uses the *same original postprocessing* on both datasets: final half of production, stored mean instantaneous Z, and a cubic spline with a finite fitted low-density intercept. This isolates the changed LJ potential from subsequent estimator improvements.

| System | Original s_ex/kB | Fixed LJ, same estimator | Cutoff-only shift |
|---|---:|---:|---:|
| 3C | -1.87985 | -1.96702 | -0.08717 |
| 4C, ε2=1.25 | -1.89091 | -1.97655 | -0.08565 |
| 4C, ε2=0.75 | -1.87730 | -1.96230 | -0.08501 |
| 5C | -1.88535 | -1.97642 | -0.09107 |

For rho≥0.20, the original LJ cutoff was already 2.5 sigma. Those reruns reproduce the original thermodynamic output, and the dense MSD arrays are numerically identical. They are **not independent transport replicas**. At rho=0.01–0.15, the corrected pressure is higher because the spurious longer-range LJ attraction is removed. Its integral raises excess free energy and consequently lowers excess entropy at all higher densities. The effect exceeds the conditional sampling errors by a large margin.

## Refined thermodynamic estimates

The primary calculation additionally uses **Z = mean(P)/(rho T_target)** with T_target=2, rather than mean[P/(rho T_instantaneous)]. The latter ratio is a different finite-sample estimator; converting to target-temperature pressure changes dense entropy by up to 0.0023 kB. This matters for small composition differences.

The integration variable is x=sqrt(rho), with F(x)=2[Z(x²)−1]/x. Debye–Hückel asymptotics give F(0)=−C, C=(4π lB Σxi qi²)^(3/2)/(12π). A cubic spline in x then integrates a finite endpoint. In the actual LAMMPS and HNC inputs the Coulomb **energy** prefactor is 1, so lB=1/T=0.5; calling the energy prefactor itself the Bjerrum length would be incorrect. The manuscript definition was clarified. See the [LAMMPS reduced-unit convention](https://docs.lammps.org/units.html).

| System | Primary s_ex/kB | 95% sampling interval | HNC s_ex/kB |
|---|---:|---:|---:|
| 3C | -1.96693 | [-1.97076, -1.96316] | -2.06350 |
| 4C, ε2=1.25 | -1.97261 | [-1.97738, -1.96814] | -2.07030 |
| 4C, ε2=0.75 | -1.96199 | [-1.96513, -1.95871] | -2.05510 |
| 5C | -1.97588 | [-1.98067, -1.97107] | -2.06770 |

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
| Original | 1.045 | 0.392 | 0.867 | 0.251 |
| Fixed LJ, original estimator | 1.228 | 0.544 | 0.585 | 0.157 |
| Fixed LJ + screened integration, stored Z | 1.227 | 0.546 | 0.677 | 0.194 |
| Fixed LJ + screened integration, target-T pressure | 1.891 | 0.692 | 0.537 | 0.099 |

The solvent-only fit improves, but it has just three points and is sensitive to entropy-estimator choices. The four-mixture fit remains poor. Fitting errors in both coordinates and obtaining independent transport uncertainties would be necessary before treating any slope as predictive. The negative descriptive slope for the ion MSD values is not evidence of a validated inverse entropy law.

## Next analysis priorities

1. Recompute species MSD and VACF in one common global reference frame, with multiple time origins and independent trajectory seeds. Verify cation/anion symmetry and species/pool additivity. The cutoff reruns do not supply new independent evidence for ion transport.
2. Extend the dilute grid and check independent thermodynamic replicas, especially for 5C, whose averaging-window sensitivity is largest.
3. Only then test whether salt addition changes association fractions and lifetimes, followed by the concentration, particle-diameter and box-size controls described in the manuscript.

The 2PT normalization/reference discrepancy is independent of this LJ cutoff correction and remains unresolved.

## Files and reproduction

Run `MPLCONFIGDIR=/tmp/entropy-mpl OPENBLAS_NUM_THREADS=1 python3 doc/analyze_fixed_lj_sweep.py` from the repository. The script validates all 60 runs and reproduces the calculations and figures. Full density results, paired-bootstrap samples, and machine-readable summaries are in `lammps_study/fixed_lj_20261030/results/`. Compact copies are in `doc/fixed_lj_dense_comparison.csv` and `doc/fixed_lj_analysis.json`.

The manuscript source, abstract, thermodynamic table, figures and conclusions have been updated. Old results and figures remain available. Build the PDFs with `bash doc/build_pdfs.sh`, which loads `texlive/20260310-GCC-13.2.0` and uses the system zlib with direct pdfLaTeX passes.
