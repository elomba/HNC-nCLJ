# HNC-guided component/entropy/diffusion study

**Status:** Completed: 35 HNC screens, eight refinements, and all 24 production allocations (100669–100692), with all 456 simulation stages verified. See [completed conclusions](results/conclusions.md), [entropy–diffusion plot](results/entropy_diffusion.pdf), and [diagnostics](results/diagnostics.json). Higher-excess-entropy solvent designs accelerate divalent ions by 18–20% relative to the original baselines; the controlled solvent-diversification effect is only about 5% and remains sensitive to transport convergence. Label-only splitting has exactly zero transport effect.

The question is whether physically adding distinct components can raise *excess* entropy and increase the self-diffusion of unchanged +2/−2 probe ions. Component count alone is not a mechanical control variable: identical labels provide an exact null comparison. We therefore separate formal ideal mixing entropy, interaction-induced excess entropy, and transport.

## Completed HNC screening

35 integer-charge/mole-fraction/dispersion candidates were solved at T*=2, rho*=0.7. The coarse grid used 2048 radial points, dr=.01, Newton dimension parameter 40. All converged. Eight selected systems were refined using the manuscript grid (4096 points, dr=.005, parameter 100). Their entropies agree with the coarse results at the reported 0.0001 kB precision. The first 20 sampled structure-factor matrices are positive for all selected fluids. This is a screening diagnostic, not a proof against demixing; MD will monitor concentration fluctuations and stationarity.

The charge-only screen included symmetric monovalent additions and mixtures of conjugate 1:3 salts. For the latter, solvent and divalent fractions, ion fraction, ionic strength, and every even charge moment are unchanged when conjugate salts are mixed. Five→seven charged-component mixtures nevertheless give slightly *lower* excess entropy (e.g. −2.1384→−2.1424 at 60% background salt). They increase ideal mixing entropy, but do not provide the desired positive excess-entropy signal. Their inputs and outputs are preserved; expensive MD is prioritized for the larger positive signals and the original monovalent-salt comparison.

## Selected compositions

Every row contains **x(+2)=x(−2)=0.03**, masses and diameters 1, primary-solvent ε=1, and ion ε=1. All charges are integers. Fractions below sum to one when the probes are included. I=½Σxi qi² is the mole-fraction ionic-strength measure; actual screening also depends on density.

| Code | Labels / distinct mechanical species | Remaining composition | I | HNC s_ex/kB |
|---|---:|---|---:|---:|
| B3 | 3 / 3 | 94% standard solvent | 0.12 | -2.0635 |
| B5 | 5 / 5 | 84% standard solvent; 5% each ±1 | 0.17 | -2.0677 |
| U4 | 4 / 4 | 64% standard + 30% ε=.25 solvent | 0.12 | -1.8434 |
| L5 | 5 / 4 | 64% standard + 15% ε=.25 + 15% ε=.25 (labels) | 0.12 | -1.8434 |
| M5 | 5 / 5 | 64% standard + 15% ε=.0625 + 15% ε=.5625 | 0.12 | -1.8102 |
| H5 | 5 / 5 | 64% standard + 15% ε=.01 + 15% ε=.81 | 0.12 | -1.7327 |
| U6 | 6 / 6 | 54% standard + 30% ε=.25; 5% each ±1 | 0.17 | -1.8478 |
| H7 | 7 / 7 | 54% standard + 15% ε=.01 + 15% ε=.81; 5% each ±1 | 0.17 | -1.7370 |

B3 and B5 reproduce the original 3C and monovalent-background 5C electrolytes, now with consistent entropy integration, larger particle populations, independent replicas and a common transport reference frame.

**Best screened candidates:** H5 without monovalent salt and H7 with the original monovalent background. HNC predicts excess-entropy increases of 0.3308 and 0.3307 kB relative to B3 and B5. These are finite-grid screening optima, not global optima over all possible models.

## What the comparisons control

1. **U4 versus L5: label-only null.** Splitting the ε=.25 solvent into two identical labels leaves forces and HNC excess entropy unchanged. If the labels are treated as distinguishable, the formal ideal mixing contribution rises by .30 ln(2)=.20794 kB. Dynamics should agree within replica uncertainty. A significant apparent difference would require a convergence/implementation audit before interpreting component-count effects.
2. **U4 versus M5/H5: physical diversification at fixed mean bare LJ amplitude.** All have 64% primary solvent, 30% secondary neutral material and 6% divalent probes. U4 has sqrt(ε)=.5 throughout the secondary solvent. M5 splits it equally into sqrt(ε)=.25 and .75; H5 into .1 and .9. Hence Σxi sqrt(εi)=.85, and the mole-fraction weighted mean pair amplitude Σij xi xj εij=(.85)²=.7225, is identical. Probe fractions, ion fraction, ionic strength, nominal diameters and masses are fixed. HNC predicts Δs_ex=+.0332 and +.1107 kB relative to U4.
3. **U6 versus H7:** repeats the uniform-versus-diverse comparison with 5% each of ±1 background ions. Solvent fraction .84, ion fraction .16, I=.17, divalent fraction .06 and the same mean bare LJ amplitude are fixed. HNC predicts Δs_ex=+.1108 kB.
4. **B3/U4/H5 and B5/U6/H7:** larger entropy ranges for testing mobility ordering. Their mean interaction amplitudes change, so these are design comparisons, not a pure component-count intervention.
5. **B3 versus B5:** an independently replicated version of the original disputed salt-addition result. It is not a fixed-ionic-strength comparison.

Mean bare interaction amplitude is not the actual mean energy, collision rate, effective excluded volume or local solvation structure. The unlike interactions and their distribution still change. Thus even a positive result supports a restricted entropy–mobility relation, not causation by the number of labels alone. H5's very weak solvent is deliberately a strong model perturbation; M5 supplies a more moderate counterpart. Every result must be checked for homogeneity.

## Prespecified predictions and decisions

Use ln(D_reduced), with D_reduced=D rho^(1/3)/sqrt(T) for all species (total density; all masses 1). At the same T and rho, diffusion ratios equal reduced-diffusion ratios. Using the *previous*, descriptive HNC solvent-series slope B≈1.4 gives provisional ratios H5/B3≈1.59, H5/U4≈1.17 and M5/U4≈1.05. These are hypotheses to test, not reliable quantitative forecasts from a three-point fit.

Primary endpoints are paired replica differences in independently measured MD excess entropy, ln(D_probe), individual ln(D+2)/ln(D−2), and ln(D_net). D_probe is the particle-averaged self-diffusion of the two divalent signs, not diffusion of a bound neutral pair. Mean self-diffusion is not conductivity or mutual diffusion.

- If physical diversification raises MD excess entropy and probe diffusion reproducibly while the label null passes, it supports an entropy–mobility association in this design family.
- If entropy rises but only the solvent-dominated mean diffusivity increases, bulk entropy fails to predict probe mobility under these interventions.
- If the label null fails, or cation/anion symmetry and estimator convergence fail, transport conclusions remain unresolved.
- If MD entropy ordering disagrees with HNC, interpret transport against MD entropy; revise the HNC screening interpretation.
- Fit the fixed-I families separately. Treat L5 as a null control, not an extra independent entropy state. Do not equate MSD linear-fit R² with entropy-scaling fit quality.

## GPU calculations

**24 production allocations: eight mixtures × three replicas**, N=4000, with 120 particles of each divalent sign. Independent seeds are 64123, 76237 and 91817; configurations and probe IDs are matched across compositions within a replica, enabling paired contrasts. Replicas within each mixture are independent. Whole density curves are paired by replica when estimating uncertainty.

Each allocation runs 18 lower-density thermodynamic states and a longer dense transport stage, supplying **19 density points**: .001, .0025, .005, .01, .02, .05, .10, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65 and .70. The transport trajectory supplies the .70 thermodynamic point, avoiding a redundant dense run. Total: **456 production simulation stages**.

Every state uses shifted LJ cutoff 2.5 with no tail correction. PPPM real-space cutoffs remain separate: r_Coul=5*(.7/rho)^(1/3), always below half the box length; accuracy 1e−5. This keeps the numerical split at a fixed fraction of the box and avoids excessively large dilute-state FFT grids. Dense benchmarks at the same accuracy favored r_Coul=5 over 2.5 and 8; the optimized dense run is about nine times faster than the initial 2.5 setting. The LJ cutoff remains exactly 2.5 in every pair coefficient. T*=2, dt=.002, NVT damping .1, Coulomb energy prefactor 1 (Bjerrum length .5). EOS stages use 50,000 equilibration and 100,000 production steps. Dense transport uses 100,000 equilibration and 500,000 production steps (1000 tau), with high-frequency velocities and ion positions during the final 100 tau. Stage-specific restart files protect completed work.

Queue gpu_v4_short, ladon28/29, one Blackwell GPU, two MPI tasks, 8 GB host RAM, 12-hour wall limit per allocation. The scheduler enforces the user's concurrent-GPU limit. Allocation scripts run stages sequentially to reduce scheduler overhead. Smoke jobs validate all compositions at low and high density before production is submitted.

## Analysis and reproducibility

`screen_hnc.py` and its input definitions reproduce all HNC screening; `hnc_screen_summary.csv`, `hnc_fine_summary.csv` and `hnc_stability_diagnostics.json` record results. `selected.json` records exact species definitions. `prepare_md.py` generates immutable inputs/configurations; each allocation manifest contains checksums and seeds. `submit.py` validates scripts with Slurm and records actual IDs in `submissions.jsonl`. Failed smoke attempts are retained separately and are not scientific results.

The common-frame online MSD subtracts the **same full-system COM displacement** for every species. NetCDF files store genuinely unwrapped coordinates. Offline MSD uses multiple time origins with several fitting windows; whole-system COM velocity is recorded separately for matched VACF corrections. All ionic velocities are retained, together with up to 200 representative particles per neutral species. Net self-diffusion from VACF must use mole-fraction weighting of species estimates. Cation/anion symmetry, species/pool additivity, time-window stability and GK cutoff sensitivity are diagnostics; no per-species COM subtraction is used.

At the end of each allocation, `analyze.py --replica <directory>` calculates EOS entropy, multi-origin MSD and species VACF integrals and updates `results/summary.json` and `results/progress.md` under a file lock. Entropy uses the configurational virial contribution and an analytic target-temperature ideal term, followed by a sqrt-density integration with the screened Debye–Hückel endpoint. For the default LAMMPS temperature DOF=3N−3, the recorded fields give Z_excess_reference=1+[mean(P)−rho*(1−1/N)*mean(T_inst)]/(rho*T_target). This removes kinetic sampling noise and the finite-N COM contribution from the ideal reference; the total-pressure estimator is retained as a diagnostic. The [LAMMPS pressure definition](https://docs.lammps.org/compute_pressure.html) specifies the DOF dependence. Finite-size changes in configurational correlations are not corrected by this identity. PCHIP and averaging-window alternatives are reported. Confidence intervals use paired differences across independent replicas with Student-t intervals, not overlapping time origins. Three replicas are an initial uncertainty estimate, not a guarantee of convergence.

The trajectory files also permit association-lifetime and local-solvation follow-up if ions decouple from bulk entropy. Demixing diagnostics use low-wavevector concentration fluctuations. If results remain inconclusive, the next stage should add independent replicas, longer association-time sampling and a box-size comparison at fixed physical state before expanding the composition grid. No extra adaptive jobs are submitted automatically.

## Commands

```bash
python3 ComponentEntropyStudy_20261030/screen_hnc.py
python3 ComponentEntropyStudy_20261030/prepare_md.py
python3 ComponentEntropyStudy_20261030/submit.py --submit
python3 ComponentEntropyStudy_20261030/analyze.py
```

Execute native HNC and Slurm commands on a configured cluster node. Submission is ledger-protected against duplicates. Preparation refuses existing replica directories. Use recorded IDs and manifests to inspect progress; do not regenerate active inputs.


## Completed 2PT rehabilitation follow-up

All 24 existing velocity records have been reanalyzed with a corrected classical 2PT implementation. Analytic and label-invariance checks pass, but composition-dependent entropy biases remain, including an overestimated diversification increment with monovalent background. TI therefore remains the primary entropy coordinate. See [the rehabilitation report](two_pt_rehabilitation/report.md), [formulation](two_pt_rehabilitation/METHODS.md), and [reproducible tests](two_pt_rehabilitation/README.md).
