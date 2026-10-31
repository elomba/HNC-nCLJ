# Completed HNC-guided MD study: conclusions

All 24 Slurm allocations (100669–100692) completed with exit code 0; all 456 stages have verification records. Eight mixtures, three independent replicas each, N=4000, T*=2 and rho*=0.7. The LJ cutoff is fixed at 2.5 throughout the entropy integration. Reported entropy is excess entropy per particle relative to the ideal mixture of the same composition; ideal mixing entropy is excluded. Diffusion coefficients use LJ units, with a common whole-system COM frame.

**Conclusion:** the selected interaction-diverse solvents raise excess entropy and divalent-ion mobility relative to the original baselines. The controlled diversification effect on the ions is smaller (about 5%) and remains provisional: its significance is sensitive to fitting window, and it is unresolved with monovalent background. Increasing the number of labels alone has exactly zero transport effect. The results support an entropy–mobility association within these model families, not a causal law based on component count.

## State averages

Composition definitions are in ../README.md. D_probe averages the +2 and −2 self-diffusion coefficients; D_net is the composition-weighted mean self-diffusion, not mutual diffusion or conductivity.

| State | MD s_ex/kB | HNC s_ex/kB | D_net | D_probe |
|---|---:|---:|---:|---:|
| B3 | -1.97101 | -2.06350 | 0.22253 | 0.20677 |
| B5 | -1.97421 | -2.06770 | 0.22279 | 0.21047 |
| U4 | -1.77209 | -1.84340 | 0.25945 | 0.23684 |
| L5 | -1.77166 | -1.84340 | 0.25945 | 0.23684 |
| M5 | -1.74170 | -1.81020 | 0.26720 | 0.24038 |
| H5 | -1.66887 | -1.73270 | 0.28401 | 0.24825 |
| U6 | -1.77603 | -1.84780 | 0.26018 | 0.23572 |
| H7 | -1.67272 | -1.73700 | 0.28196 | 0.24795 |

## Paired comparisons

Percent changes are exp(mean paired log ratio)−1; brackets are nominal 95% Student-t intervals across three independent replicas (2 degrees of freedom), without multiple-comparison correction. They do not include finite-size or estimator bias. Entropy differences are in kB per particle.

| Comparison | Delta s_ex/kB | D_net change [95% CI] | D_probe change [95% CI] |
|---|---:|---:|---:|
| L5 − U4: identical labels | +0.00042 [-0.00219, +0.00304] | +0.00% [+0.00, +0.00] | +0.00% [+0.00, +0.00] |
| M5 − U4: moderate diversification | +0.03039 [+0.02727, +0.03352] | +2.99% [+1.38, +4.62] | +1.48% [-1.57, +4.64] |
| H5 − U4: strong diversification | +0.10322 [+0.10209, +0.10434] | +9.47% [+8.96, +9.97] | +4.81% [+0.60, +9.19] |
| H5 − B3: full design change | +0.30214 [+0.29927, +0.30501] | +27.63% [+25.40, +29.90] | +20.05% [+14.30, +26.09] |
| H7 − U6: diversification with ±1 | +0.10331 [+0.09990, +0.10672] | +8.37% [+7.31, +9.44] | +5.20% [-6.14, +17.92] |
| H7 − B5: full design with ±1 | +0.30149 [+0.29823, +0.30476] | +26.56% [+25.99, +27.13] | +17.79% [+11.17, +24.81] |
| B5 − B3: add ±1 salt | -0.00320 [-0.00543, -0.00097] | +0.12% [-1.79, +2.06] | +1.79% [-1.37, +5.04] |

## Interpretation in order of evidential strength

1. **Component labels alone do not accelerate ions.** U4 and L5 have identical mechanical interactions, and their transport trajectories give exactly identical diffusion in every replica. Their small entropy difference (+0.00042, interval including zero) is consistent with thermodynamic sampling. Formal ideal mixing rises by 0.30 ln(2)=0.20794 kB if labels are counted as distinguishable; that contribution is absent from s_ex and produces no force change.
2. **HNC successfully selected higher-excess-entropy liquids.** MD reproduces the selected-state entropy ordering. HNC predicts H5−U4=+0.1107 versus MD +0.1032, and H7−U6=+0.1108 versus MD +0.1033. Absolute HNC entropies are more negative by roughly 0.064–0.093 kB; agreement of increments is stronger than agreement of absolute values.
3. **Large design changes reproducibly accelerate both ion signs.** H5 versus B3 increases pooled divalent diffusion by about 20%, and H7 versus B5 by about 18%. At the primary window, the separate +2 and −2 changes are +16.8% and +23.4% for H5/B3 and +18.0% and +17.6% for H7/B5; all four paired nominal 95% intervals exclude zero. These interventions also change mean interaction strength and local solvent friction, so they cannot isolate entropy as the cause.
4. **Controlled physical diversification has a smaller ion response than the mean response.** U4→H5 preserves composition of charged particles, masses, nominal diameters, ionic strength and the mean bare LJ pair amplitude. It raises s_ex by 0.1032 kB, mean diffusion by 9.47%, and ion diffusion by 4.81%. The primary-window ion interval barely excludes zero. With ±1 background the analogous means are +0.1033 kB, +8.37%, and +5.20%, but the ion interval includes zero. The moderate split raises entropy clearly but its ~1.5% ion response is unresolved.
5. **The original monovalent-addition effect is not established.** B3→B5 slightly lowers excess entropy (−0.00320 kB) and gives no resolved diffusion change. More charged components alone therefore do not ensure higher excess entropy or faster divalent ions.
6. **Rosenfeld scaling is descriptive, not universal here.** At common T, rho and mass the reduced diffusion D^R=D rho^(1/3)/sqrt(T) differs from D by a constant. The broad state ordering is consistent with a positive ln(D^R)–s_ex trend for both mean and probe diffusion. The controlled H5/U4 local response is B≈0.88 for the mean and B≈0.46 for the probes, substantially below the prior provisional B≈1.4 forecast; H7/U6 gives about 0.78 and 0.49. A single bulk entropy does not determine species mobility. Descriptive fits of state means give B=0.80 (mean) and 0.62 (probes) without monovalent background, and 0.78 and 0.55 with it; R² exceeds 0.99 in all four fits. These three- or four-state fits are dominated by the broad interaction changes and do not validate a universal relation. See diagnostics.json for the fit parameters.

## Convergence and limitations

### MSD fitting-window sensitivity

| Comparison | 10–50 tau | 20–100 tau (primary) | 50–200 tau |
|---|---:|---:|---:|
| H5/U4 | +5.92% [+1.87, +10.13] | +4.81% [+0.60, +9.19] | +3.94% [-8.44, +18.00] |
| H7/U6 | +5.31% [-1.17, +12.22] | +5.20% [-6.14, +17.92] | +4.46% [-5.66, +15.67] |

The positive ~4–6% ion response survives in point estimates, but H5/U4 loses nominal significance in the later window. This is evidence for a trend, not a converged precise 5% effect. VACF integrals at 5 and 10 tau support the broad mobility ordering, but the final-100-tau velocity records are too short for a precise independent small-effect determination. No reliable 20-tau plateau is claimed.

### Charge-sign symmetry

The sign-symmetric B3 Hamiltonian should have equal long-time +2 and −2 self-diffusion. The observed difference D(+2)−D(−2)=0.01039, nominal replica interval [0.00291,0.01786], flags remaining sampling/uncertainty limitations. The other states have sign-difference intervals containing zero. Inspection of generated inputs found equal masses/LJ parameters for both probes, random probe assignments and one all-particle thermostat/common COM frame; no sign-dependent input parameter was found. This does not establish that the residual asymmetry is harmless. Pooling both signs is physically appropriate but does not make the convergence flag disappear. Both signs independently speed up strongly in the full baseline contrasts.

### Entropy and homogeneity diagnostics

Across all replicas, the largest change from cubic-spline to PCHIP integration is 0.00032 kB; the largest change when retaining all/75% versus 50% of production is 0.00165 kB. These are small compared with the ~0.10 kB diversification increments. They do not test the unsampled dilute-density endpoint or finite-size bias.

Low-k concentration fluctuation diagnostics remain of order unity (individual species values approximately 0.48–1.80), without a large fluctuation growth across these runs. One box size and short density-stage trajectories cannot establish single-phase stability or exclude slower segregation. Full species diagnostics are saved in diagnostics.json.

## Most useful next checks

1. Add independent, longer matched U4/H5 and U6/H7 replicas to resolve the ~5% effect; extend B3 to test the sign-symmetry flag. Keep the same estimators and prespecified windows. More overlapping time origins are not extra independent replicas.
2. Repeat at a larger N at fixed density/composition (e.g. 8000 or 16000), matching all cutoffs and PPPM accuracy. Compare both signs and solvent mobility; the expected finite-size transport correction need not cancel when viscosity differs.
3. From existing trajectories, measure probe–solvent/ion RDFs, association populations and residence times, and species-resolved pair excess entropy. These may distinguish weaker solvation/caging from a bulk entropy proxy. A pair-entropy approximation must not be relabeled as total thermodynamic excess entropy.
4. After convergence checks, raise divalent concentration with matched uniform/diverse solvents at each concentration. Re-screen with HNC; maintain neutrality and compare within the same ionic strength at each matched pair. Higher concentration can change pairing and introduce species-specific decoupling rather than strengthen the present trend.
5. Test several densities/temperatures within each fixed composition to separate state-point Rosenfeld scaling from composition-induced trends. Varying nominal particle size is a distinct perturbation of packing and friction and should be a separate controlled series.

No new MD jobs were submitted during this analysis. All raw trajectories and original automatic results are preserved. Reproduce this report and plots with `python3 ComponentEntropyStudy_20261030/summarize_results.py`.
