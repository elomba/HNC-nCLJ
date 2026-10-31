# 2PT rehabilitation: completed tests

**Assessment: the implementation can be rehabilitated as an internally consistent classical spectral estimator, but it does not pass the prespecified quantitative excess-entropy validation for these mixtures.** It captures the broad entropy ordering and one controlled increment reasonably well. Its absolute bias is composition dependent, and it overestimates the controlled increment with monovalent background. TI remains the primary entropy estimate.

Figures: [absolute-entropy validation](entropy_validation.pdf), [entropy-increment validation](increment_validation.pdf), and [probe dynamics](probe_dynamics.pdf).

All 24 existing velocity records (eight mixtures × three independent replicas, N=4000) were analyzed. No new MD jobs were required. The primary model and practical thresholds were recorded in protocol.json before computing the baseline result. No parameter, entropy offset, partial volume, or partition exponent was fitted to TI.

State codes follow the [study compositions](../README.md#selected-compositions): B3/B5 are baselines without/with monovalent background; U4 has 30% uniform secondary solvent (epsilon=0.25); L5 is its label-only split; M5/H5 have physically split secondary solvents; U6/H7 repeat uniform/strong diversification with monovalent background. All contain 3% of each divalent sign.

## Formulation and scope

The new calculation follows the classical signed-residual R2PT construction of [Sun et al., JCP 147, 194505 (2017)](https://arxiv.org/abs/1711.02041), with partition exponent δ=1 primary and the literature δ=1.5 as an unfitted sensitivity. Species use the equal-partial-volume approximation, Vbar_i=1/rho; this is the one-fluid volume approximation discussed by [Lai et al., PCCP 14, 15206 (2012)](https://doi.org/10.1039/C2CP42011B). Equal nominal diameters do not prove equal actual partial volumes in charged, dispersion-diverse liquids. These are allocated spectral species contributions, not independently measured chemical-potential derivatives.

The density-scaled diffusivity is Δ=(8/3)(36/π²)^(1/3) D sqrt(πm/T) rho_eff^(1/3). The archived implementation instead used rho^(-1/3) and a different prefactor. Gas entropy is weighted by its fraction f. Classical oscillator weights, ideal-mixture subtraction, and the analytical Lorentzian tail are treated consistently. Quantum weights are excluded when comparing against classical TI.

The entire biased finite-record VACF is used, without subtracting each particle’s time-averaged velocity. One recorded system COM velocity is removed. Gaussian tapering preserves nonnegative total spectra; the full Nyquist band is retained. The primary taper is 10 tau, with 3, 5, 20 tau and untapered sensitivities. Spectral normalization uses measured C(0); raw target-temperature mode weights are reported, so the enforced sum rule is not mistaken for independent accuracy evidence.

The gas Lorentzian integrates analytically to 3f modes. Its entropy tail is integrated to infinity, and the residual G−Ggas remains signed as R2PT bookkeeping. It is not clipped to create a nonnegative “solid population.” Ideal mixing is added to the approximate absolute mixture entropy and subtracted in the ideal reference, cancelling exactly; arbitrary classical phase-space normalization h also cancels.

## Absolute entropy benchmark

Values are kB per particle. Parenthesized model columns are means, not independent thermodynamic determinations. Bias intervals are nominal paired 95% Student-t intervals across three replicas, df=2; no multiplicity correction or systematic error is included.

| State | TI | Primary R2PT, δ=1 | R2PT, δ=1.5 | Primary − TI [95% CI] |
|---|---:|---:|---:|---:|
| B3 | -1.97101 | -2.07445 | -1.91843 | -0.10344 [-0.10708, -0.09979] |
| B5 | -1.97421 | -2.07970 | -1.92590 | -0.10549 [-0.14993, -0.06104] |
| U4 | -1.77209 | -1.82109 | -1.68451 | -0.04900 [-0.07291, -0.02510] |
| L5 | -1.77166 | -1.82109 | -1.68451 | -0.04943 [-0.07089, -0.02796] |
| M5 | -1.74170 | -1.78826 | -1.65680 | -0.04657 [-0.08525, -0.00788] |
| H5 | -1.66887 | -1.70291 | -1.58054 | -0.03404 [-0.06069, -0.00740] |
| U6 | -1.77603 | -1.83138 | -1.69766 | -0.05535 [-0.08565, -0.02504] |
| H7 | -1.67272 | -1.69210 | -1.56681 | -0.01938 [-0.03822, -0.00054] |

The old baseline 2PT result was about −5.376, compared with corrected TI near −1.97. The corrected primary baseline is now about −2.074: the discrepancy is reduced from roughly 3.4 to 0.10 kB. This is substantial rehabilitation of the implementation, but a 0.10 kB error is as large as the controlled entropy signal. It is not an ensemble effect.

Against the prespecified absolute-bias tolerance 0.05 kB, the independent physical states that fail are B3, B5, U6. A constant entropy-offset correction would not solve the composition dependence. Changing the partition exponent moves absolute entropies substantially and does not make every mixture agree with TI. All primary absolute-bias intervals lie below zero at nominal 95% confidence, including states whose point bias meets the practical tolerance. Meeting that tolerance is not a statistical demonstration of unbiasedness.

## Entropy increments

| Comparison | TI increment | Primary increment | Primary increment bias [95% CI] | δ=1.5 increment |
|---|---:|---:|---:|---:|
| L5/U4 | +0.00042 | -0.00000 | -0.00042 [-0.00304, +0.00219] | -0.00000 |
| M5/U4 | +0.03039 | +0.03283 | +0.00244 [-0.04878, +0.05366] | +0.02771 |
| H5/U4 | +0.10322 | +0.11818 | +0.01496 [-0.00544, +0.03536] | +0.10397 |
| H5/B3 | +0.30214 | +0.37153 | +0.06939 [+0.03954, +0.09925] | +0.33789 |
| H7/U6 | +0.10331 | +0.13927 | +0.03596 [+0.01501, +0.05692] | +0.13086 |
| H7/B5 | +0.30149 | +0.38760 | +0.08611 [+0.06032, +0.11190] | +0.35910 |
| B5/B3 | -0.00320 | -0.00526 | -0.00205 [-0.04302, +0.03892] | -0.00747 |

For H5/U4, the primary increment is 0.11818 versus TI 0.10322. Its +0.01496 point bias meets the practical 0.02 kB threshold and its interval contains zero, but the wide interval does not establish equivalence within 0.02 kB. The δ=1.5 increment is 0.10397, but this good agreement in one comparison does not validate that partition generally. For H7/U6, both partitions overestimate the TI increment: the primary bias is +0.03596 kB, interval [0.01501,0.05691]; δ=1.5 gives +0.02755, interval [0.00705,0.04804]. Both fail the controlled-increment test with monovalent background. The broad baseline-to-design entropy changes are also overestimated.

An exploratory small-salt comparison further illustrates the precision limit: H5→H7 has TI entropy change −0.00385 kB (interval [−0.00669,−0.00100]), whereas primary 2PT gives +0.01081 (interval [−0.03736,+0.05898]). The 2PT estimate cannot resolve this small effect; its opposite point-estimate sign is not a significant contradiction given that wide interval. These comparisons are marked exploratory in summary.json.

The controlled primary increments also remain above TI across the 5/10/20-tau smoothing variants. See the window_sensitivity fields in summary.json for paired increments and bias intervals; changing smoothing is not used to select the closest agreement.

## Numerical and label checks

All seven analytic tests pass: the standard fluidicity polynomial; invariance of Δ to length/time units; an Ornstein–Uhlenbeck mode-count/Green–Kubo limit; the analytical Lorentzian entropy tail; the dilute ideal-gas limit; classical h cancellation; and exact synthetic label/mixing cancellation. See analytic_tests.json and test_two_pt.py.

Across all model/window/segment calculations, maximum mode-count error is 8.88e-16, reconstruction error 1.11e-16, and entropy span under h=0.1,1,10 is 7.55e-15 kB. Raw target-T mode weights for primary species range from 2.96861 to 3.01817; they are not all exactly three before equipartition normalization.

**Actual label-only null:** all three U4/L5 pairs use exactly the same sampled particle IDs and have identical hashes of the corrected velocity data for every mechanical group. Their cached VACFs are identical and entropy differences are at floating-point roundoff for every tested model, taper, and time segment. The exactly identical neutral interactions are pooled before entropy estimation; this avoids nonlinear estimator dependence on arbitrary labels. This merger is explicit, not a hidden fitting parameter.

An unmerged L5 calculation using its separately retained label samples gives a finite-sample difference of +0.00062 [-0.00428, +0.00553] kB relative to U4. Its samples and nonlinear spectral transformation differ from the matched pooled null; it is a diagnostic, not an additional independent physical state.
The signed residual has negative spectral weight of 2.29–3.89% of the total three modes across primary species. Clipping it while retaining the full gas contribution would add modes and reintroduce an arbitrary h dependence; over h=0.1 to 10 the spurious entropy span would be 0.317–0.538 kB. This is why the positive-part clipping in the archived code is not retained.

## Sampling and model sensitivity

Full 100-tau records, first/last 50-tau halves, and several taper widths were analyzed. Halves are sensitivity checks, not extra independent replicas. The table gives the span of state-mean primary entropy over taper widths 5, 10 and 20 tau, and the primary mean probe spectral-D discrepancy relative to the full-production MSD coefficient.

| State | Entropy span, tau5/10/20 | Probe spectral D error (%) | Probe fluidicity |
|---|---:|---:|---:|
| B3 | 0.01512 | +0.18 | 0.52973 |
| B5 | 0.01794 | -1.77 | 0.52934 |
| U4 | 0.01473 | -1.58 | 0.55418 |
| L5 | 0.01473 | -1.58 | 0.55418 |
| M5 | 0.01101 | -2.56 | 0.55516 |
| H5 | 0.02181 | +0.90 | 0.56918 |
| U6 | 0.02193 | -1.91 | 0.55242 |
| H7 | 0.01453 | +0.42 | 0.56790 |

Doubling/quadrupling frequency resolution changes state means by at most 0.000012 kB; doubling/quadrupling the sampling interval (a Nyquist/aliasing stress test) changes them by at most 0.002176 kB. The largest state-mean shift for either half record is 0.01199 kB. These numerical checks do not remove the partition-model bias.
Illustrative, volume-conserving ±10% changes in probe or secondary-solvent partial volume shift primary state means by at most 0.00316 kB. These perturbations are sensitivity tests, not measured partial-volume errors or a proof that equal partial volumes are adequate. Converged partial volumes would be needed for a stronger mixture-2PT claim.

The 5% spectral-D practical check fails for none of the state-mean probe estimates. This comparison uses a smoothed 100-tau velocity estimator versus a 1000-tau MSD estimator; agreement does not independently establish the entropy formula. Full replica intervals and all species diagnostics are in summary.json and replicas/*.json.

## What dynamical insight survives?

Probe fluidicity rises from approximately 0.530 in B3 to 0.569 in H5; U4→H5 increases it from about 0.554 to 0.569, and U6→H7 from about 0.552 to 0.568. This is consistent with increased low-frequency motion, but diffusivity already enters the construction, so it is not independent proof of an entropy-driven speedup. Fluidicity is not the fraction of unpaired ions.

The mean probe VACFs remain positive over the first tau in these records. Thus a pronounced early negative backscattering minimum, as previously described for the archived figures, is not supported by this new calculation. A local association or cage-escape mechanism still requires coordination and residence-time measurements.

## Recommendation

Keep TI as the entropy coordinate in the Rosenfeld analysis. Retain the corrected 2PT implementation as a transparent approximation and the spectra as supporting dynamics. Do not restore the old numerical entropy tables, call the species allocations partial molar entropies, or use the improved 2PT–diffusion correlation as an independent confirmation. For quantitative screening, validation must hold across compositions; it currently fails for the monovalent-background increment and for some absolute entropies. Extending velocity records and obtaining partial volumes can test sampling/reference assumptions, but neither should be tuned to force agreement with TI.

Reproduction: run test_two_pt.py, run_rehabilitation.py, then summarize_rehabilitation.py in this directory. The cache contains COM-corrected VACFs, particle IDs, velocity hashes and input provenance. Original MD data, TI results, and archived 2PT outputs were not overwritten.
