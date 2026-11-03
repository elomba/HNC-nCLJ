# High-divalent-concentration HNC prior study

**Status:** Complete. HNC screen, smoke jobs 101598–101605 and the 32 production allocations (jobs 101642–101673, `submissions.jsonl`, all COMPLETED with exit 0 on 2026-11-03) are finished and aggregated by `summarize_md.py` into `results/summary.json`, `results/report.md` and `results/md_contrasts.pdf`/`.png`. See *MD results* below.

This campaign repeats the 35-design HNC screen from `ComponentEntropyStudy_20261030` at
T*=2, rho*=0.7 with x(+2)=x(-2)=0.075. At N=4000 this corresponds to 300 ions
of each probe sign. The solvent fraction is reduced by 0.09 relative to each
old design; every weak-solvent and background-salt fraction and LJ epsilon is
preserved. All diameters are 1, all cutoffs are 2.5, and the electrostatic
coupling gamma is 1, so the Bjerrum length at a given temperature is unchanged.
The input records are `candidates.json`; the old additional designs are read
from the completed campaign and rescaled programmatically.

U4, L5, M5, and H5 match the mean bare LJ amplitude
sum_i x_i sqrt(epsilon_i); U6 and H7 form another matched pair. L5 splits the
U4 weak solvent into two mechanically identical labels and is an exact
thermodynamic null control. The table below lists all fractions explicitly;
unlisted species have zero fraction. Ionic strength is I = 1/2 sum_i x_i q_i^2.

| Label | Candidate | solvent | weak or A/B | background +1/-1 | probes +2/-2 | I | HNC s_ex |
|---|---|---:|---|---|---|---:|---:|
| B3 | baseline_3C | 0.850 | — | — | 0.075/0.075 | 0.300 | -2.0796 |
| B5 | mono_f010 | 0.750 | — | 0.050/0.050 | 0.075/0.075 | 0.350 | -2.0840 |
| U4 | weak_f030_e025 | 0.550 | weak_solvent=0.300 | — | 0.075/0.075 | 0.300 | -1.8598 |
| L5 | dispersion_split_d000 | 0.550 | weak_A=0.150, weak_B=0.150 | — | 0.075/0.075 | 0.300 | -1.8598 |
| M5 | dispersion_split_d025 | 0.550 | weak_A=0.150, weak_B=0.150 | — | 0.075/0.075 | 0.300 | -1.8266 |
| H5 | dispersion_split_d040 | 0.550 | weak_A=0.150, weak_B=0.150 | — | 0.075/0.075 | 0.300 | -1.7490 |
| U6 | mono_weak_uniform6C | 0.450 | weak_solvent=0.300 | 0.050/0.050 | 0.075/0.075 | 0.350 | -1.8643 |
| H7 | mono_weak_split7C | 0.450 | weak_A=0.150, weak_B=0.150 | 0.050/0.050 | 0.075/0.075 | 0.350 | -1.7534 |

## HNC outputs and interpretation

`hnc/screen/<name>/` contains every coarse attempt (2048 points, dr=0.01,
Newton dimension 40). `hnc/fine/<name>/` contains refined selected analogues
and any additional candidate beating their best coarse s_ex (4096 points,
dr=0.005, dimension 100, restarted from the coarse solution). Each solve has
its own input, solver log, thermodynamics, structure factors and `result.json`.
The summary JSON files at this directory's root collect the results, while
`results/hnc_report.md`, `results/hnc_summary.json`, and
`results/hnc_screening.pdf`/`.png` compare both concentrations and document
the TI-path and temperature-margin scans. The first 20 reciprocal-grid points
are inspected for eigenvalues of the concentration-weighted S_ij matrix and
for S_cc. A low eigenvalue is an alert, not a phase-boundary determination.

The 19-density scan follows the proposed MD thermodynamic-integration path at
T*=2. The temperature scan at rho*=0.7 continues from beta=0.5 to 1.0. It only
probes the HNC closure: failures at lower T* reflect the breakdown of HNC under
strong divalent Coulomb coupling, not a phase boundary. MD stays at T*=2.
All 152 TI-path solves converged with positive inverse compressibility and
positive sampled low-k eigenvalues. At rho*=0.7, B3, B5, U4, L5, and U6
converged through beta=1.0; M5 last converged at beta=0.9, and H5 and H7 at
beta=0.8. The later attempts failed in the HNC solver. The solver's `chi`
column is reduced inverse isothermal compressibility.

## MD campaign

`selected.json` is generated from `candidates.json` by `prepare_md.py` and
contains B3, B5, U4, L5, M5, H5, U6, and H7. Integer counts sum to N=4000,
including 300 +2 and 300 -2 ions, and total charge is zero. Four matched
replicas use seeds 41057, 52691, 67433, and 83819. The same seed gives the
same initial positions, probe IDs and velocity sample IDs across compositions.
`md/` contains 32 prepared allocations; `smoke/` contains eight. Both include
checksummed manifests and `jobs.json` or `smoke_jobs.json`. All 32 production
jobs have completed.

Each allocation follows 18 EOS densities from 0.001 to 0.65 (50,000
NVT equilibration plus 100,000 production steps), then rho=0.7 (100,000
plus 1,500,000 steps = 3000 tau production). NVT uses T*=2, dt=0.002,
Tdamp=0.1, shifted LJ cutoff 2.5 sigma at every density with no tail
correction, and PPPM tolerance 1e-5. The Coulomb real-space split is
5(0.7/rho)^(1/3) sigma; this changes the Ewald split, not the LJ Hamiltonian.
One Blackwell GPU, two MPI tasks, 16 GB, gpu_v4_short on ladon28/29 and a
12-hour allocation are requested. Dense restart files are written every
250,000 steps. **Automatic dense-stage resume is unsupported.** A rerun with
an existing checkpoint or final.restart exits before launching LAMMPS; move
the partial directory aside and regenerate a fresh allocation if needed.

The dense stage writes `thermo_avg.dat` (the 0.7 EOS point), `msd_common.dat`
(one global center-of-mass frame for every species), `positions.nc` (unwrapped,
every 1 tau), `stress_offdiag.dat` (five independent shear components every
0.01 tau), `global_vcom.dat`, and `charge_current.dat` (three components every
0.01 tau). During only the final 100 tau, `velocities.nc` records all ions
and up to 200 neutral particles per species every five steps. The per-replica
analysis writes `replica_result.json`, `entropy_density.csv`,
`msd_multi_origin.csv`, and `conductivity_running.csv`; the campaign analysis
writes `results/summary.json` and `results/report.md` and reports missing
replicas explicitly.

Conductivity uses the total charge current J=sum_i q_i v_i, with total
charge asserted zero, and sigma(t) = [3 V k_B T]^-1 integral_0^t
<J(0)·J(t')> dt'. LAMMPS `units lj` has Coulomb prefactor 1 here; q is
its dimensionless charge, k_B=1, and sigma is in reduced units
q²/(epsilon tau sigma). The Green–Kubo estimate averages the running
integral from 5 to 10 tau and also reports 2, 5, 10, and 20 tau cutoffs.
The Nernst–Einstein estimate is rho/T sum_i x_i q_i² D_i, and the reported
inverse Haven ratio is sigma/sigma_NE. Viscosity uses the same 5–10 tau
plateau window. Diffusion uses multi-origin MSD fits at 10–50, 20–100,
50–200, 100–500, and 200–1000 tau; Yeh–Hummer corrected values are included.
The TI entropy uses the configurational virial and Debye–Hückel dilute
endpoint. Paired transport contrasts use 100[exp(mean ln ratio)-1] with
Student-t confidence intervals across matched replicas.

### Validation and execution

```bash
python3 ComponentEntropyStudy_x015_20261102/prepare_md.py --smoke  # new empty smoke/ only
python3 ComponentEntropyStudy_x015_20261102/prepare_md.py          # new empty md/ only
python3 ComponentEntropyStudy_x015_20261102/test_md_analysis.py
~/.local/bin/ruff check ComponentEntropyStudy_x015_20261102/*.py
python3 ComponentEntropyStudy_x015_20261102/submit.py
python3 ComponentEntropyStudy_x015_20261102/submit.py --smoke
python3 ComponentEntropyStudy_x015_20261102/submit.py --smoke --slurm-test
python3 ComponentEntropyStudy_x015_20261102/submit.py --smoke --submit
python3 ComponentEntropyStudy_x015_20261102/check_smoke.py
python3 ComponentEntropyStudy_x015_20261102/summarize_md.py
```

`submit.py --submit` submits production and records each job ID immediately
in `submissions.jsonl` (`--smoke --submit` uses `smoke_submissions.jsonl`);
the default command only validates hashes, stage counts and shell syntax.
The Codex sandbox cannot reach the Slurm socket, so jobs were submitted from
the interactive session. Production allocations took 5.96–7.73 h each
(B5 longest), within the 12 h limit; `md/` occupies about 21 GB.

## MD results (32/32 replicas)

Paired 95% Student-t intervals over four matched replicas (`results/report.md`).

- **Entropy:** H5−U4 +0.1027 [0.1014, 0.1039] k_B, H7−U6 +0.1041, H5−B3 +0.3014,
  H7−B5 +0.3015 (x=0.03: 0.1032, 0.1033, 0.3021, 0.3015). HNC overestimates
  the H5 increments by 8–10% (0.1108, 0.3306). L5−U4 +0.0002 [−0.0012, 0.0015]
  with identical transport (exact null). B5−B3 −0.0034.
- **Probe diffusion:** H5/U4 +5.18% [3.65, 6.73], H7/U6 +6.01% [4.63, 7.41]
  (x=0.03: +5.06%, +5.90%); M5/U4 +1.70% [0.94, 2.47]; H5/B3 +19.1%,
  H7/B5 +19.5%; B5/B3 +0.12% (unresolved). Mean diffusion +8.7% and +9.1%.
- **Viscosity:** H5/U4 −9.2% [−13.3, −4.9]; H7/U6 −5.7% (unresolved).
- **Absolute shift vs x=0.03:** s_ex −0.016 k_B in every state (as HNC
  predicts), Z −0.020, u −0.078 to −0.082, D_probe −0.5 to −1.3%, η unchanged
  within uncertainty.
- **Conductivity:** the Green–Kubo running integral plateaus by 2–5 τ (2, 5,
  10, 20 τ cutoffs agree within CIs). σ/σ_NE = 0.69–0.74 per state, pooled
  0.727 [0.708, 0.746] (28 replicas, L5 excluded); no paired change in
  σ/σ_NE is resolved. B5/B3 σ +18.8% [6.6, 32.4] (NE +17.7%). Controlled
  diversification: NE predicts +5.2% (H5/U4) and +6.1% (H7/U6); measured
  σ −1.3% [−10.4, 8.6] and +4.0% [−9.7, 19.8], consistent but unresolved.

- **Rosenfeld-like fits** (`rosenfeld_fits`, `local_slopes`): B_net 0.80 (B3/U4/M5/H5)
  and 0.82 (B5/U6/H7), B_probe 0.59 in both (x=0.03: 0.62, 0.55); B_σ 0.43 and
  0.53 (noisy). Local slopes H5/U4 0.81/0.49, H7/U6 0.84/0.56 (mean/probe).
- **Sign symmetry and homogeneity:** D_+2 − D_-2 within −0.0009 to +0.0025, all
  CIs include zero; low-k concentration fluctuations 0.53–1.57.

Conclusion: every x=0.03 conclusion holds at x(±2)=0.075; ionic conductivity
follows carrier number and ion self-diffusion with a composition-independent
inverse Haven ratio of about 0.73.

## Reproduce

From the repository root:

```bash
make  # only if bin/hncr_nclj is missing
python3 ComponentEntropyStudy_x015_20261102/screen_hnc.py
python3 ComponentEntropyStudy_x015_20261102/summarize_hnc.py
~/.local/bin/ruff check ComponentEntropyStudy_x015_20261102/screen_hnc.py ComponentEntropyStudy_x015_20261102/summarize_hnc.py
```

The screen command reuses existing `result.json` files and regenerates the
same summaries. The `--x-probe` option changes the composition for an
independent run; use a separate copy of this campaign directory for that
parameter sweep to avoid mixing cached results.
