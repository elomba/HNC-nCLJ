# Dense transport follow-up (2026-11-01)

**Status:** Running. Smoke jobs 101256–101265 and 101267 (retry of a node-level MPI start failure, 101261) passed on 2026-11-01. Production submitted 2026-11-01: 30 long_N4000 jobs (101298–101327) and 20 N8000 jobs (101348–101367); IDs per replica in `submissions.jsonl`.

The completed component/entropy study left three transport questions open: the approximately 5% H5/U4 and H7/U6 divalent-probe increases vary with the MSD fitting window; the B3 charge-sign difference conflicts with the sign symmetry of its Hamiltonian; and finite-size effects have not been checked.

## Design

Both campaigns use only the dense state, rho=0.7, T=2, shifted Lennard-Jones cutoff 2.5 without a tail correction, and PPPM real-space cutoff 5.0. Each state uses the selected composition from the parent study. Within a replica, the seed, initial coordinates, and reserved probe IDs match across compositions.

| Campaign | N | States | Independent seeds | Equilibration | Production |
|---|---:|---|---|---:|---:|
| long_N4000 | 4000 | B3, U4, H5, U6, H7 | 118543, 230917, 347221, 459803, 562129, 673351 | 100,000 steps | 2,500,000 steps |
| N8000 | 8000 | B3, U4, H5, U6, H7 | 784913, 812477, 925691, 137359 | 100,000 steps | 2,500,000 steps |

The 2,500,000 production steps span 5,000 LJ time units. Positions of all atoms are recorded every 500 steps (1 tau); the last 100,000 steps record ion velocities every 5 steps. The run also records common-COM online MSD, global COM position and velocity, thermodynamic averages, and five stress components. Restart checkpoints occur every 250,000 steps. Each replica ends with `data.final` and `final.restart`.

Inputs and job scripts are under `md/<campaign>/<composition>/replica_<k>/`; smoke replicas are under `smoke/`. Manifests include the composition, seed, step counts and SHA-256 hashes. `jobs.json` and `smoke_jobs.json` index the allocations. Submitted job IDs, if any, go into separate append-only JSONL ledgers.

## Prespecified analysis

The primary endpoints are the paired H5/U4 and H7/U6 log D_probe ratios in the long_N4000 campaign at the 20–100 tau MSD fitting window, their stability across 10–50, 50–200, 100–500 and 200–1000 tau windows, the B3 D(+2)−D(−2) difference, and the dependence on N. MSD uses multiple origins and one full-system COM frame. Independent replicas provide Student-t intervals; overlapping time origins and five trajectory blocks are diagnostics. The block diffusion estimate fits 20–100 tau in each 1,000-tau block.

Additional outputs include configurational Z and u_ex/kBT from the second half of production, five thermodynamic block means, ion VACF Green–Kubo integrals at 2, 5, 10 and 20 tau, shear viscosity from FFT stress autocorrelations and a 5–10 tau plateau, Yeh–Hummer diffusion corrections, and low-k concentration fluctuations. Finite-size ratios use unpaired Welch intervals. A secondary N=4000 contrast pools the original three shorter replicas with six new replicas and is labelled accordingly.

Per-replica analysis writes `replica_result.json` and `msd_multi_origin.csv`. Locked aggregation writes `results/summary.json` and `results/report.md`, including partial completion counts. `summarize_followup.py` writes `results/followup_transport.pdf`.

## Commands

```bash
python3 prepare_followup.py --smoke
python3 prepare_followup.py
python3 submit_followup.py --smoke
python3 submit_followup.py --campaign long_N4000
python3 submit_followup.py --campaign N8000
python3 test_followup_analysis.py
python3 analyze_followup.py
python3 summarize_followup.py
```

`submit_followup.py` only submits when given `--submit`. It checks hashes, shell syntax, Slurm test scheduling and duplicate ledgers first. Smoke has 20 equilibration steps, 2,000 production steps and a final 500-step velocity window, with a 30-minute limit; it is for input validation and is excluded from scientific aggregation.
