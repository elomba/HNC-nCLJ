# 2PT rehabilitation of the completed component-entropy study

**Completed for all 24 existing velocity trajectories.** See [report.md](report.md) for results, [METHODS.md](METHODS.md) for equations, [summary.json](summary.json) for replica statistics and acceptance checks, and the comparison CSV/PDF files.

The corrected classical calculation passes analytic, normalization, reconstruction, reference-invariance, and mechanical-label tests. It does **not** pass the study-wide quantitative TI benchmark. Its composition-dependent residual biases prevent replacing TI in the Rosenfeld analysis. The numerical agreement for one diversification increment does not extend to the monovalent-background pair.

This directory is independent of the archived scripts in `../../lammps_study/two_pt_analysis/`. Original trajectories, TI results and legacy outputs remain unchanged. Cached COM-corrected VACFs make repeated model/sensitivity evaluation inexpensive. The cache is ignored by git and can be regenerated from the retained NetCDF trajectories. Sources are cited in `sources/README.md`; downloaded literature files are ignored.

Run from the repository root:

```bash
python3 ComponentEntropyStudy_20261030/two_pt_rehabilitation/test_two_pt.py > ComponentEntropyStudy_20261030/two_pt_rehabilitation/analytic_tests.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ComponentEntropyStudy_20261030/two_pt_rehabilitation/run_rehabilitation.py
python3 ComponentEntropyStudy_20261030/two_pt_rehabilitation/summarize_rehabilitation.py
```

The trajectory extraction needs NumPy, SciPy and netCDF4; plotting also needs Matplotlib. All raw-data analysis can run on a CPU. No new MD or GPU allocation is necessary. `--only <composition-name>` and `--replica <1|2|3>` select a subset; `--force-cache` recomputes cached VACFs. Change the cache version or force regeneration if modifying the particle selection or autocorrelation estimator.

Uncertainty intervals use three independent replicas. Time-record halves and overlapping correlation origins are not additional replicas. The primary model and study-specific tolerances were recorded before the first baseline evaluation, and no parameter was fitted to TI.
