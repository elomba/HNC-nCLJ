# MD all-pair excess entropy analysis

This directory contains the independent MD RDF and pair-entropy assessment for the eight selected mixtures in `../ComponentEntropyStudy_20261030/`.

The production calculation processes every saved particle position in the 24 existing dense-state trajectories. It forms finite-count normalized, species-resolved RDFs after merging labels with exactly identical physical interactions, then evaluates the Green--Wallace pair entropy. Twenty-four matched ideal-point records check finite-sample bias. The primary analysis uses the last 1000 trajectory frames; an additional fine-core extraction resolves the original short-range bin-width sensitivity. No entropy or diffusion result was fitted to the other.

## Main files

- `protocol.json`: estimator, windows, sensitivity checks and core-refinement decision.
- `extract.py`, `histogram.c`: position reader and compiled all-pair histogram kernel.
- `analyze.py`, `refine_core.py`, `summarize.py`: replica analysis, fine-core replacement and summary tables/plots.
- `hnc_validation.py`: independently integrates exported HNC RDFs and plots HNC/MD probe-solvent RDFs.
- `test_pair_entropy.py`, `validation.json`: independent checks for geometry, normalization, pair factors, label merging and numerical limits.
- `results/report.md`: scientific results, interpretation and follow-up recommendations.
- `results/summary.json`, `states.csv`, `contrasts.csv`: numerical results with replica intervals.
- `results/hnc_validation.json`, `results/hnc_md_probe_rdf.pdf`: HNC export validation and structure comparison.

The trajectories and raw histograms are not copied here. The source path and checksum for every trajectory are recorded in the cache metadata. Extraction and refinement Slurm arrays 101032 and 101057 each completed 24/24 tasks with exit code zero.

To regenerate the summary on the configured Python environment, run:

```bash
python MDPairEntropyAssessment/summarize.py
python MDPairEntropyAssessment/hnc_validation.py
python doc/generate_component_study_assets.py
```

The last command updates the pair-entropy tables and figure included by `doc/paper_pre.tex`. See `results/report.md` for estimator checks, uncertainty limits and scientific conclusions.
