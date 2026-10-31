# Manuscript data and rebuild

Both `paper_pre.tex` and `entropy_theory.tex` now use the completed eight-mixture, three-replica study in `../ComponentEntropyStudy_20261030/` as their primary transport evidence.

- Exact compositions: `selected.json`.
- State means and paired replica intervals: `results/summary.json`.
- MSD-window, VACF, symmetry, and descriptive-fit diagnostics: `results/diagnostics.json`.
- Species-resolved means: the 24 `md/*/replica_*/replica_result.json` files.
- Scientific interpretation and limitations: `results/conclusions.md`.

`component_study_methods.tex`, `component_study_results.tex`, and `component_study_conclusions.tex` are shared by both manuscripts. Six numerical table files and three `figures/component_*.pdf` figures are generated from the sources above; edit the generator rather than typing replacement numbers into those tables.

```bash
python3 doc/generate_component_study_assets.py
bash doc/build_pdfs.sh
```

The build script accepts an optional document name (`paper_pre`, `entropy_theory`, `voronoi_entropy_work_plan`, or `entropy_diffusion_prior_work`) to rebuild one PDF; without arguments it builds the two manuscripts. The build script loads `texlive/20260310-GCC-13.2.0`. Run it on a configured cluster node if the local environment is slow or lacks the module.

`voronoi_entropy_work_plan.tex` is a standalone handover plan for local Voronoi topology, chemical-neighborhood controls, comparison with TI and diffusion, and conditional cluster-entropy extrapolation. It includes the related electrolyte and graph-entropy references. Build it with `bash doc/build_pdfs.sh voronoi_entropy_work_plan`.

`entropy_diffusion_prior_work.tex` is the separate annotated literature review covering entropy scaling, ionic mixtures, high-entropy electrolytes, Voronoi/graph entropy, cluster extrapolation, and 2PT. It distinguishes published findings from our assessment and summarizes the contribution and limits of the current calculations. Build it with `bash doc/build_pdfs.sh entropy_diffusion_prior_work`.

The new study replaces the previous single-trajectory diffusion, incompatible ion-transport tables, unverified species-entropy/2PT comparisons, and old scaling fits in the research manuscript. The LJ-cutoff audit remains in its appendix. The theory report retains explicitly historical HNC benchmark tables for different compositions, corrects the interpretation of their dilute endpoint, and replaces the former two-state transport claim with the new study. Historical benchmark values are not mixed with the new replica data.

The approximately 5% controlled probe response is presented as provisional, including its fitting-window dependence and the baseline charge-sign symmetry flag. The stronger 18–20% baseline-to-design increases are distinguished from this controlled comparison.


## Subsequent 2PT validation

The completed rehabilitation tests are documented in [the 2PT report](../ComponentEntropyStudy_20261030/two_pt_rehabilitation/report.md). The corrected implementation passes internal/reference and label tests, but fails the study-wide quantitative TI benchmark. Both manuscripts now include `component_study_2pt.tex`, a generated controlled-increment table, the corrected entropy-validation figure, and a concluding assessment. The old 2PT tables remain excluded; the replicated MD/TI data are unchanged. The new assets use `two_pt_rehabilitation/summary.json` and `entropy_validation.pdf`.
