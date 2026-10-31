---
name: rebuild-manuscript
description: Regenerate the manuscript tables/figures from study results and rebuild the LaTeX PDFs in doc/. Use after study results, generators or .tex sources change.
---

Rebuild the manuscripts from the repository root (`/sif/users/e.lomba/EntropicFLuids`). Optional argument: a single document name ($ARGUMENTS), one of `paper_pre`, `entropy_theory`, `voronoi_entropy_work_plan`, `entropy_diffusion_prior_work`, `md_pair_entropy_work_plan`.

1. If MD pair-entropy results changed, regenerate their summaries first:
   `python3 MDPairEntropyAssessment/summarize.py && python3 MDPairEntropyAssessment/hnc_validation.py`
2. Regenerate generated tables and figures: `python3 doc/generate_component_study_assets.py`
3. Build: `bash doc/build_pdfs.sh $ARGUMENTS` (no argument builds `paper_pre` and `entropy_theory`).
4. On failure, read the tail of `doc/<document>.build.log` and the `.log`; report the first LaTeX error with its line.
5. Report any `undefined references/citations`, `Rerun` or `Overfull` warnings printed by the script, and `git status --short doc/` so the user sees which generated files changed.

Never fix numbers by hand-editing generated `.tex` tables; fix the generator or source data instead.
