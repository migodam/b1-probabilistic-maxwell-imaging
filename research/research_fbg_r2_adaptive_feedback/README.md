# Research FBG-R2

Start with [中文研究报告](RESEARCH_REPORT_ZH.md), [mathematics](MATHEMATICAL_THEORY.md), [algorithm](ALGORITHM.md), and [reproduction](REPRODUCE.md).

Actual frozen studies: 639 completed inference tasks, 14 fresh independent materials, 22 fresh observations, six primary methods. Fixed scalar feedback passes 8/12 clean high-contrast objects; original GlobalDual passes 3/12. This is a finite scalar-TM benchmark result, with increased average work and retained failures. Automatic qbeta passes6/12 and does not beat fixed beta.

- [All fresh digits](figures/fresh_digits_1to9.png), [all14objects](figures/fresh_all_shapes.png), [noise](figures/fresh_noise.png).
- [Results](results/RESULTS_SUMMARY.json), [case table](results/RESULTS_BY_CASE.csv), [cost ledger](audit/COST_LEDGER.csv), [independent review](audit/FINAL_SCIENTIFIC_REVIEW.md).
- `truth_after_freeze/` is the seedless disclosure AFTER prediction closure; it is an offline scorer input and is never read by `src/adaptive_solver.py`.
- Source archives and per-stage manifests distinguish frozen inference from later analysis. Some unrun extensions are explicitly NOT_RUN.

One failed Stage3qbeta run consumed60024RHS against60000; this deviation is retained, not corrected by rerunning. Nominal wall splits have measured overruns. Same-grid synthetic simulations do not establish continuous-PDE or3Dvector performance, unique GT, or posterior/root coverage.
