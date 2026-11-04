# R7 targeted secondary confirmation

Frozen before generation: 2026-11-04T09:14:01.519630+00:00

This is a secondary confirmation of two selected family/acquisition strata observed in the completed primary30 study. It cannot support a universal claim and will never be pooled with primary30. Twenty fresh independent geometries (10 curved_bar at contrast12/six full sources; 10 slit_shell at contrast18/12 half-arc sources), all complex relative RMS noise .001, seed2026110419. All20 remain in the analysis, including failures. Geometry draw distributions and order are recorded in configs/targeted_frozen.json and match generate_benchmark.geometry.

Use existing public 24x24 operators and Model, material=.1+amplitude*alpha, the original circular-complex Gaussian noise convention. Public observations contain no GT; truth is written only under sealed_truth/targeted. Online solvers never read truth. Each full-box material has576 independent real cells initialized at .1. Exactly three frozen arms: a_plus from configs/aplus_final.json; cie_dual and global_dual from configs/primary.json, with the existing full3000 auxiliary plus3000 corrector allowances. No algorithm, source, budget, acquisition, or initialization modifications.

Finite queue:20×3=60 tasks, four single-BLAS CPU workers. Anticipated3000CPU seconds; all actual generation/setup/solve costs, errors and failed attempts will be preserved. No scientific wall cap, search, replacement, or sample dropping.

All60 predictions must be sealed before any new-truth scoring. Scoring additionally waits for parent TARGETED_SEALED authorization, uses existing score_primary.material_metrics/noise/fixed-prior conventions, and makes zero new forward calls. Primary success is NRMSE<=.15 and IoU>=.7. CN99 consistency and legacy diagnostics are separate. Report paired Aplus-versus-both-baseline means and confidence intervals by object, CPU/LU/RHS checkpoints, family tables, all20×3 final images.

## Frozen baseline amendment 2026-11-04T09:17:10.901927+00:00

Before any targeted truth scoring or targeted recovery ranking, parent added cie_single from unchanged configs/primary.json on all20 existing materials. Rationale comes only from disclosed primary30: original slit_shell obj20 also recovered with cie_single. This strengthens the simple baseline comparison without changing any target, method parameter, original60 task, or selection. The extra20 tasks have the same full frozen allowance and use a separate one-worker queue. All80 results must be sealed before TARGETED_SEALED authorization. Final manifest: configs/targeted_final_manifest.json. Original manifests and results retained.
