# Saved-result tables and figure contract

`src/plot_results.py` reads persisted arrays and receipts. It does not import the physics backend, evaluate a forward model, run an optimizer, construct a new SVD, or change predictions. Run with the project Python, one BLAS thread; `--only geometry|finite|reconstruction|cost|all` selects the figure group.

## Sources and units

- `results/e1/dimensions.csv`: one row per case, linearization, source subset, material perturbation RMS, and relative data tolerance epsilon. `d_weak` counts normalized singular values satisfying `RMS*sqrt(576)*sigma <= epsilon`; `d_soft` sums the corresponding soft weights. `numerical_null_dimension` uses the separate `1e-12*sigma_max` numerical rank convention. These quantities are not interchangeable.
- `results/e1/coverage.csv`: offline squared GT-error coverage for bases constructed before GT projection. `rank` is retained material dimension; coverage is `||P e||²/||e||²`. The final figure uses saved all-rank curves in `results/coverage_cumulative.npz` and median rank requirements in `results/coverage_dimension.csv`; every basis was constructed independently of GT. The coverage panel uses A=0.4 native anchors, both clean and noisy.
- `results/e1/cases.csv`: real-material leakage geometry, exact numerical kernel threshold, and approximate dimension. The cached current dark dimension is 509; it is distinct from the material pullback kernel and from small Jacobian singular values.
- `results/e3/*.json`: each receipt contains a public-only task, full nonlinear objective top1, trace, attempts, proposal scalars, wall/CPU seconds and actual solver counts. A top1 is selected using the original nonlinear residual plus regularization, never material truth.
- `results/e3/*.npz`: saved anchor/top1 and proposal raw, filtered, before/after material arrays. `proposal_filtered` already includes scalar box-feasible ray shrink, so the figure labels it as filter + shrink. Collapse is explicitly counted at feasible proposal RMS `<1e-12`; it is retained in the denominator.
- `research_z_alg_r4_diversity/truth_after_freeze/<observation_id>.npz:m`: disclosed original GT, matched by observation ID; frozen historical estimates are the saved E3 anchor arrays, previously checked against the historical error archive. Used only offline to display GT, annotate material RMSE, and calculate the oracle candidate bestK. No truth-derived result is an online success criterion.

## Figures

1. `search_space_and_gt_coverage.png`: source-count/tolerance weak dimension, separate numerical null dimension, fixed-reference J versus anchor J versus Green nested coverage, and current-space versus real-material dimension.
2. `finite_steps_and_correction.png`: raw-to-feasible proposal sizes, box-collapse rates, strong/weak Jacobian-component retention of anchor displacement after actual TRF correction, and saved finite-step diagnostic when available. Strong/weak thresholds are those recorded by the solver, not a new rank claim.
3. `reconstruction_and_cost.png`: all 12 A=0.4 clean shapes, GT, original native, four public-objective top1 arms, and Green offline bestK. Display `(m-0.02)/0.4` uses the same fixed `[0,1]` scale in every image; values outside that range are color-clipped only for display. RMSE uses unmodified arrays. Green bestK includes anchor, all before/after proposals, and top1, scored offline against disclosed truth. It is a diversity diagnostic, not a feasible selection rule.
4. `quality_vs_cost.png`: material top1 error versus additional wall and CPU cost and top1 versus offline bestK. Counts include all available E3 tasks, clean and noisy; setup in per-task receipts is included. Historical-anchor generation cost is not charged a second time. CPU seconds and wall seconds are separate.

## Review boundaries

Low residual, approximate current darkness, a large tolerance-weak dimension, broad GT-error coverage, and improved material recovery are separate observations. Near-zero feasible rays can erase intended proposal diversity. BestK can improve even when public-objective top1 does not. Read the sealed score receipts and full cost ledgers alongside figures; plots do not establish GT-free basin access, uniqueness, or a universal reconstruction advantage.

The completed E1 supplemental `current_subset_dimensions.csv` supplies numerical real-material current-leakage kernel dimension for source counts 1/3/6/12, shown as the separate black dashed line in geometry panel B. Its relative rank threshold applies to that leakage map's own largest singular value; it does not assert equality to the data Jacobian rank. The final coverage panel also includes `fixed_reference_Green_weakest` rows from the supplemental saved fixed-reference Green projection. Original E1 cost is explicitly a lower bound after recorded native crashes; use the geometry run/supplemental summaries for provenance, not plot timing.

## Completed E3 figure cohort

The final figures use 96 persisted tasks: 24 per arm, covering the 12 shapes in clean and noisy conditions. For all 96 task anchors checked during plotting, saved anchor arrays were exactly equal to the matched historical estimates (maximum absolute difference 0).

The following medians are calculated from saved arrays and receipts; material RMSE ratios are relative to each task's frozen anchor. Oracle bestK includes the anchor, so its ratio cannot exceed 1.

| Arm | Top1 RMSE ratio | Offline bestK RMSE ratio | Additional wall seconds | CPU seconds | Collapsed feasible proposal rays |
|---|---:|---:|---:|---:|---:|
| Continue TRF | 1.000000 | 1.000000 | 0.0651 | 0.0681 | no proposals |
| Haar | 1.000000 | 0.980234 | 7.0856 | 8.6075 | 12/192 |
| Green | 1.000810 | 0.988811 | 5.0127 | 5.9799 | 12/192 |
| Soft J | 1.005430 | 0.996235 | 4.8578 | 5.5324 | 12/192 |

Ratios here are the median of individual ratios, rather than a ratio of aggregate median errors. These descriptive statistics pool clean and noisy tasks; use the formal sealed score tables for condition-specific NRMSE, IoU, stationarity, and stopping-profile comparisons. The figure generator excludes the separate E2 tight-stopping diagnostic from this E3 cohort.

Formal sealed `results/e3/SCORES.json` was cross-checked after completion: all 96 figure top1 RMSE values agree to below 1e-14, and all four arm median top1/oracle ratios agree with the table above. No physics evaluations were used for this cross-check. The plotting input now uses the released copy of the original truth used by the formal offline scorer, avoiding a dependency on the older error-analysis archive.

## Released-truth and cumulative-coverage polish

The final plotter reads the released historical labels from `research_z_alg_r4_diversity/truth_after_freeze/<observation_id>.npz:m`, matching the report package. It does not require `private_truth` or the old error-analysis archive.

Geometry panel C now reads every retained rank 1–576 from saved `results/coverage_cumulative.npz`; it does not connect the sparse 256/576 endpoints. Its inset table gives median per-case r90 and r95 from `results/coverage_dimension.csv`. This panel's cohort is **24 low-amplitude native anchors** (12 clean + 12 noisy); the report's **48 low-amplitude anchors** additionally include reciprocal TRF and can have different medians. Panel C labels its 24-anchor cohort explicitly. Other geometry panels retain their stated broader cohorts.


## Transverse retention correction

Panel C of `finite_steps_and_correction.png` uses `results/TRANSVERSE_RETENTION.csv`, not correction-step magnitudes. For each recorded strong/weak Jacobian subspace, retention is the norm of the projection of `(after - anchor)` divided by the norm of the projection of `(seed - anchor)`. The plotted bars show medians, and whiskers show 25th/75th percentiles. Rows with seed RMS `<=1e-12`, component denominator `<=1e-12`, missing ratios, or nonfinite ratios are excluded. Each arm/component retains 180 of 192 rays. Weak retention medians: Haar 0.936671, Green 0.324600, Soft J 0.106820; strong medians: 0.097694, 0.115163, 0.069324. These are measured actual-TRF outcomes and are not a theorem about correction.
