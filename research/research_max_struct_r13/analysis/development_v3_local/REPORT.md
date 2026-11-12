# R13 post-seal offline report

Campaign: `development_v3_local`. Expected tasks: 6; task issues/failures: 0.
Parent-authorized seal: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r13_worktree/research/research_max_struct_r13/audit/development_v3_local_seal.json`. Truth opened only after authorization and byte/mtime validation.

## Estimators and denominators

`primary` is frozen paid minimum-Phi selection. `secondary` uses a paid noise-feasible minimum-R candidate and falls back to primary; `secondary_only` exposes actual availability. `feasible_raw_retain` preserves an already noise-feasible raw candidate and otherwise uses primary. Raw/corrected stages are retained independently. Missing tasks are never removed from expected denominators; valid paid evidence from failed runs is labeled partial.

## Material metric definitions

RMSE uses all material cells; amplitude NRMSE divides by saved inherited amplitude for development cases, or max(cell-average truth)-host for new synthetic data. Boxrange NRMSE divides RMSE by the public allowed material box width (upper-lower); relative L2 is reported separately. Development C uses prediction Otsu256 versus saved support; other synthetic rows use fixed chi>1 support. Mean amplitude bias is mean signed material error on truth support, alongside its amplitude-normalized value. Continuous-field support IoU is descriptive. No new joint success threshold is introduced. Fresnel without trustworthy pixel GT has null material accuracy metrics.

## Paid time and costs

Prefix candidates must be legal, eligible, already audited by the requested cumulative process CPU and within total budget. Cost curves use paid physical min-Phi incumbents, never GT-selected roots. Whole external child CPU, internal audit CPU, physical LU/RHS and projector work remain separate. Missing cost counters remain null. Explicit queue partial/interrupted/failed attempt history is retained separately in QUEUE_ATTEMPT_HISTORY; unknown costs are not imputed and possible overlapping history is not added again. Complete cost records are in COST_LEDGER/COST_TOTALS. Scalar aggregate RHS and vector disjoint source/adjoint RHS use documented physical_cost_source_keys; nested deltas are never added again.

## Physics domains and statistical scope

METHOD_VIEW_SUMMARY, PAIRED_STATS and figure directories are separated by task physics_kind/cohort domain. No single pooled material accuracy mean combines classic, measured Fresnel and synthetic 3D/2D experiments. Noise/discrepancy eligibility retains the runner qualification; Fresnel/classic diagnostics are not CN noise certificates.

Paired effects resample material objects, averaging repeated conditions within an object. Fixed bootstrap seed131313,4000replicates; confidence intervals and Wilcoxon p-values are descriptive and unadjusted. Available paired counts are displayed against full case denominators. Development objects are exposed; cross-grid mismatch and prior selection can affect material error independently of data-ball feasibility.

## Artifacts and review

JSON/CSV retain all estimator, candidate, prefix, cost, failure and source rows. Paginated figures include full case IDs and shared color scale, paid CPU/LU curves, corrector objective decomposition and paired material errors. Three-dimensional pictures are central orthogonal slices, not complete volume validation. Optional --structure-gallery adds an explicitly diagnostic row-shared GT-derived physical scale, a colorbar and counts of saturated prediction cells. It never changes selection or the absolute gallery. Fake verification fixtures under delegated/evidence/generation/scorer_fixture are excluded from these campaign output roots and must never enter publication. Parent visual/scientific review and publication remain outstanding.

## Failures

No missing or failed task detected by this scorer.
