# Saved-pool selection diagnostic: limited observed result

This is POSTHOC_EXPOSED_SELECTION_DIAGNOSTIC after all earlier Austria GT scores. There were zero new forward solves or inversions, no replayed stages, and no radius/box/prior/rank changes. All original endpoints remain unchanged. Selection uses the existing relative rawdata ball delta=.12 and smallest exact inheritedR among eligible actual saved paid material/Phi pairs. Oracle calculations are explicitly offline descriptions of the same finite pools.

| Pool | Saved pairs / eligible | Registered-final NRMSE | MinimumR eligible NRMSE | Offline all-pool oracle NRMSE / residual | Eligible-pool oracle NRMSE |
|---|---:|---:|---:|---:|---:|
| direct |515/0|.986564|unavailable|.907140 / .618785|unavailable|
| Dual v1 |1027/928|1.419644|1.178598|.728915 / .711958|1.178598|
| AP |630/0|.989414|unavailable|.935826 / .569912|unavailable|
| DR |1112/0|2.131982|unavailable|1.000000 /1.000000|unavailable|
| Dual full |1011/988|1.032829|.797974|.613774 / .490351|.797974|

Within the saved Dual pools, registeredPhi selection missed better eligible material; the fixed minimumR rule selects the best eligible material in these two particular pools. Even those choices are worse than the recognizable selected auxiliary shapes, which remain outside the fixed.12 ball. Thus selection explains part of the degradation, while this pool does not contain equally good material satisfying the current ball. No general superiority follows from the coincidence between minimumR and the two eligible-pool oracles.

Direct/AP/DR have no eligible saved candidate. DR's final residual.122674 remains outside.12; it was not rounded into feasibility or granted a larger radius. Its all-pool material oracle is the initial zero-contrast background with NRMSE1, so the saved DR material pool does not hide a better feasible image lost solely by final selection. This is limited to actual recorded candidates and finite budgets, not a theorem about every iterate or possible DR run.

Coverage is real and incomplete: direct objective-call candidates; AP/DR physical audits every50projectionsteps and physical-corrector objective calls; native Dual selected auxiliary and corrector objective-call candidates. The caller did not persist the other Austria native stage arrays; no missing stage or unaudited inner material is reconstructed. Repeated paid evaluations are retained, not independent samples. Exact counts and labels appear in RESULTS.json and SELECTION_SEAL.json.

Mathematics: all stored Phi values use `Phi=.5*||Ypred−Y||²/||Y||²+R(m)`, so recorded raw relative residual is recovered by `r=sqrt(2*(Phi−R))` without a new forward. R is the unchanged32×32 mean-neighbor quadratic/smoothedTV regularizer, smooth=tv=1e-4 andepsilon=.01. Eligibility uses r<=.12 exactly; no feasibility tolerance was added. MinR ties use the first recorded candidate. Public Phi/R/ball choices were written before oracle calculations in this script, but earlier GT exposure is acknowledged and prevents a confirmatory interpretation.

The.12 threshold is a predetermined discrepancy diagnostic on noisy observations, not a known iid-CN confidence ball. Eligibility here is only a saved-data-residual condition. There is no generic feasibility, material-correctness, basin, calibration or statistical certificate, no confidence interval and no out-of-sample claim. The diagnostic used0.186379CPU seconds (RESULTS.json authoritative exact timing). Key artifacts: CONFIG.json, SELECTION_SEAL.json, SELECTED_MATERIALS.npz, RESULTS.json, per-pool metric arrays, OFFLINE_ORACLE_MATERIALS.npz, REPORT.md and diagnose_saved_pool.py.

Later separately authorized physics diagnostic: actual original observations are150²→32², and the true coarse32² material has residual.2252565 outside the fixed.12 ball. See external_physics/REPORT.md. Thus eligible-pool selection here excludes the true coarse material under that model/data mismatch; it must not be interpreted as a pure optimizer-only failure or a valid noise-only feasibility target. Results above remain unchanged. Matched32 control is separately reported in external_matched32/REPORT.md.
