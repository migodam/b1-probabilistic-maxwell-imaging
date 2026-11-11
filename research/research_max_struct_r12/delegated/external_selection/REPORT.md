# Saved Austria pool: exposed posthoc selection diagnostic

No new forward/inversion was run. Actual saved candidate arrays were paired one-to-one with recorded physical Phi values. Data residual is reconstructed only from the exact identity `r²=2(Phi−R)`, using inherited32² regularizer smooth=tv=1e-4,epsilon=.01. Eligibility is the unchanged rawdata relative ball r<=.12. The alternative rule minimizes inheritedR among eligible paid candidates, with first recorded candidate resolving ties. Registered final images remain unchanged.

This is POSTHOC_EXPOSED_SELECTION_DIAGNOSTIC after prior GT scores, not a confirmatory algorithm or parameter search. The Phi/R/ball choices were saved before this script computes material oracles; prior exposure remains explicit. Offline all-pool and eligible-pool oracles only describe available material quality and are never predictions.

| Method | Saved / eligible | Original NRMSE | minR eligible NRMSE | Offline pool oracle NRMSE (residual) | Offline eligible oracle NRMSE |
|---|---:|---:|---:|---:|---:|
|direct|515/0|0.986564|NA|0.907140 (0.618785)|NA|
|Dual|1027/928|1.419644|1.178598|0.728915 (0.711958)|1.178598|
|AP|630/0|0.989414|NA|0.935826 (0.569912)|NA|
|DR|1112/0|2.131982|NA|1.000000 (1.000000)|NA|
|Dual_full|1011/988|1.032829|0.797974|0.613774 (0.490351)|0.797974|

Coverage: every direct objective-call candidate; AP/DR physical audits every50steps plus physical-corrector objective calls; native Dual selected auxiliary and physical-corrector candidates. Unselected native-Dual auxiliary stage arrays were not saved by the caller and cannot be recovered from metadata. No missing stage/inner iterate was fabricated or replayed. Duplicate points are retained as recorded paid evaluations, not independent samples. Raw/final selected arrays are separately preserved; originals were verified present in each actual candidate pool.

A minR result or oracle is limited to its saved finite pool. A pool may contain good material outside the .12 discrepancy ball while containing no equally good eligible candidate. Conversely, a large eligible-pool/original error gap demonstrates that registered Phi selection missed a better material in that particular paid pool, not that minR universally finds it. No generic basin, noise, feasibility, calibration or statistical certificate follows.
