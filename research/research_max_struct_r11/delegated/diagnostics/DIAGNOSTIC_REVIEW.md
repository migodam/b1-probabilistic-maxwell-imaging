# Fixed-state diagnostics: observed scope and counterexamples

NRMSE in this local fixed-state report is relative L2 error ||m-GT||/||GT||. Full-run CONTRACTED_DEV_SCORES instead uses the inherited amplitude-NRMSE scoring gate.

Follow-up METRIC_SPLIT_REVIEW.md separates scalar beta contraction from Green feedback. The scalar term dominates all ten states; early cancellation must not be called Green-specific shielding.

Ten existing G19 PGP endpoints (five cases, mf3 and mf12) were used. These are post-readout endpoints. Full R10 re-readout is therefore expected to be close to its box stationarity, and this is not a neutral comparison of full versus lagged trajectories. No inverse replay, new threshold, search or rescue budget was used. All three candidates were retained. GT scoring is OFFLINE_DIAGNOSTIC after candidate generation; it does not choose candidate execution or settings.

At the four early annulus/shell/bridge/confirmation-annulus endpoints, the metric and constitutive gradients are strongly opposed: cosine ranges -0.978 to -0.945 and norm ratio ranges .771 to .958. This is evidence for cancellation in the moving-metric auxiliary gradient. It does not by itself establish that this cancellation causes every bad physical trajectory. The bridge is a direct limit: removing the metric term points against the physical descent direction at both saved states, and its finite lagged steps increase Phi.

Full box projected-gradient norms range 6.35e-8 to 4.97e-5. Lagged norms range 1.71e-5 to .0433. Raw full gradients need not be small because many material coordinates sit at box bounds. Both raw and projected alignments are retained in JSON/CSV. Denominator condition numbers range 6.15 to 38.87; no singularity is needed for the observed cancellation.

|State|metric/const|cos(const,metric)|lagged ΔPhi|lagged ΔNRMSE|ordinary ΔPhi|
|---|---:|---:|---:|---:|---:|
|r10_twod_0001_s0|0.893|-0.956|-0.0382892|-0.0223308|+0.0178758|
|r10_twod_0001_s3|0.237|-0.431|+9.87562e-07|-0.000186763|-2.1261e-05|
|r10_twod_0004_s0|0.771|-0.945|-0.276941|-0.00253172|+0.101706|
|r10_twod_0004_s3|0.681|-0.467|-0.16748|-0.0140054|+0.00125438|
|r10_twod_0022_s0|0.958|-0.978|+0.013124|-0.00836962|+0.0500073|
|r10_twod_0022_s3|0.338|0.007|+0.00309564|-0.00401732|-0.000183994|
|r10_conditional_0000_s0|0.913|-0.949|-0.12342|-0.0323207|+0.0105846|
|r10_conditional_0000_s3|0.148|-0.468|-1.71513e-07|-0.000420758|+2.0486e-06|
|r10_conditional_0003_s0|0.448|-0.832|-0.0193461|+0.0194422|+0.0184223|
|r10_conditional_0003_s3|0.375|-0.043|-0.0027231|-0.000439798|+0.00104038|

Lagged lowers Phi in 7/10 endpoints, raises it in both bridge endpoints and at the successful main-annulus mf12 endpoint. The confirmation coated object at mf3 lowers Phi but worsens NRMSE by .01944. All ten lagged steps lower raw defect norm; that alone does not guarantee a physical Phi decrease or material improvement. All lagged candidates lower their frozen-reference block objective while their moving-metric block objective rises; the objectives must not be conflated.

An initial ordinary-control wiring mistake used has_reg=False and produced the scalar seed rather than the same-prior ordinary readout. It was corrected to has_reg=True on the same ten states, without repeating full/lagged or original gradients. The initial ten seed controls and their cost are archived as SEED_ONLY_ATTEMPT_* plus *_ordinary_seed_unregularized_CANDIDATE.npz. Final totals include them: 50 physical forwards / LU, 10 data adjoints, 720 complex RHS. The final 30 comparison candidates use the requested prior. Full/lagged use the same scalar seed and 8 matrix iterations; ordinary uses its inherited at-most-60 scalar prior iterations.

Detailed evidence: STATE_INVENTORY.json, START_DIAGNOSTICS.json, FIXED_STATE_COUNTERFACTUALS.json, DIAGNOSTIC_SUMMARY.csv, DIAGNOSTIC_RECEIPT.json and ORDINARY_REPAIR_RECEIPT.json. Candidate arrays preserve every proposed material step. No same-state observed historical next material array exists.
