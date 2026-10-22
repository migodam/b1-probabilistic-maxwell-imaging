# R2 execution/probability review

Status: no unresolved fatal density or leakage error found in inspected path. One execution blocker was found and repaired externally during this review.

Review invoked 0 physics calls, 0 training updates, read no locked ground truth, and modified no source or prior results.

## E01 — BLOCKER_RESOLVED_DURING_REVIEW: Missing json import in probit anchor logger

Initial anchors.py lines1-4 lacked json; log line25 used json.dumps. Final reread finds import json at line2 inserted by another actor. Reviewer made no source changes.

Initial source could not generate finite probit anchor; no remaining missing-import blocker at final read. No physics rerun performed.

## E02 — P2: Training bypasses exclusive worker memory policy

core.py70-74 only locks non-analysis-only; tools/run_p3r2_accounted.py9 constructs analysis_only=True for arbitrary training subprocess. Live process inventory showed material training, density training and C physics concurrently.

One physics worker is enforced, but comment that physics/training retain exclusive lock and protocol paid_workers=1 do not hold for training. Controlled wall remains reserved per subprocess, so concurrency alone is not missing wall billing.

## E03 — P2: MPS path has no GPU receipt

Training CLIs accept --device mps; wrapper uses analysis_only Book, which forbids gpu deltas. Training functions report CPU and wall seconds only.

Current live commands use CPU, so no observed GPU underbilling. Future MPS use would bypass authorized GPU accounting unless separately recorded.

## E04 — CLAIM_LIMIT: Anchor ablation uses weaker caps than declared primary baseline

run_p3r2_anchors.py11 defaults48/48; run_p3r2_physics_sequence.py24-25 explicitly requests48/48. PROTOCOL primary_baseline is native96 plus incremental probit64.

48/48 anchor ablation is valid if labeled. It cannot alone establish gain against the frozen strong primary baseline.

## E05 — CLAIM_LIMIT: ROI risk must be converted to frozen full-grid metric

risk_rows computes256-cell complex ROI MSE; protocol primary_metric is576-cell full-grid complex MSE; predict_full fills known constant host outside ROI.

For legal truth and host-filled predictions full-grid MSE=(256/576)*ROI MSE. Relative reductions and paired-sign inference unchanged, absolute MSE/SEM must be labeled or converted.

## E06 — HARDENING_LIMIT: GT isolation is interface discipline, not a process access boundary

Training path guards only reject string locked; pack CLI restricts groups train/validation. Online Target.observation validates exact keys. Same filesystem process can read offline_only; analysis_only wrapper asserts child physics-free but does not enforce subprocess file/import behavior.

No actual locked access or online GT leakage found in inspected paths. Do not claim OS-enforced inaccessible truth or use filename checks as independent evidence of provenance.

## Source conclusions

- Core native inversion imports .core and uses supplied Target; runner imports p3r2.inversion, not p3l.inversion. No wrong duplicated-core Target routing found.
- Noise generation is sigma/sqrt(2) complex iid and Phi=||F-y||²/sigma²; white=sqrt(2)/sigma packing gives consistent half squared residual likelihood.
- Material supervised loss uses one global target scale and factor2 only between real-component and complex MSE; train normalization only; update0 is eligible; validation used for selection so not fresh evaluation.
- Residual inference explicitly extracts observation/residual/q0/chi0/J only; train-J geometry frame excludes GT and uses ascending/descending singular directions. Decoder outputs are point estimators, not densities.
- Conditional diag/full/MDN/flow log_prob formulas use correct Gaussian constants, mixture normalization and inverse-flow determinant sign. Full covariance diagonal is positive. Frame completion is orthonormal.
- Full66 q density is learned10(z8,eta2|obs)*prior56(z_rest)*|deta/du| wrt dz64 du2; positive .5*sum(eta²+log2pi) in logdensity is correct. Complement prior is a normalized proposal assumption, not posterior independence evidence.
- Sampling uses same fixed frame and probit map; floating-point u boundaries throw rather than silently introduce atoms. Prior unit model cancels eta density with Jacobian, recovering standard Gaussian64 times uniform2.
- Forward reservation precedes adapter physics; adapter full role becomes full_forward and consumes upfront8 credit. Receiver transfer reserves32 adjoint RHS. Ledger has exclusive update lock, stage/global caps, controlled-wall reservations, and unreconciled jobs fail closed.

## Actual-array/provenance check

Only independent parent vectors were read from actual new train/validation files: 16384/2048 unique, overlap0. R2 target config equals legacy default config at review time. No material/observation/q arrays were loaded by this check.

## Limits

Normalized proposal does not mean validated posterior. A selected validation improvement is descriptive until the frozen fresh locked evaluation. Remaining source uncertainty is recorded in REVIEW.json.

## Parent clarification received before review close

- The missing import was fixed before the first A-anchor state; parent reports no effect on earlier native calls. Reviewer verified repaired source and performed no job rerun.
- One physics plus two bounded CPU trainings is an accepted scheduling choice. Parent reports combined RSS about1.2GB and will formalize the2.5GB aggregate overlay. Existing single-worker wording needs to match this policy; no accounting loss is inferred from concurrency itself.
- Current runs are CPU-only. Rejecting future unaccounted MPS in the wrapper is planned, not reverified by this review.
- The48/48 train/validation contexts are deliberate. Frozen locked comparison will use96/64 later; context-budget shift needs explicit labeling.

No further live job polling performed after these clarifications.
