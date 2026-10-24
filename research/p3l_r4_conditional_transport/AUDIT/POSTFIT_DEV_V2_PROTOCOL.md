# Post-fit DEV and corrected-MH plan: inert

`configs/POSTFIT_DEV_V2_INERT.json` is immutable and `released:false`. Separate parent phase releases are required; do not toggle the inert protocol itself. Controller: `tools/r4_postfit_dev.py`. No new training, inference, GT read or physical action was launched during preparation. No source frozen by the running V2 fits was modified.

The plan binds `PHYSICS_TRANSPORT_PAIRED_V2_RUNTIME_RECOVERY.json`, actual `transport_v2_*` jobs, and `models/physics_transport/paired_v2`. It retains the V1 startup failure:1.936266792s,0RHS,0 optimizer updates, failed MPS float64 reduction. No budget or history is reset.

## Generator evaluation

Nine seed2026102501 methods: initial, physics-refined and supervised-continuation checkpoints for raw spline, matched raw diagonal Gaussian, and anchor spline_none. Each source must be CLOSED COMPLETE before **any** new inference; refinement sources must have completed3000 updates. Pending receipts are reported without reading live checkpoint/result contents or identities. No substituted partial/last/checkpoint is allowed.

All nine methods use exactly `inputs/anchor_conditions/dev` and DEV00000..00127. Raw-y loaders read only original observation features; anchor loaders additionally read explicit y-derived context_x. K128, source-major original observations, common sample seed20261025031, fixed index0, CPU batch<=16. Freeze methods, exact checkpoint/receipt/result identities, target, all context/observation identities, evaluation code and all settings **before inference** in `FINAL_PREINFERENCE_SPEC.json`.

The existing independent evaluator records raw full66 proper NLL; full-grid material energy; raw sample mean MSE/MAE;16 real material+4 spatial projections; phase-event scores; marginal empirical intervals and randomized inclusive rank coverage50/80/90/95; full material mean/variance; all latent draws and exact algebraic full66 logq. Primary paired comparisons are physics vs supervised material energy within each of the three families. Other scores are descriptive. Existing bootstrap uses4000 parent resamples with seed20261025099; marginal ranks and finite scores are not Full posterior certification.

Inference and offline scoring run in different child PIDs, each under its own zero-physics analysis Book. All methods/parents must finish and pass the inference gate before any new offline truth score. The controller itself never reads truth arrays. Generator failures HOLD the stage; no model, sample or parent is replaced. The panel was previously exposed and is not a fresh holdout. Raw validation uses2048 parents, anchor validation256; their training/data budgets remain separately labeled.

## Separately released Full checks

For **each of nine models**, retain sample0 for all128 parents, plus first64 samples for first8 parents (sample0 overlap counted once):632 Full evaluations/model. Total5688 evaluations,45504 RHS. Exact full66 logq comes from the sealed original128-draw artifact, with no density assigned to selected/corrected samples. Original complex observations come from `inputs/online`; keep current/field arrays. No resampling, residual selection or duplicate deletion.

The paid driver freezes the same predeclared roster after artifacts exist, then performs a bounded Full job (StageI,45504 hard RHS,2400s). Different-PID scoring reuses already paid offline truth predictions for fixed likelihood ranks and the fixed sample0 Gamma256 marginal diagnostic. Those are necessary calibration diagnostics with finite-data/IS limits; neither a pass nor NLL improvement establishes full-posterior qualification. All failed slots remain explicit.

## Separately released corrected-MH pilot

Freeze a **single validation-selected explicit checkpoint before DEV inference** using the generator release's `MH_validation_binding`. It must name one of the nine CLOSED methods and contain a nonempty validation-only rationale. There is no automated winner across raw VAL2048 and anchor VAL256, and DEV scores cannot select the checkpoint. This binding remains sealed in the final preinference spec even if MH is never released. CFM trace density is forbidden.

Pilot parents are fixed DEV00064..00067. Four common starts per parent: one online anchor plus three independent q draws, identical across all kernels.128 post-step states/chain, zero extra warmup, fixed seed20261025141. Kernels: prior pCN beta.002, q-independence, and a constant state-independent half mixture of the two individually Metropolized components. Full66 target is `exp(-Phi-.5||x||²)` in Gaussianized dx66. Independence log ratio is `Phi(x)-Phi(z)+.5(||x||²-||z||²)+logq(x)-logq(z)`; pCN acceptance uses `Phi(x)-Phi(z)` by prior reversibility. Mixture uses the selected component's correct ratio, not an incorrect mixture-density ratio.

All rejection states and endpoint references remain in the chain. No thinning, replacement, corrections/projections, GT, new Jacobians or automatic expansion. All starts are separately Full paid per kernel:4 parents×3 kernels×4 chains×129 evaluations=6192 Full evaluations=49536 RHS. StageIII,49536 hard RHS,2400s. The existing corrected runner reports same-parent finite-chain diagnostics, stuck-chain unavailability and no common Full reference. Ratios are formally correct for the explicit law; float32 model arithmetic is not a machine-exact MH certificate. No convergence or posterior qualification is assumed.

## Phase releases and commands

Inert templates exist for GENERATOR, PAID and MH under `configs/POSTFIT_DEV_V2_*_RELEASE_TEMPLATE.json`. Parent creates a reviewed release file with `released:true`, the unchanged protocol identity and matching phase. Generator additionally binds the MH method/rationale before inference. Paid/MH releases additionally bind the identity of `FINAL_PREINFERENCE_SPEC.json`. MH includes existing runner experiment/campaign IDs. Each stage refuses prior output; there is no automatic replay or blocking wait for current fits.

```text
<research-python> tools/r4_postfit_dev.py check
<research-python> tools/r4_postfit_dev.py generator --release <parent-generator-release.json>
<research-python> tools/r4_postfit_dev.py paid_probability --release <parent-Full-release.json>
<research-python> tools/r4_postfit_dev.py corrected_MH --release <parent-MH-release.json>
```

Generator artifacts: `results/postfit_dev_v2/{FINAL_PREINFERENCE_SPEC.json,online,scores}`. Paid artifacts: its `paid_probability/` subdirectory. Pilot artifacts: `mh_prepared/` (same explicit context exports) and `mh_physics/`. Pure freeze/inference/score commands are ledgered analysis jobs; paid runners own physical Books and never run inside analysis-only Books.

## Readiness evidence

`AUDIT/postfit_dev_readiness_1792878224220799000/READINESS.json`: PASS_INERT_PENDING_CLOSED_FITS, analysis-only Book COMPLETE, wall.129298s, CPU.079294s, RSS224772096B. Three initial sources were CLOSED; all six refinement sources were pending and no live checkpoint was inspected. Unreleased phases and missing preinference MH binding were rejected. Pure ratio tests give independence prior cancellation3.55e-15 and pCN cancellation2.84e-13; stuck chains do not receive ESS=N. No Target instances, physical RHS, optimizer updates, generated samples, GT reads or child launches. Exact inert protocol copy is retained beside the readiness artifact.
