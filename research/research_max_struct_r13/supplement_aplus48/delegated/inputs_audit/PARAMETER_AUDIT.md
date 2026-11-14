# A+ completion audit — FOLLOWUP_EXISTING48

This is completion of the original A+ arm on the same existing48 public observations. The user has already seen R13 results. It is **not fresh blind confirmation**, and the missing32 are not32 newly generated independent materials. No GT, scored material metrics, NPZ payloads, solver or remote resource was read/run in this audit.

## Exact reused16 and missing32

The frozen [mechanism16 configuration](../../../configs/confirmation_mechanism16.json) contains exactly16 `method=algorithm=a_plus` tasks, IDs0000,0003,...0045. All16 have original `COMPLETE` receipts and exit returncode0. The [core48 configuration](../../../configs/confirmation_core48.json) supplies the same48 original public observation paths; its protocol is literally equal to mechanism16's protocol. Every reused A+ receipt.task/protocol matches its frozen config. The shared public manifest declares576 real unknowns and initial material0.1; each observation reference and current byte/mtime metadata was recorded without opening arrays.

Exact case/task lists, all48 input references and original16 metadata/cost receipts are in [PUBLIC_INPUT_REUSE.json](PUBLIC_INPUT_REUSE.json). The32 complement indices are:

`0001 0002 0004 0005 0007 0008 0010 0011 0013 0014 0016 0017 0019 0020 0022 0023 0025 0026 0028 0029 0031 0032 0034 0035 0037 0038 0040 0041 0043 0044 0046 0047`.

Task IDs are exactly `r13_fresh_INDEX__a_plus`. Public input paths are exactly `data_public/fresh48/r13_fresh_INDEX.npz`, with no regeneration/noise copy/NN replacement. The reused16 are paid references, not16 reruns and not extra research consumption.

## Actual historical A+ arm

|Item|Frozen behavior|Supporting source|
|---|---|---|
|Dispatch|R8 `arm_config('a_plus', common, Y, E, [mf12/maxiter1500])`; common host0.1, box[0.02,20], observation-derived DPCS rank|C2 runtime ZIP `src/run_experiment.py:145–170`, captured with exact line numbers in [FROZEN_RUNTIME_EXCERPTS.json](FROZEN_RUNTIME_EXCERPTS.json); current corroborating [runner](../../../src/run_experiment.py#L194)|
|Algorithm identity|Original R7 Möbius/MM `mobius_material_readout`; `g19=False`, `mode='restricted'`, `initialize_physical_current=False`|Same C2 ZIP lines155–170; [full_current](../../../src/full_current.py#L138) installs G19 readout only if `g19`; all16 execution receipts agree|
|Transported stages|mf3:500, mf6:500, mf10:500, mf12:1500. beta_current=1 throughout; beta_material=1,0.2,0.2,0.2|[R8 arm dispatch](../../../../research_max_struct_r8/src/benchmark_solver.py#L13), then [R6 arm_config](../../../../research_max_struct_r6/src/benchmark_solver.py#L34)|
|Current update|Global PR nonlinear CG, sourcewise_cg=False; scalar current transform, onehop material transform; no independent sourcewise A+ variant|Same R6/R8 arm source; all16 AUXILIARY_EXECUTION.arm/sourcewise_cg/current_solver agree|
|Source metrics|Data source weights1/max(sum_abs_Y_squared,1e-300), state weights1/max(sum_abs_E_squared,1e-300); rho default1|[R6 arm_config](../../../../research_max_struct_r6/src/benchmark_solver.py#L64), [R8 loop](../../../../research_max_struct_r8/src/cie_fft.py#L115)|
|Möbius material update|Default2 four-color conditional-MM sweeps; step_tol1e-8; quintic stationary-root/endpoints/current candidates; aplus_arc_candidates=False and aplus_inherited_gtol=False; exact fixed-current loss fallback|[R7 readout](../../../../research_max_struct_r7/src/aplus_readout.py#L104), especially lines120–166. No config override in frozen dispatch|
|Prior|Inherited scalar regularizer smooth=1e-4, TV=1e-4, epsilon0.01; actual original prior object passed to A+|Frozen protocol, runner `aplus_regularizer=model.reg`; all16 receipt.prior equal|
|Stopping|Global coefficient_tol0, miniter50. Native maximum schedule totals3000 updates, but actual72CPU-second proposal cap can interrupt stages|[R6 frozen settings](../../../../research_max_struct_r6/src/benchmark_solver.py#L45), actual16 receipt stopping_policy/effective_stage_coefficient_tolerances|
|Acquisition/DPCS|Original per-case observation; `noise_ball_rank(S,Y,sigma)` unchanged. Existing16 receipt DPCS ranks are8/11/12/13, not a new common fixed rank|Runner dispatch; all16 DPCS_RULE.json metadata recorded in reuse ledger|

`AP_preparation=128`, `scale=6`, `omega=1` and `steps=3000` exist in the common protocol, but they do not replace A+'s inherited native stage/current/MM schedule. Generic material_inner_maxiter60/ftol1e-12/gtol1e-8 are configured by R6; the overridden Möbius readout uses its own two-sweep/step rule, with inherited-gradient guarding disabled by default. These fields must not be mistaken for new active A+ tuning.

## Budget, audits and retained estimators

Literal shared protocol: total CPU120s, proposal_fraction0.6 (72CPU-second proposal window), BLAS1, original8P-core queue, native audits plus next completed material update after each absolute2-process-CPU-second bucket; prefixes15/30/60/90/120. `audit_every50` is fallback when the CPU cadence callback is absent. Native stages and all paid physics audits/corrector remain charged; the32 follow-up must not extend native schedules or proposal fraction.

The unchanged physical corrector is L-BFGS-B: maxiter3000, maxfun12000, ftol1e-12, gtol0.01 on the objective/gradient multiplied by1e6, maxls30. It starts from the proposal's frozen incumbent, with the original material prior and box. Primary selection is minimum original paid Phi among legal eligible states; secondary is minimum originalR among legal paid noise-ball candidates, fallback primary. These are source rules, not a new material-quality gate. Frozen C2 runner excerpts include dispatch, corrector and incumbent rules. Later `min_change_paid_pool` is a separate public retention policy and does not alter this original primary identity.

## Provenance and remaining uncertainty

[CONFIRMATION_FINAL_BUNDLE](../../../audit/CONFIRMATION_FINAL_BUNDLE.json) names actual `runtime_confirmation_C2_final.zip`, shipping commit `689850d725ddcbd6b65152595ac0ddf2b6eff1bd` and algorithm freeze `3a91b35605634a18a731441693a46152eeb3363e`. Only its two necessary source members were read; exact relevant numbered lines and member byte counts are saved in the excerpts. No SHA/code-hash validation was performed.

The runtime ZIP contains the R13 runner/adapter, **not** inherited R8 benchmark_solver/cie_fft or R7 aplus_readout. The actual original execution dependency is `D:/AI/MAX_STRUCT_R11_DEV_20261107`. Local inherited sources and old execution receipts corroborate the arm, but this audit did not remotely verify dependency bytes or package those modules. Parent must reuse the unchanged original legacy dependency tree in the isolated32-task runtime; a runtime ZIP alone is not a self-contained inherited solver release.

Current public input byte/mtime records identify the existing local files; they are not cryptographic proof of historical byte equality or portable restored seals. All16 original receipt tasks point to those same original paths; original seals/provenance remain authoritative. This audit does not read pixel truth to judge any case or infer future A+ quality. The32 new run receipts, external costs, failures, sealing and follow-up authorization remain future parent work. Historical16 and new32 host/cost provenance must remain distinguishable.
