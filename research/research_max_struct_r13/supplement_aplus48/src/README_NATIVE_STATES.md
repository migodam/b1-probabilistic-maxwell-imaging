> 执行状态更新（2026-11-14）：已完成40任务、561状态、844唯一材料物理诊断，0问题。实际结果位于analysis/native_states，科学解释见TRAJECTORY_ANALYSIS；下面保留原实现交接说明，不再表示未运行。

# Native-state analyzer: source-only handoff

Implementation: `supplement_aplus48/src/analyze_native_states.py`. **Not executed.** Only source compilation has been performed. Parent chooses launch timing after terminal diagnostic queue and reviewed authorization/seal. No new GT NPZ, scoring, solver run or case choice occurred during implementation.

## Invocation

```text
<existing-python> supplement_aplus48/src/analyze_native_states.py \
  --parent-root <R13-root> --sup <supplement_aplus48-root> \
  --diag-root <downloaded-diagnostic-runtime-root> \
  --manifest <sealed-aplus48_dynamics8.json> \
  --compareauth <SUP/audit/DIAGNOSTIC_AUTHORIZATION.json> \
  --legacy-repo <exact-existing-legacy-repository> \
  --out <new-analysis-directory>
```

Optional `--no-figures` suppresses PNG/PDF creation, not computation. The code fixes BLAS threads to1 before numerical imports. Parent should capture external whole-child CPU/wall; analysis receipt timing excludes its own final write and interpreter exit.

The new authorization must retain the existing compare48 fields/seals and additionally contain:

```json
{
  "authorized_diagnostic_analysis": true,
  "diagnostic_campaign": "aplus48_dynamics8",
  "diagnostic_seal": "<absolute-local-reviewed-diagnostic-seal.json>"
}
```

The existing `COMPARISON_AUTHORIZATION.json` stays unchanged. The diagnostic seal uses existing `seal_results` schema, including parent-reviewed `authorized_scoring=true`, terminal queue status and byte/mtime inventories. These flags authorize this offline diagnostic, not primary replacement. The actual diagnostic queue must also be terminal with no running/pending tasks. Manifest must equal the seal snapshot. All guards run using stdlib **before** NumPy/model imports and public/material array access.

Original/extension prediction seals are revalidated through compare48.authorize. Each task's primary_reference_folder must be one of those sealed runs for the exact case/method. Public inputs are loaded from the original core48 manifest case reference beneath parent root, never from the Windows absolute diagnostic input_file. No observation/GT path is guessed.

## Diagnostics and scope

**Original physical objective.** For CIE-family snapshots, evaluate original beta=None Phi/gradient at m_before_readout and m, with the same declared Regularizer/box. Retain actual DeltaPhi, gPhi_before dot Delta m, raw and box-gradient-mapping cosines, exact box activity, gradient norms and unit-step projected-gradient infinity norms. Identical materials within the same case/public operator/prior reuse one paid evaluation via np.array_equal. Each unique material incurs measured actual original LU/source/adjoint RHS; retained state is reused and no Hessian/optimizer is computed.

**Counterfactual fixed-current material gradients.** At the actual fixed postcurrent J/Q, compute ordinary scalar CIE at beta_material using

`r=beta*(J-mQ)/(1+beta*m)`

`dr/dm=-beta*(Q+beta*J)/(1+beta*m)^2`

`g=rho*Re sum(conj(dr/dm)*(r Os))+g_prior`.

Also call the exact existing onehop_material_value_gradient with OneHopAbsorption(original public G), saved rho/Os, and original prior. Record before/after values, unit-step box gradient mappings, movement alignment and gradient differences. These are **counterfactual ordinary-CIE/onehop gradients even for G19**, not G19's native material derivative. Saved Os/Ls/rho are required; absent context remains unavailable. Block-context agreement with saved arrays is recorded without a new backend or root/pole sweep.

**Field defects.** For saved J/Q and each relevant m, record norms of u=J−mQ, v=Q−E−GJ, and u+m v; current and receiver discrepancy from the same paid physical J(m); auxiliary and physical residual/delta; and observed ||J−Jphys||/||u+m v|| (null at zero denominator). It is a directional amplification, **not** an operator norm. The Maxwell relation is used to organize these defects; no numerical identity residual/error-floor check is computed. One GJ and SJ is reused across before/after readout. Port difference uses already available SJ and physical prediction.

For DR, analyze c=(scale*cX,cQ), the reflected A output=(scale*aX,aQ), and the pre-iterator=(scale*preX,preQ), all against the same physical cm. Compute actual full-state memory/gap/update norms from saved arrays and retain multiplier metadata. No extra A or C projection/background factorization is constructed. Changes of cm between retained checkpoints are labeled **checkpoint changes, not one-step changes**. Post-iterator material projection stays null.

**Current-force transport (parent-specified extra diagnostic).** Rebuild only the existing NativeFourierAPCS from saved shape/mf/Vd, with no SVD. At the same postcurrent J/Q/A, use **beta_current** and saved Od/Os/rho:

`p=beta_j/(1+beta_j*m), R=p*(J-mQ)`

`gA(m)=U^*{S^*[(SJ-Y)Od]+rho[p*(R Os)-G^*((m*p)*(R Os))]}`.

Evaluate at pre/post readout m and form tau=gA_after−gA_before. Record its norm, alignment with saved native current direction, gradient rotation, normalized and raw directional derivatives before/after readout, and native **pre-current** gradient norm. The recomputed prereadout gradient is at **postcurrent** fields and is distinct from the saved pre-current PR gradient. No orthogonality identity test, PR update, material forward or parameter search is added by this component. It reuses SJ, adds one source-block S-adjoint plus two G-adjoints per snapshot, and separately counts native Fourier/DPCS actions. Context failure is reported separately without discarding existing snapshot diagnostics.

This characterizes changing current force/PR context after material readout; it is not a chaos/global-basin proof or causal prediction on the8 posthoc selected objects.

## Replay fidelity, labels and denominators

Every sealed task remains in TASK_COVERAGE, including missing/failed traces. Every requested target gets a RAW_TARGET_CHECKS row. Requested/reached/unreached flags, actual stage/global/local iteration, snapshot process CPU, aliased retention slots and stop context are preserved. DR snapshot CPU comes from its saved per-step scalar record when absent in the snapshot scalar dictionary.

At requested iterations, the full saved material vector (CIE m or DR cm) is compared with original primary run's **raw.npz m**, recording exact array equality, maximum absolute difference and RMS difference; the full difference vector is banked. An unreached/different material is not an exact replay. Even exact material equality does not establish auxiliary J/Q/PR-state equality. Original raw candidate and iteration are included so a target for a different proposal audit is not misnamed as the raw state.

Optional quality labels are copied only from already-authorized comparison48/ALL48_SCORES.csv (raw/primary rows), with comparison seal context checked. No GT array is opened, no score for a new checkpoint is computed, and labels never select snapshots/cases or modify primary prediction references.

## Outputs

Each principal table has CSV (LF line endings) plus JSON:

- TASK_COVERAGE, SNAPSHOTS, READOUT_TRANSITIONS
- FIELD_DEFECTS, CURRENT_FORCE_TRANSPORT
- RAW_TARGET_CHECKS, PHYSICAL_EVALUATIONS
- DENSE_NATIVE_HISTORY, PAID_PHYSICAL_AUDITS
- CASE_COSTS, ISSUES

GRADIENT_BANKS/<case>.npz retains physical m/gPhi/gPrior, counterfactual material gradients, coefficient current forces/tau and raw-target difference vectors. Data retention/cache is per case; physical J/LU are not banked as new predictions. Only the active model retains a LU; exact duplicate evaluations reuse stored results.

PNG/PDF per case show dense native objectives/readout steps/PR restarts; physical-readout movement/gradient alignment; observable field defects; actual paid physical audit trials; current-force transport; and sparse-as-recorded DR audit gaps. Native scalar CIE histories contain native elapsed wall, not per-update CPU; missing CPU stays null. Selected snapshot CPU is preserved. Different auxiliary losses and original Phi are never merged into one objective.

Costs include whole analysis CPU/wall, per unique physical evaluation actual counters and nested timings, per-case setup/exports/figures, explicit diagnostic G/S/DPCS actions and separately counted onehop/current-transport/basis work. Norm counters cover explicitly instrumented Frobenius contractions, not all internal library operations; whole CPU covers the rest. Nested costs are not added twice. Failed work and missing data are retained; null is not imputed as zero.

This is same-discretization post-analysis, not independent physical validation or a stationary/local-minimum/recovery certificate. Parent owns mathematical and scientific interpretation.

