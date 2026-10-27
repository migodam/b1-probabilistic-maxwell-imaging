# Track D scaffold V4 — NOT RUN

Root-selected candidate: original deterministic Q64 weights with **q KKT feature closure**, LOW `native24_probit64` (Revision1 registry), versus FULL `native96_probit64` with original actual-gradient deterministic features. Ordinary LOW actual-gradient and FULL q-closure remain explicit controls. Same NN weights, no new noise training, no Native-only matched candidate. Prior V3 sources preserved under `scaffold_v3/`; V1 under `scaffold_v1/`.

Only AST/source static review performed: no D data/physical solves/GT. Root alone freezes/executes/reallocates budget/unlocks scoring. Noise zero-shot only; source/acquisition experiments unsupported/NOT_RUN.

Root's `track_d/D_FREEZE.json` must include:

```json
{
  "status": "ROOT_FROZEN",
  "execution_authorized": true,
  "C_finalist_qualified": true,
  "groups": ["DLOW", "DHIGH"],
  "seeds": [2026102601, 2026102602, 2026102603],
  "selected_low_config": "native24_probit64",
  "anchor_max_wall_seconds": 120,
  "primary_full": "full_deterministic",
  "primary_selected": "selected_deterministic",
  "actualgradient_reference": "low_actual_deterministic",
  "distill_reference": "full_distill",
  "estimator_arms": {
    "full_deterministic": {"original": "deterministic_matched", "context": "full"},
    "selected_deterministic": {"original": "deterministic_matched", "context": "low", "feature_mode": "kkt_q"},
    "low_actual_deterministic": {"original": "deterministic_matched", "context": "low"},
    "full_kkt_deterministic": {"original": "deterministic_matched", "context": "full", "feature_mode": "kkt_q"},
    "full_distill": {"original": "distilled_w128_n32768", "context": "full"},
    "full_strong": {"original": "strongCNN", "context": "full"}
  },
  "statistics": {
    "bootstrap": 20000,
    "rng_seed": 202610270699,
    "familywise_alpha_global": 0.05,
    "conditions": 2,
    "familywise_alpha_per_condition": 0.025,
    "comparisons_per_condition": 4,
    "comparisons_global": 8,
    "per_comparison_alpha": 0.00625,
    "NI_fraction": 0.01,
    "wall_fraction": 0.8
  }
}
```

Additionally fill `manifest_identities` (DLOW/DHIGH) and `frozen_checkpoint_identities` with exact `{path:absolute,bytes,blake2b}` objects. Freeze every loaded provider/head dependency for all3 seeds. Optional `frozen_source_identities` uses the same schema and is verified each process; root should include online/common/score, shared physics, Track C anchors/frozen_heads, and kkt_closure/closure_action.py (plus needed source dependencies). Checkpoints/sources may be sourced from KKT `SELECTED_WEIGHTS_IDENTITY.json`/`CONTEXT_SOURCE.json`, but root owns completeness. Example120s is not an authorization. No source modification after root freeze except explicit recorded bug handling.

`feature_mode:kkt_q` lazily imports `track_c.kkt_closure.closure_action.Action(method='deterministic_matched',seed=seed,closure_mode='q',noise_ratio=r)`, using known r=.5/2 for the corresponding frozen noise condition. Its predict API is the original float32 measurement/context interface. This passes g=r·q0z to match original feature-training sigma with the new whitened Jacobian convention. No extra query, residual or Jacobian is synthesized; ordinary low/full controls use original inference unchanged.

After root freeze/reallocation, from R6 research directory:

```sh
../../../../.venv_nn/bin/python -B track_d/online.py generate DLOW
../../../../.venv_nn/bin/python -B track_d/online.py generate DHIGH
../../../../.venv_nn/bin/python -B track_d/online.py anchors DLOW
../../../../.venv_nn/bin/python -B track_d/online.py anchors DHIGH
../../../../.venv_nn/bin/python -B track_d/online.py predict DLOW
../../../../.venv_nn/bin/python -B track_d/online.py predict DHIGH
```

Each condition has256 fresh parents, original prior/full forward at sigma×.5/×2. Even frozen parent indices run full→low, odd indices low→full; pair order is recorded in reports/progress. Each path independently initializes from the original prior start with a fresh Target; no optimizer state/J is passed across paths. Both retain equal known-background/BLAS cache rights. Full96/64 and Low24/64 are each genuinely paid once; all six estimator arms share these paid contexts. No prediction starts until both noise groups have completed all512 full/low context pairs. After each path, transfer cache is cleared, target deleted and `gc.collect()` run. Actual cleanup seconds/object counts are recorded and included in anchor wall; Book close/admin receipt wall is separately excluded from commercial timing and still included in ledger job wall.

Completed anchor paths resume only with CLOSED=COMPLETE, GT_read=false/path_complete and matching freeze/context identity; partial jobs fail closed for root audit. Actual progress prints and appends `ONLINE_PROGRESS.jsonl`. Missing/nonfinite prediction aborts without fallback or fitting. Completed prediction seals can be revalidated without reprediction. Both conditions' predictions, phase labels, context+REPORT identities and same freeze must validate before immutable `ALL_512_PREDICTIONS_SEALED.json`; identical existing global seal is never rewritten. Online commands stop permanently after GT unlock/scoring.

After all512 predictions are sealed, root may write `GT_UNLOCK.json` with `status:"ROOT_UNLOCKED_AFTER_SEAL"` and `global_seal_identity:{path,bytes,blake2b}`, then run independently:

```sh
../../../../.venv_nn/bin/python -B track_d/score.py
```

Scorer verifies both groups' entire sealed inputs before first D truth. Parent mean of3 seed raw full-grid MSE losses is primary; ensemble losses separate. No noise-condition pooling. Global eight-comparison Bonferroni family: alpha=.05/8=.00625, one-sided confidence99.375%. Four contrasts within each condition, negative supports:

1. selected q-LOW risk −1.01 FULL originalDet risk;
2. selected wall −.8 FULL wall;
3. selected q-LOW risk −ordinary actual-gradient LOW risk (independent mechanism comparison);
4. selected q-LOW risk −FULL Distill risk (secondary, not champion selection).

`cost_GO`/overall qualification requires1+2; `mechanism_GO` additionally requires3. Failed/uncertain mechanism test is not converted into KKT superiority by cost qualification. Full-q versus Full-original remains a descriptive same-context control. No new loss or model selection.

Saved diagnostics: bootstrap parent indices+contrast draws (same parent draws for all4 comparisons within condition), compact per-parent/seed CSV, raw/fixed-passive-clip full/ROI MSE+MAE, separate ensemble actions, full-anchor phase counts/per-cell risk/contributions, physical full/low point baselines, optimizer statuses/fallback and phi. Phases never use GT. No full-posterior/probability claim. One batch64 inference timing per seed does not estimate repeated hardware variability; root may require additional no-GT timing qualification before deployment claims.
