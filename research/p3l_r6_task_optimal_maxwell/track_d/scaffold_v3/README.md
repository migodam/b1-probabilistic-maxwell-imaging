# Track D scaffold v2 — NOT RUN

Root-selected candidate: **original frozen deterministic_matched**, three original seeds, low `native48_probit64` (Revision1 registry) versus the **same estimator** on full `native96_probit64`. Failed Native-only matched heads are excluded. No new model training. Prior scaffold sources are preserved unchanged under `scaffold_v1/`.

Only AST static checks were run. No D data/anchors/predictions/GT or shared budgets were changed. Root owns qualification, budget reallocation, execution and later GT unlock. `source_subset` remains NOT_RUN/unsupported.

Root creates `track_d/D_FREEZE.json`:

```json
{
  "status": "ROOT_FROZEN",
  "execution_authorized": true,
  "C_finalist_qualified": true,
  "groups": ["DLOW", "DHIGH"],
  "seeds": [2026102601, 2026102602, 2026102603],
  "selected_low_config": "native48_probit64",
  "anchor_max_wall_seconds": 120,
  "primary_full": "full_deterministic",
  "primary_selected": "selected_deterministic",
  "distill_reference": "full_distill",
  "estimator_arms": {
    "full_deterministic": {"original": "deterministic_matched", "context": "full"},
    "selected_deterministic": {"original": "deterministic_matched", "context": "low"},
    "full_strong": {"original": "strongCNN", "context": "full"},
    "full_distill": {"original": "distilled_w128_n32768", "context": "full"}
  },
  "statistics": {
    "bootstrap": 20000,
    "rng_seed": 202610270699,
    "familywise_alpha_global": 0.05,
    "conditions": 2,
    "familywise_alpha_per_condition": 0.025,
    "comparisons_per_condition": 3,
    "comparisons_global": 6,
    "per_comparison_alpha": 0.008333333333333333,
    "NI_fraction": 0.01,
    "wall_fraction": 0.8
  }
}
```

Root must additionally fill `manifest_identities` for DLOW/DHIGH and `frozen_checkpoint_identities` (a list) with actual `{path:absolute,bytes,blake2b}` identities. Freeze **all loaded original dependencies**, including Gaussian provider, deterministic heads, Strong heads and Distill provider/head checkpoints across three seeds. The guard rechecks supplied identities before commands; root owns confirming this list covers the historical loader's exact dependencies. The120s wall above is an example, not an authorization. No scaffold creates this freeze or changes resource caps.

After root freeze/reallocation, from R6 research directory:

```sh
../../../../.venv_nn/bin/python -B track_d/online.py generate DLOW
../../../../.venv_nn/bin/python -B track_d/online.py generate DHIGH
../../../../.venv_nn/bin/python -B track_d/online.py anchors DLOW
../../../../.venv_nn/bin/python -B track_d/online.py anchors DHIGH
../../../../.venv_nn/bin/python -B track_d/online.py predict DLOW
../../../../.venv_nn/bin/python -B track_d/online.py predict DHIGH
```

Each group has256 parents, original prior and full forward, sigma×.5 or×2. Both full96/64 and low48/64 paths are genuinely paid in the singleton stage-D ledger; no cached Native-only candidate or free counterfactual low path. Original `track_c.anchors.run` accepts the merged original/Revision1 AnchorConfig registry.

Anchor progress prints after each completed parent/config and appends `ONLINE_PROGRESS.jsonl`. Completed contexts resume only if Book CLOSED=COMPLETE, report GT_read=false/path_complete, exact freeze identity and context identity all match. Partial jobs/directories fail closed for root audit; no implicit replay. Generation uses the separate dedicated subprocess and refuses existing observations. Completed prediction seals can be revalidated and global-sealed without repredicting; partial prediction directories require root audit.

Every configured original action runs all three frozen seeds. Full-anchor field defines phase labels (never truth). Both complete group seals, every prediction/phase identity, exact parent order and same freeze are verified before `ALL_512_PREDICTIONS_SEALED.json` is written.

Only then may root create `GT_UNLOCK.json` with `status:"ROOT_UNLOCKED_AFTER_SEAL"` and `global_seal_identity:{path,bytes,blake2b}` for the exact global seal, then run a separate process:

```sh
../../../../.venv_nn/bin/python -B track_d/score.py
```

The scorer verifies all512 predictions before first truth access. Parent mean of three seed **losses** is primary; no seed resampling or prediction ensemble substitution. Two noise conditions share one global six-comparison Bonferroni family (global alpha=.05, per-condition=.025, per-comparison=.008333333333333333, one-sided confidence99.1666667%). Per noise condition:20,000 paired parent draws; negative adjusted upper bounds support:

1. selectedDet −1.01 fullDet risk noninferiority;
2. selected wall −.8 full wall;
3. selectedDet −fullDistill zero-shot risk difference, **secondary**, never used for champion selection.

Overall reduced-anchor qualification uses1+2 only. FullStrong/fullDistill raw risks and full-anchor phase contributions remain descriptive references. No pooling of the noise conditions; global6 multiplicity correction is explicit. Zero-shot noise robustness only, no acquisition/source experiment and no posterior/probability claim. Bootstrap parent indices and contrast draws are saved for exact statistical replay; identical parent draws serve all three comparisons within condition.

Timing includes actual independently paid anchor wall plus measured batch64 original inference per parent; checkpoint loading is separately reported. Inference is one pass per seed, not a repeated hardware-variability experiment. Runtime is NOT_RUN until root authorizes; root may require additional no-GT timing qualification before a deployment speed claim. Missing/changed assets stop without fallback or extra model fitting.

Scoring additionally saves raw/fixed-passive-clip MAE and ROI MSE/MAE, separate prediction-ensemble metrics, phase cell counts/per-cell risk/contributions from full-anchor phase, full/low physical-point material baselines, optimizer status/fallback and phi summaries, and compact parent×seed CSVs. These are secondary diagnostics, never model or loss selection. Before first D truth, sealed full/low context bytes and REPORT identities are verified for all512 parents.
