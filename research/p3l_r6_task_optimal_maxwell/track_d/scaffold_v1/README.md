# Track D scaffold — NOT RUN

Only root may freeze/execute. No D data, anchors, predictions, scores, models or budgets were changed by scaffold creation. AST syntax checks only; no Target was constructed. `source_subset` is explicitly unsupported/NOT_RUN.

Root must create `track_d/D_FREEZE.json` with:

- `status: "ROOT_FROZEN"`, `execution_authorized: true`, `C_finalist_qualified: true`;
- `groups: ["DLOW","DHIGH"]`, `seeds: [2026102601,2026102602,2026102603]`;
- `manifest_identities`: DLOW/DHIGH exact `{path:absolute,bytes,blake2b}` for the existing config manifests;
- `selected_low_config`: actual existing Track C CONFIGS name (e.g. `native48_probit32`; unused when selected context is native);
- `anchor_max_wall_seconds`: root-assigned per-anchor Book wall cap;
- `matched_arms`: required `matched_full` with `context:"full"`, required `matched_selected` with `context:"native"` or `"low"`, optional additional matched arms. Every arm supplies `checkpoints:{"2026102601":{path,bytes,blake2b}, ...}` from frozen matched TRAIN512 heads. No new training occurs;
- `statistics` exactly `{ "bootstrap":20000, "rng_seed":202610270699, "familywise_alpha_per_condition":0.05, "comparisons_per_condition":3, "NI_fraction":0.01, "wall_fraction":0.8 }`.

Checkpoint identities are rechecked before any command. Per condition: 256 fresh parents, original Gaussian64×Uniform2 prior, full original forward at sigma×.5 or ×2. Native96+Probit64 always paid once as full baseline. Native96 context/state/J is captured **before** Probit using its already paid J. A selected Native candidate reuses this prefix without another solve; its cold standalone cost remains the original paid Native phase plus target/context-persistence cutoff, not zero. A genuinely different low path is separately paid. Original frozen Strong/Distill/Deterministic and all root-frozen matched arms run all three seeds.

From the R6 research directory, after root freeze only:

```sh
../../../../.venv_nn/bin/python -B track_d/online.py generate DLOW
../../../../.venv_nn/bin/python -B track_d/online.py generate DHIGH
../../../../.venv_nn/bin/python -B track_d/online.py anchors DLOW
../../../../.venv_nn/bin/python -B track_d/online.py anchors DHIGH
../../../../.venv_nn/bin/python -B track_d/online.py predict DLOW
../../../../.venv_nn/bin/python -B track_d/online.py predict DHIGH
```

Each condition seals every prediction and full-anchor phase labels with BLAKE2b. Both conditions must share the same root freeze. `ALL_512_PREDICTIONS_SEALED.json` appears only after both group seals exist. No score is run by the online command.

After all512 predictions are sealed, root may create `track_d/GT_UNLOCK.json` with `status:"ROOT_UNLOCKED_AFTER_SEAL"` and `global_seal_identity:{path,bytes,blake2b}` matching that global seal. Then, in a separate process:

```sh
../../../../.venv_nn/bin/python -B track_d/score.py
```

The evaluator verifies every artifact of both conditions before first D truth access. Primary estimates average three seed **losses** within parent; parent is bootstrap unit. Each noise condition has its own20000-draw/three-comparison Bonferroni family: selected-minus1.01full risk NI, selected-minus historicalDistill risk, selected-minus.8full wall. No pooled-noise claim. Raw baseline MSE and full-anchor host/transition/pure contributions are preserved; phases never come from truth. There is no probability/posterior claim.

Limitations requiring root review before execution:

1. Static syntax pass is not runtime validation; the controller uses Track C private trace/hook helpers and requires the original Native paid endpoint J before Probit. Missing assets stop without extra solve or fallback.
2. Native cold cost uses the measured Native phase of the paid full run. It is explicitly a counterfactual standalone estimate, not a separately measured deployment run. `standalone_speedup_confirmed=false` remains for Native-selected outputs even if its arithmetic wall contrast passes.
3. Inference wall is measured once per method/seed at batch64 and apportioned to parents. Parent bootstrap does not establish repeated timing variability; checkpoint load is separate. A deployment speed claim needs the root's repeat-timing qualification.
4. The native capture transiently persists current/field but no LU; actual shared campaign costs counted once. Cold standalone cost is separately annotated, never billed again as new RHS.
5. No implicit partial resume or overwrite. Interrupted Book jobs, partial context directories, prediction directories or scores require root audit.
6. Existing `src/generate.py` assigns DLOW/DHIGH to stage D correctly. Scaffolding does not alter shared authorization caps or guarantee that the remaining budget suffices for all256+256 full plus optional low paths.
