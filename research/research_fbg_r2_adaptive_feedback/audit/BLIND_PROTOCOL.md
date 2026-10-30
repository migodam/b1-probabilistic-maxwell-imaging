# FBG R2 finite experiment and blind boundary

## Scientific scope

The experiment uses the existing 24×24, 576-cell scalar-TM nonlinear volume-integral model with 48 receivers. Stage 1, Stage 2, and Stage 3 use previously disclosed R4 truth and are **DEVELOPMENT**. Stage 4 is a fresh, separately generated holdout. A low physical residual does not by itself establish material recovery, a unique inverse, a new information channel, or generalization outside these finite cases.

The frozen R4 material gate is NRMSE ≤ 0.15 and reconstruction-only Otsu IoU ≥ 0.70. The data gate is full nonlinear relative residual ≤ 0.001 for clean observations and ≤ 0.003 for noisy observations. Joint success requires both. Report RMSE, NRMSE, MAE, Otsu IoU, full physical residual, LU, actual complex RHS, wall, stop reason, failures, and `NOT_RUN`. Do not silently drop missing or nonfinite endpoints. The known-midpoint threshold is offline oracle information and must never be used in selection.

## Frozen finite stages

1. Stage 1 has five clean disclosed R4 observations in `data_contract/stage1_cases.json`: low-contrast curved bar, moderate ring, high-contrast digit 8, high ring, and high curved bar. Pair frozen baseline and adaptive feedback on each. Preserve trigger, accepted updates, full physical residual, endpoint, and work counters. Any refinement chosen from these cases is development. `stage1_mechanism.py` writes only a small diagnostic case/arm plan; the parent's frozen `data_contract/stage1_remote.json` is the actual 64-task execution manifest and supersedes that diagnostic plan.
2. Stage 2's parent-owned actual execution plan has 360 jobs across six arms: native, reciprocal, data, GD all, β=0.2, and new qβ. It uses the same public initial states: K=1,4,8 prefixes, with per-start nested K=8 results and fixed-total K=4/8 runs. The eight starts are constants 0.1, 1, 4; two smooth random fields; two correlated fields; and one block field, all clipped to [0.02, 20]. Their public seed is fixed before inversion. Each candidate is an independent optimization run. The selected endpoint minimizes the finite public nonlinear physics ranking score across returned starts, with no GT access. Report **per-start** budgets, whose total allowance rises K-fold, separately from **fixed-total** budgets split across K. Report actual measured aggregate costs in both; do not call K independent samples or call a solver trajectory multiple roots. The older `stage2_multistart.py` planner is diagnostic and is not the parent execution manifest.
3. Stage 3 uses two predeclared hard high-contrast disclosed cases, digit 8 and ring. It constructs p = 3, 6, 12, 24 physical disk-pulse illuminations with 48 unchanged receivers. Uniform p=3/6 are subsets of the original 12 positions and uniform p=24 adds angular midpoints. A seeded jittered random 12-position geometry has nested p=3/6 prefixes; p=24 adds 12 new physical positions. The incident field formula is exactly `0.25j*hankel1(0,k*distance)`. One model LU handles all p source RHS for a fixed material. `per_source` gives every source amplitude 1, so total incident amplitude-squared grows with p. `fixed_total` would scale every source by sqrt(12/p), holding the sum of source amplitude-squared at the p=12 reference. Because this global scale multiplies both E0 and Y, the normalized data and normalized dual objectives are analytically unchanged; record this equivalence without duplicate optimization. The 16 clean observations arise from 16 physical forwards and only eight acquisition operators. Avoid interpreting source count or energy changes as an isolated algorithm gain.
4. Stage 4 is generated **only after** final model, selection rule, budgets, and the method/start plan have been frozen by the parent. It contains 14 new objects: nine deformed raster digits 1-9 in DejaVu Sans/Serif Bold (different from R4 Arial/Georgia), a deformed ring, two islands, and curved bar, all A=12; plus a ring A=0.4 and curved bar A=2.5 as contrast controls. Digits receive secret rotation, translation, and thickness changes; geometry parameters are secret. All 14 have a clean observation. Only predeclared hidden shape slots 0, 3, 6, 9 have two independently drawn 10^-3 norm-noisy observations each. This is 22 observations from 14 objects, not 22 independent objects; no other noisy observations are planned or counted `NOT_RUN`. Uniform and seeded jittered random physical 12-source layouts alternate by hidden shape slot. They use the same 576-cell G and 48-receiver S and two public operator NPZs. Each public record has its own `operator_path`; the solver sees source positions and Y but no slot, digit decoder, or shape template. Layout and shape effects are confounded, while each method is paired on the same observation. Case IDs and order are randomized. Fresh truth, labels, amplitudes, backgrounds, and generator entropy live only in a private local directory, outside the public payload or solver upload. No solver process, adaptive tuning, or parent interpretation reads these until all planned prediction dispositions are frozen.

## Source mixing

For a fixed p, invertible right mixing maps `E0 -> E0 M` and `Y -> Y M`. It spans the same incident subspace and adds no independent information. If original source noise has covariance σ²I, the mixed source covariance is σ²MᴴM. Use the corresponding generalized least squares rule; right-multiplying a mixed residual by `M^{-1}` restores the original norm and is a valid whitening implementation. Compare mixed and unmixed records as a numerical/control check only. A rank-deficient mix loses information and is excluded from the invertible control.

## Blind close, score, and leakage rules

- Public generation writes only `Y`, opaque ID, two source-only public operator files, and a manifest of IDs/noise/`operator_path`. `private_out` must be a separate directory. Never upload `mapping_PRIVATE.json` or private NPZ to GPU/solver hosts. Do not pass private entropy, amplitude, label, hidden slot, or private path to inversion code. The public source-layout seed describes acquisition only.
- Freeze a JSON mapping of `method -> planned_start_count` before generation and pass it with `--method-starts`. Save exactly one aggregate prediction JSON per public observation ID. Format: `{"observation_id":"obs_...","methods":{"method_name":{"starts":[{"start_id":"s0","status":"OK","result":{"m":[576 real values],"ranking_objective":0.0,"relative_residual":0.0,"counts":{...}}}, ...]}}}`. A direct method payload with `m` is accepted only as its first start; all additional planned starts become `NOT_RUN`. The parent may retain richer raw receipts separately. Adaptive choices and online start selection use only public observations and physical ranking scores.
- After all planned tasks finish or are explicitly closed as failed/`NOT_RUN`, run `stage4_blind_imaging.py seal --close`. It copies every planned method/start material endpoint, public score, and missing/failure disposition into one write-once JSON before truth access. Unplanned methods/starts are preserved as protocol deviations and excluded from planned scoring. Preserve raw receipts. Do not rerun or retune a failed case after seeing any score.
- Run `score` only against this closed seal. It performs no physics and writes per-method online top-1 gates and separately labeled offline best@K from the finite planned starts. Offline best@K uses GT and is not a deployable selector or independent-root evidence. Clean/noisy are paired variants of the same object, not independent objects. Summarize equal object groups and disclose 14 objects/22 observations. Preserve costs and any protocol deviation separately.

## Commands and interfaces

From the new study directory, using an existing Python with NumPy/SciPy:

```sh
python experiments/stage1_mechanism.py
python experiments/stage2_multistart.py plan --output data_contract/stage2_plan
python experiments/stage3_illumination.py --out data_contract/stage3_dev --go
```

Stage 3 reads disclosed R4 private material only in its local generation process. Its `public_data/*.npz` contains `Y` and a relative `operators_path`; each manifest record also exposes a task-level `operators_path` relative to the manifest directory. The referenced operator NPZ contains matching G/E0/S and actual source positions. The inference adapter must support arbitrary p in E0 and Y, while preserving the same scalar-TM physics and shared 576-cell material. The new adapter must accept explicit initial materials for Stage 2; the archived frozen R4 solver enforces 12 sources and a fixed 0.1 start and is not silently repurposed.

Stage 4 is deliberately **not generated by preparation**. After selection/freeze:

```sh
python experiments/stage4_blind_imaging.py generate \
  --public-out /separate/local/public_data \
  --private-out /separate/local/private_truth \
  --method-starts /separate/local/method_starts_FROZEN.json --go
python experiments/stage4_blind_imaging.py seal \
  --manifest /separate/local/public_data/manifest.json \
  --predictions /separate/local/predictions \
  --output /separate/local/audit/ALL_PREDICTIONS_SEALED.json --close
python experiments/stage4_blind_imaging.py score \
  --seal /separate/local/audit/ALL_PREDICTIONS_SEALED.json \
  --private-dir /separate/local/private_truth \
  --output /separate/local/RESULTS_SUMMARY.json
```

Those placeholder paths must be replaced with actual isolated local paths. `seal` is a closure action by the parent, not permission to inspect truth early. The score artifact is the first truth-disclosing object.
