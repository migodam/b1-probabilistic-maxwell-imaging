# Confirmation scoring preparation V1

Status: **PREPARED_NOT_EXECUTED**. No confirmation seal, score, private mapping, material, clean observation, or noise payload was opened by this preparation. Source files were read and syntax parsed only. The later development phase-pi/grid32 authorizations were handled separately; their old13-row score table was preserved.

## Entry points

- `experiments/offline_score_common.py`: only NumPy metrics, serialization, metadata and known/unknown summation. Exact frozen R2 Otsu256; RMSE/12; public sigma/99% noise radius; unchanged physical prior and box. No model/algorithm import.
- `experiments/score_confirmation.py`: separate `seal` and `score` operations; immutable outputs. No inverse, generator, or physics import.
- `experiments/score_development_addenda.py`: separate write-once phase-pi/grid32 development scores, not a replacement of the original13 table.

Run from this worktree with the project's existing scientific Python. The scorer has **not** been executed for confirmation.

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r4/experiments/score_confirmation.py seal
```

This public-only operation expects downloads rooted at `results/confirmation_remote/`, with `results/confirmation/<task_id>.json` and `predictions/<observation_id>/<method>.npz`. It reads the frozen local public manifest/config/execution manifest and all18 paired outer/inner dispositions. A missing receipt or RUNNING disposition blocks closure. An explicit terminal failed/unknown disposition is retained in the fixed denominator, never silently omitted.

For missing killed/lost tasks, an optional `--dispositions /ABSOLUTE/PATH.json` may explicitly close them. Shape: `{"records":[{"task_id":"<frozen ID>","status":"UNKNOWN","closed":true,"reason":"documented terminal outcome"}]}`. Parent must supply the actual disposition; the scorer invents none. Failed/unknown output metrics remain unavailable, not zero. Partial NPZ metadata is retained without being treated as a scientific prediction.

The seal snapshots metadata for receipts/predictions/public observations, runtime scoring sources, protocol, shared generation receipt and available runner summary/affinity. It recomputes original Phi from saved m/Yprediction and public Y, with no GT or forward. Secondary selection is sealed here: strict minimum original Phi across the3 endpoints, with all exact ties retained, no added near-tie tolerance and no GT tiebreak. If an arm is unavailable, the unknown full-pool minimum and available endpoints are both retained.

Only after all18 terminal dispositions have been sealed and the parent explicitly authorizes scoring:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r4/experiments/score_confirmation.py score --authorize-confirmation-truth
```

The script verifies sealed path/size/mtime metadata, then opens only the local confirmation mapping/truth files. It writes `results/CONFIRMATION_SCORES.{json,csv,md}` once, prints each object and method/group counts, keeps all6 original objects per method and secondary denominator, and makes no significance claim. Existing outputs require a new versioned `--output` prefix. No hashes are computed. No scoring command was run during preparation.

## Required/optional download fields and remaining gaps

| Evidence / field | Requirement / disposition |
|---|---|
| All18 `<task_id>.json` | Currently download root absent. Need every outer transport receipt and inner `result.status`; outer OK alone is insufficient. Missing tasks require explicit closed failure/UNKNOWN records. |
| Successful NPZ m and Y_prediction | `(576,)` and `(48,12)`, finite. Missing or invalid saved predictions become UNKNOWN in the fixed denominator. Partial arrays are never called predictions. |
| Inner method/configuration/physical objective | Returned by the frozen adapter. Method identity and exact frozen configuration agreement are recorded; original Phi is independently recomputed from saved outputs/public data before GT. Missing receipt-Phi deltas are null. |
| Top-level `result.physical_model_cost` | Required for complete physical accounting. Use it once; never add nested setup/auxiliary/ranking/corrector/finalization physical subtotals again. Actual source count12; source/adjoint RHS remain separate from small-block RHS. |
| `result.auxiliary_operation_counts` | Preserve small4x4 factor/inverse RHS/matrix apply or onehop divisions as separate local work. An exception may lack completed counters; mark unavailable/incomplete, not zero. |
| `total_worker_s`, adapter wall/CPU | Recorded task wall includes public read and serialization; adapter wall is nested and must not be added again. Summed worker wall is work time, not concurrent elapsed. |
| `worker_affinity.json` | Download it. Establishes actual launched workers/pids/core masks. No complete worker import total is claimed without matching charged-once receipts for every worker. |
| `runner_summary.json` | Download it. Legacy `cuda_workers=2` is fixed runner metadata; actual manifest/task backends are all CPU, actual GPU tasks/workers0. |
| Shared generation/setup | Existing local `audit/confirmation/GENERATION_RECEIPT.json`:6 LU/72 source RHS, charged once across all18 tasks. Local model setup is already included in each task total. No historical operator assembly rerun. |
| Worker solver imports | Sum only marked `solver_import_charge.charged_once_per_worker` receipts once per worker. Do not charge the same shared worker import per task. |
| Earlier runner/bootstrap imports | Not timed by existing runner; unavailable/null. Do not infer them from adapter or task timer. |
| Overall remote elapsed | Not supplied by unchanged runner summary. Keep null unless parent supplies a separate measured execution receipt. |
| Peak memory | Adapter reports worker/process lifetime high-water mark, not an isolated per-job peak. Missing values remain null. |
| Currency/model cost | Not supplied by execution receipts; unknown. No zero-cost assertion. |
| Near-tie tolerance | Frozen secondary contract gives no numeric tolerance. This implementation adds none; strict original-Phi ties retain every endpoint. |

The JSON cost table contains each task, primary method task sums, each3-arm secondary pool and actual all18+shared-once campaign total. The secondary summary is the same already executed costed pool, not an extra execution to add. Unknown physical totals keep known subtotals and null full totals. Small6 results remain descriptive counts.

## Authorized development addenda already scored

After separate parent authorization and pre-GT seals:

- `results/DEVELOPMENT_PHASE_PI_SCORES_V1.{json,csv,md}`: phase-pi RMSE/12 .558703269, IoU .2578125, original joint FAIL.
- `results/DEVELOPMENT_GRID32_SCORES_V1.{json,csv,md}`: onehop RMSE/12 .053067749/IoU .935483871, global .053102210/.935483871; both original joint PASS.
- `audit/DEVELOPMENT_PHASE_PI_SEALED_V1.json` and `DEVELOPMENT_GRID32_SEALED_V1.json` precede their corresponding authorized development truth reads.

These are disclosed development controls, not fresh blind observations; confirmation GT remains unopened. Scoring adds0 LU/0 RHS. Final interpretation belongs to parent Codex.
