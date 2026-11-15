# R15 finite review export — implemented, not run

Entry point: `src/export_review.py`. Implementation-only validation performed AST parsing, bytecode compilation and prohibited-call inspection. No module import/execution, source arrays, ZIP generation, physics or GT access occurred. Static receipt: `audit/REVIEW_EXPORT_STATIC_CHECK.json`.

After the actual main264 **and** fork56 queues are terminal:

```powershell
python "$r15Root/src/export_review.py" --root "$r15Root"
```

The script depends on the actual R15-local `src/package_raw.py` terminal validation/boundary helpers and `src/queue_interventions.py` V2 `source_checkpoint`/`complete_checkpoint` APIs. Those files must be deployed alongside it. It invokes neither queue main program. The shared terminal validator checks actual `FROZEN_MANIFEST.json` for main264 and `FROZEN_TASKS.json` for fork56, exact terminal task IDs, per-task final records/logs, and permits real failed/missing-result tasks with retained records.

Outputs, all exclusive-create without overwriting an existing attempt:

- `archives/R15_REVIEW_EXPORT_V1.zip`
- `archives/R15_REVIEW_EXPORT_V1_MANIFEST.json`
- `archives/R15_REVIEW_EXPORT_V1_FILES.csv`
- `archives/R15_REVIEW_EXPORT_V1_README.md`
- `archives/R15_REVIEW_EXPORT_V1_RECEIPT.json` on completion; `..._ERROR.json` on a failed partial ZIP attempt.

This is expressly **a reduced review copy, not a full raw release or a startup/source package**. All archived file paths remain original R15-root-relative paths. It does not delete or rewrite any source. `package_raw.py` supplies the complete original j/Q, histories/checkpoints/caches and RUNNING status snapshots omitted here. Public inputs, source/runtime and rank-basis setup are not copied because the review machine already has them.

Main264: copies actual small JSON/JSONL/CSV/log artifacts, including RECEIPT/HISTORY/CYCLES/ACTUAL_CONFIG/TWOFOLD_SETUP/WRAPPER_PROVENANCE and actual other emitted metadata. Original `CANDIDATES.npz` and all five view JSON files are copied byte-for-byte. Each existing five-view NPZ is transformed to **only its actual `m` and `Y_prediction`**, preserving original shape/dtype and real-material m. No other field is loaded to perform this transform, and no missing Y is fabricated. Readback opens the actual outer-ZIP member and performs `np.array_equal` plus exact shape/dtype and retained-key checks against original arrays. Failure/missing slots remain explicit rather than becoming material-only fallbacks.

Common-hop family = twofold_hop / twofold_dr / twofold_ap / twofold_full_dr / twofold_reset. For all24 cases it retains first4 original CYCLE BEFORE/AFTER full NPZs. The correct denominator is **24×5×4=480 pairs, hence960 expected NPZ files**. Every existing file is copied, and every missing slot is retained in the manifest. Main NATIVE_HISTORY is retained only for0009/0020/0021/0023 ×these5 methods; the remaining long native histories stay in raw assets.

For the same4 hop tasks, V2 prefers complete `CHECKPOINT_00512.npz`, otherwise latest complete numeric checkpoint<=512, otherwise earliest complete numeric checkpoint>512. It excludes CYCLE/STAGE/FINAL as source choices and uses metadata completeness only, with actual phase/beta recorded. The returned original checkpoint, `FIXED_OPERATOR_CACHE.npz` and any real fixed/background dependencies are copied. This source selection is a readout-analysis input choice, not a new prediction or modification of the source recorded by fork56. Other full cycle snapshots may require raw-release caches for solver restoration; they can be read as full saved-state evidence with the existing codec.

Fork56: every admitted original file under RESULTS/fullstate_interventions_v1 is copied, including STATE/FORK/CANDIDATES/HISTORY/NATIVE/CYCLES/receipts/cache files. Both queue exports contain final RECEIPT, actual FROZEN file, MAIN_TERMINAL_SNAPSHOT where emitted, RESOURCE/EVENTS, TERMINAL_JOBS and logs; all RUNNING STATUS snapshots are omitted. No oldendpoint, other audit/analysis, truth, credential/private, root startup ZIP or raw-release archive tree is traversed.

Text uses DEFLATE1; NPZ and other binaries use outer ZIP_STORED. Original copied bytes remain untouched; main view transformations and omitted original fields are recorded per file alongside original filename/bytes, output bytes, roles and array checks. ZIP member paths/sizes and original size/mtime are verified, without SHA. Actual archive bytes, source/admitted bytes, exported bytes, missing/failure lists, IO/readback CPU and wall time are recorded. Existing inverse execution CPU remains a separate queue provenance snapshot, neither charged again nor refunded. No physics/optimization/GT/selection operation exists, and no scientific acceptance judgment is made.
