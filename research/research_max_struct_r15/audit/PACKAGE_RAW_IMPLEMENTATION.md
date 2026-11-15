# R15 raw archive implementation — not executed

Implemented `src/package_raw.py` after reading the actual R15 queue sources. This task performed AST parsing and bytecode compilation only: the module was not imported/executed, no queue/result data was opened, and no ZIP or physics/GT operation ran. Static receipt: `audit/PACKAGE_RAW_STATIC_CHECK.json`.

Run on the actual remote R15 root only after both queues have terminated. Here `$r15Root` is the existing R15 deployment directory, not the R14 root or an extracted startup ZIP:

```powershell
python "$r15Root/src/package_raw.py" plan --root "$r15Root"
python "$r15Root/src/package_raw.py" create --root "$r15Root"
```

`plan` writes only manifests to a new `archives/r15_raw_plan_v1`. `create` independently verifies the current frozen sources and writes actual ZIPs to a new `archives/r15_raw_v1`. Neither directory can already exist. `--output` accepts a different **new directory beneath the same R15 root/archives** for an explicitly recorded later attempt; failed attempts and source files remain untouched.

The three separately named groups are:

- `main264_raw`: all admitted `RESULTS/twofold_exposed24_v1` and `queues/twofold_exposed24_v1` files. Exact source seam: `queue_independent.py` writes `FROZEN_MANIFEST.json`, `TERMINAL_JOBS/<task_id>.json`, `logs/<task_id>.log`, terminal `RECEIPT.json` and intermediate `STATUS_*` files.
- `fork56_raw`: all admitted `RESULTS/fullstate_interventions_v1` and `queues/fullstate_interventions_v1` files. Exact source seam: `queue_interventions.py` writes **`FROZEN_TASKS.json`**, `MAIN_TERMINAL_SNAPSHOT.json`, the same per-task records/logs and terminal receipt. Failed/preflight/launch tasks may legitimately lack result outputs; the manifest records those missing artifacts and preserves their original terminal records/logs instead of fabricating result files.
- `source_public`: `src`, `runtime`, `configs`, `data_public`, the complete original `audit/rank_all24` real/complex basis caches and metadata, plus the exact top-level files `EXPERIMENT_CONTRACT.md`, `THEORY_DERIVATIONS.md`, `ORIGINAL_REQUEST_ZH.txt` and optional `.gitignore`. Root remote-input ZIPs are outside this whitelist. No other audit tree, analysis, local old-endpoint GT-derived arrays, release tree, credentials or private operations are traversed/admitted.

Both frozen task inventories must be exact (264 and 56), final receipts must be `COMPLETED` or `COMPLETED_WITH_ERRORS` with no pending/running work, and every frozen task must have its matching terminal record and log. All other admitted files—including states, candidates, caches, failures and intermediate queue statuses—remain in their original root-relative paths. Source membership, bytes and mtime are rechecked after actual archival. Scientific successes are not required; failed jobs remain included.

Volumes target **600,000,000 bytes**, with a conservative bound on raw DEFLATE growth plus ZIP/name/Zip64 overhead. A singleton exceeding the target may be admitted only within the **650,000,000-byte hard maximum**. Every actual ZIP size is checked after closure. An individual source file whose conservative bound exceeds650MB causes an explicit recorded failure; no implicit fragmentation, truncation, or numerical reduction is implemented. This limitation has not been exercised against the remote terminal data.

NPZ and other binary files use `ZIP_STORED`; JSON/JSONL/CSV/log and other text use lossless `ZIP_DEFLATED` level1. NPZ files are never loaded. Generated archives are exclusive-create; no source deletion, SHA, publication, upload, scorer, physics or optimization call exists. Forbidden credential/private path components and symlinks are excluded before opening their contents, with exclusions recorded. A conservative literal scanner blocks possible credential-bearing admitted text without logging matched contents.

Each invocation writes `PACKAGE_PLAN.json`, three `*_MANIFEST.json` files, `FILE_INVENTORY.csv`, and `RESTORE.md`. Final **`PACKAGE_RECEIPT.json`** distinguishes `PLANNED_NO_ARCHIVES` from `ARCHIVED_ORIGINAL_BYTES_NOT_PUBLISHED` and records actual part sizes, files, per-part/process CPU and wall time. Failure writes **`PACKAGE_ERROR.json`**, retaining completed/partial archives and their measured cost. Source queue CPU totals are retained as nested provenance snapshots, not added to archive CPU. Normal ZIP directory/member-size/compression verification is performed without SHA or numerical transformations.
