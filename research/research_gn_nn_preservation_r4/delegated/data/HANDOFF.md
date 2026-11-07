# Data and GN preparation handoff

Generation is complete: 690 independent mothers, 1380 public observations; TRAIN/VAL/TEST mothers 504/78/108. TEST truth is sealed and has never been scored. `audit/data_contract.json` and `data/GENERATION_RECEIPT.json` record contracts and actual generation costs.

The strict primary GN queue is running through a live local SSH process (`exec` session 54327). Its remote root is `D:/AI/GN_NN_R4_20261107`; only this root is writable. Four workers use P cores 4,6,12,14 (affinity mask 20560), BLAS=1. Historical mask85 was avoided because unrelated session tasks occupy cores0/2. RTX4060 was occupied by another session and is not allocated to R4.

Local fetch monitor (`exec` session 9338, process `delegated/data/watch_preparation.py`) downloads complete NPZ+receipt pairs every ~120s. It writes `audit/data_progress.json`. All primary states must have `config.discrepancy_stop=false`, strict 1e-12 tolerances, 45 seconds/500 LU/30000 RHS budgets and complete 8→12→24 parameter release. Primary filenames are `anchors/<obs_id>.npz`.

The monitor publishes `audit/DATA_READY.json` only after all 1008 TRAIN and 156 VAL primary anchors exist and all are full24; training must use those exact obs IDs from the manifest, never a wildcard over all suffixes. The network worker knows this release interface. No marker is published for a partial dataset.

After all1380 primary anchors are fetched, the monitor starts an independent four-worker VAL/TEST discrepancy and weakreg1e-5 queue. Primary strict anchors already are the fullTRF control; no duplicate fullTRF run. Discrepancy thresholds are1e-6 clean/.0011 noisy, checked only after full24 analyticJ. Suffixed controls are `_discrepancy` and `_weakreg1e-5`, each with its own frozen45s/500LU/30000RHS budget. `audit/data_baseline_progress.json` must reach744/744. After final fetch the monitor exits; CIE CPU work can then use these four cores.

The initial17 TRAIN discrepancy attempts were mistakenly placed on primary paths, corrected before training, preserved with `_discrepancy` suffix, and fully charged in `audit/ANCHOR_PROTOCOL_CORRECTION.json`. Up to4 interrupted in-flight attempts have unknown costs, not zero. TRAIN seen_image/seen_parent flags were corrected True; VAL/TEST remain False. The remote package contains no truth or evaluator manifest.

Manual recovery (study working directory and established Python):

```
python src/remote_jobs.py status
python src/remote_jobs.py fetch --split train
python src/remote_jobs.py fetch --split val
python src/remote_jobs.py fetch --split test
```

If the live SSH queue fails, first inspect `delegated/data/remote_anchor_strict_session.log`, remote `anchors.log` and receipts, then rerun `python src/remote_jobs.py run` to resume missing primary IDs only. The primary runner skips complete NPZ+receipt pairs. Baselines may be resumed using `python src/remote_jobs.py run --kind baselines`. Never use detached Windows Start-Process: Windows OpenSSH tears down its descendants when the SSH session closes. Preserve partial evidence and do not expand scientific budgets to rescue a case.

To update remote source for the later CIE adapter, `python src/remote_jobs.py upload-code` only uploads source, not observations or GT. Additional CIE runner invocation is owned by the parent. No TEST GT metrics or plots have been produced by this worker.

## Authorized M4 acceleration

A second finite controller now handles all156 VAL observations followed by the final256 TRAIN observations in fixed online-manifest order. The exact public IDs are `delegated/data/local_val_ids.json` and `local_train_ids.json`; allocation and first legitimate throughput measurement are in `audit/LOCAL_PREPARATION_SPLIT.json`. TRAIN remainder for i9 is752 (1008−256), and all216 TEST observations remain i9 work. Two M4 BLAS1 CPU workers retain the strict per-case45s/500LU/30000RHS budgets and full24 release.

Local controller session58619 (`delegated/data/run_local_preparation.py`) uploads complete local NPZ+receipt pairs every30s. Each destination is created exclusively, with JSON committed last; an already complete remote canonical pair wins and is never overwritten. Progress is `audit/LOCAL_PREPARATION_PROGRESS.json` and `LOCAL_PREPARATION_UPLOADS.json`. The first VAL case `obs_bb500efe50f35ed6a2c5` used2.527s solve/100LU/5376RHS. Subsequent noisy examples can use the full45s, so this is a measured case, not a universal speed claim.

Fresh remote worker subprocesses now recheck completed pairs, returning `REUSED_COMPLETE_PAIR` with zero new physics. If a reserved local ID has not arrived, they return `DEFERRED_TO_LOCAL` with zero new physics and let the FIFO advance. Public-only reservation file `anchors/LOCAL_RESERVED_IDS.json` was installed remotely before starting the M4 batch. It remains evidence of ownership; once all412 are completed/uploaded, its IDs may be cleared if an explicitly authorized missing-case recovery is later needed. No existing remote physics was stopped.

`audit/data_canonical_costs.json` is refreshed by each new fetch. It sums each primary canonical observation receipt once by actual M4/i9 origin, separates misconfigured attempts, and records reuse/defer actions from downloaded completed batch receipts. The old already-running remote coordinator cannot hot-reload its aggregate logic: its legacy aggregate may include imported M4 case costs, so never add it to M4 batch costs. Use the canonical ledger and preserve actual action counts; partial batch action coverage is explicitly marked until remote batches finish.
