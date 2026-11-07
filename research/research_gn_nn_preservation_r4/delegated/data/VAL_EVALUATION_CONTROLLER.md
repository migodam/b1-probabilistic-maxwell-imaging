# Finite VAL public-correction controller

`watch_val_evaluation.py` is a finite script for this R4 run. It reads only the public online manifest, prediction/evaluation receipts, and preparation/request status. It never reads truth or the evaluator manifest and never invokes scoring.

The initial owner-controlled request is `audit/VAL_EVALUATION_REQUEST.json`: the five explicitly registered seed7 cascade2 runs, `closed=false`. The parent owns future edits to `runs`/`closed`; the controller never writes the request. Write updates atomically where possible. Only a run explicitly registered at launch is eligible, and all156 VAL prediction NPZ/JSON pairs must exist with valid JSON.

It waits until the existing M4 controller reports `LOCAL_PHYSICS_COMPLETE`, `LOCAL_COMPLETE_UPLOADED`, or `LOCAL_COMPLETE_UPLOAD_INCOMPLETE`, with no active preparation process PID. It then starts exactly two existing `experiments/run_evaluation.py` subprocesses per run, VAL-only with fixed partition count2 and indices0/1. BLAS/OMP thread limits are1. Existing completed evaluation items are skipped by the existing runner.

After all156 evaluation NPZ/JSON pairs complete, it sums canonical per-case charged costs once, records process returns/wall time, then rechecks the parent's latest request. Case charged costs may include proposal preparation per the existing evaluation contract; do not add that preparation again. `audit/VAL_EVALUATION_CONTROLLER.json` and `delegated/data/watch_val_evaluation.log` record status. Per-partition logs are under `delegated/data/val_evaluation_<run>_P<index>_<timestamp>.log`.

It exits when `closed=true` and every currently requested run is complete, or at the fixed absolute deadline2026-11-07T10:44:27Z. At that deadline it stops only its own two process groups, preserves completed artifacts, and marks in-flight costs unknown rather than zero. It never stops preparation, training, i9, or baseline tasks. A failed/incomplete run causes a HOLD and exit, with no automatic scientific retry. The parent can inspect and explicitly resume the same script; existing completed cases are retained.

Started local exec session78255. Initial observed state: `WAIT_M4_PREPARATION_END`; no numerical subprocesses launched.
