# Authorized R9 finite queue, attempt02

36 frozen tasks:16 baseline noise-ball,12 landscape continuation and8 historical fixed-.001 controls. Maximum4 simultaneous local children, each1 BLAS thread. The first4 jobs use distinct case/method pairs. No automatic retry, scoring or GT read is performed. All task commands, frozen manifest copies, start/end events, output paths and errors are in PLAN.json, EVENTS.jsonl and RECEIPT.json; stdout/stderr is retained per task under logs/.

Reproduce from R9 repository root with a fresh receipt directory and fresh scientific outputs:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" \
 research/research_max_struct_r9/experiments/run_queue.py \
 --manifests research/research_max_struct_r9/configs/baseline_development.json \
 research/research_max_struct_r9/configs/landscape_development.json \
 research/research_max_struct_r9/configs/baseline_historical_fixed.json \
 --workers 4 --receipt-dir research/research_max_struct_r9/delegated/PARENT_FRESH_QUEUE
```

The executable path intentionally preserves the virtual-environment symlink entry rather than resolving it to base Python. Existing result directories are never overwritten; they are reported as non-launched errors. Parent signals stop dispatch only and never terminate children or unrelated processes.

Attempt01 mistakenly resolved that symlink, causing28 import failures with zero physics. Its complete queue logs remain at ../development_queue/ and12 failed landscape outputs were moved into audit/queue_attempt_01_import_failure/, with an archive map retained. A separate shell-launch typo also failed before any queue process and its log is retained. These are engineering attempts, not new scientific task denominators. Attempt02 changes no scientific algorithm, budgets or frozen input.
