# Reproduce the frozen two-dimensional pipeline

Run from a complete repository checkout with the existing NumPy/SciPy/psutil/threadpoolctl environment. Keep the inherited `research_max_struct_r2/r4/r5/r6/r7/r8`, `research_fbg_r3_experiments`, `research_z_alg_r4_diversity`, and `research_z_alg_r4_basin` Python modules available. No new backend is installed. On this Mac use the existing `Gaussian/.venv_nn/bin/python` symlink without resolving it to the base interpreter. On XINAN use `D:/python/python.exe`. The frozen scientific contract is [CONTRACT.md](CONTRACT.md).

Generate once into fresh output directories (refuses existing public/truth/receipt outputs):

```sh
python research/research_max_struct_r10/experiments/twod_generate.py --execute
```

Generator-only inputs are `configs/twod_generator.json` and the inherited public operator. Generated inverse input consists of `data_public/twod/manifest.json`, anonymous observations, `operators.npz`, `configs/twod_tasks.json` and `configs/twod_solver.json`. Do not transfer `sealed_truth/twod`, generator metadata or scoring labels to inverse workers. No scoring command is provided or authorized yet.

Run a single complete public task, keeping its existing mature allowance:

```sh
python research/research_max_struct_r10/experiments/twod_run.py --task-id r10_twod_0000__cie_single__noise_ball_99
```

The finite Windows queue validates the complete manifest and fresh outputs, samples live CPU load, saves resource/launch receipts, and runs up to four single-thread BELOW_NORMAL children. A preparation invocation performs no inversion; give it a separate fresh receipt directory:

```powershell
D:/python/python.exe research/research_max_struct_r10/experiments/twod_queue.py --receipt-dir research/research_max_struct_r10/delegated/twod/queue_plan --vendor-path vendor
D:/python/python.exe research/research_max_struct_r10/experiments/twod_queue.py --execute --receipt-dir research/research_max_struct_r10/delegated/twod/queue --vendor-path vendor --cpus 0 2 4 6
```

`vendor/threadpoolctl.py`, when used, is a copy of the already installed pure-Python module, not an installation. Keep SSH connected for the finite foreground Windows job. Monitor `delegated/twod/queue/RECEIPT.json` and `EVENTS.jsonl`; task logs are `queue/logs/<task_id>.log`. To stop dispatching, create this task's own `queue/STOP_DISPATCH` file: active children finish, and unstarted tasks are recorded. It does not terminate another process. Existing scientific outputs are never overwritten or implicitly reused.

For another host, use a fresh repository location/output tree and explicit authorized CPU numbers. The source public manifests use repository-relative paths and are portable; per-run frozen adapter/receipt files deliberately retain actual provenance paths. Hardware adaptations do not alter physics, precision, frozen scientific allowances or candidate selection. No claim of bitwise identical optimization trajectories across BLAS platforms is made.

The actual inverse-only packaging command is `python research/research_max_struct_r10/experiments/twod_package.py`; it excludes fresh generator code/configuration and truth. Use `--output-dir <fresh_owned_path> --manifest research/research_max_struct_r10/configs/denominator_confirmation_tasks.json --extra-file research/research_max_struct_r10/src/lowrank_denominator.py` for the separately root-frozen candidate bundle. Packaging uses only existing repository modules and copies installed threadpoolctl3.7.0 into vendor. Transfer/unpack commands and their return codes are retained in each transport/TRANSPORT_RECEIPT.json.

Actual baseline remote root: `D:/AI/MAX_STRUCT_R10_2D_20261107`. Candidate root: `D:/AI/MAX_STRUCT_R10_DENOM_20261107`. The candidate queue passes its explicit manifest and `--cpus 8 10`; its results must stay separate by `material_readout.label`, despite canonical dispatch arm=a_plus. Existing baseline files/processes are not changed by the candidate bundle.

Local foreground launch/status commands, after a fresh bundle is transferred/unpacked and its preparation receipt is saved:

```sh
python research/research_max_struct_r10/experiments/twod_remote.py launch
python research/research_max_struct_r10/experiments/twod_remote.py status
python research/research_max_struct_r10/experiments/twod_remote.py launch --tag DENOMINATOR --remote-root D:/AI/MAX_STRUCT_R10_DENOM_20261107 --manifest-name denominator_confirmation_tasks.json --cpus 8 10
python research/research_max_struct_r10/experiments/twod_remote.py status --tag DENOMINATOR --remote-root D:/AI/MAX_STRUCT_R10_DENOM_20261107
python research/research_max_struct_r10/experiments/twod_fetch.py --cohort baseline
python research/research_max_struct_r10/experiments/twod_fetch.py --cohort denominator
python research/research_max_struct_r10/experiments/twod_costs.py
python research/research_max_struct_r10/experiments/twod_costs.py --cohort denominator
```

Live local state is REMOTE_ACTIVE/REMOTE_QUEUE_RECEIPT or DENOMINATOR_ACTIVE/DENOMINATOR_QUEUE_RECEIPT under delegated/twod. Cost compilation reads JSON only and never scores arrays; generation is additive only in the baseline cohort, reused without a second charge in the candidate cohort. Unknown earlier metadata-query/transfer CPU is explicitly unmeasured, not zero. The two foreground SSH invocations each retain a finite process until their queue ends.

Fetch is permitted only after its queue is terminal. It archives only RESULTS/twod and queue metadata, records remote archive CPU and local controller/SSH CPU, rejects existing local scientific files, and never copies sealed truth or invokes scoring. Separate archive paths map candidate and baseline queue metadata to fetched_denominator_queue/fetched_baseline_queue. Preserve a failed fetch receipt rather than silently retrying or deleting it.
