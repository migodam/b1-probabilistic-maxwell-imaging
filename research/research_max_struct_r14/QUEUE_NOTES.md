# R14 campaign preparation and finite queue

Implementation only. No public inputs were copied, manifests prepared, predictions launched, truth read, or scoring run by this implementation task.

## Preparation interface

`src/prepare_campaign.py --root <R14-study-root> --history-root <historical-R13-study-root>`

The history argument is the directory containing R13's `data_public/fresh48` and `configs/confirmation_core48.json`, not its repository root. Current local historical root:

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r13_worktree/research/research_max_struct_r13`

Current new study root:

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r14_worktree/research/research_max_struct_r14`

Preparation copies the48 original public NPZ files unchanged into `data_public/exposed48/`, preserving source file metadata with copy2. It reads the historical core48 JSON's g19 task observation contracts, without reading any private file or scientific array. It writes three immutable new manifests and `audit/CAMPAIGN_PREPARATION.json`:

|Manifest/campaign|Cases|Tasks|
|---|---:|---:|
|dualmetric48_v1|48|192|
|dualmetric_pilot12_v1|12|48|
|dualmetric_remaining36_v1|36|144|

Pilot object indices are fixed at2,5,7,16,17,24,29,37,38,42,44,45. The remaining manifest is the exact complement in ascending case order. All three share the same frozen protocol; task campaign directories differ. Run the main manifest **or** the pilot+remaining partition; running all three duplicates experiments.

Task IDs are `<r13-case-id>__pair00` etc; task method/algorithm are pair00/pair01/pair10/pair11; `pair` is the two-character string00/01/10/11. Explicit current_metric/material_metric fields carry M0/M1. Every arm uses global PR, cold original affine-current initialization, box[.02,20], original source normalizations, original material prior,120 whole-process CPU seconds,72 proposal seconds, two-CPU-second audits and common original corrector.

The stage schedule is mf3/6/10/12, maxiter500/500/500/1500. Every stage includes b0=1 and b1=1/.2/.2/.2 (also beta_current=1 and historical beta_material=b1). Per the R14 contract, the runner must use b0=1 for an M0 material readout in00/10; historical beta_material must not override it. M1 readouts in01/11 carry the exact R7 final Arc options: mobius, arc_candidates=True, inherited_gtol=True, sweeps2, step_tol1e-8. These options appear in protocol/readout_options and task/readout_options. The runner owns scientific interpretation of the fields.

Preparation refuses any existing target public directory/manifests/receipt. A partially failed copy leaves a partial new output, which is not automatically retried or overwritten. Historical source_commit is explicitly renamed historical_source_commit in the new protocol; the new implementation is not mislabeled as the old commit.

## Queue interface and scheduling

`src/run_queue.py --root <R14-study-root> --legacy-repo <inherited-repository> --manifest <one-prepared-manifest> [--cores 0,2,4,6,8,10,12,14]`

The queue runs in the foreground until its finite task set completes. It groups arbitrary manifest tasks by case_id, retaining case first-occurrence order and each case's task order even if input tasks were noncontiguous. A lane takes one entire case. Its four pair tasks run sequentially on the same logical CPU; after the fourth finishes, the lane takes the next pending case immediately without waiting for other lanes. An arm failure does not trigger a retry and does not move the remainder of that case to another core. An individual launch failure is recorded and the finite queue continues.

Default eight lanes are logical IDs0,2,4,6,8,10,12,14. These are the frozen Windows target's declared P-core choices, not a portable topology detector or a claim of exclusive cores. A later explicit E-core diagnosis can pass another distinct nonnegative core list; actual IDs are recorded. The queue validates all requested IDs against its own allowed affinity and rejects a host without process affinity (including the default macOS psutil interface). No process outside this queue is modified.

The queue sets six common BLAS/vector-library thread limits to1. Each child is pinned and given Windows NORMAL priority before importing the numerical runner. Its command interface is the parent's `src/run_experiment.py --root --legacy-repo --manifest --task-id --cpu`. PYTHONPATH contains the declared legacy repository and its vendor directory. A read-only resource probe records other Python processes; it does not stop, reprioritize or inspect their data. CPU load remains a shared-host observation.

A new `queues/<campaign>` directory is required. The queue freezes the manifest there and gives children that immutable copy. It refuses any existing expected RESULTS/task directory before starting. Logs and individual TERMINAL_JOBS receipts are new files. RUNNING statuses use monotonically numbered STATUS_*.json files; only the final terminal status uses RECEIPT.json, all opened exclusively. EVENTS.jsonl is the sole append-only event stream. There is no mutable status-alias rename, retry daemon, GT step, service, detached launcher, or automatic scientific continuation.

## Cost and error accounting

Windows whole-child lifetime CPU is captured by the actual inherited `research/research_max_struct_r11/experiments/twod_queue.py:external_cpu` (`GetProcessTimes`, user+kernel), after child termination while its process handle remains available. It includes interpreter startup/imports and final serialization. A missing/error capture remains null and makes a Windows task fail queue completion; it is never replaced by an internal runner timer. On non-Windows affinity-capable hosts the inherited helper reports UNKNOWN_EXTERNAL_NONWINDOWS/null, so actual whole-child CPU remains unavailable.

Per-task external wall time, queue wall time, controller CPU, external child CPU and nested child internal CPU have separate fields. Never add wall time to process CPU or add child internal timing to external child timing. The aggregate finished_child_CPU_seconds is null if any measurement is missing; the separately labeled measured subtotal and missing task IDs remain available. Controller CPU measured from dispatch setup is a disjoint queue-controller interval, not the runner's120-second scientific budget.

A nonzero child exit, missing/malformed child receipt, child FAILED status, or Windows CPU capture gap is reported as failure; other finite tasks still proceed. There are no retries. A launch failure records no-child CPU0 explicitly. Interrupted controllers are not resumable through this script: existing immutable queue/output paths prevent silent relaunch, and another process is never killed for cleanup.

## Verification and limits

Both source files compile in memory. A source-only helper check verified noncontiguous manifest tasks are grouped into the correct ordered four-arm case groups, duplicate task IDs are rejected, and invalid/duplicate/negative CPU lists are rejected. No campaign was prepared and no subprocess/solver/score was executed. Windows affinity, priority, whole-child CPU capture, actual parent runner schema and end-to-end execution remain untested here. Source provenance is the inspected R13 finite_queue and R11 external_cpu implementation; no SHA checks were performed.
