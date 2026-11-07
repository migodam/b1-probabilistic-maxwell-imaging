# Fixed TEST host partition and finite execution

Implementation only. No upload, numerical execution, physical check, scoring, or truth/evaluator access was performed for this task. Four Python entry files and their embedded remote Python blocks passed syntax compilation. PowerShell transport and the new orchestration have not been exercised end to end.

## Public allocation and comparison scope

`audit/TEST_RESOURCE_PARTITION.json` was created from the five-field online manifest only. Enumerate TEST parent IDs by first occurrence: parent index modulo3 equal0 goes to M4, all others to i9. It contains no class, family, font, support, truth, or material metadata.

- `local`:36 mothers,72 observations; at most2 BLAS1 numeric children.
- `remote`:72 mothers,144 observations; at most4 BLAS1 numeric children, inherited own P-core mask20560 (cores4/6/12/14).
- Both clean/noisy observations of each mother remain on one host. Every method, including cold CIE, is postprocessed on that same assigned host.
- `global_online_indices` maps all216 TEST observations to their indices in the complete original online manifest. NN/continue/independent method order rotates using that index. A single-observation subprocess therefore retains its real index. Existing VAL rotation/calls retain their prior filtered-index behavior.

This is a resource allocation. Hardware differences and shorter queue wall time are not algorithmic advantages. A completed execution partition is coverage, not physical or GT success. Actual per-case FAIL and budget-stop statuses remain in their receipts.

## Required public frozen files

All three are bundled by source upload:

1. `audit/METHOD_FREEZE.json`: `test_runs`, `nn_runs`, `baselines`, `test_observations=216`, `test_mothers=108`.
2. `audit/FINAL_RUNS.json`: `plans`, `runs`, `nn_runs`. `runs` must exactly match plans in order, and match the freeze. Each plan has a unique `run`; `kind` is exactly `NN`, `continue_trf`, `independent_baseline`, or `cie`. CIE plans additionally supply `method`, for example `{ "run":"cie_dual_cold55", "kind":"cie", "method":"cie_dual" }`. Checkpoints may remain relative to the study; weights are not uploaded to the CPU host.
3. `audit/TEST_RESOURCE_PARTITION.json`: declared `local`/`remote` parent and observation ID lists plus original global indices.

There are no alternate-schema guesses. The parent writes the first two after VAL selection. Missing or inconsistent freezes hold execution rather than selecting another method or sample. Full CIE all108 scope still has216 observation-method combinations for the one frozen CIE-dual method. Other-host and unfinished combinations remain explicit `NOT_RUN`; they are never dropped from the global plan.

## Actual entry points, for later parent execution

Use the established study Python. Commands below are documented, not executed by this implementation task.

```
python src/remote_jobs.py upload-code
python src/remote_jobs.py upload-predictions
```

Source upload includes vendor Python/provenance and the three available public frozen files, without truth, evaluator manifest, checkpoints, or weights. Existing remote frozen JSON must agree before replacement is allowed. `upload-predictions` reads `FINAL_RUNS.nn_runs` and transfers only each registered NN run's144 remote-assigned TEST NPZ/JSON pairs. It requires all expected pairs, preserves existing canonical remote pairs, and holds partial destinations. It uploads no local-assigned predictions and no weights.

Wait for the preexisting744 final-only discrepancy/weakreg baseline pairs and release of the same four own i9 P cores. New remote queue entry points check both conditions and do not stop another session. If the separately declared all-scales control is registered, run it independently before its verification inputs are needed:

```
python src/remote_jobs.py run --kind discrepancy-all-scales
python src/remote_jobs.py fetch --split val
python src/remote_jobs.py fetch --split test
```

That new queue calls the existing all-scales baseline for all372 VAL/TEST observations with its unchanged45s/500LU/30000RHS case budget. It neither substitutes for nor erases the old744 computations. Existing anchor fetch remains verbatim. Original anchors/baselines command branches retain their original arguments.

The two host entry points are:

```
python experiments/run_test_partition.py --execution-partition local --workers 2 --execute
python src/remote_jobs.py run --kind test
```

The remote command invokes `run_test_partition.py --execution-partition remote --workers 4 --execute`. Each host first runs its frozen CIE subset through `run_cie_reference.py --scope all108 --execution-partition <host>`, then starts a finite ready-observation FIFO. When an observation finishes, the next starts immediately. One observation subprocess calls the existing `run_evaluation.py` for all frozen non-CIE methods. No daemon, service, new physics solver, or optimizer is introduced.

`run_evaluation.py` now supports partition counts1/2/4, records platform/machine and global method-order index, and maps `discrepancy_allscales45_anchor` to `_discrepancy_allscales`. Its public TEST calls require the declared host and registered methods. VAL commands keep their existing arguments and solver paths. NN correction, continuation, independent verification and cold CIE retain their existing10s/3s8nfev/55s and LU/RHS contracts; no scientific budget is expanded.

Fetch completed remote evaluation pairs and queue receipts with:

```
python src/remote_jobs.py fetch-evaluation
```

Only frozen run names and remote-assigned TEST IDs are fetched. Complete canonical pairs are preserved; partial existing pairs and missing/mismatched assigned-host provenance are held. Queue receipts are also copied into timestamped immutable snapshots, so earlier progress/error versions remain available. If a historical TEST case lacks the new assigned-host provenance, execution preserves it and reports a HOLD; it does not infer hardware or overwrite the old case.

## Receipts, time and interruption costs

- Per-case scientific costs and actual statuses remain under `evaluation/<run>/<obs_id>.json` and `.npz`, once per canonical method/observation.
- `delegated/data/test_queues/PARTITION_<host>_<stamp>.json`: whole partition setup, phase/queue wall time, child outcomes, canonical costs and all planned missing cases.
- `delegated/data/test_queues/POST_FIFO_*.json`: refill-on-completion FIFO history; per-observation logs are adjacent.
- `delegated/cie/PLAN_*.json`, `QUEUE_*.json`, `SUMMARY_*.json`: full global CIE denominator, assigned subset, explicit NOT_RUN cases, queue time and costs.
- `audit/REMOTE_TRANSFER_*.json`: source/prediction packaging, transfer/import/fetch and foreground queue wall times, exceptions and retained interruptions.
- `audit/REMOTE_evaluation_<stamp>.zip` and `REMOTE_EVALUATION_SNAPSHOTS/<stamp>/`: immutable transferred versions of public outputs/queue receipts.

Saved case charged totals can include proposal preparation or independent baseline construction under the existing evaluation contract. Do not add that preparation again, and do not sum repeated canonical/batch summaries as fresh physics. Reused completed cases incur no new physics; failed or interrupted in-flight work is `UNKNOWN_NOT_ZERO`. No automatic scientific retries occur in the new TEST FIFO.

The absolute deadline is2026-11-07T10:44:27Z. New queues stop dispatching and terminate only their own tracked process groups/trees at that deadline. Other preparation, training, sessions and GPU work are not stopped. Unstarted methods/cases remain NOT_RUN; completed artifacts, failures, raw logs and unknown interruption costs are preserved. A lost SSH connection is not treated as proof that remote work stopped, and it is never automatically retried.
