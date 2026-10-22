# Memory interruption recovery review

Read-only source/receipt audit, 2026-10-22. No physics, truth files, training, network, or large arrays opened. No journal scan performed. Source findings describe the version inspected; subsequent parent edits require re-review.

**Verdict: reservations and interruption receipts are retained; recovery implementation needs corrections before an equal-cap continuation claim.** The bounded native cache fixes a concrete historical strong-reference retention mechanism. It does not establish a process-wide RSS bound or enforce one paid worker.

## Verified from source and small receipts

- `src/p3l/core.py:41–53,156–157`: campaign adapter transfer cache is an `OrderedDict` capped at four entries. Eviction drops only the adapter reference; live `Point` objects retain their own state/T/J. Upstream `receiver_transfer` uses ordinary membership/read access, so this is insertion-order eviction rather than read-promoted LRU. This distinction does not defeat capacity four.
- `src/p3l/inversion.py:80–115`: the optimizer also owns a per-invocation dictionary of evaluated Points. Its state-count cap bounds that lifetime; transfer-cache eviction cannot free states still owned by that dictionary. The optimizer dictionaries become collectible after return. Context runner clears adapter transfers and calls GC after each parent. This limits cross-parent accumulation; no runtime peak-memory acceptance was performed here.
- `Book.reserve_state` precharges eight RHS, one state attempt and one LU; callback credits prevent double charging. Receiver solves reserve their own RHS. Recovery invokes `_update` with controlled wall and no negative RHS/states/LU adjustments. Failed and incomplete endpoints remain paid.
- Six interruption closures match every corresponding recovery-receipt field. Totals: **109,488 paid RHS; 3,408 state/LU attempts; 3,405 completed endpoint files; 3,408 attempt files**. Three missing endpoints remain charged. Recorded physical wall is **965.848711 s**. Reconciled controlled wall is **3,623.999538 s**, versus **1,900.077719 s** observed lower bound. Receipts explicitly disclose possible tail overcharge and zero new recovery physics.
- The small ledger snapshot inspected has `active_wall: {}`. This is a snapshot observation, not proof of OS process exit or independent reconstruction of the full journal.
- Existing completed contexts keep their NPZ and source-job metadata; recovery audits the saved cache and checks q/prediction agreement. Missing context caches incur a paid forward/J rebuild. Incomplete optimizer runs restart coarse TRF from the saved best q; optimizer trust-region state is not resumed. Previous paid attempts are counted from that parent's trace.
- Online context/recovery paths read online observations, paid caches and public metadata; no truth import/access found in those paths. The packing tool deliberately accesses only train/validation latent q after all public contexts preflight. It does not open locked/development truth or import a physical backend/trainer. Its output guard refuses any existing final artifact. Closure receipts are explicitly not OS-liveness proof.

## Actionable defects and limits

1. **Per-parent cap can be exceeded** — `tools/run_p3l_contexts.py:63–66`. The restore forward is outside `remaining`, and `max(1, evaluations-previous_attempts)` grants fresh optimizer work even when the original cap is exhausted. Restoring a best point also restarts the coarse schedule. All added physics is ledger-billed, but the metadata's original `context_evaluation_cap` does not alone describe total fresh forwards. Retain distinct original optimization cap, recovery overhead and actual totals; do not claim an unchanged total cap unless enforced.

2. **False zero-new-physics recovery field** — line 52. `new_physics_for_validation=not cache_dest.exists()` is evaluated after `persist_cache` created a missing file, making the paid missing-cache branch report false. Capture existence before branching and report actual extra RHS/restores. Existing ledgers remain authoritative for paid cost.

3. **Multi-file interruption window blocks reuse** — lines 72–74 write context NPZ before cache and metadata; lines 44–47 assume metadata exists whenever NPZ exists. A stop in between produces an orphan NPZ that later resume cannot handle. Per-file atomic replacement does not provide a three-file transaction. Preserve the orphan, recover provenance from paid job records, validate it, then publish a completion receipt last.

4. **Single-worker requirement is not enforced** — `Book` locks one job directory and holds `budget.lock` only during ledger updates. Distinct jobs can run concurrently. `paid_workers=1` is a written constant, not measured concurrency. A parent-owned exclusive campaign worker lock/process inventory is still necessary; launching only one process operationally also satisfies the requirement while held.

5. **RSS guard has gaps** — lines 53 and 76–79. Completed-context reuse immediately continues and bypasses RSS checking; normal generation checks only after a full parent finishes. Neither prevents a transient within-parent peak. No additional paid work should be started from this audit.

6. **Recovery is not crash-idempotent** — reconciliation appends/releases a job before writing its closure, and publishes the aggregate receipt only after all jobs. A stop in either gap can leave reconciliation without a closure or an incomplete aggregate. It also selects active jobs from the snapshot despite the journal being authoritative. Do not blindly rerun after another interruption; reconcile the bounded journal tail and per-job records first. The six current closures agree with the aggregate, so no such partial state is demonstrated in this audit.

7. **Wall upper bound is conditional** — `recover_p3l_interruption.py:32` clamps elapsed wall to the original reservation. This is conservative only when the real-time quota did stop the process within its reservation. A deferred Python signal during a long native call, suspension or wall-clock changes can break that assumption. The six receipts have lower bounds below their 7,200 s reservations; no violation is evidenced here. Recovery wall includes waiting until reconciliation, as disclosed, rather than measured process CPU time.

## Other provenance/overwrite observations

Normal context runs refuse an existing context and `reconstruct` refuses an existing trace. The registered resume token permits rewriting `.recovery.json`, and missing-cache restoration can overwrite an orphan cache; no atomic parent-level claim prevents two distinct jobs targeting the same parent. Enforcing one paid process and disjoint ranges is necessary until a parent lock/transaction is added.

Packing validates every interrupted-source parent separately, but recovered source/validator closures are omitted from its `context_job_closures` summary because they are not inserted into `closures`. Validation is still executed; provenance reporting should include these validated pairs. Packer checks public frame metadata and orthogonality but does not independently prove frame construction or source-parent eligibility. Those are trusted generation provenance, not established anew by packing.

No substantive scientific conclusion, leak-free execution claim, memory stress-test pass, or authorization to restart is supplied by this review.

## Follow-up after corrective source edits

Read-only follow-up on 2026-10-22; this section supersedes the initial findings where explicitly marked fixed. No jobs launched, array/truth files opened, physics performed, or original receipts altered by the reviewer.

**Observed outcome:** `contexts_train_0000_16_memory_resume1/CLOSED.json` reports `COMPLETE`. Its final sampled `MEMORY_CURRENT.json`, at `SBI_train_0127`, reports **442.8125 MiB RSS**, zero native transfer entries and one paid worker. A next-shard sample at `SBI_train_0232` reports **445.484375 MiB**. These are post-parent samples, not peak-memory measurements or a full-campaign completion claim.

| Initial finding | Follow-up status |
| --- | --- |
| Recovery restore excluded from residual cap; forced minimum one optimizer attempt | Corrected for known partial cases: `remaining = evaluations - previous_attempts - int(initial is not None)`; explicit stop if below one. See remaining ordering issue below. |
| Missing-cache validation falsely reported zero new physics | Fixed: `cache_was_missing` captured before restoration and reused in the recovery receipt. Ledger remains authoritative for exact RHS. |
| Context NPZ without metadata could not resume | Recovery path added: derive a unique original interrupted job from paid best-point provenance, record reconstructed metadata, then equation-audit cache and compare saved q/prediction. Multi-file writes remain nontransactional, but the identified missing-sidecar state is now handled. |
| One paid worker only asserted | Fixed for `Book` users: nonblocking campaign-wide `paid_worker.lock` acquired before authorization/job setup and held until `Book.close`. A second participating paid job cannot start while the first holds it. |
| Reuse branch skipped RSS checks | Fixed: both reused and freshly generated parents invoke the shared memory-check function, which clears transfers and collects garbage. Within-parent peak is still unmeasured. |
| Recovered source/validator closures omitted from packing manifest | Fixed: each recovered parent now receives a distinct `job::parent` entry with original interrupted status, recovery receipt path and successfully completed validation-job closure. |

**Actual cap evidence for the already completed old-formula case:** `SBI_train_0078.json` records one previous optimizer attempt, a paid restore, and old nominal `remaining_optimizer_cap=15`. Its resumed `result.json` records **14 actual fresh forward states**, status `OPTIMIZER_RETURNED`, no failure. Therefore **1 + 1 + 14 = 16 actual paid forward attempts** for this parent. The old formula permitted an overrun in principle, but **no realized overrun is demonstrated for this parent**. Do not infer an actual 17-state cost from its nominal cap. Other resumed partials with three original attempts, one restore and twelve permitted new attempts satisfy the 16-state bound if executed under the inspected corrected formula.

**Remaining hard edge:** `t.evaluate(...restore_paid_best_after_memory_stop)` still precedes the `remaining < 1` check. If a future recovered parent had already consumed all sixteen attempts, the script would pay a seventeenth restore before refusing optimizer work. Check the complete remaining allowance before restoring, or use a previously persisted full-state cache without a fresh forward. This edge does not change the verified 16-state count for `SBI_train_0078`; no exhausted-parent execution was run here.

**Remaining operational limits:** recovery utility remains single-use and non-idempotent after a crash between journal reconciliation, closure publication and aggregate publication; current six completed interruption receipts match exactly and require no replay. Reservation-clamped wall remains a conditional conservative bound. Global lock protection assumes all paid work participates through `Book`; post-parent RSS sampling is not a guaranteed within-operation memory ceiling. New metadata-recovery uses unique original-job provenance, rather than restoring the complete optimizer state, and original NPZ/context evidence remains unchanged.
