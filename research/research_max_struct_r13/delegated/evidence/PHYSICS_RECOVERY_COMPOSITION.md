# Physics dependency recovery composition (source only)

`src/compose_physics_recovery.py` prepares a new isolated local evidence root. It has not been invoked on study results. No truth, prediction arrays or quality values are loaded. Original queue/results/ZIP stay unchanged.

Parent invocation after both terminal queues exist:

```text
python src/compose_physics_recovery.py --root ORIGINAL_R13_ROOT --original-manifest ORIGINAL_R13_ROOT/configs/physics_v1.json --recovery-manifest ORIGINAL_R13_ROOT/configs/physics_recovery_v1.json --out NEW_ISOLATED_EVIDENCE_ROOT --archive-name physics_v1_recovered_full_results_v2.zip
```

The tool requires exactly all60 original attempts,18 original technical successes and42 failures. The18 original successes must be the classic Austria cohort. Recovery must be exactly those42 failed IDs and complete technically; no quality/alpha threshold is applied. Case/method/input reference/budget and all other task fields, plus the entire numerical protocol, must match after slash normalization; only campaign/explicit parent-attempt metadata (including narrowly named parent_campaign and technical_recovery) is ignored. Both manifests must refer to the same original public observations. No NPZ content comparison or reconstructed numerical certificate is claimed. Result RECEIPT.task/protocol must also match. Input files are not copied: retain the original public input source or restore separately.

New evidence layout:

- RESULTS/physics_v1:18 original successful copies +42 recovery copies; recovered receipt task.campaign remains provenance of its actual attempt.
- queues/physics_v1/attempts/initial_RECEIPT.json: byte-preserved initial terminal receipt.
- audit/physics_failed_initial/ID: byte copies of original failed/empty directories where they existed; absence does not become invented evidence.
- queues/physics_recovery_v1 and RESULTS/physics_recovery_v1: original retry evidence for nested accounting references.
- queues/physics_v1/RECEIPT.json: COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY,60jobs,42failed_attempts, recovery_source pointer.

CPU semantics: finished_child_CPU_seconds is only the60 final primary jobs; failed_attempt_child_CPU_seconds is the additional42 original failed attempts, separate from that total. Missing any term gives null, never zero. Controller CPU is the original plus recovery controller, counted once in merged parent. Wall is the sum of the two controller active windows, explicitly excluding intervention gap; never child-wall sum. The cost collector now consumes failed_attempts as separate online_failed_attempt rows and excludes nested recovery jobs/controller through recovery_source. Original failed physical counters come only from their copied receipts, otherwise null.

The merged receipt declares archive_filename=physics_v1_recovered_full_results_v2.zip (or parent explicit distinct version); package_artifacts respects that field for final evidence and additionally requires/publishes immutable initial physics_v1_full_results.zip as its own failed-initial asset linked to attempts/initial_RECEIPT.json. Missing either original ZIP or receipt blocks prepare, keeping initial42 failure logs accessible. Composer creates no archive and copies no original archive. New output paths require new seals and parent scoring authorization. Source helper checks used synthetic metadata only: valid60/42 mapping and Windows slashes accepted; changed budget/input, wrong retry cohort and active queue rejected. Syntax checks passed. No composition/scoring/solver execution performed.
