# Conditional4 scorer preparation

Status: PREPARED_NOT_EXECUTED. New conditional truth reads: 0. No conditional score, reconstruction figure or material result has been generated.

The scorer is `experiments/score_conditional_confirmation.py`. It reuses the original frozen `material_metrics` function: RMSE divided by sealed amplitude `max(true chi)-.1`, strict chi>1 support IoU, primary success NRMSE<=.15 and IoU>=.7 with legal material, plus both explicitly separated legacy Otsu metrics and region RMSE. It uses saved paid predictions only: no new forward or inverse calculation.

The public freeze is `delegated/twod/conditional_confirmation/FREEZE.json`, seed2026110711, four prespecified high objects (annulus12, slotted_shell18, bridge18, coated_inclusion18). The unchanged CIE prefix triggers BOTH cold fallbacks only if its saved physical data residual divided by public sigma strictly exceeds the public single-case CN99 radius. Each conditional deployment retains minimum original Phi between the prefix and its own fallback, with exact ties favoring the prefix. Invalid prefixes, failed fallbacks and untriggered NOT_RUN opportunities remain explicit.

Before any new conditional truth read, the scorer requires a separate authorization token and note, checks all public terminal queues/actual prediction receipts, verifies the public discrepancy trigger and exact fallback identities, independently verifies both min-Phi selectors, and compares copied selected prediction arrays with their source. It writes `ALL_ACTUAL_AND_SELECTOR_SEAL.json` before opening conditional score metadata or true material arrays. The intended output directory is write-once `delegated/conditional_confirmation_score/`.

Preparation inspected only the freeze, public manifest and task definitions, generation cost receipt, producer source and existing original24/right3 public setup formats. Reading producer source establishes the private schema without opening those private files. Neither the conditional scorer nor its main function has been executed. Syntax compilation passed.

## Prepared outputs

- `DEPLOYMENT_VIEWS.csv` and `DEPLOYMENT_SUMMARY.csv/.json`: CIE baseline, conditional CIE→Dual, conditional CIE→G19, all four objects, material/noise/Phi/stop/CPU/LU/RHS.
- `ALL_ATTEMPTS_WITH_NOT_RUN.csv`: every actual prefix/fallback and both untriggered fallback opportunities; NOT_RUN never becomes an inferred score.
- `ACTUAL_ATTEMPT_STAGES.csv`: saved native physical stages and Phi-selected corrector start, without extra forward work.
- `ACTUAL_ATTEMPT_COSTS.csv`: unique scientific attempts; shared prefix charged once. Its startup-repair field retains a documented zero-physics startup failure once.
- `RAW_WORKER_PROCESS_COSTS_NOT_ADDITIVE.csv`: each original/repaired process separately, as an alternative provenance view. It must not be added to the unique-attempt table.
- `CAMPAIGN_COSTS.json`: unique prefix+all actual fallback CPU, own queue processes, own pipeline segments, generation once. Nested inverse timers and alternate deployment views are not summed. Unknown values remain null. Physical576 and small19 LU remain separate.
- `ALL4_DEPLOYMENT_MATERIALS.png/.pdf`: four rows × GT/CIE/CIE→Dual/CIE→G19, common material scale, no case selection.
- `REVIEW.md`, seal and scoring receipt: separate cohort, no modifications to main24 or exposedright3 results.

## Producer interface reconciliation

The producer confirmed local public pipeline paths `delegated/twod/conditional_confirmation/{RECEIPT,TRIGGER_DECISIONS,FALLBACK_TASKS}.json`, with prediction/selector trees `RESULTS/conditional_confirmation/{prefix,fallback,deployed}`. Queue-own GetProcessTimes CPU is in `queues[].external_queue_process_cost`; the pipeline own controller time is `pipeline_controller_CPU_seconds`. Local SSH/transport wrapper cost remains separately itemized engineering cost, not silently allocated across deployment alternatives.

One original prefix worker for conditional0001 rejected an output-path guard before physics (reported0.125s). The producer performs one explicitly documented repair of this previously unexecuted solve, preserving its original failure and original pipeline receipt. The scorer accepts only this exact zero-physics exception: one prior failed job with no physical cost and the recorded output-guard message, replaced by the finished `startup_repair_queue` job. It preserves both process costs, adds the startup cost to that prefix's deployment cost, and requires the sealed selector to have the same prefix cost. It does not treat this as a second scientific prefix or a numerical retry. Unexpected duplicate task IDs fail closed.

The final producer reconciliation and root's explicit postseal authorization remain pending. Preparation is not a claim that those final gates have already passed.

## Later invocation, only after root authorization

Use the existing project Python with `experiments/score_conditional_confirmation.py --authorize-postseal ROOT_EXPLICIT_CONDITIONAL4_POSTSEAL_AUTHORIZATION --authorization-note '<actual separate root authorization and terminal reconciliation>'`. Do not use this preparation note as authorization.

Preparation CPU/wall were not measured as a complete process window and remain unknown; physicalLU=0, RHS=0, newinverse=0, newconditionalGTreads=0. No algorithm, budget, threshold, denominator, main24 score or right3 score changed.
