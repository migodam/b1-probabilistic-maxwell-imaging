# Offline scoring design — paired R12 supplement

`src/score_and_plot.py` is offline-only: no Maxwell model, solver, optimizer, or training imports. It must not be run until root explicitly authorizes offline scoring. Compilation/import inspection is safe and does not read truth.

## Gate and invocation

Root seal JSON must contain status `ALL44_PREDICTIONS_SEALED`, `root_authorized_offline_scoring: true`, and `task_ids` equal to all44 frozen config task IDs. Optional `files` entries contain `path`, `size`, and `mtime_ns`, checked before private access. All44 terminal receipts are required, including failures; online truth declarations must be NONE. Candidate/history lengths and paired snapshot metadata are checked. Failed tasks remain NA in every statistical denominator, while partial snapshots remain auditable.

After root authorizes:

```
python supplement_paired/src/score_and_plot.py \
  --seal supplement_paired/audit/ALL44_PREDICTIONS_SEALED.json \
  --controller supplement_paired/queues/main/RECEIPT.json
```

No output is written before the gate passes. Then `RESULTS/offline_v1/TRUTH_ACCESS_BEGIN.json` records first authorized private access. Existing offline_v1 is never overwritten.

## Material metrics and truth provenance

A: use actual study sealed_truth/basin if present, otherwise original sibling r10_worktree sealed_truth/twod/<fullID>.npz, matching the original R12 basin scorer. B: delegated/confirmation/sealed_truth/confirm24/<fullID>.npz. C: original sibling r8_worktree sealed_truth/warm_test/<fullID>.npz and score_metadata amplitude, as in delegated/neural/score_nn_interface.py. Actual paths are recorded in TRUTH_SOURCES.

A/B amplitude NRMSE = (||m−GT||₂/24)/(max(GT)−0.1), relativeL2 = ||m−GT||₂/||GT||₂. IoU uses fixed χ>1 sets. C preserves original R8/R7 amplitude metadata and prediction-only 256-bin Otsu segmentation vs saved mask; its IoU is explicitly labeled distinct. C float32 bounds retain original one-ULP box tolerance; A/B use inherited 1e−12. Old descriptive success is legal box + NRMSE≤.15 + IoU≥.7. No new success cutoff, no GT selection, no pooled IoU across different definitions.

## Pairing and denominator

A three exposed bad endpoints × four methods. B all24 original G19 primary endpoints × four deployment strategies; only four public CN99-triggered cases run rescue. Other20 retain identical original materials, physical metadata and parent acquisition cost. C four identical frozen NN starts × four methods. Failure/missing selection always retains the object in denominator. Means identify their available count; inherited descriptive successes use full denominators. Paired-vs-direct table includes paired availability. This is exposed posthoc descriptive evidence, not fresh confirmation or novelty judgment.

Store initial/raw/corrected/primary/secondary scores separately. Primary is runner minPhi among paid eligible candidates. Secondary remains unavailable if no legal CN99 candidate exists; never silently impute primary. Prefix choice is recomputed strictly from public eligible CPU/Phi history and checked against PREFIXES, then scored. All source histories and candidate arrays stay untouched. Direct raw is the shared initial incumbent by runner contract, explicitly distinct from an auxiliary proposal.

## Costs

External controller `jobs[*].external.child_process_CPU_seconds` is authoritative full-exit CPU. `jobs[*].cpu` is affinity, **not CPU seconds**. Internal pre-exit CPU and nested physics/projection/G19 counters are itemized but never added. Parent acquisition is recorded separately, once per case; B full strategies charge original parent plus added rescue. C shared preprocessing/import/model load total is charged once, not once per arm or added again to its nested feature/inference values. Per-case warm costs are separate deployment views, not extra overhead on the shared total.

## Figures and outputs

FIGURES receives PNG+PDF full-ID mapping; A/C all-object GT/start/four-method raw and corrected grids; B all24 original/four-strategy material and signed-error grids (8 objects/page); every-object paid CPU curves for NRMSE/Phi; corrector arrows for D/R/Phi/NRMSE/CN99 ratio. All unavailable slots remain visible. Continuous metrics are primary visual evidence; old success is a descriptive table field.

Historical DR24 and vector3D saved decomposition arrays are copied into a separately labeled context JSON and plotted separately. Their objectives/budgets differ, so no pooled mean or equal-cost inference is made. No new forward solves or receiver operators are assembled.

Outputs: SCORES, PAID_PREFIX_SCORES, CORRECTOR_DECOMPOSITION, PAIRED_VS_DIRECT, SUMMARY, CHILD_COST_LEDGER, PARENT_COST_LEDGER, STRATEGY_COST_TOTALS and TRUTH_SOURCES (JSON+CSV); COST_ACCOUNTING, HISTORICAL_CONTEXT_SEPARATE, SCORING_RECEIPT JSON. Final scientific interpretation belongs to root.

## Secondary deployment and G19 comparisons

Original secondary rows retain conditional availability and conditional means. The derived `secondary_deploy` view chooses secondary whenever an available sealed secondary candidate exists, otherwise primary, and records `secondary_available`, `fallback_to_primary`, and `deployment_source_view`. This adds no online candidate, selection step, or GT criterion. It includes all31 IDs under each method (A3/B24/C4), and all20 retained B cases deploy their unchanged original primary as both primary and secondary. `SECONDARY_DEPLOYMENT` provides the explicit full-denominator table; SUMMARY includes that view and fallback counts.

`PAIRED_VS_G19` complements PAIRED_VS_DIRECT for raw/corrected/primary/secondary/secondary_deploy. Each frozen case/method difference retains paired availability; failed or unavailable pairs have null differences. `G19_INITIAL_PROJECTION_DEFECT` extracts the saved inherited LS projection diagnostics from G19 receipts, with explicit missingness and no new current evaluation.

## Plot presentation and actual-cost comparison

A/C mapping says FIXED PAIRED; trigger/retain labels apply only to B. Material grids share a labeled χ colorbar fixed at [0.02,20]; error grids share labeled signed χ error at ±20. Corrector arrows use consistent G19/DR/AP/direct colors and a legend; direct raw is explicitly the initial incumbent with no auxiliary step. Positive Phi prefix axes use log scale. B paid-prefix panels are paginated at eight objects/page.

PERFORMANCE_VS_ACTUAL_MEAN_CPU gives primary and secondary-deployment cohort means, full object and available counts, and separate actual mean parent, added rescue and total deployment CPU. B divides four rescue-child costs by all24 objects, retaining original parent cost for all24. This is descriptive observed cost/performance, with no GT-selected points or implied equal-spend comparison.

Historical context now includes residual/δ where saved data supports it. Original DR24 ratios derive from already saved relative residual and unchanged public Y/CN99 radius. Vector3D uses raw.noise.residual_over_CN99_delta and corrected.residual_over_CN99_delta exactly; missing baseline raw physical residuals stay NA (bounds are not silently treated as exact ratios).

## Cost correction and portable offline reproduction

Version2 prioritizes `parent_original_actual_cost.completed_external_child_CPU_seconds` over the nested pre-exit `whole_child_CPU_seconds`. B's authoritative original-parent sum is 4290.1875 CPU seconds; v1 retained the internal 4286.734375 sum. Each parent ledger row now identifies its CPU scope. If A has no completed external measurement, its original receipt scope is explicitly retained rather than relabeled external. This correction affects costs only, not selected materials or physical/scientific metrics.

`--output-dir <new-directory>` allows root to preserve offline_v1 and write offline_v2. The directory must not already exist. Material/prefix selections and online run outputs are unchanged.

After seal validation and authorized truth loading, the scorer exports each of31 cases to `offline_truth/{A,B,C}/<ID>.npz` with m, support and amplitude, plus provenance JSON. These are offline-only reproduction inputs; the online runner does not read this directory. Future postseal scoring prefers these local portable exports, so sibling R8/R10 worktrees are unnecessary for material truth loading. Exports and provenance are never overwritten.

B grids label the inherited column “original G19 primary” and use “TRIGGERED” for rescue attempts. Historical vector arrows now identify methods by color and legend with DR emphasized. Direct raw remains the initial incumbent (no auxiliary proposal). Current frozen-NN supplement CPU and historical NN-interface CPU were measured on different hardware/configurations; do not infer speedups across those runs.
