# Normal-memory Maxwell mechanism experiment: offline evaluation

## Scope

Three already exposed physical instances; exact published field/current projectors; shared states after 3000 AP steps; ordinary, fractional-event and integer-event DR continuations of 3000 steps. No corrector. No truth is used in initialization, threshold calculation, trajectories or minimum-Phi selection. This file independently scores the frozen material arrays after the runs.

The original `offline_scores.json` is preserved. Its 12 scores agree with recomputation to 1e-12. Integer-event scores are newly appended to the combined artifact.

## Material results

|Case|Original material|AP shadow|Ordinary DR|Fractional event|Integer event|
|---|---:|---:|---:|---:|---:|
|A: slotted shell|0.515137|0.521084|0.187079|0.177757|0.186115|
|B: persistent failure|0.365795|0.368014|0.365795|0.365795|0.365795|
|C: exposed NN case|0.273914|0.334856|0.334856|0.133181|0.128761|

The B continuations do not improve the original incumbent. Its best iteration is −1; therefore the selected B material is the original input, not the AP endpoint or a successful reflected reconstruction.

## Finite-event prediction

|Case|Predicted normal t*|Predicted ordinary index|Observed large box switch|Relative AP-pair defect|
|---|---:|---:|---:|---:|
|A: slotted shell|5.731722|5|5|0.000947753|
|B: persistent failure|12.495138|12|12|0.000906983|
|C: exposed NN case|14.079063|14|14|0.00598649|

The ordinary branch starts at normal coordinate t=1, so the first strict switch is predicted at index floor(t*). Small interior motion before that event is not a full box switch. The three AP pairs are approximate: the prediction is checked against actual DR, not asserted as an exact trajectory theorem for these finite-precision states. The B negative example shows that a correctly predicted switch need not lead to recovery.

## Cost and limits

Measured probe CPU plus all nine timed continuation branches: **241.448473 s**. This is a measured subtotal, **not** an end-to-end CPU total. The three probe times already include 3000 AP steps, both 500-step exploratory branches, setup and their physical audits. Six later process imports, projector setup, initial physical scores, serialization and process exit are outside continuation timers and cannot be reconstructed as seconds from the JSON. Do not label the saved field `shared_AP_cpu_seconds...` as AP-only cost.

Exact physical work derived from receipts and source: **624 LU / 6240 complex RHS**, including six initial forward solves omitted from branch counters. The projector was constructed nine times, hence nine H Cholesky factorizations. Full projection-action counts are preserved per case in `combined_offline_scores.json`.

This study does not establish an end-to-end speedup, improved generalization, or superiority over published G19/DR. It starts after a paid 3000-step AP preparation and uses only three exposed objects. In C, published DR from the original NN material reaches NRMSE 0.077988; the fractional and integer event continuations here reach 0.133181 and 0.128761. In A, ordinary and event continuations all recover part of the shell, without the published physical corrector. B remains a negative result.

## Outputs

`combined_offline_scores.json`, `comparison.csv`, `finite_thresholds.csv`, `physical_traces.csv`, and the three PNG/PDF figure pairs under `figures/`. The complete plotting/scoring program is `work/plot_normal_memory.py`.
