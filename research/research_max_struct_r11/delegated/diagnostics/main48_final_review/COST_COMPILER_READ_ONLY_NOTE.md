# Cost-compiler read-only note: review elapsed time is not CPU

Bounded source reviewed: only `experiments/summarize_completed_costs.py`; it was not run or modified, and neither COST_LEDGER file was run/updated.

The completed independent review receipt records elapsed wall232.443981s and CPU=null. The elapsed value measures workflow time including reasoning/tool waits, not process CPU. Text-only report/template correction likewise has no measured CPU; keep null. Do not use elapsed/wall as a replacement, do not impute zero, and do not add232.44 to known CPU totals.

Current compiler mechanisms already support this correctly: `add` lines25–35 retains `cpu=None` and emits an unknown-cost reason; totals lines213–215 sum only non-null CPU; Markdown formatting at229–230 prints null. There is currently no explicit main48_final_review entry. At final consolidation root can add one completed editorial/review entry sourced from this directory's RECEIPT/TEXT_CORRECTION_RECEIPT, with CPU=None and zero physical operations, preserving elapsed only as nonadditive wall context. Source capture at31–32 recognizes wall_seconds but not review_elapsed_wall_seconds, so the latter needs an explicit separate context field or note; never rename it to a CPU key.

Suggested record semantics: status=COMPLETE_CPU_UNMEASURED; CPU_known_lower_bound_seconds=null; CPU_measurement_scope='Review/editorial CPU not instrumented'; physical solve counts0; reported_review_elapsed_wall_seconds=232.443981, nonadditive. If an external CPU measurement is later found for a genuinely disjoint mechanical child, add only that recorded measurement once; do not substitute the compiler's own current process_time for the historical review.

No source change was made. Root will refresh the total ledger after scientific work completes. This note is engineering accounting only.

A later, disjoint mechanical mixed-phase summary has a real external wait4 measurement:0.104296CPU seconds, source `mixed_phase_summary/EXTERNAL_PROCESS_COST.json`. It can be added once at final consolidation, replacing its nested self timer. This does not turn review elapsed232.443981s intoCPU, and neither total ledger nor compiler was executed/modified here.
