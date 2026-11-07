# R11 cost ledger

Generated 2026-11-07T01:58:19.854539+00:00. Metadata-only compilation; no science, forward, optimizer, GT read or remote command.

Known nonoverlapping recorded CPU lower bound: **35211.857694 seconds**. This mixes documented exact child measurements and partial timers as a lower bound; it is not an exact total or a cross-hardware speed comparison. Ongoing48/3D work and unmeasured setup/import/exit costs prevent a complete total.

|Physical matrix dimension|Known LU|Known complex RHS|Known source RHS|Known adjoint RHS|
|---:|---:|---:|---:|---:|
|576|382015|9158172|2430552|2424444|
|5184|4272|470880|236760|234120|
|17496|0|120|120|0|

Counters are lower bounds when producer splits are missing. 17496 uses two fine-grid iterative forward batches, not LU. Small denominator LU/RHS remain auxiliary counters in JSON.

|Entry|Status|Recorded additive CPU s|Scope/source|
|---|---|---:|---|
|twod_remote_DEV|COMPLETED_KNOWN_CPU_LOWER_BOUND|1165.234148|Remote 2D DEV attempts/setup/transport/failed metadata; `delegated/twod/TWOD_DEV_COSTS.json`|
|twod_independent48|PARTIAL|18990.300679|Remote 48x6 confirmation completed lower bound, generation and transport; `delegated/twod/independent48/COSTS.json`|
|contracted_DEV3|COMPLETE|169.950853|Three actual local exposed-DEV inversions; `delegated/diagnostics/CONTRACTED_DEV_COST_LEDGER.json`|
|contracted_DEV3_scoring|COMPLETE|0.006708|Standalone cached postseal scoring; `delegated/diagnostics/CONTRACTED_DEV_COST_LEDGER.json`|
|BothBlocks_DEV3|COMPLETE|225.008589|Three actual local exposed-DEV inversions; `delegated/diagnostics/BOTHBLOCKS_DEV_COST_LEDGER.json`|
|BothBlocks_DEV3_controller|COMPLETED|0.067796|Distinct local queue controller; `delegated/diagnostics/BOTHBLOCKS_QUEUE.json`|
|BothBlocks_DEV3_scoring|COMPLETE|0.013115|Standalone cached postseal scoring; `delegated/diagnostics/BOTHBLOCKS_DEV_COST_LEDGER.json`|
|current_only_DEV3|COMPLETE|128.297977|Three actual local exposed-DEV inversions; `delegated/diagnostics/current_only_dev/CURRENT_ONLY_DEV_COST_LEDGER.json`|
|current_only_DEV3_controller|COMPLETED|0.055809|Distinct local queue controller; `delegated/diagnostics/current_only_dev/CURRENT_ONLY_QUEUE.json`|
|current_only_DEV3_scoring|COMPLETE|0.009891|Standalone cached postseal scoring; `delegated/diagnostics/current_only_dev/CURRENT_ONLY_DEV_COST_LEDGER.json`|
|BothBlocks16_inversions|COMPLETE|4281.798299|80 actual local confirmation children; `delegated/diagnostics/bothblocks_confirmation/COST_LEDGER.json`|
|BothBlocks16_generation|COMPLETE|1.372750|16 generation children/exports; `delegated/diagnostics/bothblocks_confirmation/COST_LEDGER.json`|
|BothBlocks16_controller|COMPLETE|0.448456|Distinct controller, no science; `delegated/diagnostics/bothblocks_confirmation/COST_LEDGER.json`|
|BothBlocks16_postprocess|COMPLETE|77.417023|Main scorer/frontier/two plot renders incl layout refinement; `delegated/diagnostics/bothblocks_confirmation/POSTSEAL_EXTERNAL_PROCESS_COSTS.json`|
|BothBlocks16_prediction_seal|COMPLETE|0.056023|Public-only selection/seal process; `delegated/diagnostics/bothblocks_confirmation/PREDICTION_SEAL_RECEIPT.json`|
|BothBlocks16_common_LU_selection|SEALED_PUBLIC_COST_AND_PHI_BEFORE_GT|2.041168|Public-only selection/seal process; `delegated/diagnostics/bothblocks_confirmation/COMMON_PHYSICS_CAP_SELECTION_BEFORE_GT.json`|
|fixed_state10_counterfactuals|COMPLETE|1.511641|10 exposed endpoints; final30 candidates plus10 archived seed-control attempts; `delegated/diagnostics/DIAGNOSTIC_RECEIPT.json`|
|metric_split|COMPLETE|0.042238|Independent zero-FW diagnostic / cached scoring; `delegated/diagnostics/METRIC_SPLIT_RECEIPT.json`|
|BothBlocks_wiring|PASS|0.053630|Independent zero-FW diagnostic / cached scoring; `delegated/diagnostics/BOTHBLOCKS_WIRING_CHECK.json`|
|current_curvature|COMPLETE|4.273310|Independent zero-FW diagnostic / cached scoring; `RESULTS/current_curvature_exposed/RECEIPT.json`|
|current_remainder|COMPLETE|0.047847|Independent zero-FW diagnostic / cached scoring; `RESULTS/current_remainder_exposed/RECEIPT.json`|
|exposed_DEV_scoring|COMPLETE|0.014224|Independent zero-FW diagnostic / cached scoring; `RESULTS/exposed_dev_scores/lagged/RECEIPT.json`|
|stage_compare|COMPLETE|0.080438|Independent zero-FW diagnostic / cached scoring; `delegated/diagnostics/BOTHBLOCKS_STAGE_COMPARISON_RECEIPT.json`|
|contracted_wiring|PASS|null|Standalone contracted material closure check; `delegated/diagnostics/CONTRACTED_WIRING_CHECK.json`|
|exposed_figures|COMPLETE|9.688470|Two distinct exposed figure render processes; `FIGURES/exposed_development/PLOT_PRIOR_RENDER_RECEIPT.json`|
|finite_noise_prior_audit|COMPLETE|1.261185|All16x5 cached finite-noise/prior audit; `delegated/diagnostics/bothblocks_confirmation/finite_noise_prior_audit/EXTERNAL_PROCESS_COST.json`|
|vector_top32_audit|PUBLIC_SAVED_FIELDS_ONLY|1.517860|Standalone vector diagnostic/setup; `delegated/vector/top32_audit/RECEIPT.json`|
|vector_source_basis|PUBLIC_BACKGROUND_SOURCE_BASIS|0.802634|Standalone vector diagnostic/setup; `delegated/vector/source_basis/RECEIPT.json`|
|vector_source_basis_audit|PUBLIC_SAVED_FIELDS_ONLY|1.191356|Standalone vector diagnostic/setup; `delegated/vector/source_basis_audit/RECEIPT.json`|
|vector_derivative_check|PASS|0.118662|Standalone vector diagnostic/setup; `delegated/vector/derivative_check/RECEIPT.json`|
|vector_z_derivative_check|PASS|0.178673|Standalone vector diagnostic/setup; `delegated/vector/z_derivative_check/RECEIPT.json`|
|vector_fixed_state_cuda|COMPLETE|1.875000|Standalone vector diagnostic/setup; `delegated/vector/fixed_state_cuda/RECEIPT.json`|
|vector_fixed_state_cuda_gate_helper|COMPLETE|0.031250|Separate cooperative-gate helper process; `delegated/vector/fixed_state_cuda/RECEIPT.json`|
|vector_source_gn_cuda|COMPLETE|1.953125|Standalone vector diagnostic/setup; `delegated/vector/source_gn_cuda/RECEIPT.json`|
|vector_source_gn_cuda_gate_helper|COMPLETE|0.046875|Separate cooperative-gate helper process; `delegated/vector/source_gn_cuda/RECEIPT.json`|
|vector_fixed_state_transport|COMPLETE|0.439798|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/fixed_state_transport/TRANSPORT_COST.json`|
|vector_source_gn_transport|COMPLETE|0.900750|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/source_gn_transport/TRANSPORT_COST.json`|
|vector_dev_ordinary_transport|COMPLETE|7.126572|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/dev_ordinary_transport/TRANSPORT_COST.json`|
|vector_remaining_transport|RUNNING_NO_FINAL_RECEIPT|null|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/remaining_transport/TRANSPORT_COST.json`|
|vector_MEASUREMENT_OVERLAY_COST|COMPLETE_KNOWN_COMPONENTS|0.049431|Separate local source/measurement overlay transport; `delegated/vector/remaining_transport/MEASUREMENT_OVERLAY_COST.json`|
|vector_LINALG_PROTOCOL_OVERLAY_COST|COMPLETE_KNOWN_COMPONENTS|0.042362|Separate local source/measurement overlay transport; `delegated/vector/remaining_transport/LINALG_PROTOCOL_OVERLAY_COST.json`|
|vector_LINALG_TIMING_AMENDMENT_COST|COMPLETE_KNOWN_COMPONENTS|0.037423|Separate local source/measurement overlay transport; `delegated/vector/remaining_transport/LINALG_TIMING_AMENDMENT_COST.json`|
|vector_ordinary_DEV6|COMPLETE|4088.453125|Six ordinary DEV arms in one shared scientific process; `delegated/vector/dev_ordinary_process/PROCESS_COST.json`|
|vector_ordinary_DEV6_gate_helper|COMPLETE|50.890625|Distinct gate helper of ordinary shared process; `delegated/vector/dev_ordinary_process/PROCESS_COST.json`|
|vector_confirmation_generation|COMPLETE|119.206795|10 independent vector observations generated, no inverse score; `delegated/vector/confirmation_generation/PROCESS_COST.json`|
|vector_dev_denominator|COMPLETE|4292.234375|Shared vector scientific phase including setup/probe and task arms; `delegated/vector/remaining_controller/PHASE_COSTS.json`|
|vector_dev_denominator_gate_helper|COMPLETE|48.859375|Distinct phase gate helper; `delegated/vector/dev_denominator_process/PROCESS_COST.json`|
|vector_confirmation|RUNNING_OR_QUEUED_NOT_SCORED|671.750000|Running/queued vector phase; completed nested producer lower bound only; `delegated/vector/confirmation_process/PROCESS_COST.json`|
|vector_remaining_controller|RUNNING|null|Still-running controller, own final CPU not yet received; `delegated/vector/remaining_controller/PROCESS_COST.json`|
|vector_DEV_early_postseal_controller|COMPLETE|1.179619|Independent local waiting/export controller; `delegated/vector/dev_early_postseal/PROCESS_COST.json`|
|vector_DEV_early_postseal_transport_children|COMPLETE|6.036545|Distinct local pack/download transport children; `delegated/vector/dev_early_postseal/PROCESS_COST.json`|
|vector_DEV_early_scorer|COMPLETE|14.214806|Old exposed DEV postseal scorer only; `delegated/vector/dev_early_postseal/PROCESS_COST.json`|
|vector_DEV_early_remote_pack|COMPLETE|4.234375|Independent remote archive producer, no physics; `delegated/vector/dev_early_postseal/REMOTE_PACK_COST.json`|
|first_CN99_hit_posthoc|COMPLETE_POSTHOC_OFFLINE_DIAGNOSTIC|2.941704|Two standalone zero-FW public selection then GT-scoring processes; `delegated/diagnostics/bothblocks_confirmation/first_cn99_hit/EXTERNAL_PROCESS_COSTS.json`|
|posthoc_budget_control_children|COMPLETE_POSTHOC_BUDGET_CONTROL|62.857739|Only triggered restarted original-Maxwell corrector children; original80 costs excluded; `delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control/COST_LEDGER.json`|
|posthoc_budget_control_controller|COMPLETE_POSTHOC_BUDGET_CONTROL|0.044582|Distinct finite continuation controller; `delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control/COST_LEDGER.json`|
|posthoc_budget_control_postprocess|COMPLETE|1.042965|Public32pair seal then separate postseal scorer; `delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control/POSTPROCESS_EXTERNAL_COSTS.json`|
|oracle_trace_boundary|COMPLETE|0.478958|Root saved-scored158626node CSV audit, no new GT arrays; `delegated/diagnostics/bothblocks_confirmation/oracle_trace_boundary/RECEIPT.json`|
|mixed_low_phase_plots|COMPLETE|4.124288|Only new all8 mixed low-scale plate and paired low-RMSE/mean chart; `FIGURES/bothblocks_confirmation/MIXED_LOW_PHASE_EXTERNAL_COST.json`|
|exposed_endpoint_path_diagnostic|COMPLETE|0.998588|Root fixed3x9 legal endpoint-segment actual physical evaluations; `RESULTS/exposed_endpoint_paths/EXTERNAL_PROCESS_COST.json`|
|exposed_endpoint_path_plot|COMPLETE|1.824294|Only plotting saved27 rows, no new physics or GT; `FIGURES/exposed_development/EXPOSED_ENDPOINT_PATH_PLOT_EXTERNAL_COST.json`|
|exposed_face_Hessian|COMPLETE_WITH_EXTERNAL_CPU_UNMEASURED|0.923438|Root exact face-Hessian/JVP wiring diagnostic at3 exposed endpoint materials; `RESULTS/exposed_face_curvature/RECEIPT.json`|
|absorptive_DEV_wiring|COMPLETE|0.870579|One public-state complex-adjoint and real-material derivative wiring process for K19/K0; `delegated/diagnostics/absorptive_dev/WIRING_EXTERNAL_COST.json`|
|absorptive_DEV8_children|COMPLETED|759.281010|Fixed4 exposed DEV x2 complex-shift controls, one local Mac lane; `delegated/diagnostics/absorptive_dev/QUEUE.json`|
|absorptive_DEV8_controller|COMPLETED|0.141694|Own finite sequential local queue; `delegated/diagnostics/absorptive_dev/QUEUE.json`|
|absorptive_DEV8_postprocess|COMPLETE|2.121386|Public all8 prediction seal then separately authorized offline scorer; `delegated/diagnostics/absorptive_dev/POSTPROCESS_EXTERNAL_COSTS.json`|
|vector_mac_confirmation_controller_own|RUNNING|null|Active finite Mac queue, final own CPU not yet written; `delegated/vector/mac_confirmation_controller/CONTROLLER_COST.json`|
|vector_host_split_engineering|COMPLETE_KNOWN_COMPONENTS|0.043008|Distinct local host dispatch/overlay transport only; `delegated/vector/remaining_transport/HOST_SPLIT_ENGINEERING_COST.json`|
|absorptive_DEV_material_plot|COMPLETE|1.699283|Saved4x6 exposed-material plate and source table only; `FIGURES/absorptive_dev/EXTERNAL_PROCESS_COST.json`|
|local_public_document_reference_check|PASS|0.040500|Local link and saved panel-source existence check only; `delegated/diagnostics/absorptive_dev/LOCAL_DOCUMENT_REFERENCE_CHECK.json`|
|unmeasured_setup_and_failures|KNOWN_GAPS|null (view only)|Recorded failed environment starts, source preparation and uninstrumented remote extraction; `delegated/vector/IMPLEMENTATION_NOTES.md`|
|historical_R10_reference|REFERENCE_ONLY_NOT_CHARGED|0.000000 (view only)|Read-only reused R10 states/operators/scores, not rerun R10 campaign; `no new execution source`|
|portfolio_views|NONADDITIVE_VIEW|0.000000 (view only)|Overlapping deployment/secondary views, no new solver run; `delegated/diagnostics/bothblocks_confirmation/PORTFOLIO_RULE_BEFORE_DATA_SCORE.json`|

## Running and incomplete coverage

48 progress snapshot: RUNNING / 128 terminal / 156 pending.
3D: RUNNING_NOT_SCORED. Whole running-phase/helper/transport CPU is null until final receipt; completed engineering probe is only a replaceable lower bound.

## Deduplication and gaps

- Never add nested arm/corrector/aux timers to their whole process.
- 576,5184 LU and17496 iterative solves remain separate; never infer matrix order from legacy accumulated full_factorization_dimension.
- Count source/primal and adjoint aliases only once.
- Fixed-state merged50FW already includes ordinary repair10FW.
- Portfolio scenarios overlap underlying attempts; never add them to campaign total.
- GPU stream interval, wall time, active kernel time and CPU are not interchangeable.
- Mac scientific wait4 is inclusive of its child RSS helpers; Windows GetProcessTimes is own-process and its independent helper is separate.
- Mac spool/write/flush/close telemetry is nested; phase2 inclusive parent wait4 is an overlapping view, never added to per-object scientific children.
- Missing/running/unmeasured values stay null; displayed totals are lower bounds, not complete exact project cost.

- Two early CUDA audit receipts contained a template suggestion of external OS CPU measurement. The producer confirms none was captured: post-import CPU is retained as a lower bound, imports/exit stay unknown.
- DEV/48 transport scopes are divided by their authoritative upstream producer ledgers; this compiler never adds their raw transport logs again.
- Engineering6LU probe is included in its future denominator shared process, and six ordinary vector arm receipts are inside one ordinary process.
- CPU for ordinary repair is a separate process timer, while its10FW are already merged in the fixed-state50FW count.
- Plot layout-refinement attempts are included; mathematical/scientific experiments are not rerun to fill accounting gaps.

Recompile after refreshing local producer ledgers:

```
python research/research_max_struct_r11/experiments/twod_compile_costs.py
python research/research_max_struct_r11/experiments/twod_independent_costs.py
python research/research_max_struct_r11/experiments/summarize_completed_costs.py
```

The compiler does not fetch remote results or stop/restart running queues. For exact scope, null reasons, raw auxiliary counts and source byte/mtime provenance, use COST_LEDGER.json.
