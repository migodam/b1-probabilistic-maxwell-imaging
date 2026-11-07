# R11 cost ledger

Generated 2026-11-07T08:24:59.402255+00:00. Metadata-only compilation; no science, forward, optimizer, GT read or remote command.

Known nonoverlapping recorded CPU lower bound: **133106.212714 seconds**. This mixes documented exact child measurements and partial timers as a lower bound; it is not an exact total or a cross-hardware speed comparison. All frozen scientific tasks are terminal; unmeasured setup/import/exit/editorial costs still prevent an exact campaign total.

|Physical matrix dimension|Known LU|Known complex RHS|Known source RHS|Known adjoint RHS|
|---:|---:|---:|---:|---:|
|576|802447|19236588|3780732|3771432|
|2304|8|96|96|0|
|5184|20507|2449020|1230420|1218600|
|17496|0|120|120|0|

Counters are lower bounds when producer splits are missing. 17496 uses two fine-grid iterative forward batches, not LU. Small denominator LU/RHS remain auxiliary counters in JSON.

|Entry|Status|Recorded additive CPU s|Scope/source|
|---|---|---:|---|
|twod_remote_DEV|COMPLETED_KNOWN_CPU_LOWER_BOUND|1165.234148|Remote 2D DEV attempts/setup/transport/failed metadata; `delegated/twod/TWOD_DEV_COSTS.json`|
|twod_independent48|ALL288_DESIGNATED_TERMINAL_END_TO_END_COST_INCOMPLETE|51125.855796|Remote 48x6 confirmation completed lower bound, generation and transport; `delegated/twod/independent48/COSTS.json`|
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
|vector_top32_audit|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|1.517860 (view only)|Standalone vector diagnostic/setup; `delegated/vector/top32_audit/RECEIPT.json`|
|vector_source_basis|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.802634 (view only)|Standalone vector diagnostic/setup; `delegated/vector/source_basis/RECEIPT.json`|
|vector_source_basis_audit|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|1.191356 (view only)|Standalone vector diagnostic/setup; `delegated/vector/source_basis_audit/RECEIPT.json`|
|vector_derivative_check|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.118662 (view only)|Standalone vector diagnostic/setup; `delegated/vector/derivative_check/RECEIPT.json`|
|vector_z_derivative_check|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.178673 (view only)|Standalone vector diagnostic/setup; `delegated/vector/z_derivative_check/RECEIPT.json`|
|vector_fixed_state_cuda|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|1.875000 (view only)|Standalone vector diagnostic/setup; `delegated/vector/fixed_state_cuda/RECEIPT.json`|
|vector_fixed_state_cuda_gate_helper|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.031250 (view only)|Separate cooperative-gate helper process; `delegated/vector/fixed_state_cuda/RECEIPT.json`|
|vector_source_gn_cuda|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|1.953125 (view only)|Standalone vector diagnostic/setup; `delegated/vector/source_gn_cuda/RECEIPT.json`|
|vector_source_gn_cuda_gate_helper|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.046875 (view only)|Separate cooperative-gate helper process; `delegated/vector/source_gn_cuda/RECEIPT.json`|
|vector_fixed_state_transport|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.439798 (view only)|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/fixed_state_transport/TRANSPORT_COST.json`|
|vector_source_gn_transport|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.900750 (view only)|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/source_gn_transport/TRANSPORT_COST.json`|
|vector_dev_ordinary_transport|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|7.126572 (view only)|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/dev_ordinary_transport/TRANSPORT_COST.json`|
|vector_remaining_transport|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|23.478604 (view only)|Local transport/controller only; remote scientific CPU excluded; `delegated/vector/remaining_transport/TRANSPORT_COST.json`|
|vector_MEASUREMENT_OVERLAY_COST|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.049431 (view only)|Separate local source/measurement overlay transport; `delegated/vector/remaining_transport/MEASUREMENT_OVERLAY_COST.json`|
|vector_LINALG_PROTOCOL_OVERLAY_COST|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.042362 (view only)|Separate local source/measurement overlay transport; `delegated/vector/remaining_transport/LINALG_PROTOCOL_OVERLAY_COST.json`|
|vector_LINALG_TIMING_AMENDMENT_COST|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.037423 (view only)|Separate local source/measurement overlay transport; `delegated/vector/remaining_transport/LINALG_TIMING_AMENDMENT_COST.json`|
|vector_ordinary_DEV6|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|4088.453125 (view only)|Six ordinary DEV arms in one shared scientific process; `delegated/vector/dev_ordinary_process/PROCESS_COST.json`|
|vector_ordinary_DEV6_gate_helper|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|50.890625 (view only)|Distinct gate helper of ordinary shared process; `delegated/vector/dev_ordinary_process/PROCESS_COST.json`|
|vector_confirmation_generation|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|119.206795 (view only)|10 independent vector observations generated, no inverse score; `delegated/vector/confirmation_generation/PROCESS_COST.json`|
|vector_dev_denominator|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|4292.234375 (view only)|Shared vector scientific phase including setup/probe and task arms; `delegated/vector/remaining_controller/PHASE_COSTS.json`|
|vector_dev_denominator_gate_helper|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|48.859375 (view only)|Distinct phase gate helper; `delegated/vector/dev_denominator_process/PROCESS_COST.json`|
|vector_confirmation|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|11018.875000 (view only)|Shared vector scientific phase including setup/probe and task arms; `delegated/vector/remaining_controller/PHASE_COSTS.json`|
|vector_confirmation_gate_helper|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|127.281250 (view only)|Distinct phase gate helper; `delegated/vector/confirmation_process/PROCESS_COST.json`|
|vector_remaining_controller|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|4.187500 (view only)|Controller own lifetime, no science; `delegated/vector/remaining_controller/PROCESS_COST.json`|
|vector_DEV_early_postseal_controller|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|1.179619 (view only)|Independent local waiting/export controller; `delegated/vector/dev_early_postseal/PROCESS_COST.json`|
|vector_DEV_early_postseal_transport_children|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|6.036545 (view only)|Distinct local pack/download transport children; `delegated/vector/dev_early_postseal/PROCESS_COST.json`|
|vector_DEV_early_scorer|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|14.214806 (view only)|Old exposed DEV postseal scorer only; `delegated/vector/dev_early_postseal/PROCESS_COST.json`|
|vector_DEV_early_remote_pack|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|4.234375 (view only)|Independent remote archive producer, no physics; `delegated/vector/dev_early_postseal/REMOTE_PACK_COST.json`|
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
|vector_MAC_object_c11_matched12|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|5969.315371 (view only)|Five scientific arms in one actual Mac object child; `delegated/vector/mac_confirmation_controller/PROCESS_RECEIPTS.json`|
|vector_mac_confirmation_controller_own|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|4.438037 (view only)|Own finite Mac queue CPU plus distinct monitor helpers; `delegated/vector/mac_confirmation_controller/CONTROLLER_COST.json`|
|vector_MAC_object_c15_matched12|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|6190.119862 (view only)|Five scientific arms in one actual Mac object child; `delegated/vector/mac_confirmation_spooled_controller/PROCESS_RECEIPTS.json`|
|vector_MAC_object_c14_matched12|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|6315.360375 (view only)|Five scientific arms in one actual Mac object child; `delegated/vector/mac_confirmation_spooled_controller/PROCESS_RECEIPTS.json`|
|vector_mac_confirmation_spooled_controller_own|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|5.404159 (view only)|Own finite Mac queue CPU plus distinct monitor helpers; `delegated/vector/mac_confirmation_spooled_controller/CONTROLLER_COST.json`|
|vector_mac_spool_handoff_own|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.852774 (view only)|Finite parent handoff helper, no physical worker interruption; `delegated/vector/mac_spool_controller_handoff/PROCESS_COST.json`|
|vector_mac_phase2_inclusive_view|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|12510.889336 (view only)|Overlapping phase2 parent+all descendants measurement; `delegated/vector/mac_spool_controller_handoff/PHASE2_INCLUSIVE_WAIT4_VIEW.json`|
|vector_host_split_engineering|SUPERSEDED_NONADDITIVE_VECTOR_VIEW|0.043008 (view only)|Distinct local host dispatch/overlay transport only; `delegated/vector/remaining_transport/HOST_SPLIT_ENGINEERING_COST.json`|
|absorptive_DEV_material_plot|COMPLETE|1.699283|Saved4x6 exposed-material plate and source table only; `FIGURES/absorptive_dev/EXTERNAL_PROCESS_COST.json`|
|local_public_document_reference_check|PASS|0.040500|Local link and saved panel-source existence check only; `delegated/diagnostics/absorptive_dev/LOCAL_DOCUMENT_REFERENCE_CHECK.json`|
|material_column_energy_fixed10|COMPLETE|0.793096|Ten original exposed saved G19/ordinary same-state material-column-energy comparisons; `delegated/diagnostics/material_column_energy/EXTERNAL_PROCESS_COST.json`|
|paired_grid8_generation|COMPLETE|8.816213|Eight independent fine geometries with matched24/cross48to24 paired observations,16FW total; `delegated/diagnostics/crossgrid8_protocol/GENERATION_V2_EXTERNAL_COST.json`|
|paired_grid8_public_package|PUBLIC_ONLY_PACKAGE_PREPARED_NO_UPLOAD|0.141961|Local public-only inverse input archive, no remote upload or inverse; `delegated/diagnostics/crossgrid8_protocol/PUBLIC_PACKAGE_RECEIPT.json`|
|paired_grid8_inverse_and_postprocess|ALL64_SCORED_AND_PLOTTED_POST_PUBLIC_SEAL|10950.631673|64 scientific children plus their disjoint controller/transport/seal/scoring/figure producers; `delegated/twod/crossgrid8_paired/COSTS.json`|
|paired_grid8_root_public_seal_review|ROOT_EXPLICIT_POST64SEAL_SCORE|0.188900|Root-only verification of64 terminal identities and576 public predictions before GT authorization; `audit/PAIRED_GRID8_ROOT_POSTSEAL_AUTHORIZATION.json`|
|main48_endpoint_table_analysis|COMPLETE_KNOWN_CHILD_COST|0.242041|Saved public endpoint/authorized score join, including two preserved engineering failures; `delegated/diagnostics/main48_endpoint_analysis/COST_LEDGER.json`|
|main48_scaled_projection_correction|COMPLETE|0.279910|Distinct saved-table arithmetic process; `delegated/diagnostics/main48_endpoint_analysis/scaled_projection_correction/EXTERNAL_PROCESS_COST.json`|
|main48_mixed_phase_table|COMPLETE|0.104296|Distinct saved-table arithmetic process; `delegated/diagnostics/main48_final_review/mixed_phase_summary/EXTERNAL_PROCESS_COST.json`|
|main48_independent_review_editorial|COMPLETE_CPU_UNMEASURED|null|Bounded completed review and text-only report corrections; `delegated/diagnostics/main48_final_review/RECEIPT.json`|
|Mac_postrun_package_inventory|POSTRUN_READONLY_PACKAGE_METADATA_NOT_LAUNCH_RECEIPT|0.107974|Read-only package metadata observation after runs, not a launch receipt; `audit/MAC_POSTRUN_PACKAGE_METADATA.json`|
|root_vector50_public_seal_review|ROOT_AUTHORIZED_FRESH3D_OFFLINE_SCORING_AFTER_ALL50_PUBLIC_SEALS|0.715675|Root review of all50 public material/data/prefix selections before freshGT; `audit/VECTOR50_ROOT_POSTSEAL_AUTHORIZATION.json`|
|lossless_publication_archives|COMPLETE_KNOWN_COMPONENTS|351.140912|Publication only: archive or repository/release transport; zero new science; `audit/PACKAGING_RECEIPT.json`|
|publication_control|COMPLETE_EVIDENCE_RELEASE_PUBLISHED|42.375180|Publication only: archive or repository/release transport; zero new science; `audit/ROOT_RELEASE_PROCESS_COST.json`|
|publication_snapshot_supplement|COMPLETE|0.195328|Publication only: archive or repository/release transport; zero new science; `audit/PACKAGING_SUPPLEMENT_RECEIPT.json`|
|publication_supplement_transport|COMPLETE|0.464214|Publication only: archive or repository/release transport; zero new science; `audit/ROOT_RELEASE_SUPPLEMENT_COST.json`|
|unmeasured_setup_and_failures|KNOWN_GAPS|null (view only)|Recorded failed environment starts, source preparation and uninstrumented remote extraction; `delegated/vector/IMPLEMENTATION_NOTES.md`|
|historical_R10_reference|REFERENCE_ONLY_NOT_CHARGED|0.000000 (view only)|Read-only reused R10 states/operators/scores, not rerun R10 campaign; `no new execution source`|
|portfolio_views|NONADDITIVE_VIEW|0.000000 (view only)|Overlapping deployment/secondary views, no new solver run; `delegated/diagnostics/bothblocks_confirmation/PORTFOLIO_RULE_BEFORE_DATA_SCORE.json`|
|vector_unified_top32_audit|PUBLIC_SAVED_FIELDS_ONLY|1.517860|Independent measured diagnostic/build scope; zero new physical Maxwell solves; `delegated/vector/COSTS.json`|
|vector_unified_source_basis|PUBLIC_BACKGROUND_SOURCE_BASIS|0.802634|Independent measured diagnostic/build scope; zero new physical Maxwell solves; `delegated/vector/COSTS.json`|
|vector_unified_source_basis_audit|PUBLIC_SAVED_FIELDS_ONLY|1.191356|Independent measured diagnostic/build scope; zero new physical Maxwell solves; `delegated/vector/COSTS.json`|
|vector_unified_derivative_check|PASS|0.118662|Independent measured diagnostic/build scope; zero new physical Maxwell solves; `delegated/vector/COSTS.json`|
|vector_unified_z_derivative_check|PASS|0.178673|Independent measured diagnostic/build scope; zero new physical Maxwell solves; `delegated/vector/COSTS.json`|
|vector_unified_fixed_state_cuda|COMPLETE|1.875000|Independent measured diagnostic/build scope; zero new physical Maxwell solves; `delegated/vector/COSTS.json`|
|vector_unified_fixed_state_cuda_query_helpers|COMPLETE|0.031250|Distinct helper CPU,not in Windows parent CPU; `delegated/vector/COSTS.json`|
|vector_unified_source_gn_cuda|COMPLETE|1.953125|Independent measured diagnostic/build scope; zero new physical Maxwell solves; `delegated/vector/COSTS.json`|
|vector_unified_source_gn_cuda_query_helpers|COMPLETE|0.046875|Distinct helper CPU,not in Windows parent CPU; `delegated/vector/COSTS.json`|
|vector_unified_fixed_state_transport|COMPLETE|0.439798|Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work; `delegated/vector/COSTS.json`|
|vector_unified_source_gn_transport|COMPLETE|0.900750|Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work; `delegated/vector/COSTS.json`|
|vector_unified_dev_ordinary_transport|COMPLETE|7.126572|Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work; `delegated/vector/COSTS.json`|
|vector_unified_remaining_transport|FAILED_SCIENCE_EVIDENCE_RETRIEVED|23.478604|Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work; `delegated/vector/COSTS.json`|
|vector_unified_c18_append_transport|CONTROLLED_CANCEL_WAITING_ZERO_FW_RAW_FAILURE_RETAINED|0.210509|Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work; `delegated/vector/COSTS.json`|
|vector_unified_c18_append_transport_remote_pack|COMPLETE|0.031250|Distinct remote archive producer,zero physical solve;not local SSH CPU; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt1_transport|CONTROLLED_NEXT_ARM_MKDIR_DRAIN_RAW_FAILURE_RETAINED|1.133095|Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt1_transport_remote_pack|COMPLETE|0.343750|Distinct remote archive producer,zero physical solve;not local SSH CPU; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt2_transport|COMPLETE|11.176869|Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt2_transport_remote_pack|COMPLETE|8.531250|Distinct remote archive producer,zero physical solve;not local SSH CPU; `delegated/vector/COSTS.json`|
|vector_unified_overlay_DISPATCH_COUNT_METADATA_COST|COMPLETE|0.054102|Local engineering overlay transport; `delegated/vector/COSTS.json`|
|vector_unified_overlay_HOST_SPLIT_ENGINEERING_COST|COMPLETE|0.043008|Host manifest/CPU guards/archive-lifetime overlay transport,no new science/probe; `delegated/vector/COSTS.json`|
|vector_unified_overlay_LINALG_PROTOCOL_OVERLAY_COST|COMPLETE|0.042362|Local overlay transport only; remote extractionCPU not separately instrumented, no science/GPU yet; `delegated/vector/COSTS.json`|
|vector_unified_overlay_LINALG_TIMING_AMENDMENT_COST|COMPLETE|0.037423|Local amendment transport; no GPU/science executed by this command; `delegated/vector/COSTS.json`|
|vector_unified_overlay_MATERIAL_STORAGE_LIFETIME_COST|COMPLETE|0.041770|Future-source overlay transport only, no physical solve/GT; `delegated/vector/COSTS.json`|
|vector_unified_overlay_MEASUREMENT_OVERLAY_COST|COMPLETE|0.049431|Overlay transport childCPU not separately instrumented; NOT_MEASURED, no new science/GPU; `delegated/vector/COSTS.json`|
|vector_unified_overlay_ROOT_35_15_HOST_AMENDMENT_COST|COMPLETE|0.033249|Local engineering overlay transport; `delegated/vector/COSTS.json`|
|vector_unified_ordinary_DEV6|COMPLETE|4088.453125|Whole ordinary6 shared process timer;replaces all6nestedarmCPU; `delegated/vector/COSTS.json`|
|vector_unified_ordinary_DEV6_query_helpers|COMPLETE|50.890625|Distinct nvidia-smi helpers once; `delegated/vector/COSTS.json`|
|vector_unified_new_DEV6_including_one6LU_library_probe|COMPLETE|4292.234375|GetProcessTimes whole scientific child replaces process/arm CPU;originalfailedphase includes15sealed+failedarm,not15armplusphase; `delegated/vector/COSTS.json`|
|vector_unified_new_DEV6_including_one6LU_library_probe_query_helpers|COMPLETE|48.859375|Windows GetProcessTimes parent excludes nvidia-smi child CPU;add once; `delegated/vector/COSTS.json`|
|vector_unified_original_GPU25_failed_after15seals_plus315accepted|FAILED|11018.875000|GetProcessTimes whole scientific child replaces process/arm CPU;originalfailedphase includes15sealed+failedarm,not15armplusphase; `delegated/vector/COSTS.json`|
|vector_unified_original_GPU25_failed_after15seals_plus315accepted_query_helpers|COMPLETE|127.281250|Windows GetProcessTimes parent excludes nvidia-smi child CPU;add once; `delegated/vector/COSTS.json`|
|vector_unified_remaining_controller_own|COMPLETE|4.187500|Controller own only;phase children above separately counted; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt1_scientific|CONTROLLED_ZERO_FW_NEXT_ARM_DRAIN|651.671875|Additional whole external child CPU,not replacement of originalfailedphase;replaces ownnestedarmCPU; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt1_query_helpers|COMPLETE|10.640625|All successful and retried query-helper CPU returned by wrapper;raw timeout logs overlap and are not added again; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt1_controller_own|COMPLETE|0.093750|Own finite controller only;scientific child separate; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt1_preflight_helper|COMPLETE|0.015625|One standalone public resource snapshot helper,not scientific gatehelper; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt2_scientific|COMPLETE|5682.343750|Additional whole external child CPU,not replacement of originalfailedphase;replaces ownnestedarmCPU; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt2_query_helpers|COMPLETE|5.937500|All successful and retried query-helper CPU returned by wrapper;raw timeout logs overlap and are not added again; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt2_controller_own|COMPLETE|1.296875|Own finite controller only;scientific child separate; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_attempt2_preflight_helper|COMPLETE|0.015625|One standalone public resource snapshot helper,not scientific gatehelper; `delegated/vector/COSTS.json`|
|vector_unified_confirmation_generation|COMPLETE|119.206795|Whole measured generation process;private writer only,not scoring or inverse; `delegated/vector/COSTS.json`|
|vector_unified_mac_object_c11_matched12|COMPLETE|5969.315371|Macwait4 complete scientific process and ALL children;replacesnestedprocess,all5armCPU,RSShelperCPU; `delegated/vector/COSTS.json`|
|vector_unified_mac_confirmation_controller_own|CONTROLLED_IDLE_PARENT_DRAIN|4.438037|Own controller+its ps/vm_stat helpers only;scientificwait4 children separate; `delegated/vector/COSTS.json`|
|vector_unified_mac_object_c15_matched12|COMPLETE|6190.119862|Macwait4 complete scientific process and ALL children;replacesnestedprocess,all5armCPU,RSShelperCPU; `delegated/vector/COSTS.json`|
|vector_unified_mac_object_c14_matched12|COMPLETE|6315.360375|Macwait4 complete scientific process and ALL children;replacesnestedprocess,all5armCPU,RSShelperCPU; `delegated/vector/COSTS.json`|
|vector_unified_mac_confirmation_spooled_controller_own|COMPLETE|5.404159|Own controller+its ps/vm_stat helpers only;scientificwait4 children separate; `delegated/vector/COSTS.json`|
|vector_unified_mac_object_c18_cross18to12|COMPLETE|9480.803012|Macwait4 complete scientific process and ALL children;replacesnestedprocess,all5armCPU,RSShelperCPU; `delegated/vector/COSTS.json`|
|vector_unified_mac_object_c18_matched12|COMPLETE|9497.191632|Macwait4 complete scientific process and ALL children;replacesnestedprocess,all5armCPU,RSShelperCPU; `delegated/vector/COSTS.json`|
|vector_unified_mac_c18_return_controller_own|COMPLETE|11.235425|Own controller+its ps/vm_stat helpers only;scientificwait4 children separate; `delegated/vector/COSTS.json`|
|vector_unified_mac_handoff_helper_own|COMPLETE|0.852774|Handoff helper own+querychildren only;PHASE2_INCLUSIVE_WAIT4_VIEW overlaps allspooledcontroller andobjectcosts,neveradd; `delegated/vector/COSTS.json`|
|vector_unified_cancelled_waiting_c18_GPU_controller|CONTROLLED_CANCEL_WAITING_APPEND_ZERO_SCIENCE|0.234375|Whole terminated OWN waiting controller,0scientific child,measured GetProcessTimes;duplicatecopy in c18_gpu_append_controller excluded; `delegated/vector/COSTS.json`|
|vector_unified_c18_cancel_action|COMPLETE|0.640625|Separate cancellation action own+its powershell child,not targetcontroller; `delegated/vector/COSTS.json`|
|vector_unified_c18_return_mac_engineering_local_action|COMPLETE|0.124898|Own local action+SSH/SCPchildren,not remote actionCPU; `delegated/vector/COSTS.json`|
|vector_unified_gpu_recovery_controlled_drain_local_action|COMPLETE|0.089491|Own local action+SSH/SCPchildren,not remote actionCPU; `delegated/vector/COSTS.json`|
|vector_unified_remote_controlled_drain_sentinel_action|COMPLETE|0.125000|Distinct own metadataaction;not scientificworker;nextarm0FW provenfinallycounts; `delegated/vector/COSTS.json`|
|vector_unified_DEV_early_postseal_controller_own|COMPLETE|1.179619|Own finite local controller only; `delegated/vector/COSTS.json`|
|vector_unified_DEV_early_postseal_transport_children|COMPLETE|6.036545|Local SSH/SCP childCPU,not remote pack or scorer; `delegated/vector/COSTS.json`|
|vector_unified_DEV_score_and_all3Dplots|COMPLETE|14.214806|Externalwait4 scorer includes allplots/exports;replacesnesteddev_postseal_scores/SCORING_COST; `delegated/vector/COSTS.json`|
|vector_unified_DEV_early_remote_pack|COMPLETE|4.234375|Distinct remote archive CPU; `delegated/vector/COSTS.json`|
|vector_unified_DEV_scope_metadata_correction|COMPLETE|0.008083|Only postscoremetadatafield correction,no new GT/FW; `delegated/vector/COSTS.json`|
|vector_unified_fresh50_public_merge_monitor|COMPLETE|8.336681|Publiconly collectionwait/shape+box/copy validation,noGT,walloverlaps science; `delegated/vector/COSTS.json`|
|vector_unified_two_stopped_public_observers|COMPLETE|0.700000|pscentisecond ownCPU observations before controlledstop,lowerbound;onlymetadata monitors,0FW; `delegated/vector/COSTS.json`|
|vector_unified_public50_selection_and_registry_finalization|COMPLETE|0.476263|Onlypublicselections/attemptregistry/flattenedmetadata;originalvalues+receipts retained,zeroGT/FW; `delegated/vector/COSTS.json`|
|vector_unified_fresh50_score_and_all3Dplots|COMPLETE|39.323735|Externalwait4 scorer includingallplots;replacesnested SCORING_COST; `delegated/vector/COSTS.json`|
|vector_unified_fresh50_score_controller_own|COMPLETE|0.060576|Own finite postseal controller only,not childscorer CPU; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_confirmation_process|COMPLETE|4.640625 (view only)|Cold sharedsetup/basis/libraryprobe view already nested in scientificparent,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_dev_denominator_process|COMPLETE|21.140625 (view only)|Cold sharedsetup/basis/libraryprobe view already nested in scientificparent,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_gpu_recovery_attempt1_process|COMPLETE|7.953125 (view only)|Cold sharedsetup/basis/libraryprobe view already nested in scientificparent,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_gpu_recovery_attempt2_process|COMPLETE|7.328125 (view only)|Cold sharedsetup/basis/libraryprobe view already nested in scientificparent,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_c11_matched12|COMPLETE|3.334380 (view only)|Mac cold sharedsetup/basis is inside objectwait4,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_c14_matched12|COMPLETE|2.003413 (view only)|Mac cold sharedsetup/basis is inside objectwait4,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_c15_matched12|COMPLETE|1.903663 (view only)|Mac cold sharedsetup/basis is inside objectwait4,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_c18_cross18to12|COMPLETE|1.931523 (view only)|Mac cold sharedsetup/basis is inside objectwait4,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_nested_setup_c18_matched12|COMPLETE|1.866517 (view only)|Mac cold sharedsetup/basis is inside objectwait4,NOT additive; `delegated/vector/COSTS.json`|
|vector_unified_cost_compiler_all_recorded_runs|COMPLETE|1.918037|Each metadata-only compiler run ownCPU exactly once;history costs also appear here so entries sum matches platform totals; `delegated/vector/COSTS.json`|
|vector_unified_completed_CSV_analysis_and_REPORT_all_runs|COMPLETE|0.131020|Own CSVsummary/REPORTmd+HTML/COMPLETED_REPRODUCE generation;allrecordedruns once; `delegated/vector/COSTS.json`|
|vector_unified_completed_cost_finalizer_own_all_runs|COMPLETE|0.069276|Own wrapper only;child compiler already counted separately; `delegated/vector/COSTS.json`|

## Completion and remaining measurement gaps

48 progress snapshot: COMPLETED / 288 terminal / 0 pending.
3D cost producer: COMPLETED_VECTOR_ALL50_SCORED_AND_REPORTED_KNOWN_COSTS_WITH_EXPLICIT_NULL_GAPS; 50 / 50 synchronized terminal predictions. Platform-specific entry sums are reconciled with this producer. Historical CURRENT_RUN snapshots are retained in JSON and are not completion authority.

## Deduplication and gaps

- Never add nested arm/corrector/aux timers to their whole process.
- The unified vector producer replaces every legacy vector view for CPU and physical counts; its platform subtotals are reconciled during compilation.
- 576,2304,5184 LU and17496 iterative solves remain separate; never infer matrix order from legacy accumulated full_factorization_dimension.
- Count source/primal and adjoint aliases only once.
- Fixed-state merged50FW already includes ordinary repair10FW.
- Portfolio scenarios overlap underlying attempts; never add them to campaign total.
- GPU stream interval, wall time, active kernel time and CPU are not interchangeable.
- Mac scientific wait4 is inclusive of its child RSS helpers; Windows GetProcessTimes is own-process and its independent helper is separate.
- Mac spool/write/flush/close telemetry is nested; phase2 inclusive parent wait4 is an overlapping view, never added to per-object scientific children.
- Missing/running/unmeasured values stay null; displayed totals are lower bounds, not complete exact project cost.

- Two early CUDA audit receipts contained a template suggestion of external OS CPU measurement. The producer confirms none was captured: post-import CPU is retained as a lower bound, imports/exit stay unknown.
- DEV/48 transport scopes are divided by their authoritative upstream producer ledgers; this compiler never adds their raw transport logs again.
- Engineering6LU probe is included in the denominator shared process, and six ordinary vector arm receipts are inside one ordinary process.
- CPU for ordinary repair is a separate process timer, while its10FW are already merged in the fixed-state50FW count.
- Plot layout-refinement attempts are included; mathematical/scientific experiments are not rerun to fill accounting gaps.

Recompile after refreshing local producer ledgers:

```
python research/research_max_struct_r11/experiments/twod_compile_costs.py
python research/research_max_struct_r11/experiments/twod_independent_costs.py
python research/research_max_struct_r11/experiments/vector_compile_costs.py
python research/research_max_struct_r11/experiments/summarize_completed_costs.py
```

The compiler does not fetch remote results or stop/restart running queues. For exact scope, null reasons, raw auxiliary counts and source byte/mtime provenance, use COST_LEDGER.json.
