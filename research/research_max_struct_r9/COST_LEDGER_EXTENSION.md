# R9 extension cost union

Status: **FINAL_KNOWN_COST_LOWER_BOUND**. Generated 2026-11-06T06:57:06.989862+00:00. Read-only metadata snapshot; no GT/large-array access or new scientific solves. Existing COST_LEDGER.md and raw receipts remain unchanged.

CPU sums are retained consumption lower bounds by host platform. They are not cross-platform speedup estimates. Unknown CPU is null, never inferred as zero. Task-wall sums, overlapping queues, paused processes and GPU event times do not become a project makespan.

|Platform|Known CPU lower bound (s)|Charge records|Unmetered charges|
|---|---:|---:|---:|
|Mac_AppleSilicon|16755.944084|49|3|
|XINAN_Windows_i9|37353.515625|17|0|
|UNKNOWN_PLATFORM|null|1|1|

|Component|Authority|State|Known CPU lower bound (s)|
|---|---|---|---:|
|phaseI|canonical|COMPLETE_RECORDED_WORK|4353.102408|
|small_vector3D|canonical|COMPLETE_RECORDED_WORK|5631.278697|
|small_sphere_IGT|canonical|COMPLETE_RECORDED_WORK|2.603262|
|known_interface_Mie16|canonical_two_processes|COMPLETE_RECORDED_WORK|5.646354|
|classic2D|canonical|COMPLETE_RECORDED_WORK|36294.245629|
|original_electrical_size3D|producer_CM_union_plus_disjoint_IGT_processes|COMPLETE_RECORDED_WORK|6221.913316|
|original_residual_transport|single_metadata_process|COMPLETE_RECORDED_WORK|0.008087|
|original3D_physical_refinement|canonical_disjoint_process_rows|COMPLETE_RECORDED_WORK|1560.180511|
|publication_metadata|unmetered|COMPLETE_RECORDED_WORK|40.289916|
|extension_cost_compiler|own_invocation_history|COMPLETE_RECORDED_WORK|0.191529|

## Disjoint charge scopes

### phaseI

- `research/research_max_struct_r9/audit/ALL_COSTS.json` — Mac_AppleSilicon, CPU 4353.102408s; 113 Mac inverse receipts + diagnostics/generation/analysis/measured failures/compiler invocations, once; historical R8 excluded.
- **Not added again:** groups/inverse_receipts/compiler_invocations/standaloneweak/portfolio recipes; Already in canonical total or attributed standalone method recipes, not additional executed work.
- **Gap/limit:** {"kind": "import_failures_unmeasured_baseline_children", "attempts": 16, "CPU_seconds": null, "wall_seconds": null, "physical_calls": 0}
- **Gap/limit:** {"kind": "queue_launch_path_typo", "attempts": 1, "CPU_seconds": null, "wall_seconds": null, "physical_calls": 0, "source": "research/research_max_struct_r9/delegated/development_queue_launch_path_typo_FAILED.log"}
- **Gap/limit:** {"kind": "first_scoring_schema_error", "attempts": 1, "CPU_seconds": null, "wall_seconds": null, "physical_calls": 0, "source": "research/research_max_struct_r9/delegated/scoring/FAILED_SCORING_ATTEMPT_01.json"}
- **Gap/limit:** {"kind": "unmetered_orchestration_literature_seal_metadata_exports_plot_imports", "CPU_seconds": null, "wall_seconds": null, "note": "No complete CPU/LLM/token/network/power/money metering; never infer zero"}

### small_vector3D

- `research/research_max_struct_r9/RESULTS/vector_extension/PHASE2_COMPLETE.json` — Mac_AppleSilicon, CPU 3441.388072s; Producer deduplicated platform lower bound; includes whole workers and own measured ancillary costs/migration attempts.
- `research/research_max_struct_r9/RESULTS/vector_extension/PHASE2_COMPLETE.json` — XINAN_Windows_i9, CPU 2189.890625s; Producer deduplicated platform lower bound; includes whole workers and own measured ancillary costs/migration attempts.
- **Not added again:** scientific_workers_by_hardware, Mac_incremental_ancillary, XINAN_ancillary, rows, nested generation/Mie receipts; All already included in retained platform total; do not add nested setup/source costs again.
- **Gap/limit:** {"CPU_seconds": null, "reason": "Original already-running parent CPU has an explicit polling bound, not exact postmortem measurement."}
- **Gap/limit:** {"CPU_seconds": null, "reason": "Killed Mac control/metadata processes use immediately pre-signal observed CPU lower bounds; terminal unobserved CPU is not falsely zero."}
- **Gap/limit:** {"CPU_seconds": null, "reason": "Small metadata ps/upload/remote-extract helpers were not separately CPU-metered; transport record gives wall times. These are unmeasured overhead, not scientific zero-cost operations."}

### small_sphere_IGT

- `research/research_max_struct_r9/RESULTS/vector_extension/cell_integrated_mie16/COST.json` — XINAN_Windows_i9, CPU 2.218750s; Whole remote diagnostic process, checks/setup/oneLU/exports.
- `research/research_max_struct_r9/RESULTS/vector_extension/cell_integrated_mie16/COST.json` — Mac_AppleSilicon, CPU 0.384512s; Local metadata/figure/report only.
- **Not added again:** measured_new_CPU_seconds and nested kernel timers; Sum/subsets of the two platform charges.
- **Gap/limit:** {"CPU_seconds": null, "reason": "separate upload/extract/pack/scp helpers not CPUmetered; knownoneLU total includesallnumeric checks/setup/imports/export"}

### known_interface_Mie16

- `research/research_max_struct_r9/RESULTS/vector_extension/interface_mie16/COST.json` — XINAN_Windows_i9, CPU 2.468750s; One separate known-interface forward process; exactly1LU/1RHS; worker setup is included.
- `research/research_max_struct_r9/RESULTS/vector_extension/interface_mie16/COST.json` — XINAN_Windows_i9, CPU 2.375000s; One separate known-interface forward process; exactly1LU/1RHS; worker setup is included.
- `research/research_max_struct_r9/RESULTS/vector_extension/interface_mie16/COST.json` — XINAN_Windows_i9, CPU 0.468750s; Parentonlyimports/dispatch/seal/offlineMieassessment; workers separatelycharged.
- `research/research_max_struct_r9/RESULTS/vector_extension/interface_mie16/COST.json` — Mac_AppleSilicon, CPU 0.333854s; Local offline report/figure only.
- **Not added again:** known_incremental_CPU_seconds and per-worker operator setup/LUT timers; Sum/subsets of disjoint process and reporting charges.
- **Gap/limit:** {"CPU_seconds": null, "reason": "Small transfer/packing helpers not CPU-metered; not zero"}

### classic2D

- `research/research_max_struct_r9/RESULTS/classic_extension_xinan/PHASE2_COST_LEDGER.json` — Mac_AppleSilicon, CPU 6209.105004s; Canonical whole classic platform union; categories are nested, not added again.
- `research/research_max_struct_r9/RESULTS/classic_extension_xinan/PHASE2_COST_LEDGER.json` — XINAN_Windows_i9, CPU 30085.140625s; Canonical whole classic platform union; categories are nested, not added again.
- **Not added again:** actual_new_known_CPU_seconds plus all original classic task/partial/benchmark/transport receipts; Canonical Mac+Windows platform totals already cover all of them; GPU event is a nonadditive device diagnostic.
- **Gap/limit:** "interactive code authoring/inspection and py_compile CPU not instrumented; no physical solves"
- **Gap/limit:** "Mac final2 worker startup/output tails and stopped dispatcher shutdown exactCPU unknown"
- **Gap/limit:** "first remote mkdir/unpack/status helper CPU and shared SSH server overhead not instrumented"
- **Gap/limit:** "engineering syntax and optional psutil probe failures have0physics but unknown CPU"
- **Gap/limit:** "final report/canonical output serialization tail unharvested"

### original_electrical_size3D

- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.054789s; mainbody.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 2.151182s; whole probe.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 2.062030s; whole probe.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 475.440874s; mainbody; incomplete-source GMRES callbacks unknown.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 2224.278452s; mainbody;52qualified source states reused.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 106.921875s; mainbody including both timed solves/setup; GPUstate reused.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 35.984375s; whole process;53resumed verification +7new.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 76.562500s; whole process; no grid sweep.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 616.000000s; whole process; common precompute+wiring already included; possible in-flight Green remainder unknown.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 10.171875s; whole process; inner8.859CPU probe is included, not added again.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 87.484375s; Whole process; child timers included.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 8.968750s; Whole process; child timers included.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 2570.593750s; Whole process including imports/data prep/shared G/E/S/Vr/two Y fits/all9arms/endpoint validation/cache/I/O.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.220459s; Offline material scoring from stored predictions after all9seals; CPU includes imports.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 2.523994s; First offline rendering, hostCPU-axis convergence; retained after main axis clarified to wall; zero physics.
- `research/research_max_struct_r9/RESULTS/original_vector3d/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 2.494036s; Offline plots from all9sealed predictions/scores; common0..4scale; fixedz±.7 and max projection.
- **Not added again:** recorded_CPU_seconds_deduplicated, benchmark subtimers, shared/wiring inside616CPUpartial, layout mainbody, resumed old52source costs; None.
- **Gap/limit:** {"CPU_seconds": null, "reason": "First uncropped interrupted-source GMRES callback count unknown; its actual Green calls and all CPU retained"}
- **Gap/limit:** {"CPU_seconds": null, "reason": "Old strided4 inverse lacked an exported total Green counter;230complete auxiliary steps and120physical RHS observed; possible trailing partial-step Green count unknown, all process CPU charged"}

### original_residual_transport

- `research/research_max_struct_r9/RESULTS/original_residual_transport/COST.json` — Mac_AppleSilicon, CPU 0.008087s; One metadata-only residual-transport analysis, including imports and exports except receipt tail.

### original3D_physical_refinement

- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 1557.953125s; Wholeprocess imports/sharedCM G/E/S/airFDwiring/all3optimizationarms/rejectedtrials/cache/I/O; no GT.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 0.296875s; Own remote helper only; GPU child process CPU excluded.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.583233s; Local imports/sealguard/offline savedarray scoring/two figures/report; zero solves.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.366371s; Plot-only after allthree sealed/scored; no new metric or solve; attempt1 retained.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.528221s; Plot-only after allthree sealed/scored; no new metric or solve; attempt1 retained.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.005726s; Metadata-only final report interpretation and cost notes; no solve or rescore.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — XINAN_Windows_i9, CPU 0.015625s; Remote isolated public-code/data archive extraction only.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.056849s; Own local_controller_CPU_seconds only; remote worker/helper CPU excluded.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.095875s; Own local_transport_child_CPU_seconds only; remote worker/helper CPU excluded.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.031050s; Own local_controller_CPU_seconds only; remote worker/helper CPU excluded.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.157924s; Own local_transport_child_CPU_seconds only; remote worker/helper CPU excluded.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU null (unmeasured); Syntax-only AST checks and first connection CUDA availability/affinity probe before transport, whole CPU unmeasured; no scientific RHS.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU null (unmeasured); First failed offline guard plus dependent ledger call; no truth read, no new physical RHS.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU null (unmeasured); System python3 stdin encoding error before executing report append; corrected by explicit existing venv.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.010675s; Independent metadata compiler invocation 1; excludes physics/score/transport.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.019360s; Independent metadata compiler invocation 2; excludes physics/score/transport.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.009497s; Independent metadata compiler invocation 3; excludes physics/score/transport.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.021783s; Independent metadata compiler invocation 4; excludes physics/score/transport.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.009462s; Independent metadata compiler invocation 5; excludes physics/score/transport.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.008903s; Independent metadata compiler invocation 6; excludes physics/score/transport.
- `research/research_max_struct_r9/RESULTS/original_vector3d/physical_refinement/ACCOUNTING_FINAL.json` — Mac_AppleSilicon, CPU 0.009957s; Independent metadata compiler invocation 7; excludes physics/score/transport.
- **Not added again:** RESULTS/original_vector3d/physical_refinement/air/RECEIPT.json; Producer nested/attribution view; included_in=RUN_PROCESS_COST.json.
- **Not added again:** RESULTS/original_vector3d/physical_refinement/GlobalDual_stage2/RECEIPT.json; Producer nested/attribution view; included_in=RUN_PROCESS_COST.json.
- **Not added again:** RESULTS/original_vector3d/physical_refinement/AplusGlobal_stage2/RECEIPT.json; Producer nested/attribution view; included_in=RUN_PROCESS_COST.json.
- **Not added again:** RESULTS/original_vector3d/physical_refinement/GRADIENT_WIRING.json; Producer nested/attribution view; included_in=RUN_PROCESS_COST.json.
- **Gap/limit:** "Unmeasured connection-probe CPU and pre-export process/output tails remain unknown, not zero; no claim of total GPU elapsed from host CPU/wall."

### publication_metadata

- `research/research_max_struct_r9/audit/INTERIM_PUBLICATION.json` — UNKNOWN_PLATFORM, CPU null (unmeasured); Interim/final commit, review, readback and site helpers unmetered; measured archive/upload excluded; not zero.
- `research/research_max_struct_r9/audit/PACKAGING_RECEIPT.json` — Mac_AppleSilicon, CPU 34.977978s; Measured final lossless archive packaging process; no scientific solve; manifest/export tail not separately measured.
- `research/research_max_struct_r9/audit/RELEASE_UPLOAD_RECEIPT.json` — Mac_AppleSilicon, CPU 5.311938s; Upload controller plus gh child; archive creation excluded.
- **Gap/limit:** {"CPU_seconds": null, "reason": "Publication/web/network/LLM/energy/money not separately metered"}

### extension_cost_compiler

- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.012819s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.013822s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.011229s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.012741s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.011751s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.012340s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.011590s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.011897s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.013165s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.012006s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.013149s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.013256s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.012999s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.014084s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.
- `research/research_max_struct_r9/experiments/compile_extension_costs.py` — Mac_AppleSilicon, CPU 0.014681s; This metadata compiler pre-export CPU only; export/final-tail CPU not separately timed, no new solve.

## Counts and completion

Physical counters remain in dimension-tagged source views in audit/EXTENSION_COSTS.json. Scalar576, scalar4096, vector1536/3000/3408 and vector30³/60³ are never collapsed into one LU-equivalent speed metric. FFT calls, current Cholesky, source/adjoint RHS, zero-air RHS and partial GMRES counters retain producer semantics.

Canonical classic and original ledgers replace their fallback rows. A completed task overrides all of its earlier active CPU snapshots. Small3D29 platform totals already include ancillary/migration costs; IGT sphere control is additional. Original partial616CPU already contains its shared precompute/wiring; the whole layout receipt covers its inner probe. Parent CPU and child CPU are added only when explicitly disjoint.

Rebuild: `python experiments/compile_extension_costs.py` from this study directory (or use its repository-relative path). When all writers have ended, root can use `--final`; it still reports IN_PROGRESS if required canonical/diagnostic work is missing. The flag does not convert unknown costs into exact totals.
