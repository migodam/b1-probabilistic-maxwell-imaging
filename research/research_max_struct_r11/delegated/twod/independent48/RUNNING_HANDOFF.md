# Running independent48 campaign

Root owns final science/GT scoring/publication/git. This worker only executes/records frozen2D study.

- Worktree: /Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r11_worktree.
- Remote: jayzh@100.121.97.39, D:/AI/MAX_STRUCT_R11_DEV_20261107, Python D:/python/python.exe, existing vendor/threadpoolctl.py. Key path is in twod_dev_transport.py; never print key bytes.
- 48 objects generated once, status GENERATED_SEPARATE_TRUTH, new seed202611071107,48LU, generator absoluteCPU1.4375s. All public arrays/tasks/frozen solver synced. No numeric confirmationGT access; private arrays remain remote.
- Six methods: CIE/Dual/Global/fullG19/freezeA/contractedT, same mature3000aux/max3000corrector, whole-child480CPU cap, checkpoints60/120/240/480, no padding. No SHA, new backend, driver install, or algorithm changes.
- Original queue baseline-first was safely STOP_DISPATCHed.5 completed tasks retained, no rerun. Its receipt: delegated/twod/independent48/queue/RECEIPT.json.
-283 unstarted tasks reordered by within-case method rotation case index modulo6. Engineering receipt DISPATCH_REORDER_RECEIPT.json; source configs/twod_independent48_interleaved_remaining_tasks.json. Scientific task dicts unchanged.
- Current finite queue: delegated/twod/independent48/interleaved_queue/RECEIPT.json, P8/P10 only. Foreground local tool session60743 keeps SSH connected; old session67078 ended with expected code1 from deliberate STOP_DISPATCH.
- Finite monitor local session84530 runs experiments/twod_independent_monitor.py, polls every60s, syncs12 new terminal tasks and compiles costs. Live local snapshot: delegated/twod/dev_transport/INDEPENDENT48_STATUS.json. Full trace/costs are preserved remotely and batch-synced.
- Monitor terminal sync originally includes sealed-byte-copy argument. Sync script now deliberately downgrades to PUBLIC ONLY unless ROOT_POSTSEAL_SCORER_APPROVED.json exists; root requested no GT read pending scorer review. No such approval file has been created. Do not decode private arrays or run scorer automatically.
- Scorer experiments/twod_independent_score.py: all288 unique terminal jobs gate, accept original stopped queue + current completed queue, manifest must be originalall288 once. Separate root authorization required. Includes χ8 region, NRMSE/IoU, fixedCPU and original-Phi LU frontier, stratified paired bootstrap10000 seed202611071108, final paid Y_prediction vs publicY relres/normσ/CN99, no new forward. CPU absolute/import/body starts explicit. Compile + synthetic metrics/data checks PASS. Not run on confirmationGT.
- Plot/report source experiments/twod_independent_plots.py and configs/twod_independent48_plots.json prepared, not executed. Eight6×7 pages, all objects/failures,0..20, reliable local fonts, paired metrics, actualCPU and fixedCPU/LU display curves; factual MD/HTML only. Root reviews interpretation.
- Resource policy: do not touch P0/2/4/6 other sessions. Latest probe had them98–100%; P8/P10 own2worker100%; P12/P14 reserved vector, no GPU use. Root says notify if at least2 of P0/2/4/6 naturally end/unbind and stay idle. Only then may future finite queue be adjusted to4 independent Pcores, preserving inflight and480/BLAS1/science. Never use Ecores/HT. Such a future phase needs aggregate priorqueue job handling and monitor/scorer queue-list update; current monitor hardcodes two queues.
- DEV6 attempts (onefull replay+5freezeA) all sealed, local RESULTS/exposed_dev. Full replay exactly equals R10publicprediction m/j/Q/Y. DEV actualchildCPU1157.796875s; ledger TWOD_DEV_COSTS.*. Root independently scored DEV; worker did not readGT.

## Later root instructions implemented

- Existing Mac scientific interpreter is exactly `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`; preserve its venv symlink entry, never resolve to base, no installs. Runtime config is configs/twod_local_runtime.json.
- Root requested R10 conditional policy secondary extension. Frozen registration is configs/twod_independent48_secondary_policy.json, UTC2026-11-07T00:25:14.390693Z. It explicitly follows root16-object scoring and was registered during ongoing48 primary execution, before any48GT read; it is NOT initialR11 primary preregistration. R10 original pipeline snapshot is under secondary/.
- experiments/twod_secondary_seal.py must run only after all288 unique terminal jobs, PUBLIC ONLY,48 CIE-first triggers and96 CIE_then_G19_PGP/CIE_then_DualCIE choices. Strict stored public noise_radius inequality; source originalPhi minimum; exact tie CIE. Reuses cap-eligible cold attempts; no new solves, no GT. Costs are actualCIE plus only triggered fallback full actualCPU/LU/RHS; alternative views shareCIE and are nonadditive, not sameCPU as one CIE480 attempt.
- Final sync source now invokes this PUBLIC seal using the existing science interpreter, with exact RUSAGE_CHILDREN local child cost. It does not decode/fetch private GT. Primary scorer now hard-requires the96 public seal before it can read48GT, and writes UTC public gate/start/finish timestamps.
- experiments/twod_secondary_score.py never readsGT: after primary scoring, it reindexes already-authorized primary CSV by the previously sealed source-task choices, checks seal timestamp precedes primaryGT gate, and reports factual metrics/deployment costs. Not executed.
- Primary scorer includes final sealed Y_prediction vs publicY relative residual, residualnorm/sigma, CN99pass with rawCN scale. Synthetic checksPASS, no confirmationGT access.
- Prepared plot script remains NOT_EXECUTED; root owns final interpretation. Root's new16 finite-noise audit must not alter48 methods or root selection.

## Engineeringv2 four-Pcore phase (latest)

- Other sessions' P0/2/4/6-bound scientific Python processes naturally exited. Sustained probes RESOURCE_UPDATE_4/5/6 retained; availableRAM19.55GB. P6 retained about18% background load, so root explicitly chose4 lanes and said not to pursue6 or migrate again.
- Own interleaved queue safely STOP_DISPATCHed and all inflight finished. Exactly68 primary tasks retained,220 unstarted. Sourcev2 twod_reorder_remaining_v2.py aggregates both prior disjoint queue receipts and preserves cyclic case/method order for unstarted only.
- Current actual queue: delegated/twod/independent48/expanded_queue/RECEIPT.json, CPU logical[0,2,8,10], BLAS1,480CPU perchild. Source manifest configs/twod_independent48_expanded_remaining_tasks.json. EXPANDED_DISPATCH_RECEIPT.json is definitive. All48 stillWindows, no E/HT/P12/P14/GPU.
- Foreground SSH tool session25565 keeps the current finite220 queue connected. Both old SSH sessions67078/60743 ended with expected STOP_DISPATCH return1, not scientific failure.
- Old LOCAL METADATA monitorPID21753 was stopped with SIGINT after exact executable/argv verification; no remote physics process was signaled. New monitorPID45819 / tool session97976 reads all THREE queues. Source monitor84530 is obsolete. This is bookkeeping reload only.
- Sync/compiler/scorer now preserve source_queue_receipt/cpu_logical; jobs remain unique across allthree phases. Source v1 snapshots are under engineering_v1/. Do not sum STATUS and COSTS totals; STATUS may lead snapshot. COSTS uses terminal GetProcessTimes and queue controller once each.
- Public secondaryseal automatically receives allthree queue receipts; primaryGT gate still requires96 choices sealed first. No48GT read.

### Definitive monitor correction

PID reuse/wrap made the first local monitor reload stop the newer own metadata instance instead of the older one. Both were metadata-only, no remote worker or scientific process was signaled. The remaining old metadata monitor was then explicitly stopped and a single fresh three-queue monitor started. Definitive current local monitor PID25290 / tool session46250 is recorded in ACTIVE_METADATA_MONITOR.json; use that file, not earlier PID/session notes. Current physics SSH remains25565. Last verified72/288 terminal,0 failures, four independent Pcore workers active. Old monitor sessions84530/97976 are closed; theirKeyboardInterrupts are bookkeeping only.

## Actual-end CPU fairness secondary (registered after root request)

The original primary has fixed60/120/240/480CPU caps and matched-min physicalLU, but had no pair-specific common actual-end spend. Root explicitly requested this, so configs/twod_independent48_actual_cpu_fairness.json was registered UTC2026-11-07T01:53:09.447871Z while48 stillGT-unread. It is a late-registered secondary, not original primary.

-8 pair types: CIE-single or Dual versus fullG19/freezeA/contractedT/Global,48 objects each,384 pairs. Public B=min(two actual full child GetProcessTimes). No new solve/FW and no GT budget/state choice.
- experiments/twod_actual_cpu_fairness_seal.py chooses originalPhi minimum among recorded cap-eligible paid incumbents atCPU≤B and copies their m into one public NPZ. Keeps actualT_end, actual last physical record, B-lastphysical gap, chosen point time/LU/RHS, nominal480 unused/overrun. A gap after last exported physical point can contain paid aux/computation/serialization; it is NOT free/unspentCPU. No claim both runs actually stopped/spent exactlyB.
- Final public sync automatically seals96 conditional views then384 actual-CPU pairs using the root-provided science interpreter; each standalone child has wait4/RUSAGE cost once. Main GT gate now hard-requires both public seals plus all288 terminal before any48GT read.
- Main authorized scorer grades the already-sealed fairness m snapshots after the same GT gate, separately outputs ACTUAL_CPU_FAIRNESS_SCORES.csv and ACTUAL_CPU_FAIRNESS_PAIRED_BOOTSTRAP.csv (10000 high/mixed stratified paired draws, frozen seed202611071108). Original six-method primary480 tables stay unchanged. ArbitraryB Y_prediction usually not retained; no newFW to fill it.
- None of the newGT scoring or seals has been executed early. Private package remains held until root authorization after allpublic seals.

## Append-only48 Global conditional selections (latest root addition)

Root explicitly requested a Global-Dual fallback control after the exposed16 attribution showed Global could explain much of G19-versus-CIE improvement. Original96 rule/config/master outputs remain unchanged. configs/twod_independent48_secondary_global_appendix.json registered UTC2026-11-07T01:58:54.651542Z is clearly a late-added secondary before any48GT, not primary preregistration. No other policies are added.

- New experiments/twod_secondary_global_seal.py reuses the EXACT48 public CIE trigger decisions from original96, and same attempt/minPhi/tie-CIE helpers; adds48 CIE_then_GlobalDual choices from the already-paid cold Global runs.
- Total conditional choices now144 (original96 + separate48 appendix); primary remains288 and actual-CPU fairness remains384 pairs. Final sync seals96, then48 Global appendix, then384 fairness. All postprocess children have disjoint wait4 cost receipts.
- Main GT gate and private-byte-copy hold now require all288 +144 +384 public choices before any48GT. No GT or new FW has been used.
- Secondary scorer reindexes all144 public choices from the eventual authorized primary scoredCSV; no extraGT decoding or root selection. Three views shareCIE prefixes and must not be summed as campaign cost, and conditional two-attempt480 allowances are not same actual cost as CIE alone.
