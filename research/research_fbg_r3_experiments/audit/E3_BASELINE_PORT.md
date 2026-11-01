# E3 mature-baseline engineering port

Status on preparation: **PROTOTYPE WRITTEN**. Subsequently the parent explicitly authorized exactly one four-arm real576 diagnostic on the historically disclosed public observation `obs_9cffd133ec9f435d` (random12): native TRF, fixed beta .2 feedback, single CIE and dual CIE, each45s, lanes2, public m0=.1. Its status is **LIMITED_MATURE_BASELINE_DIAGNOSTIC; N1 E3 SCIENTIFIC GATE NOT_PASSED**. This does not authorize another batch or upgrade E1's scientific gate.

Owned new files are only `src/cie_tsom_cuda.py`, `experiments/e3_fullwave_imaging.py`, and this audit. Frozen pipeline/R2 sources remain read-only.

## Reused implementation and scope

The port first read `research_z_alg_cuda_pipeline/OPTIMIZED_OPERATOR_CATALOG.md`, its `REPRODUCE.md`, `src/device.py`, `src/executor.py`, and R2 `src/adaptive_solver.py`. It reuses existing `SharedOperators` and `ObjectiveLane`; no new Maxwell operator, CUDA factorization, Fourier propagation backend, precision, driver, daemon or scheduler is introduced. The unchanged CPU CIE readout and Fourier basis construction are imported from `src/cie_tsom_baseline.py`.

```python
result = run_cie_cuda(shared, Y, m0, config)
# m, candidates, history, cost, last_joint_m, status, receiver_rank,
# shared_setup_cost, host_result_wait_seconds
```

The current solver accepts the same explicit stage/rank/material/prior choices as the CPU solver. It deliberately retains sourcewise PR+ and exact quadratic line steps, the original projected Fourier coefficient coordinates, zero-padding continuation, direct shared real material readout, box-before-inverse mapping, fixed positive source weights and correctly transported single/dual beta. It does not replace the mature method with R2's 100-step prototype.

## Device/host partition

Public receiver SVD and projected Fourier columns are computed once on CPU per immutable operator/rank/mask schedule; their setup is separately timed and counted. The same roundoff-only singular-action cleanup as the CPU solver is applied in the original Fourier coordinates. The prepared buffers contain no Y, initializer or label. `G/E0/S`, projected columns, `G@basis`, `S@basis`, DPCS current, APCS coefficients, current gradients and PR steps stay on GPU. Each task owns an existing private ObjectiveLane stream and completion event; all tasks share one CUDA process/context and read-only operator upload.

Each alternating step downloads **one packed J/Q/diagnostic buffer**, invokes unchanged `_readout` on CPU, and uploads the complete real material. The common material prior is exactly the R4 smooth/TV prior with width 19.5, weights 1e-4/1e-4 and TV epsilon .01 unless the parent supplies the same changed public weights consistently. Its gradient is accumulated over the same horizontal and vertical edges. Full regularized material-block derivatives include the off-manifold `DP e` term because `_readout` is reused. This is a faithful hybrid partition, not an assertion that every part is on GPU.

For n=576 and p=12, the planned J/Q packed transfer is `4*n*p*8+16 = 221200` bytes per current step, plus `576*8=4608` material bytes returned to GPU. The real scalar diagnostic and data term are packed into the same buffer. Per-step transfers, original physical endpoint return sizes, source-weight download, pinned-result waits and CPU material evaluations are explicitly recorded. These are predicted transfer sizes, **not measured throughput**.

The port uses the original Fourier range with explicit matrix columns, not FFT acceleration. The requested masks mf=3→6→10 generate **36→144→400 Fourier coefficients per source** before rank deficiency from DPCS projection. These are coefficient counts, not guaranteed independent APCS ranks. Actual retained rank and real current capacity are recorded in stage histories; the material always has 576 independent real values.

## Full-physics candidate scoring and common corrector

Initial and stage endpoints are evaluated through the existing complex128 ObjectiveLane, using every original source and receiver and the same R4 prior. Only these physical scores select the current incumbent. The runner then calls unchanged R2 `solve(..., initial_pixels=...)` with a **full 24×24 material stage** and original-data corrector. No support or outside material is fixed.

The parent manifest controls whether that common corrector is data L-BFGS or native Full-J TRF. The current-proposal fraction defaults to .55, inherited as an engineering split from R2; all remaining time in the 45-second envelope goes to the corrector. The full cold setup is separately charged. A stage or atomic inner/physical operation may complete just beyond its polling cap; actual wall and completed operations remain visible, not truncated to the nominal 45 seconds. Incomplete current schedules must be labeled as budget-limited, not mature convergence.

For the authorized four-arm diagnostic the CIE corrector is frozen to **data L-BFGS**, matching R2's final corrector. There is no stronger TRF corrector granted to CIE. The runner now reproduces R2 main's second reference LU exactly and installs `shared.feedback_reference_j/q/lu/pivots`; total shared reference setup is **2 LU and24 source RHS**, not one and not free.

Task caps **max_lu=2000, max_rhs=60000** propagate into the current phase, common corrector and shared task context. All current endpoint validations count. A runner-local subclass of the frozen FeedbackLane reserves room for a final physical closure and calls the base physical closure with budget checks active, rather than using R2's budget-disabled `close` wrapper. This changes only cap enforcement, not objectives or updates. The conservative reserve is one final LU/source solve plus one next complete evaluation; no validation cap exemption is used. Own process affinity21845 and single CPU/BLAS threads are applied; other process affinity and priorities remain untouched.

The parent-requested schedule is beta=1→.2→.2, with optional later `beta_current=1` for the dual arm; it is an explicit local heuristic, not borrowed 3D published values. No additional beta or material-driven schedule tuning has been introduced. Default adequate-duration settings are coefficient_tol=1e-8, miniter=50, maxiter=500 per stage; the parent may freeze a different common budget point explicitly.

## Finite execution and data isolation

`e3_fullwave_imaging.py MANIFEST` reads metadata and prints **NOT_RUN**. It executes only with explicit `--execute` after the parent sets `resource_idle_confirmed` based on a fresh resource snapshot. `--inspect-resources` is read-only. A finite FIFO queue keeps at most 2 or 4 futures active; one completion immediately permits the next submission, and only the main thread writes receipts. There are no callbacks writing files or device-wide hot-loop barriers. An optional existing hybrid SVD slot may be enabled by `svd_backend: "hybrid"`; the port does not change its implementation.

Manifest essentials (parent supplies paths, finite task list and frozen configuration):

```json
{
  "pipeline_root": "D:/AI/ZALG_CUDA_PIPELINE_20261030",
  "r2_root": "D:/AI/FBG_R2_20261030",
  "operators_path": "PARENT_PUBLIC_OPERATORS.npz",
  "output_dir": "FRESH_PARENT_RESULT_DIRECTORY",
  "vendor_path": "EXISTING_THREADPOOLCTL_VENDOR_DIRECTORY",
  "lanes": 2,
  "resource_idle_confirmed": false,
  "config": {"wall_seconds": 45, "proposal_fraction": 0.55},
  "tasks": [
    {"task_id": "PARENT_ANONYMOUS_ID", "method": "cie_single", "public_path": "PARENT_PUBLIC_OBSERVATION.npz"}
  ]
}
```

Public operator loading opens only `G`, `E0`, `S`; observation loading opens only `Y`. Known historical labels, `m_disclosed`, private material arrays and GT metrics never enter solver selection or stopping. The original twelve R2 targets are already historically disclosed, so this can only be a disclosed-target diagnostic. Any later evaluator and scientific decision remain with the parent.

Physical counts are aggregated from the actual existing lane objects, including partial failed-task operations; source, adjoint, receiver and tangent RHS are summed into actual complex RHS. Shared setup, basis acquisition, CPU readout, H2D/D2H, failed attempts and wall time are separate. A process failure within a GPU primitive can still leave an interrupted operation incompletely counted; the failure receipt must not be interpreted as a successful zero-cost attempt.

## Checks performed and outstanding gate

- Python syntax compilation passed for both new modules.
- Runner `--help` executed locally; no CUDA context was initialized by that check.
- One read-only remote snapshot returned GPU utilization **0%**, memory **0/8188 MiB**, and no Python process entries at that instant. This is an instantaneous resource observation, not a reservation; check again before an authorized run.
- No E3 numerical benchmark, current trajectory, full-wave imaging task or GT evaluation was started during preparation. The later explicitly authorized run is isolated under `D:/AI/FBG_R3E_20261101`; complete receipts and exceptions are retained, and a fresh GPU/Python idle check precedes launch. No additional batch is authorized.

No speed recommendation between this hybrid port and a CPU-resident basis implementation can be made without the authorized one-whole-DEV-case measurement. The likely competing costs are CPU regularized material readout, packed J/Q transfer and medium-width complex products; no bottleneck is claimed measured. If the port cannot complete the parent-authorized DEV case faithfully within the bounded engineering window, the outcome is **NOT_RUN / implementation incomplete**, not a weak R2 substitute. Runtime API, CPU/CUDA trajectory agreement, end-to-end wall and final original physical residual remain the concrete acceptance items.

## Authorized launch record

The isolated deployment is `D:/AI/FBG_R3E_20261101`, reusing the frozen pipeline and R2 paths above. The first detached Windows `Start-Process` launch reported PID17772 but that child exited on SSH-session closure before producing any receipt. `attempt1.stdout.log` and `attempt1.stderr.log` were both zero-length, and no live Python process or task JSON remained. Its startup cost is **NOT_MEASURED**, not zero; no scientific endpoint is claimed. The identical manifest is then executed with the SSH connection retained, preserving separate attempt2 logs. This transport correction does not alter algorithms, data, starts, budgets or the authorized task count.

The source has been inspected at formula level. As specifically directed by the parent, no tiny CPU/GPU error-floor or separate surrogate microbenchmark is added; the authorized full576 original Maxwell endpoint validation and closed cost receipts supply the actual numerical execution evidence. This does not imply bitwise CPU/CUDA trajectory equality.

## Closed authorized diagnostic: actual results

**COMPLETED: 4/4 task receipts OK, retained remote Python exited0; no CUDA-port runtime exception and no scientific-setting repair.** The only failed launch was the empty detached transport attempt described above. No next batch was started or authorized. The scope remains **LIMITED_MATURE_BASELINE_DIAGNOSTIC**, with **N1 E3 scientific gate NOT_PASSED**.

All four tasks used the same disclosed R2 observation, public initializer .1, [.02,20] box, 12 complex sources,48 receivers, complex128 scalar-TM operators, common prior, nominal45-second task cap, max_lu2000/max_rhs60000, and two FIFO lanes. Native/feedback preserve frozen8→12→24 material stages; CIE current/material unknowns are full576 throughout and their corrector uses full24 only. These represent different algorithms under a shared cap, not matched actual costs.

| Method | Actual task wall/s | Selected original relative data residual | Original data+prior score | Actual LU | Actual complex RHS |
|---|---:|---:|---:|---:|---:|
| Native Full-J TRF | 9.88204 | 1.00901360 | .5090542391 |23 |1,380 |
| Fixed beta=.2 feedback + data corrector |39.27340 |.0002569203 |.000005495137 |1,999 |47,940 |
| CIE single beta + same data corrector |45.00696 |.0003107143 |.000005399318 |1,548 |37,104 |
| CIE transported dual beta + same data corrector |45.01383 |.0002956448 |.000005386313 |1,896 |45,444 |

The selected point minimizes **original data plus the common prior**, not the data residual alone. For example, a phase may have a smaller raw data residual yet lose on the public regularized score. These observations do not establish GT material accuracy or a reliable method ordering. Native stopped by its optimizer's gtol/ftol rules before consuming45s; no extra start was added to fill its unused time. Fixed feedback stopped at the paid-closure LU reserve before its wall cap; its highest actual LU count is1,999, within the2,000 cap. CIE stopped on the wall envelope during its original-data corrector; its measured few-millisecond overrun is disclosed rather than clipped.

### CIE current phase and subsequent physical correction

Both CIE arms traversed all three Fourier stages with shared prepared APCS ranks36→144→400. The DPCS rank and full native coefficient history are retained in raw receipts. A spatial complex mode has independent coefficients for12 sources; actual real current capacities are864→3,456→9,600, while material capacity remains576 real cells.

| Current arm | PR steps | Current-phase wall/s | Raw final original physical data residual | Current stage stop conditions |
|---|---:|---:|---:|---|
| Single beta |1,176 |22.59242 |.01650149 |MAXITER; COEFFICIENT_TOL; COEFFICIENT_TOL |
| Dual beta |1,284 |21.11957 |.02530480 |MAXITER; MAXITER; COEFFICIENT_TOL |

Thus the current-only endpoints are materially less data-compatible than the endpoints after the shared full-material data corrector. The original Maxwell solve, not an auxiliary current residual, certifies each listed physical residual. `SCHEDULE_COMPLETE` means every prescribed stage ran to its explicit cap or stop; the MAXITER stages are not silently relabeled converged. The final corrector yields roughly3e-4 relative data residual for both CIE variants. No claim is made that this proves CIE alone reconstructs the target or that the current stage dominates the common corrector.

The single/dual arms downloaded260,131,200 /284,020,800 current-buffer bytes and uploaded5,419,008 /5,916,672 material bytes, with1,176 /1,284 packed current returns. Material inner evaluations were5,880 /6,420. These are actual recorded operation/transfer counts; no isolated CPU-versus-GPU bottleneck or SM occupancy was measured.

### Shared setup and complete queue cost

The run's **shared setup wall was1.532859s**: existing shared operators and the two paid reference solves took1.119773s; receiver/Fourier preparation took.371149s; remaining measured setup includes observation preloading and wrapper work. Shared basis setup counted one receiver SVD, three APCS SVDs,580 explicit Fourier columns,580 G-column and580 S-column actions, and5,515,144 basis-upload bytes. The two shared LU solves used24 complex source RHS. These are in addition to task costs.

The two-lane finite queue completed in **86.115981s** from setup start through task receipt writing. It is not the sum of overlapping task wall times and excludes Python imports, SSH/SCP transfer and the failed detached launch's unmeasured startup. Summed completed physical cost including shared references is **5,468 LU /131,892 complex RHS**. LU/RHS caps were respected; wall polling overruns are reported above. Shared cold setup is explicitly separate and may be amortized only once across this four-task run.

### Evidence files and limits

Complete local evidence is under `results/e3_limited_diagnostic_20261101/`:

- `limited_digit2_attempt1/SUMMARY.json` and four task JSON files: full phases, all current histories, candidate material vectors, cost and stop reasons.
- `manifest_limited_digit2_attempt1.json`: exact authorized finite queue and public paths.
- `attempt1.*.log`: empty failed detached-launch logs; `attempt2.*.log`: successful retained-SSH run, no stderr.
- `executed_source_snapshot/`: copies of the exact owned current/runner/CPU-baseline sources deployed for this run.

No GT/material/support arrays, NRMSE or IoU were opened, and no target label entered the solver. The public observation was already historically disclosed; this is not a new blind or generalization result. N1 was not run here. Data-level success on this single target does not pass the N1 expansion gate, establish independent N1 gain, imply a posterior sample, or establish cross-object superiority. Another batch requires a new parent authorization.

## Bounded completion authorization after the pilot

After reviewing the pilot, the parent explicitly authorized **44 additional tasks**: the same four arms on all remaining11 historical amplitude12/clean objects, excluding the already-closed pilot. The completed pilot supplies the twelfth object and is not repeated. The purpose is to close the original task's missing mature-baseline comparison across all12 previously disclosed structures, rather than selecting only favorable pilot cases. **N1 E3 scientific gate remains NOT_PASSED; no N1 full-size run or new-method superiority claim is authorized.**

Selection uses only the already-disclosed `truth_after_freeze/mapping_PRIVATE.json` metadata fields `observation_id`, `amplitude`, `noise`, joined with `stage4_public/manifest.json` public paths. No truth file or GT material/support array is opened; shape labels never enter GPU tasks. There are exactly11 distinct new observation IDs and44 tasks. The new manifest deep-copies the pilot configuration per arm, changing only public observation/operator paths, anonymous task IDs and the new output directory. All budgets, stops, stage masks, betas, priors, initializers, proposal fraction, caps, lane count, precision and code remain fixed. No GT- or residual-based tuning occurred.

The finite queue uses the known retained foreground SSH transport after another idle-resource check, isolated output `D:/AI/FBG_R3E_20261101/results/e3_mature_baseline_remaining`, and separate `remaining.stdout.log` / `remaining.stderr.log`. The earlier pilot and failed-launch costs remain separate. This authorization ends after those44 tasks; no third beta, extra iteration budget, additional initials or subsequent expansion is permitted without a new parent decision. Final status and cost will be appended after all receipts close.

## Bounded completion: closed48-endpoint delivery

**COMPLETED: all44 new tasks returned OK; the retained process exited0. Combined with the four reused pilot endpoints, there are48 unique task receipts over12 distinct historical clean/amplitude12 observations and four methods.** No pilot repeat, runtime repair, new initializer, beta, extra stage or subsequent batch occurred. Remaining stderr is zero-length. The worker has stopped all numerical execution for this authorization.

Final evidence is under `results/e3_mature_baseline_remaining_20261101/`: complete `e3_mature_baseline_remaining/` task JSONs and SUMMARY, exact `manifest_mature_baseline_remaining.json`, and both remaining-run logs. The unchanged executed-source snapshot resides with the pilot evidence above; the same remote sources were reused without edits or another upload. The parent receives full576-value selected material matrices, all current candidates and all phase matrices in the closed receipts, for its independent scorer/images. No GT material/support array was opened by this worker.

### Combined data-level and cost summary

The following threshold is **only the original-data criterion** `relative_residual<=.001`. It is not NRMSE/IoU, not joint image success, and not a new scientific gate. All48 selected matrices are finite576-vectors in [.02,20]. The largest recorded selected original Maxwell state residual is approximately3.16e-15. All task LU/RHS caps hold; task wall polling overruns remain recorded.

| Method | Data criterion count/12 | Median original relative data residual | Worst original relative data residual | Mean actual wall/s | Mean actual LU | Mean actual complex RHS |
|---|---:|---:|---:|---:|---:|---:|
| Native Full-J TRF |0 |.8596483 |1.0090136 |24.7401 |99.333 |5,960 |
| Fixed beta=.2 + original-data corrector |8 |.0003638201 |1.0055630 |35.1143 |1,699.250 |40,752 |
| Single-beta CIE + same data corrector |9 |.0003671625 |.1811433 |45.0129 |1,327.167 |31,805 |
| Dual-beta CIE + same data corrector |10 |.0004149848 |.1674014 |45.0195 |985.417 |23,606 |

The unequal actual time/work and physical-prior ranking prevent interpreting this table as a clean universal speed or material-recovery ordering. In particular, the two CIE data-count improvements are not independent N1 gains; N1 was never run at full size here. All failure objects and poor-residual endpoints remain in the delivered cohort. The parent must apply the complete frozen material/data joint gate to the selected matrices before stating imaging success.

### Required maturity limitation: truncated current schedules

The fixed45-second envelope did **not** provide full mature convergence across objects. Current phase outcomes are:

| CIE arm | Fully executed prescribed schedule | Current phase WALL_BUDGET | Reached mf=10 | Total PR steps | Stage MAXITER events |
|---|---:|---:|---:|---:|---:|
| Single beta |3/12 |9/12 |12/12 |14,446 |21 |
| Dual beta |1/12 |11/12 |9/12 |13,294 |21 |

Even a fully executed schedule can include500-iteration MAXITER stages; it is not automatically converged. Three dual-beta cases did not reach the last400-coefficient Fourier stage within the shared proposal allocation. Thus this is a faithful mature-structure **budget point**, not the mature CIE performance ceiling or an adequate complete-schedule convergence study for every case. No extra budget was used to rescue this limitation. Claims that depend on mature convergence must remain **HOLD / NOT_ESTABLISHED** from this diagnostic.

Recorded current transfers across all24 CIE trajectories are6,136,088,000 packed J/Q/diagnostic D2H bytes and127,825,920 material H2D bytes, in27,740 alternating steps. Shared-basis/physical-validation transfers remain separately counted in each receipt. No separate CPU-resident baseline timing was run, so no CPU-versus-CUDA speedup or isolated bottleneck claim follows.

### Complete paid costs and stop

The remaining44-task two-lane queue took **844.163249s**, including **1.943975s shared setup**. It created separate uniform12 and random12 shared operator caches: each paid2 LU/24 source RHS, one receiver SVD, three APCS SVDs,580 explicit Fourier columns,580 G and580 S basis-column actions. Both caches' setup and transfers are retained; they are not amortized away against the separate pilot run.

Pilot plus remaining queue wall is **930.279230s** (about15.5min), including **3.476834s measured shared setup**. Sum of actual task costs is49,334 LU/1,225,476 complex RHS; adding all three paid operator-cache reference setups yields **49,340 LU /1,225,548 complex RHS**. These totals exclude SSH/import/transfer overhead and the first detached launch's explicitly unmeasured startup, and do not sum overlapping task walls as elapsed queue time.

The exact bounded run is now closed. **N1 E3 scientific gate remains NOT_PASSED.** No new-method superiority, global convergence, blind generalization, posterior or material-accuracy claim is licensed by these data-only results. Any next numerical action needs a new parent authorization.

## Final post-closure score/coverage and comparison caveats

The later independent scorer `results/e3_all_scores.json` is now `POST_CLOSURE_SCORED` and records GT access only after the finite prediction queues closed. This supersedes the earlier data-only limitation **only for finite-set material metrics**; it does not change any online method, endpoint, gate or budget. This worker's final check read the closed score table and scorer formula, not raw GT arrays, and performed no new solve.

The frozen joint gate is NRMSE≤.15 normalized by this cohort's amplitude12, reconstruction-derived Otsu-IoU≥.7 and original physical relative data residual≤.001. Exact selected-incumbent results are:

- **Dual CIE10/12 joint and10/12 data:** the nine specific digit1–digit9 objects plus two_islands pass; ring and curved_bar fail.
- **Single CIE8/12 joint and9/12 data:** digit1,2,3,4,6,8,9 and two_islands pass. Digit5 has data/NRMSE pass but Otsu-IoU.67826<.7; digit7, ring and curved_bar fail. Do not report9/12 joint.
- **Fixed beta=.2 feedback8/12 joint/data:** digit2,3,5,7,8,9, ring and two_islands pass; digit1,4,6 and curved_bar fail.
- **Native Full-J TRF0/12 joint/data.**

These are the same48 online-selected material endpoints, not oracle best@K or newly optimized images. All576 final material variables remain free, with no support/shape restriction, target-label input or known outside material. Native/feedback earlier coarse stages are public initialization only; all twelve of each reach and select stage24. CIE material and full24 corrector expose all576 cells throughout. Stored GT support appears only in independent closed evaluation, and the GT/Otsu failure of single-CIE digit5 has not been repaired by changing its gate.

The corrector is original-data L-BFGS from the frozen R2 path, with common prior, not a stronger TRF corrector granted to CIE. The strong-baseline result belongs to the entire current-plus-corrector pipeline. Most current schedules remain wall-truncated (single9/12; dual11/12), and dual reaches mf10 only9/12. Therefore these exact nine-digit results do **not** establish mature CIE convergence or a performance ceiling, nor independently validate N1. The scientific N1 E3 gate is still NOT_PASSED.

The directly read2019 primary source supports the implemented receiver-SVD / projected Fourier decomposition, normalized residuals, fixed-material current gradient and one-step PR update mathematics, material alternation, staged Fourier release and the later-stage different-beta variant. Single and dual beta must remain separate arms: dual uses two state metrics after transporting through the same real m, not a single-beta objective with a reused numeric contrast. PR+ and a non-descent reset, restricted real material, the common-prior material inner solve, the local beta/mask schedule and original-data correction are disclosed adaptations. Explicit dense IDFT columns and GPU matrix products preserve the current subspace; they do **not** confer a literal FFT/CG-FFT complexity claim or reproduce2019's3D vector/complex-material acquisition.

No policy/provider/model setting, physics backend, precision, acquisition, prior, initial material or frozen solver schedule changed between the pilot and completion batch. The runner-local cap wrapper only enforces paid closure within LU/RHS caps; it is disclosed rather than claimed bitwise equivalent to the older cap-exempt wrapper. There is no hidden GT restriction. There are, however, resource-comparison limits: two FIFO lanes do not allocate fixed SM shares; co-running methods, non-randomized arm/order, warm-up and CPU inner solves affect wall-limited step counts. Early optimizer termination was not replenished with new starts. Equal45-second ceilings are not equal actual wall, LU/RHS, current capacity or exploration opportunities. Shared setup/caches are paid and reported separately, not treated as free per-arm compute. No isolated throughput, occupancy, randomized resource-allocation or CPU-versus-GPU causal speedup study was run.

Two current scorer metadata gaps must not erase the raw provenance: its generic NRMSE-definition string says “8 or18” although E3 rows and formula use amplitude12; and summaries list `stopping_counts=NOT_RECORDED`, `physical_caps_recorded_runs=0` / zero aggregate cap stops despite budget stops in raw phase receipts. These are reporting-adapter omissions, not evidence that caps were absent or unused. Use raw receipts and this audit for stop/cap conclusions. The parent's separate native-TSOM batch remains separately sourced until its own closure/scoring record is incorporated; no in-progress result enters this four-arm table.
