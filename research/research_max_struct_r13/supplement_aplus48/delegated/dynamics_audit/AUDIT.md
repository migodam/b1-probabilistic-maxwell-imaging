# A+ supplement dynamics: source and artifact audit

Scope: source/metadata inspection only. No new GT, material-quality scoring, candidate selection, physical solve, replay, numerical experiment, or algorithm modification. No frozen src/RESULTS file was edited. Parent's clarified plan is authoritative: primary A+48 uses the exact copied C2 runner/full_current/runtime with **no observer overhead**; any trajectory replay is a separately paid, post-score OFFLINE_DIAGNOSTIC, with parent-selected cases. This audit selects no cases.

## 1. Frozen provenance and actual A+ path

The local archive is `runtime_confirmation_C2_final.zip` (109 entries). Direct source-text comparison, without hashes, found:

- `zip::src/full_current.py` (311 lines) exactly equals current `src/full_current.py`.
- `zip::runtime/field_current.py` (534 lines) exactly equals current runtime file.
- `zip::src/run_experiment.py` is the original **306-line** runner, not current 355-line runner. Current adds later GlobalDual/TRF/online controls and Mac/path amendments. The A+ supplement primary must use the C2 copy, not substitute current runner merely because its A+ branch appears similar.
- The C2 archive does **not** embed the inherited R7 A+ or R8 cie_fft modules. They are imported through `--legacy-repo`. The parent identifies original C2 deployment as `D:/AI/MAX_STRUCT_R13_CONFIRM_20261112/src`; this audit did not query remote files. Preserve the exact inherited source tree/environment in the new deployment receipt. Local inherited source is inspected below, not falsely certified as a newly compared remote binary.

Frozen C2 runner source references:

| Frozen entry | Relevant behavior |
|---|---|
| `zip::src/run_experiment.py:145–170` | scalar method routing, public noise-derived DPCS, R8 arm_config, explicit A+ Möbius readout, cold coefficient start, restricted family, original paid audits |
| `:151–166` | A+ uses R8 `arm_config('a_plus',...,[mf12,maxiter1500])`, original prior and onehop readout; `initialize_physical_current=False` |
| `:167–170` | original absolute process-CPU audit callback; original `run_full_current(...g19=False)` |
| `:247–253` | actual DR update, audit using pre-update constitutive cm, sparse projection-history append, final X/Q only after complete loop |
| `:268–276` | saved physical snapshots, full paid material pool, optional final auxiliary endpoint, physical/audit/event metadata |
| `:277–300` | original prefix selection and whole-process receipt |

Frozen mechanism protocol: total120 process-CPU seconds, absolute proposal fraction.6, common corrector remaining at most48 seconds, physical audit CPU interval2 seconds, fallback every50 updates, 3000 native stage updates, scale6, original box[.02,20], prior smooth=tv=1e-4, epsilon=.01, single BLAS thread. Existing source-qualified native A+ schedule is mf3/6/10/12 with500/500/500/1500 updates; beta_current=1 throughout, beta_material=1/.2/.2/.2; global PR+, scalar-current CIE transform, onehop material transform, two Möbius sweeps by default. Do not enable inherited_gtol, arc candidates, projected readout, sourcewise PR, or physical warm start. Sources: R6 `src/benchmark_solver.py:34–80`; R8 `src/benchmark_solver.py:13–18,29–33`; R7 `src/aplus_readout.py:103–166`.

A+ minimizes its fixed-current onehop material loss with the actual normalized prior. It is neither the constitutive projector PC nor a common physical-Maxwell descent update. Global coefficient PR belongs to the current step; the material step is four-color conditional MM with a scalar seed and exact-loss fallback. The dual-beta outer loop has no established single common physical objective.

## 2. Existing evidence that can be reused without replay

Read only JSON metadata, file inventories and NPZ **headers** (stdlib zip/npy-header inspection, no array payloads).

| Saved campaign/method | Existing coverage | What it cannot establish |
|---|---|---|
| confirmation_mechanism16 / A+ | 16 runs, all termination.proposal=PROPOSAL_PROCESS_CPU; all retain HISTORY, CANDIDATES.m, physical initial/raw/corrected/primary and sometimes secondary, DPCS file, aggregate auxiliary receipt | No per-update auxiliary J/Q/A, no native auxiliary history or PR gradients/directions; FINAL_AUXILIARY_STATE absent16/16. Cannot infer native schedule completion or recreate full continuation. |
| representation_v1 / A+ | 4 exposed runs, all proposal CPU interrupted; same artifact categories | No auxiliary endpoint in4/4; no saved PR/native history. Development protocol differs from confirmation audit cadence. |
| confirmation_core48 / DR | 48 runs, all FIXED_STEPS_COMPLETE; final auxiliary X/Q present48/48, paid material history and sparse projection audit records | No every-step iterator/shadows/directions; final X/Q is not the physical raw/primary state. A final state alone cannot reconstruct earlier dynamics. |
| development_v3_local / DR | One exposed DR run with final X/Q and sparse audit-gap history | Exposed development only; no dense trajectory from sparse records. |

For example, header-only inventory showed A+ CANDIDATES.m shape(806,576), physical snapshots j/Q shape(576,12), and no auxiliary endpoint; a core DR example had CANDIDATES.m(817,576) and final X/Q(576,12). These are shape examples, not case-quality selection.

A+ physical HISTORY proposal-extra keys are phase/global_iteration/stage/iteration/adapter_CPU/relative_residual/physical_cost/objective. **They do not contain A+ inner-loss/PR rows.** DR proposal extras are iteration/shadow_gap. PROJECTION_HISTORY stores iteration/CPU/shadow_gap/projector counters at paid audits only. Its saved shadow_gap is ||PA(2c−z)−c||, not ||z−c|| or the approximate mutual-AP defect.

Reusable all-case material-history diagnostics include the time-indexed incumbent, paid trial versus accepted/unknown status, material displacement between actual audited candidates, legal/bound activity, raw/corrector separation, original D/R/Phi/residual and prefix candidate indices. CANDIDATES rows map to HISTORY.candidate. Corrector candidates include line-search **trials**, not just optimizer accepted iterates. An A+ proposal candidate's m is its actual post-material-readout value, but the file is sampled by physical audits, not every outer update. No unsaved proposal states should be interpolated or called observed.

Physical snapshot j/Q values are Maxwell fields at the snapshot m. They are **not** the A+ auxiliary J/Q used to produce that m. Do not replace auxiliary states with those physical fields during failure explanations.

The inherited R8 loop returns native history/stage_states on normal return (`R8/src/cie_fft.py:350–353`), but the R13 runner does not serialize those return fields. Moreover, a BudgetEnd escapes before returning them. Retaining `result['history']` at the end alone would therefore miss all16 existing interrupted A+ runs.

## 3. Best instrumentation seams, no arithmetic change

### A. Existing public hooks: minimum implementation, partial coverage

R8 `cie_fft.py:320–322` already invokes iteration_callback with copies of post-update m, auxiliary J, coefficients A, stage/iteration and counters. R13 full_current captures caller_iteration at:127 and forwards after the periodic physical audit at:183–190. Existing stage_callback receives m/J/Q/A at stage end (R8:331–337).

Use these unchanged callbacks if only post-readout m/J/A and completed stage states are needed. They do not expose Q at every step, native history[-1], grad, direction, old_grad/old_direction or pre-readout m. Recomputing Q from GJ in the callback adds unneeded arithmetic/actions. An audit or budget exception in full_current occurs **before** caller_iteration and may suppress a completed update's snapshot. Stage callback follows physical stage validation, which may likewise be interrupted.

An A+ material_readout wrapper is another existing seam: input blocks/J/Q/m/beta_m/Os/Ls/rho/bounds/reg/config/cost and returned new_m/inner expose exactly fixed J/Q, pre/post m and inner loss/fallback metadata (R8:283–286; R7:103–166). Call the original readout **exactly once**, return its original objects, and never alter arrays/costs. Copy only into an independent observer buffer. This captures useful frozen-state material directions without new Maxwell/Green actions. It still does not expose coefficient PR arrays.

### B. Recommended complete outer-update hook in a diagnostic-only loop copy

Copy the inherited R8 _run_cie into a **new supplement diagnostic module**, keeping its update expressions and operation order. Add optional observer calls that merely copy already-existing arrays/metadata, with observer absent by default. Use the existing private function/globals binding pattern in full_current:36–47,253–256 to supply the diagnostic loop to a diagnostic adapter; do not monkeypatch shared modules or edit primary legacy files.

Recommended emission points:

1. **Static setup/stage entry:** after receiver JD/Vd construction(:169–170) and after stage J/Q(:213–214). Save actual loop JD/Vd, source metrics Od/Os/Ls, bounds/rho, stage betas, frequency ordering and A. Saving actual loop basis avoids assuming a separately saved DPCS SVD has identical phase conventions.
2. **Completed outer update:** immediately after native history append(:306–319), **before existing iteration_callback**(:320) and before the existing increments/checks. Save current A/J/Q, m_before_readout/m, grad/direction, alpha/pr, restart flag, history[-1], counters. At this point old_grad=grad and old_direction=direction(:296); no duplicate copies are needed. Tag that original callback/audit and PR step-counter increment have not yet returned. Do not increment frozen counters earlier.
3. **Stage boundary:** preserve actual old_coefficients frequency transport and PR reset(:195–196,221,338–339). Emit stage-end after existing code reaches it; stage reset is meaningful, not a solver failure.
4. **Interruption:** outer diagnostic owner flushes already-buffered snapshots in finally and records original exception/last emitted position. It must not synthesize a post-interruption update, finish a readout after cancellation, resume a previous PR chunk, or force a physical audit.

Minimal patch shape is in MINIMAL_PATCH_SHAPE.md. Parent owns implementation and decides any diagnostic update/CPU stopping protocol. No callback is proposed for the primary120CPU results.

### C. A+ inner-Möbius hook, only if needed for a fixed question

The input/output wrapper suffices for most outer failure dynamics. If parent needs why an MM cell moved, add an optional copy-only hook in a diagnostic copy of the readout:

- _seed: scalar seed, old/seed objective, chosen start (R7:13–23).
- Each color after the **existing** update(:144–150): sweep/color/cells, old and new cell material, q/c/b/h/s already computed(:138–143), newz, and objective candidate policy.
- If branch identity is needed, expose existing roots/legal/candidate/objective arrays and chosen index just before original return in _conditional_candidates(:89–100). Do not call eigvals again or change root tolerance/tie ordering.
- Final: completed sweeps, fallback flag, old/seed/final block loss and returned material(:154–166).

Do not turn on aplus_inherited_gtol merely to obtain a gradient. It invokes new gradient work and changes stopping(:121–129,152). Exact onehop value/gradient diagnostics can be computed later from a saved frozen block using the existing R4 derivative(:101–120), as separately paid analysis, not substituted into the original update.

## 4. Necessary arrays and exact time semantics

| Record | Required payload | Why |
|---|---|---|
| Static per run | public input/source paths and runtime environment; shape/ns/complex128; s/bounds/prior; actual JD/Vd; Od/Os/Ls/rho; actual blocks.indices/diag_g/goff; stage frequencies/beta schedule | Identifies operator, source metric, current family and actual material loss; no new SVD/Green construction needed solely for recording |
| A+ outer step | stage/local/global index; pre/post m; post-current J/Q; post-step A; grad/direction used; native PR_plus/line_alpha/restart/delta; native history row; operation counts and clocks | Separates current decrease, fixed-current material move, PR memory, physical audit and corrector |
| Optional current substep | A_before or actual change=direction*alpha already computed(:269), incoming old_grad/old_direction if snapshot is before assignment | Exact checkpoint phase; for a post-update checkpoint A+grad+direction+m encode next-step PR recurrence within stage |
| A+ inner color | scalar seed, accepted/fallback start; sweep/color; actual cells and before/new; existing q/c/b/h/s; optional already-computed root/candidate arrays | Identifies conditional MM surrogate/tie choice, not just aggregate objective |
| DR step | **pre** X/Q, cX/cQ/cm from existing PC; reflected A-shadow aX/aQ; **post** X/Q; last_A_info, indices/clocks/counters | Constitutive shadow and propagated reflector are distinct; next iterator differs from both |
| Failure/final | last completed snapshot position; original stop exception; stage/iteration; count of complete/incomplete records; full buffer/write cost and bytes | Interrupted state is explicit; missing states stay missing |

Retain physical candidate IDs as separate references, never relabel a J/Q auxiliary snapshot as a paid physical candidate. At DR iteration k, cm is from **pre-update** z_k; its physical audit is after z_(k+1) was computed. Save both indices. The final saved iterator is z_K, but the last usual audited shadow is PC(z_(K−1)); no extra PC call is warranted just to make names match.

DR's actual C projection and last_A_info are in unchanged runtime `field_current.py:388–425,382–385`; the multiplier key is **multiplier**, with pseudoinverse_limit and minimum_reachable_residual, not guessed mu/lambda. Capture immediately after the corresponding affine_ball call; later projections overwrite it. Do not switch to the separate runtime run_projection implementation(:428ff), whose timing/audit contract differs from the C2 inline runner.

## 5. Frozen-state constitutive directions: save inputs, compute later

For an A+ outer update, save actual fixed auxiliary J/Q plus incoming/outgoing m and exact inner-loss metadata. The observed material displacement m_after−m_before is its onehop/MM update direction. It is **not** PC(J/s,Q)−m_before and does not inherit the normal-ray event theorem.

For DR, actual pre-state X/Q and its already-computed cX/cQ/cm allow d=z−c, per-cell v/w, energies/correlation and angular competition diagnostics without another propagation solve. Existing `src/material_events.py:13–47` specifies the finite-threshold inputs/outputs and amplitude-normal/stationarity defects. Store source-coordinate complex arrays rather than only scalar norms so direction alignment is inspectable later. Material scale s must be applied consistently (A+ J is unscaled; DR X=J/s).

A generic DR iterator is not a point in A, so its PC-normal does not automatically define an approximate mutual-AP event pair. If a diagnostic needs a=PA(c) and c'=PC(a), that is an **additional declared projection** with its own paid cost. Do not silently add it to the actual DR update or substitute it for an observed state. Likewise, any receiver-adjoint or physical-gradient direction requires separately paid actions; copy-only recording does not magically provide them.

Onehop denominator/surrogate directions are not original Maxwell descent or recovery certificates. Parent owns the comparison math. No spectral parameter sweep is needed or recommended.

## 6. Cost and low-overhead storage

Copy-only hooks preserve numerical update expressions and ordering, **not zero runtime cost** or bitwise whole-trajectory equivalence under a process-CPU deadline. Allocations, serialization, cache/memory pressure and clocks change CPU buckets, audit eligibility and stopping. Therefore parent’s no-observer primary extension is essential; replay endpoints must not replace it.

Recommended diagnostic storage:

- append a compact native scalar/metadata row each completed update;
- keep full state checkpoints at a parent-frozen iteration schedule, plus last-completed snapshot and stage boundaries; if every outer state is required, state that explicitly and budget disk/memory accordingly;
- write complex128/float64 arrays in bounded uncompressed chunks or preallocated memmap rather than thousands of compressed NPZ calls in the hot loop; no float32 cast, new dependency or asynchronous observer queue;
- copy arrays when accepted into a buffer: J/Q/A/current change later, and readout current changes in place; no shallow references masquerading as history;
- flush in diagnostic finally without calling model.state, P.material, P.affine_ball, value(), gradient(), SVD or eigvals; offline derived quantities are separate;
- charge imports/setup, original arithmetic, native physical audits, copy/trace CPU, serialization CPU/wall, output bytes and peak-buffer memory. Observer time is **nested within** whole replay CPU; do not add it twice or subtract it to grant more primary computation.

Scale estimate, not a benchmark: N576, ns12 means one complex J+Q checkpoint is221,184 bytes. At3000 updates that pair alone is663,552,000 bytes. Across the full configured mf3/6/10/12 schedule, three coefficient streams A/grad/direction add664,704,000 bytes; pre/post m adds27,648,000 bytes. Total is about **1.356GB decimal per full dense3000-step replay**, before DR shadows, inner-color snapshots, metadata or overhead. ns6 roughly halves the complex part. This motivates bounded chunks or sparse declared checkpoints, not changing numerical precision.

A diagnostic copy may still be cut off by its own budget. Record completed update count and truncation; do not force3000 iterations silently. A longer separately authorized diagnostic run is a different compute contract, even if update arithmetic is identical.

## 7. Recommended handoff

1. Complete A+48 primary with exact C2 copied source/params/selection and original120CPU contract; existing16 subset runs are not automatically independent new objects or free replacement results.
2. Reuse paid material histories for all48 without new forward calls, clearly naming physical audit sampling and corrector trials.
3. Only after the parent fixes diagnostic cases/protocol, run separately billed source-copy observers; attach original primary provenance and never reselect primary from replay candidates.
4. Add complete-array/native-history coverage labels to the supplement, preserving the old archive gap. Offline replays do not retroactively fill the frozen original trajectory.
5. No source edit, run, new GT read or scientific category was performed by this worker.

Review CPU/wall for the source-and-metadata audit were not measured as a whole and remain unknown, not zero.

