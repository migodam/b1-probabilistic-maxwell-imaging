# R12 paired supplement — implementation review

Read-only code audit, 2026-11-12. Only this review was written by the audit worker. No inversion, optimization, GT access, source changes, or scientific scoring was performed. The parent owns scientific interpretation and the final supplement contract.

## Scope and approved contracts

- A: three already saved bad endpoints × G19/DR/AP/direct; each child has 240 whole-process CPU seconds. Proposal methods use at most 144 proposal + 96 corrector; direct uses 240.
- B: retain original primary G19 result for all 24 confirmation cases. Restart only endpoints outside the public global CN99 noise ball, using all four methods with 120 whole-process CPU seconds, split 72/48 for proposal methods. Retain already-inside endpoints unchanged. This audit did not inspect which cases are inside/outside.
- C: four existing NN materials × the same four methods, budgets as A. Reuse frozen public NN initializations; inference provenance/cost stays separately attributable.
- Every method includes its common material m0 as a physical incumbent. Online selection has no GT access. Save full paid physical prefixes, proposal/raw and corrected selectors, and the registered secondary selector.

## Exact reusable APIs

Paths below are relative to the repository worktree.

| Purpose | Exact implementation | Required wiring |
|---|---|---|
| Full scalar-TM physical model | `research/research_fbg_r3_experiments/src/physics.py:16`, `MaxwellModel(public)` | Public G/E0/S/Y; original full source count and grid. `state(m)` returns j/Q/Y and caches exactly the last material. `physical_score(m,Y)` is original normalized data loss plus original regularizer. `data_value_gradient(m,Y)` includes original prior and one all-source adjoint. |
| Corrector | `research/research_max_struct_r8/src/benchmark_solver.py:11` reexports R6 `refine_observed`; body at `research/research_max_struct_r6/src/benchmark_solver.py:164` | Original L-BFGS-B, scale 1e6, box [.02,20], maxfun=4×maximum, ftol=1e-12, gtol=.01, maxls=30. Use maximum=3000 for proposal-method correctors; historical direct allowance is 6000. Runtime cap remains authoritative. |
| Paid physical observation | `research/research_max_struct_r6/src/benchmark_solver.py:99`, `PaidStateObserver` | `state` captures auxiliary validations; `data_value_gradient` copies model._current after the paid adjoint, including line-search trial states. Preserve them all, not just accepted iterations or incumbent improvements. |
| Mature native CIE loop | `research/research_max_struct_r8/src/cie_fft.py:21`, `run_cie_fft`, `_run_cie:71` | Same PR current steps, source metrics, DPCS rule, native Fourier representation, readout and stages. A custom readout is used only through the onehop/a_plus branch at lines 284–286. |
| G19 configuration | `research/research_max_struct_r8/src/benchmark_solver.py:13`, `arm_config('a_plus',...)` | The confirmation manifest uses **arm=a_plus**, algorithm_label=fullG19, not global_dual. `research/research_max_struct_r12/delegated/confirmation/data_public/confirm24/baseline_tasks.json:73` registers the exact custom hooks. |
| G19 readout | `research/research_max_struct_r11/src/material_readouts.py:97`, `prepare_full`; line 101 alias `full_material_readout=matrix_material_readout` | `prepare_full` delegates to R10 `prepare_green_strong`; readout is R10 `matrix_material_readout` with exact moving-denominator derivative, not lagged/freeze-A/freeze-T. |
| G19 fixed denominator setup | `research/research_max_struct_r10/src/lowrank_denominator.py:113`, `_prepare`, `prepare_green_strong:128` | Fixed public G/S basis, threshold .1, green_strong, Galerkin projection, matrix_inner_maxiter=8. Setup writes DENOMINATOR_SETUP.json and DENOMINATOR_BASIS.npz next to `receipt_path`; that parent directory must exist. |
| AP/DR | `research/research_max_struct_r12/src/field_current.py:125`, `DenseFieldCurrentProjectors`; `run_projection:428` | scale=6, bounds [.02,20], components=1, grouped_ports=True, full G/E/S/Y, raw global CN99 delta; omega=1; existing 3000-step/50-audit schedule unless parent explicitly freezes another schedule. |
| Whole-child cap mechanics | `research/research_max_struct_r11/src/twod_budget.py`, `WholeChildBudget` | Uses absolute `time.process_time()`; imports/setup are paid. Its observer gives a good template for completed-state eligibility and atomic overruns, but needs the supplement's generic per-phase deadlines. |
| Public-only input validation | `research/research_max_struct_r8/experiments/run_case.py:44`, `load_public` | Rejects private/sealed_truth paths and truth-like keys, checks conflicting public arrays, maps E to E0. Do not reuse unrestricted `arrays()` for arbitrary manifest inputs. |
| DPCS/noise contract | `research/research_max_struct_r10/src/twod_adapter.py:77,110`, `require_iid_raw_noise`, `noise_ball_rank` | Use original raw iid-complex sigma and public frozen noise-ball rank. Replacing dpcs_rtol with resulting dpcs_rank reproduces the confirmation adapter. |

### Exact G19 numerical configuration

Base config: `research/research_max_struct_r12/delegated/confirmation/DEV6_BASELINE_SOLVER.json`. Its common settings are host=.1 and box=[.02,20], with public noise-ball DPCS rank replacing its fallback dpcs_rtol=.001. Pass continuation mf12/maxiter1500 into `arm_config('a_plus',...)`.

Inherited stages are mf3/beta_m1/beta_j1/maxiter500, mf6 and mf10/beta_m.2/beta_j1/maxiter500, then mf12/beta_m.2/beta_j1/maxiter1500. Current PR is global/non-sourcewise. `arm_config` enforces coefficient_tol=0, miniter=50, scalar material inner maxiter60 and unchanged source normalization. `prepare_full` subsequently sets the G19 matrix-inner maxiter to 8.

Set `fft_geometry={coords:public coords or centers, radius, k}`, `regularizer=lambda m:model.reg.evaluate(m)[:2]`, `aplus_regularizer=model.reg`, `physical_sigma=norm(Y)`, and an output-owned `receipt_path`. Apply `prepare_full`, install `full_material_readout`, and call `run_cie_fft` directly. Calling generic `run_auxiliary` afterward would replace the custom readout with default A+ Mobius/projected logic.

## Same-material G19 restart and initialization defect

The inherited R8 code already supports `initialize_physical_current=True`. At its initial physical validation it retains j0=J(m0). Its receiver SVD computes JD from the fixed DPCS fit, then the first mf3 stage solves

`A0 = lstsq(UU, j0-JD, rcond=1e-12)`, where `UU=NativeFourierAPCS(shape,mf,Vd).apply(I)`.

Thus the initial auxiliary current is `Jaux=JD+UU A0`, in precisely the unchanged historical affine family. Later nested masks transport the same coefficients. Material m0 itself remains the first paid physical candidate even when its full physical current is not representable. This is a legitimate inherited R8 projected-current restart; label the supplement's choice explicitly because the old confirmation tasks defaulted this flag to false.

The mismatch is generally nonzero for two reasons: fixed DPCS is set by Y rather than the possibly outside-ball j0, and the first Fourier complement is truncated. Even a full Fourier mask cannot remove the fixed DPCS mismatch. Record `norm(UU A0-(j0-JD))/norm(j0)`, absolute defect, initial_projection_rank and columns. Do not claim that a material start implies exact auxiliary physical-field retention.

Prefer capture at the already existing least-squares call rather than a second receiver SVD or a second basis fit. An observation-only, module-local NumPy/linalg proxy can capture this one `lstsq` result and its residual without changing any numerical return. Limit its scope to the inherited cie_fft module and restore it in `finally`; avoid replacing global numpy.linalg.lstsq. Any residual calculation/copies are paid. If using an R12-owned copied solver instead, verify numerical identity apart from the observation callback and record that provenance.

Exclude `current_anchor='physical_complement'`: the R8 source itself calls it development-only and it translates the affine family. It is not the same inherited G19 control.

## Budget and candidate hazards

1. **Do not clone basin_v1 timing literally.** `experiments/basin_escape.py` starts its timer after imports and input loads, uses the same proposal/corrector split for direct, performs unconditional final audits and gradient evaluations after the cap, and optimizes only callback/last materials into its saved incumbent. It is an implementation reference, not the supplement's runner.
2. **Do not clone neural timing literally.** `delegated/neural/run_nn_interface.py` prepares the four intended public starts and records all corrector function evaluations, but its historical proposals are uncapped and its 120-second cap is corrector-only. Its `execute` also starts its timer after imports. Reuse data/provenance and physical hooks, not those budgets.
3. **Do not call confirmation run_baseline as an unmodified generic child.** It hardcodes 144/96/240 repeatedly, restricts output ownership to confirmation RESULTS, and assumes m_init for initialization files. A generic R12-owned adapter should parameterize T/P/C and support the basin m key explicitly; do not modify archived code.
4. Use absolute child CPU deadlines P=144/72. After proposal ends at CPU s, corrector deadline is `min(T,s+C)`, so an early proposal finish does not enlarge the correction allowance. For direct, skip projector/G19/DPCS setup and run original physical L-BFGS-B up to T; use no artificial 144/96 or 72/48 restart of direct optimization.
5. `run_projection` audits m0 unconditionally even when remaining proposal budget is zero. It checks budgets only at iteration boundaries, so projector construction and a scheduled physical audit may overrun. Check absolute deadline before constructing it; pass `max(0,P-time.process_time())` after setup. If setup consumed P, skip proposal and continue from the already retained m0 with the remaining corrector allowance.
6. G19 prepare/readouts are atomic relative to outer cap hooks. Check before/after preparation and at native iteration callbacks. If an atomic step crosses P, retain accounting but exclude newly completed over-cap proposal states. An incomplete stage has no newly paid physical endpoint; do not silently validate its last auxiliary m for free. Original stage physical validations are the legitimate physical candidate schedule.
7. The old observer updates its unrestricted best before eligibility is known. Use separate all-paid bookkeeping and eligible selector state; never seed a corrector from an over-P state, and never let an over-T state replace final predictions. Preserve excluded rows with reason/deadline and charge the overrun.
8. Every paid full physical state, including all direct/corrector line-search trials, belongs in the paid ledger with m, Phi, prior R, raw residual, phase, completed absolute CPU/wall, costs, eligibility and provenance. Native G19 auxiliary iterations are not paid physical states unless a full solve was performed. Save auxiliary raw material/iteration records separately without assigning a physical residual to their nonphysical J/Q.
9. The physical model cache avoids an additional LU for repeated m0 auditing, but repeated recording, regularizer evaluation and serialization still cost CPU. Prefer one common physical m0 record, then cache-backed initial callbacks where required by inherited code. Export retained arrays; no finalization solve, post-cap gradient, or selected-state reevaluation.
10. Secondary CN99 selection uses an already computed full-source residual, minimizes original R and breaks exact ties by original Phi; it never feeds back into iterations. Keep proposal/raw and corrected pools explicit and retain m0 in both relevant selectors. No feasible candidate yields an explicit absent status, not a fallback mislabeled feasible result.
11. B's inside/outside partition must use the original primary G19 material and raw residual with the same public Y/sigma, not the earlier secondary selection or GT. All methods start from the same saved material on an outside case. Link original cost separately from additional restart cost; do not overwrite the original prediction. For an inside case retain the original exact arrays and mark restart not requested.

## Minimal runner wiring

Use one generic child per (cohort,case,method), set one BLAS thread before numeric imports, enforce new-output-only ownership, and load only allowlisted public operators/Y and the explicitly supplied m or m_init initializer. Start absolute CPU accounting at process birth. Build the unchanged MaxwellModel and record the shared full physical m0 before any proposal.

Create a supplement observer that copies all already paid states and stores separate unrestricted, phase-eligible, proposal/raw, corrected and secondary selectors. For DR/AP, return `PhysicalAudit` from the same observer; use checkpoint_callback for incremental metadata/material persistence because an outer budget exception can prevent run_projection from returning its ProjectionRun. For G19, use the exact a_plus/custom-readout dispatch above, install a native budget check and an accepted-iteration observation callback, and retain stage validations through observer.state. The paid candidate list—not the solver's unfiltered best—chooses the corrector seed.

For each proposal method, hand the best eligible original-Phi proposal state (including m0) to `refine_observed`. Record all corrector evaluations; catch a dedicated budget control-flow exception outside the inherited solver. A BaseException subtype, as in WholeChildBudget, prevents broad inherited `except Exception` blocks from mislabeling the intentional global cap. If using a normal Exception for proposal handoff, reclassify its inherited FAILED execution receipt explicitly, as run_baseline does.

At exit, serialize already retained raw/corrected/secondary arrays, the full paid-material prefix and metrics, excluded atomic states, projector or G19 operation counters, phase CPU deadlines and final whole-child costs. Serialization and failed setup are part of actual CPU/wall costs; do not add nested timers to the whole-child timer. Empty/missing eligible-state outcomes require an explicit receipt. No GT-bearing path is needed in this child.

## Required verification only

No numerical run was performed for this audit. The following narrowly targeted implementation checks are justified before launch; no new error-floor study, parameter search, or GT test is required.

- Verify manifest invariants without loading truth: exactly A=12 and C=16 tasks; B covers all24 original primary endpoints, four children per outside endpoint and no children per inside endpoint; method budgets and initializer/public paths agree pairwise. All children preserve full E/S/Y dimensions and original prior parameters.
- Test budget/selector state transitions with synthetic paid rows and a fake clock: setup exhausts P, proposal atomic overrun, corrector atomic overrun, early proposal finish, no feasible secondary, inside-ball B identity retention. Confirm direct has the entire T and m0 remains selectable. This tests control flow, not Maxwell accuracy.
- On the first already-authorized real supplement tasks, inspect instrumentation invariants rather than launching a separate physics experiment: each objective/gradient row matches the material in model._current; all output predictions originate in eligible paid rows; no final forward is issued; G19 readout is full_material_readout, arm a_plus, current sourcewise false, and initialize_physical_current true; projection defect is measured from the existing LS call; physical_complement is absent.
- Reconcile serialized prefix lengths/timestamps and selected indices with the ledger after each completed child. Raw/corrected/secondary records must be reconstructible from paid m/Phi/R/residual/CPU alone. Review one normal completion and any actual cap-overrun receipt before batch continuation.

## Uncertainties left to parent

- The parent has proposed the legal projected-current G19 flag true; this audit agrees with its implementation legitimacy but does not claim it preserves complete j0 or proves fairer scientific performance.
- The original B primary residual partition, endpoint arrays, and actual scientific outcomes were deliberately not inspected here.
- Raw selector nomenclature must distinguish native auxiliary material snapshots from full physical paid candidates. If the desired raw pool means a new physical validation cadence inside G19, that is an additional explicitly charged algorithmic protocol and must be frozen before execution.

## Actual runner pre-dispatch review

Reviewed `supplement_paired/src/run_paired.py` read-only after implementation. No optimization, physical forward, GT read, or runner edits were performed by this reviewer. No decisive numerical wiring blocker was found in the reviewed source.

- Full public E0/S/Y are loaded unchanged into MaxwellModel; the inherited model handles all source RHS. The original regularizer and physical normalization are retained.
- G19 calls `arm_config('a_plus', ...)`, so PR is global/non-sourcewise; its stages/betas remain the inherited dual-beta 3/6/10/12 schedule. It applies `prepare_full` and `full_material_readout` before calling `run_cie_fft` directly. The correct custom readout is not overwritten.
- Projected initialization uses the existing LS call with flag true and no physical_complement anchor. The wrapper returns the exact existing LS result, records the affine-family mismatch and additional diagnostic G/S RHS, and restores lstsq in finally. The diagnostic products consume the same paid CPU allowance.
- The paid observer records all direct/corrector value-gradient evaluations and auxiliary full physical validations; completed over-phase states are logged but excluded. The m0 row remains selectable. Corrector deadline is min(total, actual proposal end + fixed allowance); direct is one full-cap solve.
- AP/DR `checkpoint` accesses `AuditRecord.gap` and `.costs`, which match the supplied runtime dataclass. Proposal cutoffs retain the already paid material ledger, even when no ProjectionRun returns.
- `HISTORY.json` plus `CANDIDATES.npz` are sufficient to reconstruct primary, secondary, raw/proposal and corrected selectors. `PREFIXES.json` itself currently contains primary indices only; secondary/raw/corrected prefix indices can be generated offline from the full ledger without new physics.

One receipt correction is required: `PhaseEnd` derives from BaseException, while inherited `run_cie_fft` catches only Exception. Its finally block consequently returns a G19 execution receipt still labeled `RUNNING` on an intentional proposal cap. At the runner's G19 PhaseEnd catch (reviewed lines 131–132), reclassify the captured execution status to `PROPOSAL_CPU_PREFIX_HANDOFF` and retain the stop reason/CPU. This does not alter trajectories or candidate eligibility, but prevents a completed capped experiment from carrying a falsely live execution receipt.

The TASKS.json manifest was still being prepared during this read; task counts and B inside/outside dispatch must be checked once that file exists. The runner's public input loader relies on the preparation allowlist; the preparation manifest records all 31 inputs and full source dimensions.

## Sealed numerical/source review after authorized offline scoring

Reviewed the actual `RESULTS/offline_v1/` products and sealed online outputs after the parent explicitly authorized offline GT reads. This review ran only saved-array/metadata arithmetic: no forward, adjoint, optimization, training, new algorithm, or output-data modification. The parent retains scientific interpretation. All 672 online file size/mtime stamps still matched the seal at review time.

### Consistency findings

- All44 child receipts are COMPLETE, the finite queue is COMPLETED, and the public acceptance has no contract errors or warnings. There are 18 completed atomic states excluded by their phase/whole-child CPU deadline; independent reconstruction confirmed the same eligibility decisions, primary minimum-Phi selectors, secondary minimum-(R,Phi) feasible selectors, and paid-prefix indices. Excluded states did not enter the scored predictions.
- Each saved selected material matched its indexed CANDIDATES.npz row exactly. Every scored selected state's candidate/Phi/D/R/phase/CPU metadata matched its saved snapshot JSON. Independent recomputation from the saved GT/material arrays found no RMSE, NRMSE, relative-L2, amplitude or IoU discrepancies in available score rows, including secondary-deployment sources. This is offline verification of frozen outputs, not online selection by GT.
- The four methods have exactly identical saved initial m, j, Q and predicted Y for each of the 11 actually restarted cases. Thus their original full physical state and source acquisition agree before auxiliary representation changes. A/C caps are 240, B caps120; proposal/corrector allocations remain144/96 and72/48. Direct receives the complete cap.
- Table counts reconcile: 744 score rows =31 cases ×4 methods ×6 views; 400 paid-prefix rows; 124 corrector decompositions; 124 secondary deployments; 11 G19 initial-projection records; 44 child-cost rows; 31 GT-source rows. Twenty untriggered B cases are retained under each method/view and remain in the all24 denominator.
- The seal timestamp is `2026-11-11T20:37:16.986741+00:00`; offline TRUTH_ACCESS_BEGIN is later, `20:37:42.471276+00:00`. Online adapters import the public physical solver and public inputs, all44 online receipts declare truth_access NONE, and offline scoring imports no forward solver. This establishes the declared code/artifact separation; it is not an operating-system audit of every file access.

### Important limits when presenting the new numbers

1. **B has three physical/material improvements, two inherited descriptive successes.** DR reaches CN99 and improves material recovery on0001/0003/0022;0010 remains outside with NRMSE .365604. Corrected NRMSE for the three improvements is .107162/.119076/.028289. Case0001's IoU .659091 is below the inherited .7 descriptive threshold, so “3/4 success” would conflate noise feasibility/material improvement with that joint label. The joint count is2/4 among triggered cases, or all24 success9→11. All24 DR mean NRMSE is .110176 raw→.106261 corrected, whereas the four triggered cases alone are .178524→.155033. Keep these denominators explicit.
2. **The corrector contributes materially.** A DR mean NRMSE changes .200448 raw→.116873 corrected. In B0003, raw DR is .167742/IoU .584, corrected .119076/IoU .728033. In B0022, raw .063057→.028289, but IoU slightly declines .928571→.914798. Report the raw search and original-physics correction separately instead of attributing all final material quality to reflection alone.
3. **A's ring counterexample remains.** Ring0001 G19 primary NRMSE .114232 beats DR .140651. On shell0004, DR .112699 beats G19 .515137; on bridge0022, DR .097269 beats G19 .371690. The A mean improvement does not establish universal DR superiority.
4. **C's average comparison is close and selector-dependent.** Primary means G19 .0867682 and DR .0860935 differ only .00067466. The registered secondary means reverse their order: G19 .0801039 versus DR .0829008, with descriptive successes4/4 versus3/4. AP's poor primary mean .187543 is dominated by0001 (.487912); its other three cases are .080356/.062564/.119340. Direct also has a severe C0022 outcome (.289829). Four exposed NN cases cannot support a broad generalization claim.
5. **Lower Phi after correction does not guarantee lower material error.** C DR improves its mean .100916→.086094 but worsens NRMSE on0004 by .004787 and0022 by .011402. C G19 worsens0004 by .022675. A G19 worsens bridge NRMSE by .034687 while optimizing the original public objective. Preserve these per-object counterexamples.
6. **G19 starts from the same physical material but then loses substantial current information.** Its legal original-family LS projection loss is A=.818/.956/.845; B=.806/.798/.908/.980; C ranges .437–.716. This is the inherited admissible affine DPCS/Fourier projection, not an exact retention of the initial full physical current. It qualifies conclusions about what the control comparison isolates; it does not invalidate the recorded same-m/full-physical-start equality.
7. **Equal CPU caps are not equal actual expenditure.** For B0001/0003/0010/0022, direct used45.921875/.734375/15/.84375 whole-child CPU seconds, while DR used91.4375/91.4375/61.921875/92.078125. Direct stopped naturally at the unchanged local optimizer termination and had the full120 cap available. The comparison is between the approved algorithms under the same maximum cap, not matched actual CPU or an exhaustive direct multistart control. Use the paid-prefix and actual-cost tables for efficiency statements.
8. **Secondary availability must retain fallbacks.** A G19 has1/3 feasible secondary states, DR3/3, AP/direct0/3; B has20/24 for G19/AP/direct and23/24 for DR (twenty are retained originals); C has4/4 for G19/DR and3/4 for AP/direct. The feasible-only secondary mean omits unavailable cases. `secondary_deploy` correctly restores the full denominator through its explicit primary fallback, and that fallback must not be described as noise feasible.
9. **offline_v1 has a known parent-cost accounting defect; its scientific scores are unaffected.** Independent summation of the original B parent external completion receipts gives4290.1875 CPU seconds; offline_v1 PARENT_COST_LEDGER sums4286.734375, an underestimate3.453125. Root assigned a separate scoring worker to produce corrected v2 cost tables. Treat v1 parent/strategy totals as superseded; this reviewer did not alter them or infer a corrected scientific result from that cost correction.

No additional numerical/source mismatch was found. These findings validate the saved-output consistency and delimit presentation; they do not constitute a final novelty, mechanism or publication judgment.

## Final report fact-check against offline_v2

Read the completed supplement RESEARCH_REPORT_ZH.md, COST_LEDGER.md and the parent's PAPER_DRAFT.md paired-follow-up section against offline_v2 and the actual receipts. No new science, optimization, forward evaluation or data edits were performed.

No numerical or denominator blocker was found. Independent comparison confirms zero changes across the744 score rows' RMSE/NRMSE/IoU/relativeL2/Phi/D/R/candidate fields from v1 to v2. B's corrected parent external sum is4290.1875, DR additional336.875 and combined4627.0625 CPU seconds; the increment is7.8522209%. The A/B/C tables, raw/corrected means, B3/4 material-and-CN99 improvements versus2/4 inherited descriptive successes, twenty B bypassed originals, G19 ring win, C near tie/secondary reversal and G19 initial-current losses match the saved evidence.

All22 AP/DR children completed3000 A and C projections each; their66000/66000 totals and reported H setup/iteration RHS match aggregate receipts. Reported47716 physical LU,529224 source RHS,513390 adjoint RHS and1042614 total complex RHS also agree. C percentages recompute to31.1290% less child CPU than G19,39.6358% lower primary mean NRMSE than direct and45.4268% more child CPU than direct. The upstream NN preprocessing record is .763496 CPU seconds; the report correctly separates that old-platform shared cost from the new same-platform child comparison.

The report and paired manuscript section retain the important maximum-cap versus actual-expenditure, G19 PR-reset/affine-projection, same-start, exposed/posthoc and corrector counterexample qualifications. Scientific claims and final publication judgment remain with the parent.

Minor wording ambiguity only: report section7 says “之后一次授权读取GT”. Both v1 and v2 performed authorized offline scoring after the same immutable prediction seal, so “封存后经授权进行离线评分” would avoid implying there was only one actual GT-read pass. The report already discloses both offline versions and the cost correction; this wording does not affect numbers or the seal-before-truth chronology.
