# MAX-STRUCT-R4 source reuse and overlap audit

2026-11-03. Read-only audit of inherited source and already published/frozen receipts. Worktree HEAD: `1bd07e492c3d07a702817fa0429d6993de1a4cac`. The five inherited study directories inspected here had no local Git modifications at audit time. No inversion, derivative check, parameter sweep, GT-array load, commit, or push was performed. Only this new audit file was written. Historical reports and score summaries are treated as inherited evidence, not new R4 results.

## 1. Decisive findings

1. **Exact material Hessian/HVP already exists.** R4-BASIN implements the residual-Hessian term, exact smoothed-TV curvature, fixed-face eigenanalysis, and finite negative-curvature ray probes. Reusing it is preferable to claiming or implementing a new Hessian method.
2. **Finite rational exploration already exists at two scales.** R4-BASIN has a full dense Schur transfer on any material direction and scans its complete box-feasible interval. FBG-R3 has exact 4/6-cell finite response. MAX-STRUCT-R3 has rank-1 analytic and rank-2 polynomial-root paths, all-pair exchanges, and an observable full-wave interaction graph.
3. **Observable/cooperative pair selection is also inherited.** MAX-STRUCT-R3 selected pairs by `||O_i C_ij Q_j + O_j C_ji Q_i||_F²`; it already combines source excitation, dressed propagation, receiver visibility, and the interference cross term. Calling any observable feedback graph new would overlook this implementation.
4. **Current-space Krylov and Schur-induced material direction selection are inherited.** FBG-R2 has common-core Galerkin current exchange; FBG-R3 has residual-seeded Green-Krylov pools, joint-GN current elimination, material-task scores, exact physical JVP finalist checks, and signed finite current seeds.
5. **The inherited `MaxwellModel.regularizer()` matrix is IRLS curvature, not the exact regularizer Hessian.** Exact-curvature work must explicitly replace the returned weights; see Section 2. Historical joint-GN results remain historical surrogate models and must not be rewritten as exact Hessian results.
6. **No audited historical result establishes a general new basin/recovery advantage for these mechanisms.** Existing positive recovery belongs primarily to mature CIE; finite N1 has limited small-grid positive evidence, while exact current profiling, rational rays, and rank-2 cooperative paths retain their recorded negative/limited outcomes.

These are implementation/overlap findings. Selection of the new R4 mathematical direction and final novelty judgment remain with the parent Codex scientist.

## 2. Regularizer: exact Hessian versus the returned IRLS matrix

The real inheritance chain is:

- `research/research_fbg_r3_experiments/src/physics.py:12-13` imports `Regularizer` from `research/research_z_alg_r4_diversity/src/physics.py`.
- R4-diversity `Regularizer.evaluate`, lines 209-222, sets `d=D@m`, `q=sqrt(d*d+epsilon**2)` and returns
  `R=smooth*mean(d*d)+tv*mean(q-epsilon)`, the exact gradient, and
  `w_IRLS=(2*smooth+tv/q)/edges`.
- FBG-R3 `MaxwellModel.regularizer`, lines 57-60, constructs `D.T @ diag(w_IRLS) @ D` and names it `hessian`.
- FBG-R3 `joint_current_model.regularizer`, lines 59-75, passes that matrix through. `JointCurrentModel.residual/linearize` stores it as `Lambda`; `profile_schur` then uses `B.T(I-QQ.T)B + Lambda + lm*I`.

For the exact Hessian of the stated regularizer, the existing R4-BASIN implementation uses

```text
w_exact = (2*smooth + tv*epsilon**2/q**3)/edges
H_R h = D.T @ (w_exact * (D @ h))
```

This is explicitly implemented at `research/research_z_alg_r4_basin/src/hessian_geometry.py:33-35,47-48,147-156`. Its older `research/research_z_alg_r3_digits/src/multiscale.py:42-55` has the same objective/gradient/IRLS convention.

The two matrices coincide at an edge with `d=0`; otherwise the TV part of IRLS overweights exact TV curvature by `tv*d²/(edges*q³)`. Consequently, an exact-physics Hessian plus the inherited IRLS matrix is still **not** the exact Hessian of the inherited objective. It can hide negative curvature or change the smallest eigenvalues. This statement follows directly from the inspected formulas; no new numerical formula test was run.

`MaxwellModel.regularizer_residual()` (lines 62-71) is a separate residual/Jacobian representation for TRF. Its `J_reg.T J_reg` is another GN curvature and should also not be labeled the full scalar-objective Hessian. The objective and exact gradient used by the original L-BFGS corrector remain consistent with the stated regularizer.

## 3. Reusable interfaces and their boundaries

| Existing source/API | What is actually implemented | Boundary relevant to R4 |
|---|---|---|
| `research_z_alg_r4_basin/src/hessian_geometry.py`, `run_case` local `hvp(h)` (26-94) | One retained Maxwell LU; one base adjoint; each HVP solves a tangent and incremental adjoint, then adds exact regularizer HVP. Data contribution includes residual curvature. | Hardwired historical input/report workflow; reuse the local operator/action, not execute its historical driver. Eigenproblem is on coordinates farther than `1e-6` from bounds, not the full critical cone. |
| Same file, `dense_freeface_followup()` (132-171) | Assemble a small fixed-face Hessian by coordinate HVP and call dense symmetric `eigh`. Historical cutoff was 80 free coordinates. | The historical function loads fixed cases and has a hardcoded holdout-7 observation in the executed low-dimensional branch. It is not a general standalone library API; adapt explicitly to R4's public Y and endpoint. |
| `research_z_alg_r4_basin/src/rational_rays.py`, `RationalRay(model,state,direction)` (14-37), `scan_ray(...)` (40-78) | Exact full dense complex-Schur resolvent for an arbitrary real material direction; all source columns; complete box-feasible scalar interval; sampled valleys then bounded scalar refinement; full prior scored at every candidate. | Dense setup pays `n+ns` retained-LU RHS plus `n×n` complex Schur. No reduced-order approximation or complete global-minimum certificate. All separated **sampled** valleys, not all mathematical minima, are refined. |
| `research_z_alg_r4_basin/src/basin_solver.py`, `direction(state,name)` (175-194) | `amplitude`, `shape`, `offset`, normalized original-objective gradient, GlobalDual/small-Gram gradient, and damped GN/IRLS direction. | A generic feedback/dual gradient followed by exact rational scanning already overlaps. `bridge` gives distinct valleys corrector slots; `ray_best_only` is an existing ablation. |
| `research_fbg_r3_experiments/src/finite_response.py`, `ExactFiniteResponse(model,m0,support,state=None,ledger=None)` (63-139) | `F(d)=F0+O(I-diag(d)C)^-1 diag(d)E_A`, exact fixed-outside physical-cell response; analytic local Jacobian; source/receiver channel selection. Retained anchor LU used for support RHS. | Constructor deliberately accepts exactly 4 or 6 distinct cells. It is not an arbitrary-rank or arbitrary global-direction class. The complement is frozen. |
| Same file, `full_data_jacobian(model,state,ledger)` (142-147) | Analytic full-material Jacobian from retained LU and receiver-transpose RHS; correct source/receiver/real-imag packing. | Useful interface/packing reference, not a Hessian. |
| `research_fbg_r2_adaptive_feedback/src/schur_current.py`, `current_schur(L,core,candidate,rhs)` (10-22) | Retained/enriched Galerkin tangents with exact finite **subspace** Schur exchange `Cbar Gamma^-1 innovation`. | It is exact exchange between the specified subspaces, not an exact full Maxwell tangent unless their span suffices. |
| Same file, `krylov_spaces(G,Q,core_rank=None,candidate_rank=None)` (25-32), `projector_tail(...)` (35-46) | `U=orth(Q)`, candidates from projected `[GQ,G²Q]`, and an explicit exact-versus-approximate tangent projector tail diagnostic. | Current-space bank; source covariance and a small tail do not prove material recovery or a basin bound. |
| `research_fbg_r2_adaptive_feedback/src/schur_feedback_lane.py`, `make_schur_feedback(lane,x,kind,cfg)` (29-82) | CUDA common-core `U+2p` Schur current construction for seven public spatial material directions, source-normalized ridge feedback fit, one paid exact tangent tail, whole-box safety. | Exploratory `schur_diag` arm, not full A17 exchange strategy. Use only with the existing compatible resident CUDA lane; audit does not authorize GPU work. |
| `research_fbg_r3_experiments/src/joint_current_model.py`, `JointCurrentModel` (122-185) | Exact nonlinear joint residual and realified derivative. For material derivative it keeps `-w Q + dw e`, including nonphysical DPCS-anchor defect. Fixed positive source whitening is explicit. | `linearize` is joint GN plus inherited regularizer surrogate; it is not exact nonlinear joint Hessian or exact profile Hessian. |
| `research_fbg_r3_experiments/src/profile_schur.py`, `profile_schur(...)`, `release_candidate(...)`, `replacement_steps(...)` (45-124) | Current residual QR elimination; constrained material QP; soft/hard release in QR coordinates; equal real-current-DOF replacement. | Standard joint-GN/Schur algebra. Exact for that local quadratic model, not exact finite nonlinear elimination. Unconstrained bound needs its stated assumptions and is not automatically a constrained-step certificate. |
| `research_fbg_r3_experiments/src/current_mode_selection.py`, `build_common_pool(...)`, `evaluate_modes(...)` (174-404) | Common Fourier/receiver/Green-Krylov/random pool; material-step, energy, normal and fixed-family controls; proposal screening; paid exact physical-JVP finalist ranking. | Default proposal map is explicitly Born with **no error bound** unless `physical_J` supplied; finalists get exact JVP. That distinction must survive any reuse. |
| `research_fbg_r3_experiments/src/controlled_release.py`, `finite_normal_seed(...)` (23-134) | Material response `-H^-1 D.T w` induced by a current residual direction; both signs; box QP; current/defect caps; common bounded amplitude chosen by at most eight feasibility trials. | Finite **joint-current** seeds can be nonphysical before readout/corrector; not the same as exact physical finite-material rays. No proof of maximal feasible amplitude or monotonic feasibility. |
| `research_max_struct_r2/src/profile_current.py`, `run_profile(...)` (12-45) | Exact full-current elimination by one dense Cholesky per material evaluation; envelope gradient includes derivative of `beta/(1+beta*m)`; L-BFGS-B material solve. | Uses total data/incident energy scalar normalizations to share one factor across sources. It is not native sourcewise CIE's exact same metric. No exact reduced Hessian implementation here. |
| `research_max_struct_r2/src/cie_global.py`, `run_cie(...)`, `transport_source_metric(...)` (19-98) | Full fixed Hermitian SPD data/state source metrics, transported under mixing; global PR; source-whitened shared-real readout. | Optimizer covariance stated for unitary mixing. General metric transport identity does not establish arbitrary nonunitary PR/stop invariance. Tiny exact-current Kronecker/SVD diagnostic is capped, not the scalable main arm. |
| `research_fbg_r3_experiments/src/full_maxwell_corrector.py`, `corrector(model,Y,m0,config=None)` (8-88) | Mature bounded analytic Full-J TRF / reciprocal TRF; original or feedback analytic-adjoint L-BFGS-B; archive/incumbent selected by unchanged original physical score. | Return chooses best physical-score evaluated incumbent, not merely last auxiliary minimizer. No strong convergence certificate from generic stop text alone; preserve `max_nfev`, wall, and actual counters. |

Paths in the table are below repository `research/`. The source fields differ across adapters: FBG `MaxwellModel.state` provides `Y/y`, `Q/E`, `j`, `lu`; its backing `_physics` owns the retained `factor` and operation ledger. R4-BASIN `RationalRay` expects the backing-model `_matching(state)`, `factor`, and `operation` interface. Do not pass the wrapper itself without an explicit adapter. State/LU ownership must remain coherent when another material is evaluated.

### Existing exact HVP action

The old code's symbols are `L=I-diag(m)G`, `Q=E0+Gj`, `R=Sj-Y`, `s=||Y||²`, and `a=L^-H S*R`. It performs:

```text
U    = L^-1 (h Q)
dQ   = G U
da   = L^-H [G* (h a) + S* S U]
Hh   = Re sum_sources(conj(da) Q + conj(a) dQ)/s + H_R h
```

This transcribes `hessian_geometry.py:41-49`; it is inherited mathematics, not a new derivation. Each HVP pays `2*ns` complex RHS after the base forward/adjoint. `GN data Rayleigh + exact reg Rayleigh` is separately reported to expose the residual-Hessian contribution.

### Existing feedback directional curvature

`research/research_fbg_r2_adaptive_feedback/src/objective_geometry.py:10-38` implements `directional_geometry(Q,R,S,G,B,T1,T2)` for a **fixed** B, with transformed anchor V supplied as its `Q` argument and `T2=2 L^-1 diag(h)G T1`. It separates positive weighted-J, physical-residual second, Gram-cross, and Gram-second signed terms. This is one exact directional second derivative of the feedback objective; it is not a material HVP routine, and does not include the material regularizer.

`feedback_fit.py:443-482` supplies `feedback_value_gradient(shared,state,feedback,factor)` and `feedback_jvp(...)`. The main gradient requires a Torch lane state, retains the physical LU, uses one source Gram and one adjoint, and requires reference factor and B fixed during the inner solve. There is no continuously varying B(m) derivative here.

## 4. Cooperative material directions and exact finite paths already tried

The closest overlap to a proposed `O C^k Q` direction is MAX-STRUCT-R3, not only the older R4-BASIN rays:

- `experiments/finite_coordinate_paths.py:28-34` supplies `resolvent(model,m,cost)`: one full inverse from retained LU, `C=G L^-1`, `O=S L^-1`, with all `n` inverse RHS charged. `onecell_scan` forms exact rank-1 Möbius response and quadratic data-stationarity roots, then evaluates the unchanged prior.
- `experiments/paired_material_rays.py:23-53` supplies `coefficients(...)`, `scaled_derivative(...)`, and numerical degree-at-most-six roots for directions `(1,1)` and `(1,-1)` on two physical cells. `run` uses the common scalar interval and full prior; only actual original-objective improvement can invoke the common corrector. Data stationary roots plus endpoints do **not** exhaust stationary points of the regularized total objective.
- `experiments/dressed_neighbor_rays.py:16-46` supplies `select(model,C,O,Q,BO,BQ,public,cost,directory)`. Per cell it selects one partner maximizing
  `||O_i C_ij Q_j + O_j C_ji Q_i||_F²`, evaluated by source/receiver Gram contractions including the interference term. The union contained 567 unordered pairs, including 338 nonlocal pairs. The score is an unsigned response norm, not signed reduction of the observed residual.
- `experiments/all_pair_swaps.py` evaluated all 165600 coordinate exchanges at the same 576-material endpoint. This is a finite exchange set, not all two-dimensional material changes.

The reusable distinction is therefore **unsigned observable pair interaction versus a proposed signed, genuinely collective direction criterion and its finite outcome**. A higher-support direction could still be useful, but size/rank alone is not novelty, a proof of escape, or permission to expand a failed search. Both exact rational search and current/Schur material sensitivity already have nearby implementations.

## 5. Source metrics and corrector parity

Mature native CIE (`research_fbg_r3_experiments/src/cie_tsom_baseline.py:155-159`) fixes `wd_s=1/||Y_s||²`, `ws_s=1/||E_s||²`, with the same real material shared across all sources. FBG-R3 native-affine N1 uses equivalent source whitening (`experiments/e1_current_release.py:105-109`) while explicitly setting `physical_data_scale=||Y||`; source-weighted joint score and original physical score differ.

R2 full-current `run_profile` deliberately uses total-energy normalization. R3's fixed-material diagnostic **does** implement exact elimination in the same native per-source metric (`experiments/fixed_material_diagnostics.py:35,56-65`): it factors `wd_s S*S + ws_s H*H` separately for every source. The driver opens already published GT and must **not** be executed in a no-GT R4 inference process. Its algebra can be extracted from source with the public data and candidate material supplied explicitly. General distinct full data/state source matrices need coupled/Sylvester treatment rather than an assumed single common Cholesky.

For R4 comparisons preserve: the actual original objective, material box/prior, the affine DPCS offset, independent APCS rank, actual real source-current DOF, source metrics, and the physical corrector allowance. Similar DOF or a lower auxiliary profile is insufficient. R3 `src/cie_strong_release.py:11-33` is already an ordinary final-stage DPCS release control; it monkeypatches/restores the mature Fourier hook and is intentionally sequential/not thread-safe.

The R3 scaled-corrector result also matters: `results/development/scaled_original_restart/RECEIPT.json` records objective scaling `1e6`, 745 iterations/748 evaluations, and termination by projected-gradient tolerance equivalent to `1e-8` in original units. It reached `PGinf=9.915385845715718e-9`, `Phi=0.005096102692209315`, residual `0.10073368671742994`, at 748 LU/17964 complex RHS. This is the qualified first-order endpoint for a new curvature question; it is not a strict local-minimum certificate.

## 6. Frozen evidence and what it permits

| Study / exact evidence path | Frozen or published result | Permitted interpretation |
|---|---|---|
| R4-BASIN `diagnostics/hessian_geometry.md` and `.json` | h7 and h9 ARPACK stopped at 250 HVP without eigenpairs. h8: smallest exact fixed-face eigenvalue `-5.65881e-6`; corresponding GN+exact-reg Rayleigh `8.643652e-5`; PGinf about `2.001e-4`. Positive feasible ray improved objective by `2.800224e-4` at RMS `0.057502`. | Exact residual-Hessian negative curvature exists at that old endpoint. Nonzero first-order gradient and boundary truncation prevent stationary-saddle/nonlocal-escape claims. Eight probes collapsed to three unique feasible alpha values; all costs remain paid. |
| Same, dense h7 followup | 31-dimensional fixed-face eigenvalues start `1.760459765e-8, 2.411887526e-8, 3.546555008e-8`; added 31 HVP, 1 LU, 768 RHS. h9 free dimension 200 was NOT_RUN. | Nonnegative tested fixed-face spectrum; not full box local/global optimality. Earlier capped ARPACK work remains inconclusive and paid. |
| R4-BASIN `audit/V3_SOURCE_SNAPSHOT/FREEZE.json`, `audit/SCIENCE_FREEZE.json` | Source freeze `2026-10-30T03:11:58Z`; final scientific freeze `03:33:52Z`; 18 fresh predictions, zero joint successes; `NO_FRESH_GT_IMPROVEMENT; PAPER_HOLD`. | Existing rational/continuation/basin machinery did not establish fresh material recovery. Historical geometry-specific positives are separate. |
| FBG-R2 `audit/STAGE4_FREEZE.json`, `RESEARCH_REPORT_ZH.md:3,87,109-142` | Six methods frozen before new generation; 14 independent materials/12 high-contrast, 14 clean plus 8 noisy observations. Fixed beta .2: 8/12 high-contrast clean joint successes; automatic beta 6/12; original GlobalDual 3/12. Schur imaging probe 2/5 versus 2/5. | Fixed feedback has scoped positive evidence; fitted B/Schur did not independently outperform. Missing mature CIE fresh strong comparison in that study; no universal advantage. |
| FBG-R3 `audit/SCIENTIFIC_FREEZE.json`, `RESEARCH_REPORT_ZH.md:7-24,185-222` | Numerical closure `2026-11-01T07:07:10Z`; mature dual-beta CIE+corrector 10/12 on disclosed historical 576-material objects, including all nine digits. Native-affine finite-current development gave generic finite coverage 7/12 versus 6/12; new confirmation did not reproduce an overall joint gain. One clean 36-material task recovery had NRMSE≈.000310 and residual≈8.9e-5, but original IoU gate remained FAIL. | Preserve limited N1 positive material/data evidence; no task-exclusive finite advantage, no overall joint gain, and N1/N2 576-superiority both NOT_RUN. |
| FBG-R3 `RESEARCH_REPORT_ZH.md:226-242`, `results/e2_scores.json` | N2 algebraic conditional roots gave no independent complete-image gain; algebraic and block multistart each 2/40 joint successes, algebraic successes from anchor fallback. | Exact fixed-outside conditional roots are not new complete imaging roots or an algorithm advantage. |
| MAX-STRUCT-R2 `results/BLIND_SCORES.csv`, `RESEARCH_REPORT_ZH.md:28-57` | Native/global 2/3, full-current profile 1/3, direct physical 0/3. Post-hoc mf12 bar: residual `.100735→.097693`, NRMSE `.510695→.558744`. | No new superiority over mature CIE; lower physical residual alone is not improved material recovery; post-hoc capacity test outside blind denominator. |
| MAX-STRUCT-R3 `results/development/paired_neighbor_rays/RECEIPT.json` | 2208 sign/pair cases, 1223 nonzero-interval root problems, 2470 finite candidates; no accepted move from scaled endpoint. | The specific exact paired-ray family failed, not all two-dimensional variation. |
| MAX-STRUCT-R3 `results/development/dressed_neighbor_rays/RECEIPT.json` | 567 observable-interaction pairs; final `Phi=0.005096102692103583`, versus seed `0.005096102692209315` (about `1.06e-13` decrease); final residual `.10073368671518108`. | No substantive escape from this dressed full-wave pair graph. An accepted machine-scale decrease must not be described as recovery. |
| MAX-STRUCT-R3 `RESEARCH_REPORT_ZH.md:27-45`, `results/DEVELOPMENT_SCORES.md` | All 11 developed controls failed original joint imaging conditions; passive path reduced residual to `.051017` but worsened material NRMSE to `.817495`. | The current failure remains solver/reachability evidence with explicit scope, not finite-data nonidentifiability or a proof that all higher-rank paths fail. |

This audit did not load any sealed truth or prediction NPZ arrays. Numeric findings above were read from existing human-readable reports and JSON receipts. No historical score was recomputed.

## 7. Concrete R4 reuse boundary

- Reuse the exact R4-BASIN HVP action with explicit exact regularizer weights and a current public endpoint; preserve free-face versus feasible-cone scope and count tangent/incremental-adjoint RHS.
- Reuse the existing finite-transfer identity/Schur evaluator for whichever specific direction the parent qualifies. Maintain a scalar box-feasible interval; independent coordinate clipping changes the intended ray. Any sampled search remains sampled and numerical roots remain uncertified unless separately justified.
- If the new direction comes from current release, state whether the score is exact joint GN, a Born proposal, exact physical JVP, feedback-objective curvature, or original-physical exact Hessian. These are distinct objects already represented by different modules.
- Use existing `full_maxwell_corrector.corrector` and the same declared physical objective for handoff. Record solver termination and actual physical/auxiliary work separately; do not infer local minimality from an iteration cap or ordinary ftol stop.
- Treat the mechanism difference from `dressed_neighbor_rays`, GlobalDual rays, and finite-current Schur seeds as the unresolved scientific issue. This audit finds no existing general signed collective-material direction selector with an established finite-escape guarantee in these inspected modules; that is a bounded source-inventory observation, not a novelty claim.

## Audit provenance

The local memory registry was used only to locate inherited R2/R4-BASIN scope and avoid conflating publication with recovery; all scientific/source facts used above were checked against the current worktree. Memory pointers: `/Users/migodam/.codex/memories/MEMORY.md:49-68,85-106`, prior chats `01a19a85-dae6-7e33-be05-7bd6cb7693d3` and `01a18851-dc9b-7a53-bf0e-7ef44f639189`. No memory file was changed.
