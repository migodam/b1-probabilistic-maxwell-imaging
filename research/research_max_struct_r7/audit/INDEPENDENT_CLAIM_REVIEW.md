# R7 claim review against source and retained evidence

Date: 2026-11-04. Status: **PROVISIONAL CLAIM REVIEW; final sealed statistics NOT_REVIEWED**.

This review reads `THEORY_DERIVATIONS.md` and `PAPER_DRAFT.md` against the implemented scalar A+ and rational actions, the complete disclosed development records, vector-DDA source/receipts and reports, and Fresnel adapter/inversion/audit. It executes no numerical experiment, forward solve, optimizer, new derivative test or scoring operation. No raw/sealed material truth, scorer-only arrays or new R7 primary method rankings were opened. Previously reported 3D aggregate material scores were read in `BENCHMARK_3D.md`; they were not independently recomputed from truth here. No external literature/DOI verification or exhaustive novelty review is claimed.

**Reviewer independence limitation:** this is a separately tasked, source-based claim review by the same native Codex workstream that implemented A+. It is not an external or author-blind review. The rational implementation has a separately authored implementation audit. Parent Codex retains final theory validity, novelty, interpretation and scientific acceptance.

The two explicit manuscript placeholders `[FINAL_PRIMARY_TABLE]` and `[FINAL_CONCLUSION]` are recognized as pending sealed-statistics work, not graded as failed final results.

## Overall assessment

The conditional scalar Möbius algebra, unchanged-prior MM acceptance and actual four-color implementation are consistent under their finite-model assumptions. I found no fatal algebraic discrepancy in that chain. The finite-action resolvent and port-Gram formulas also match the reviewed source. The draft generally respects the key distinction between evaluating an update correctly and improving material recovery.

Before final scientific sign-off, tighten the complete-surrogate notation, make the sourcewise/global PR comparison confound explicit, and carry the Fresnel host-selection recording limitation into the manuscript. Make clear that the 3D controls are a different vector-adapted experiment, not execution of scalar A+ in 3D. Final new-material success, cost frontiers, uncertainty and initialization conclusions remain pending the parent-authorized sealed-statistics review.

## Reporting findings to resolve

| Priority | Location at reviewed draft | Finding and required clarification |
|---|---|---|
| P2 | `PAPER_DRAFT.md:39` | The proposition calls Q a **prior** majorizer but then writes `f(new)<=Q(new|old)`. If f includes the state loss, Q must mean the **complete** surrogate. Use `q_R` for the prior bound and `M=f_state+q_R`, then write `f(new)<=M(new|old)<=M(old|old)=f(old)`. This is a notation gap, not a failure of the actual implementation. |
| P2 | `PAPER_DRAFT.md:27,55–63` | Explicitly identify the same-current control: A and Global-Dual use global Frobenius PR; Dual-CIE uses sourcewise PR. A versus Dual changes both material readout/metric and PR coupling, so their difference cannot isolate the A material mechanism. A versus Global-Dual addresses that material-position contrast, while A+ versus original A addresses the readout implementation. Source metrics/beta schedules and sourcewise/global choices should be visible together. |
| P2 | `PAPER_DRAFT.md:79–81` | The Fresnel paragraph reports identical endpoints and no advantage, correctly, but omits why this is weaker than three successful auxiliary inversions. All first-corrector objectives equal the host objective exactly; original auxiliary candidate material arrays were not retained. The evidence is **consistent with retaining the host**, not a strict proof that the chosen m equals host, and does not show any auxiliary method improved the measured trajectory. Carry that recording limit into the manuscript. |
| P2 | `PAPER_DRAFT.md:73–77` | Explicitly state that the vector experiment does **not** execute the scalar four-color A+ proposal or the mature 2D PR/Fourier algorithm. Its joint current solve is an ordinary full-space quadratic Cholesky update; it uses a quadratic-edge plus ridge prior and its declared vector normalization/box. These are valid adapted controls, but their agreement after correction is not validation of the scalar A+ theorem in 3D. `BENCHMARK_3D.md` already gives the appropriate distinctions. |
| P3 | `PAPER_DRAFT.md:45–49` | State the nonsingularity assumptions alongside the finite-response formula: invertible anchor and perturbed full operators, equivalently legal nonsingular reduced solves. `THEORY_DERIVATIONS.md:55` states them already. Avoid reading “exact finite response” as a guarantee at singular candidate systems. |
| HOLD | `PAPER_DRAFT.md:57–63,95` | Main-table counts, uncertainty, matched-budget frontiers, initialization/historical results and the final conclusion require the completed seal/scoring evidence. No pass/fail assessment of those missing final numbers is made here. Claims that all main predictions were sealed must be supported by the actual final seal receipts when the text is finalized. |

Suggested concise manuscript additions, without new experiments:

> A-material and Global-Dual share the global Frobenius PR current update; Dual-CIE instead uses sourcewise PR. Comparisons against Dual-CIE therefore do not by themselves isolate the material metric/readout.

> In the measured pilot, all first-corrector objectives equal the common host objective. Since auxiliary candidate material arrays were not retained, this is consistent with host retention and provides no evidence that the auxiliary methods improved the measured trajectory.

## 1. Conditional state, prior-MM and colored implementation

### 1.1 Exact state structure, not exact regularized minimization

For fixed J,Q, source metric and beta, set `b=beta-Gii`, `z=m/(1+b*m)` and `H=Q+b*J`. The implemented state residual is exactly `beta*(I+diag(z)*Goff)*(J-diag(z)*H)`. Goff has zero diagonal and is confined to the actual disjoint 2x2 blocks; hence every residual row affected by one z_i is affine in that coordinate. This does not make the residual jointly affine.

After whitening with the actual supplied Cholesky Ls, the code constructs the entire block coefficient:

- own row: `B_i=beta*(-H_i+(Goff@V)_i)`;
- each other row j: `B_j=-beta*z_j*Goff[j,i]*H_i`.

`src/aplus_readout.py:114–148` includes all affected rows and sources. With `q=rho*sum(abs(B)^2)` and `c=rho*sum(conj(B)*A)`, its rational coefficients are `n1=Re(c)`, `n2=q/2+Re(c*b)`. The stationary numerator `n1+2*n2*m+(n2*d1-n1*d2)*m^2` in `THEORY_DERIVATIONS.md:43` matches the code.

The final arc branch (`src/aplus_readout.py:58–79`) proposes only feasible state-stationary roots, endpoints and the old coordinate. It **does not** find the minimizer of the regularized conditional objective. Degenerate roots, roundoff filtering or an incomplete proposal set do not invalidate the acceptance argument: global minimization is unnecessary when a correctly evaluated complete surrogate and old-point candidate are retained. They can affect progress, not the exact-arithmetic nonincrease implication.

### 1.2 Actual prior and per-color proof

The inherited scalar prior uses normalized edge differences `(m_j-m_i)/19.5`, and both penalty terms are means over all original edges. `_prior_context` extracts actual sparse D values and checks two-cell edges and absence of a same-color prior edge (`src/aplus_readout.py:26–47`). It does not silently replace that prior with raw differences or a penalty on transformed material.

For an incident edge, tangency of the square-root bound gives coefficient

`w_ei = [2*smooth + tv/sqrt((D_e*m_old)^2+epsilon^2)] * D_ei^2 / edges`.

The code forms `h_i=sum_e w_ei`, `s_i=sum_e w_ei*m_neighbor`; therefore `.5*h_i*m_i^2-s_i*m_i+const` is the correct conditional prior surrogate (`:140–143`). Zero-filled absent neighbor slots carry zero coefficient, so they do not impose an artificial boundary penalty.

Let C be one color and freeze all other coordinates at the current point. Every 2x2 propagation block contains exactly one C coordinate, and no prior edge connects two C coordinates. The complete material objective and its majorizer therefore decompose into independent coordinate terms plus a constant. Choosing each candidate by its **state plus prior** surrogate yields `M_C(new)<=M_C(old)`. Tangency and majorization then give `f(new)<=M_C(new)<=M_C(old)=f(old)` for that color. Recompute the tangent weights at the next color and repeat. The code performs exactly this sequence: batch each color, update the cached full block residual, update V, then proceed to the next color. It is not an unjustified simultaneous update of mutually coupled coordinates.

The initial scalar seed is also accepted by the original complete A loss (`:13–23`), so it can be prepended to this chain. The final original-loss comparison/fallback (`:154–157`) protects the retained readout from observed floating-point increases. The mathematical statement remains conditional on finite pole-free legal coordinates and correct exact-arithmetic surrogate evaluation; the measured numerical nonincrease record is separate evidence.

### 1.3 dP, inherited tolerance and meaning of success

The candidate construction keeps the material-dependent transform's derivative implicitly through the complete exact coordinate dependence. Its stopping test calls the inherited exact gradient (`src/aplus_readout.py:122–128`). The inherited function contains both `du` and `da*Goff*u` contributions (`research_max_struct_r4/src/cie_onehop.py:105–119`); dP is not frozen away from a constitutive root.

The final four-object arc development executed 12,000 readouts, 24,000 sweeps and 36,000 inherited exact-gradient checks: **all reached the two-sweep budget**. Do not explain the final structural speed improvement as fewer sweeps or as a demonstrated gradient-tolerance acceleration. The substantive algorithm change is removal of quintic companion eigensolves. Earlier guarded versus unguarded timings also vary with shared-machine conditions; the recorded small sweep-count changes do not explain all CPU differences.

`material_inner_success=True` in the readout means the declared update procedure completed; it is not proof that the material subproblem reached first-order tolerance. The paper presently describes a fixed two-sweep update and does not claim inner convergence; preserve that meaning.

### 1.4 What this does not prove

The fixed-current A state objective differs from the scalar residual used in the inherited current update. Thus block descent is not descent of one common complete alternating objective. It also does not establish original physical-Phi monotonicity, identical PR history, retention of a basin, uniqueness or GT recovery. The manuscript correctly states these limits. Projected8's b20d public-Phi failure despite retained block-loss nonincrease is direct development evidence against equating these properties.

## 2. Development evidence and same-current attribution

The 20 full disclosed runs are four objects under five readout variants, not 20 independent materials. `delegated/aplus/FULL_TRAJECTORY_COMPARISON.md`, `DEVELOPMENT_COST_SUMMARY.json`, `FINAL_COUNTERS.json` and `APLUS_FREEZE.json` support the draft's stated values:

| Observation | Evidence and permissible interpretation |
|---|---|
| Final arc total CPU197.781556s versus original A284.209923s | Measured four-object complete cost reduction; one old bar remains about4.8% slower. It is not a universal speed result or a new-material population estimate. |
| Final arc quintic roots/eigensolves0 | Structural operation removal; 27,648,000 quadratic stationary-root proposals remain, plus full surrogate/prior and gradient work. |
| All60,000 retained material readouts show zero logged block-loss increases | Numerical implementation record over all20runs, not an original-Phi or material-recovery theorem. |
| Projected8 b20d Phi.0058324 versus original A6.4736e-6 | About901-fold public-objective degradation after the complete declared trajectory; rejects that replacement on this disclosed case without needing truth. |
| Four final arc objects retain the disclosed public physical states | Both old bars' public Phi/residual scales remain; slit shell remains at approximately.087 residual. Actual material recovery retention awaits parent independent scoring. |

The current solver identity must be precise. `research_max_struct_r6/src/benchmark_solver.py:64–80` sets `sourcewise_cg=True` for CIE/Dual and False for Global-Dual/A. `src/cie_fft_r7.py:222–249` computes sourcewise or global PR ratios accordingly. A and Global-Dual share the global-current algorithm and beta schedule; they do not have identical numerical J trajectories after different material updates. “Same current update” means the **same conditional algorithm**, not equal currents throughout a run. A+ versus A likewise shares that algorithm while changing the material-update mapping.

No new-main-material conclusion is inferred from these disclosed development objects. The arbitrary1.01relative-Phi gate was proposed then removed; all four final controls completed. It is not a scientific acceptance threshold and must not reappear as a hidden selection gate.

## 3. Rational finite actions and Hessian interpretation

`src/rational_search.py:11–24,37–42,68–86` uses the correct anchor order `I-X0*G`, selected columns of `G*A0^-1`, and response `(I-D*K)^-1*D`. The general-q cross term `Re(vdot(Z,T))` and quadratic trace use the correct conjugation. The rank2 four-feature Gram also has the correct determinant/numerator signs. The local prior increment uses original affected edge differences and divides by the original global edge count (`:48–66`). These agree with the theory, under the stated nonsingularity and fixed scalar observation-normalization assumptions.

The candidate online solve is small **after paid setup**. Coherent rank32 actually obtains576selected resolvent columns, in addition to receiver-transpose RHS and full physical endpoint checks. Two group amplitudes do not make that action state-rank2. The source selects support/sign using weak-Hessian/interface/gradient directions and then uses box-distance weights (`experiments/run_rational_coherent.py:45–53`); it does not follow the original eigenvector component ratios. The draft's wording already respects these limits.

The true physical Hessian includes GN, residual feedback, and the true original prior Hessian. The coherent source uses the correct TV second-derivative weight `epsilon^2/(d^2+epsilon^2)^(3/2)` rather than the IRLS/MM weight. For the shell, the reported smallest computed free-face eigenvalue4.341386e-9 is restricted to50selected interior variables. It is not a feasible/critical-cone test at all active bounds, not proof of stationarity, and not a box-constrained local-minimum or basin certificate. The paper's “diagnostic only” qualification is correct.

Public rational receipts show selected full-versus-reduced endpoint Phi gaps of order1e-17. They validate executed endpoint wiring, not every candidate point, continuum accuracy or escape capability. On ring, anchor continuation reaches.004929460949 while coherent rank32 reaches.004929485902; on42cb, continuation reaches.002182490593 versus coherent.002183314073. These support “no additional basin capability demonstrated”; they do not prove globally that only extra budget matters or that no better action family exists. The shell's finite/corrected changes remain tiny, so its positive free-face curvature must not be used to upgrade this bounded negative search into an impossibility theorem.

The serial rank2/rank8 comparator budget arrays include earlier unrelated correctors. The draft correctly excludes treating them as independently matched-CPU curves. `RationalAnchor.setup_cpu` also excludes anchor and later per-action Gram setup; only complete process accounting supports study-total costs. The separate `audit/INDEPENDENT_IMPLEMENTATION_REVIEW.md` documents these restrictions in detail.

## 4. Vector-DDA evidence

The actual source has125/216real isotropic material variables and375/648vector state components. It uses full dyadic propagation with scalar isotropic material repeated over three vector components (`src/vector_inverse.py:14–31,40–52`). The A material derivative includes the derivative of the transform (`:61–75`). The stored off-root derivative receipt activates constitutive residual.2048 and reports A relative directional error8.6841e-9; a root-only check would indeed have been inadequate.

The joint current update is a full ordinary data/state normal-equation Cholesky solve (`:54–60`), shared by the adapted ordinary/local-resolvent/A material variants. The material prior is quadratic nearest-neighbor smoothing plus ridge (`:36–39`), and the vector physical objective has its own normalization. These are legitimate disclosed vector controls but not a faithful mature3D CIE benchmark, and not the scalar four-color A+ algorithm.

The published aggregate continuation table in `BENCHMARK_3D.md` supports the draft's quoted post-continuation A NRMSE values and small differences among methods. It retains300-cap exits for the two-sphere/bracket joints; similar endpoints are not proof of the same exact optimum or convergence. I did not reopen material truth or recompute those scores in this review.

The actual-grid Mie validation is an independent forward comparison: binary sphere error6.584%→27.588% fails monotone-refinement behavior; the single fixed volume-fraction approximation gives8.396%→8.199% and still substantial error. It is not an inversion-data repair, not a calibrated cut-cell material model, and not evidence of continuum-qualified3D imaging. The manuscript's negative physical-transfer conclusion is appropriately limited.

## 5. Measured Fresnel evidence

The adapter (`src/fresnel_adapter.py:14–52`) subtracts target-absent incident E from total E, conjugates both conventions consistently, keeps missing channels absent, and builds the training split from actual measured channels. The cylindrical source gain uses only target-absent on-axis incident data. Held-out scattered E is evaluated only after training-based selection (`experiments/fresnel_invert.py:37–61`); it does not set gains or select an iterate in this source.

The six receipts/audit support234training and1530withheld channels per target, train residuals.141931/.100624 and withheld residuals.267555/.234442. These are relative complex fit errors, not noise-normalized statistical acceptance: no calibrated covariance is supplied. Horn field mismatch, limited aperture, the inherited disk-pulse parameterization and contrast-cap contact restrict absolute material interpretation.

`RESULTS/fresnel_auxiliary_selection_audit.json` shows every first corrector objective equals host exactly. Together with identical saved final material endpoints, this supports a common physical-corrector explanation and no demonstrated auxiliary contribution. Because stage candidate arrays were not retained, do not turn that consistency into a strict recovered-material identity at auxiliary selection. Twin correctors hit the declared cap; do not call them converged. The draft's no-pixel-GT/no-absolute-permittivity/no-A-advantage conclusion is supported, with the additional host-recording limitation requested above.

This review does not reopen raw measurements, certify the download/redistribution license, or verify the eventual release filter. Those are separate provenance/release tasks.

## 6. Required final-statistics follow-up

After parent notifies that complete predictions and final scoring are sealed, review the filled primary table and conclusion against their exact artifacts. In particular:

1. Keep material NRMSE/IoU separate from Phi, physical residual and declared discrepancy.
2. Keep30new objects separate from repeated initializations, historical digits, disclosed development, Fresnel targets and3Ddevelopment objects.
3. Verify paired object-level uncertainty and the common CPU/LU/RHS support/checkpoint resolution; nominal update caps do not imply equal cost.
4. Distinguish same-grid finite data from cross-grid data with shared physical backend, and noise-only discrepancy from model mismatch.
5. Treat actual old-bar material recovery retention as a sealed-score result, not a consequence of public Phi retention.
6. Update abstract/conclusion only after these results are available, without upgrading development/transfer limits to a universal recovery or novelty claim.

**Current verdict:** scoped finite algebra and mechanism descriptions are consistent; the draft requires the reporting clarifications above and the separate final sealed-statistics review. No new experiment, physics call, GT access or numerical code change is requested by this review.
