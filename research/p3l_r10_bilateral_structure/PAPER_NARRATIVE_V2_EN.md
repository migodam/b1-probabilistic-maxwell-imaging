# Fixed Bilateral Scattering Operators: Amplitude Reuse and the Limits of Material Possibility Coverage

Research manuscript, 30 October 2026. This version incorporates sealed Fresh10 scores, independent score recomputation, the DEV compression scan, and the separately frozen heldout-shape kernel probe. It supports a scoped computational benefit and a negative candidate-coverage result; it does not establish publication-scale qualification.

Alternative titles: **Residual-Guarded Amplitude Reuse for Fixed-Shape Inverse Scattering**; **Bilateral Structure and Material Possibility in Nonlinear Inverse Scattering**. The first title emphasizes the supported computational result, with coverage retained as a separate negative study.

## Abstract

Nonlinear electromagnetic inversion can admit materially different objects that explain the same noisy measurements. Useful material possibility sets therefore require evidence of coverage, alongside computational efficiency. We examine what fixed bilateral scattering structure contributes to these two objectives. Exact host elimination separates illumination, propagation, and observation from a material-dependent resolvent. At fixed shape, complex material amplitude defines a rational transfer family. Building on established eigenpermittivity and realization theory, we combine standard block Krylov reduction with a full-state residual and a coercivity-based numerical error guard. An order-64 model, selected on eight development shapes, passes both output and likelihood-potential gates at all 320 amplitude nodes on ten heldout anchor shapes without fallback. For the executed workload of 32 queries per shape, aggregate cold kernel time across all ten shapes is 0.1150 seconds, compared with 0.3173 seconds for full-order Schur reuse and 0.2387 seconds for direct factorization, reductions of 63.76% and 51.82%. Separately, a sealed ten-object experiment compares full randomized maximum-likelihood targets, shape-only targets, and a fixed bilateral covariance variant using eight proposals and a common full-coordinate solver. Material coverage at RMS tolerance 0.075 is 3/10, 3/10, and 2/10, respectively, against 2/10 for the anchor; the bilateral variant provides no supported action-error advantage. Three conditional amplitude references also satisfy their numerical task gates, without qualifying a full-material posterior. An executed 836-node fixed-shape probability replay saves 59.22% with the original nodes and weights, while reusing historical grid construction. These findings support repeated fixed-shape computation while limiting the tested target heuristic. The reduced order is empirical, the guards are numerical discrete bounds, and no new scattering identity, globally minimal realization, calibrated posterior, or universal coverage claim is established.

## 1. Introduction

Inverse scattering is often presented as the recovery of one material image from measured fields. A single image is useful when a downstream action demands a point estimate, but it can conceal alternatives that remain compatible with the same measurement. The relevant scientific question is therefore not merely whether an optimizer can find multiple solutions. It is whether a limited collection of legal material configurations represents the possibilities that matter, including the generating object or a task value, with a useful extent and a defensible computational cost.

This question does not require an accurate posterior as an entry condition. For a declared distribution of objects and noise, empirical coverage of neighborhoods around a finite candidate set is a meaningful imaging objective. Posterior probabilities, credible intervals, and importance weights have stronger requirements and should be evaluated separately. Conversely, omitting a posterior claim does not excuse a cloud that is wide but misses the generating objects.

Our preceding R9 study exposed precisely this distinction. With sixteen direction proposals per object, a fixed weak-subspace method produced many legal, data-compatible alternatives, but covered the generating material within RMS 0.075 for only 6 of 32 objects, compared with 3 of 32 for a shared anchor. Its mean reconstruction did not improve upon the anchor. Post hoc geometry implicated omitted directions and finite candidate placement; merely preserving more points along the same trajectories did not add primary-scale hits. These observations motivate a change of proposal objectives and allocation rather than another identity-accuracy test or a larger repetition of the same weak-space search.

There is substantial prior art. SOM/TSOM organize induced-current optimization using fixed propagation structure; RTO and rMAP/RML generate candidates by randomized optimization; LIS and goal-oriented reduction balance data, prior, and task information. T-matrix decompositions and classical realization theory describe input–output behavior. Our contribution cannot be the rediscovery of these components. The narrower questions are whether this structure lowers repeated amplitude-query cost and whether a specific bilateral construction improves randomized material targets under the original nonlinear model and a finite proposal budget. See [the evidence ledger](PRIOR_ART_AND_NOVELTY.md) for verified identities, support levels, and conflicts.

The closest precedent for separating fixed shape from material amplitude is more specific: [material-independent modes](https://doi.org/10.1103/PhysRevB.94.201406) and [eigenpermittivity-based Green-tensor expansion](https://doi.org/10.1103/PhysRevApplied.11.044018) already provide reusable geometry-dependent representations. We therefore treat rational amplitude separation as established structure and evaluate the practical value of our particular solver, error guard, and target-distribution choices.

## 2. Scientific object and claims

Given observation y, let \(\mathcal C_K(y)\) be the finite set of distinct candidate materials that belong to the original decoder image and satisfy the original noise-based data test. Include the anchor only if it qualifies. At material tolerance epsilon, define
\[
\mathcal V_{K,\epsilon}(y)=\bigcup_{m\in\mathcal C_K(y)}B_M(m,\epsilon).
\]
The primary material quantity is the repeated-object event \(M_{true}\in\mathcal V_{K,\epsilon}(Y)\). This is not a posterior mass statement for one observation, nor a Cartesian product of independently chosen pixel intervals.

Three contributions can be stated at their appropriate levels:

1. A measured fixed-shape computational mechanism: standard block Krylov/Galerkin reduction, complex Schur amplitude reuse, and full-state residual guards reduce cold repeated-query cost on the declared discrete model. The exact identities and reduction components are established mathematics.
2. A reproducible negative comparison of three randomized target distributions under one full-coordinate nonlinear solver: the tested shape and bilateral variants do not improve primary coverage over full RML. Block elimination is shared and standard.
3. A qualification framework that separates empirical useful order from exact minimal realization, finite material coverage from probability calibration, and numerical fixed-shape references from Full66 inference. All attempted costs, exposure status, and conditional scopes remain explicit.

## 3. Model and bilateral information

After exact fixed-host elimination,
\[
F(D)=F_h+O(I-DC)^{-1}DE.
\]
O, C, E and the host response are shared. D contains the spatially varying material, and every arm uses the same exact discrete ROI solver. The derivative is
\[
\mathrm dF=O(I-DC)^{-1}\operatorname{diag}(\mathrm dd)(I-CD)^{-1}E.
\]
Consequently the host Born operator has Gram \((O^*O)\circ\overline{EE^*}\). This fixed Gram is useful proposal information; it is not the exact nonlinear Fisher matrix throughout material space.

For \(D=\delta A(\theta)\), \(P^2=A\) gives
\[
F-F_h=\delta(OP)[I-\delta(PCP)]^{-1}(PE).
\]
This equality allows repeated amplitude evaluations at one shape to share a Schur factorization. It does not make all shapes share a factorization. The response is rational in amplitude, and a full nonlinear amplitude profile remains different from an ordinary linearized block step. Its minimal realization is the reachable/observable quotient, equivalently the sufficiently large exact block-Hankel rank. The subsequent reduction study tests orders 32, 64, 96, and 128. Order 64 is the smallest tested order passing the numerical guards on all development nodes and then passes the heldout-shape probe. This is an empirically useful order, not a measured or proved exact minimum realization.

## 4. Actual algorithm

All three new arms minimize randomized objectives of the form
\[
\Psi_{\eta,\xi}(x)=\tfrac12\|r_y(x)-\eta\|^2+\tfrac12\|x-\xi\|^2,
\]
where \(r_y\) uses the original observation and whitening, \(\eta\) is a common 512-dimensional standard Gaussian draw, and x is the original Gaussianized 66-coordinate representation. The last two coordinates map to the original bounded amplitude variables through the normal CDF.

`FULL66_RML` uses a standard Gaussian target in all 66 coordinates. `SHAPE_TARGET_RML_BLOCK` sets the two amplitude-target perturbations to zero. `BILATERAL_TARGET_RML_BLOCK` additionally transforms the first 64 target coordinates using a fixed covariance constructed at x=0 from the host bilateral Gram and decoder pullback. **Every arm still optimizes all 66 coordinates**, including amplitudes; zero amplitude-target perturbation is not a frozen-amplitude constraint.

For positive trace, the bilateral construction scales the 64-dimensional pulled-back curvature H by \(64/\operatorname{tr}H\), takes \((I+64H/\operatorname{tr}H)^{-1}\), clips its eigenvalues to [0.25,4], and rescales the trace to 64. This retains all 64 shape directions. Clipping and trace normalization are declared heuristics; they do not derive an exact posterior covariance. The altered target distribution is a proposal mechanism, while the coordinate mapping and quadratic regularization metric remain the original ones.

The solver, initial anchor, random draws, full-coordinate freedom, and caps of twelve forward trials and eight Jacobian states are matched across the three new arms. Each uses the same ordinary amplitude-block Schur solution of the joint GN system. True coupled nonlinear amplitude profiling is `NOT_RUN`. The endpoint minimizes its own randomized objective among visited states, including the anchor. It is then tested against the **original, unperturbed observation** and deduplicated in material space. Rejected trials and duplicated endpoints remain charged.

### Residual-guarded amplitude kernel

For each fixed original occupancy a, form P=diag(sqrt(a)), B=PCP, L=OP, and R=PE. Construct an orthonormal block Krylov basis from R, BR, and subsequent propagation, stopping at the preselected order. Project B and use a general complex Schur factorization; B need not be Hermitian or normal. For each amplitude, reconstruct an approximate full state u_hat and evaluate the full 256-dimensional residual
\[
e=R-(\delta^{-1}I-B)\widehat u.
\]
With \(\beta=\lambda_{min}[(B-B^*)/(2i)]\) and \(b_\delta=\operatorname{Im}\delta/|\delta|^2+\beta>0\), use
\[
\epsilon_w=\frac{\sqrt2\|L\|_F\|e\|_F}{\sigma b_\delta},\qquad
|\Delta\Phi|\le\|\widetilde r_y\|\epsilon_w+\tfrac12\epsilon_w^2.
\]
Acceptance requires epsilon_w ≤ 0.01 and the potential bound ≤ 0.1; otherwise the frozen policy returns an exact full-order result and charges its setup and evaluation. The actual reference error is diagnostic only. Negative beta is retained. The positive material-loss contribution supplies coercivity in the tested rectangle; neither normality nor a positive-semidefinite discrete Im(B) is assumed. Floating-point eigenvalues and residuals are not outward-rounded certificates. This kernel has not been integrated into the material candidate solver.

## 5. Completed experiments

### A. Fixed-shape cost benefit

The DEV scan uses eight exposed shapes and 30 amplitude nodes per shape. Order 32 passes 0/240 guards and therefore is assigned the cost of the already paid exact reference under the fallback policy (no separate fallback solve was launched); orders 64, 96, and 128 each pass 240/240. Order 64 saves 69.78% against full Schur at the 30-query cold workload, including the conservative DEV basis-construction cost. Failure of the order-32 guard does not prove that no 32-state approximation can be accurate: the guard is conservative and the basis family is fixed.

The subsequent probe freezes order 64 before constructing bases for ten heldout anchor shapes, then evaluates 32 new prior-amplitude nodes per shape. All 320 nodes pass both gates, with zero fallback. Maximum actual white output error is 9.13e-11; the maximum computed white bound is 8.01e-7, and the maximum potential bound is 8.32e-5. These diagnostic errors are reported to audit the guard, not as an identity-precision contribution.

| Executed workload: ten shapes × 32 amplitudes | Cold kernel seconds | Saving by order 64 |
|---|---:|---:|
| Direct LU256 | 0.238701 | 51.82% |
| Exact complex Schur256 reuse | 0.317346 | 63.76% |
| Order 64 + full residual + both gates | 0.115014 | — |

Each arm includes its required shape formation, decomposition, transforms, node solves, and guards; the reduced arm constructs only the 64-dimensional basis. Common host-cache loading and artifact writing are included in the separately recorded 1.111748-second job, not in the cold kernel table. The actual probe charges 5,120 RHS256, 2,560 RHS64, and 320 LU256. Cross-order RHS counts are kept separate. The historical host setup cost is 264 RHS576; no new RHS576 was used.

Only M=32 is an executed full workload in this probe. Extrapolations from the measured average node times show direct LU faster at M=1 and M=4, with the reduced setup amortized by M=16 among the tabulated counts. These are timing estimates, not additional executed experiments. The ten shapes are held out relative to the DEV compression scan; their material truth had already been opened elsewhere before this kernel protocol. No truth enters this probe, but it is **not a new blinded material-generalization test**. Sources: [kernel report](REDUCED_RESPONSE_FRESH10.md), [results](results/reduced_response_fresh10/RESULTS.json), and [independent audit](INDEPENDENT_QA.md).

### B. No new candidate-coverage advantage

DEV8 uses the first eight exposed R9 objects and K=16. Fresh10 uses a new frozen seed namespace, ten new objects, and K=8 for all six arms. The fresh protocol and statistics were frozen after DEV scores were opened but before fresh generation; the three new mechanisms and their solver caps were already frozen before DEV scoring. Fresh predictions were sealed before their truth was opened. Independent recomputation of 240 parent/arm/prefix rows and 463 saved endpoints found no scoring discrepancy.

| Arm, K=8 | Coverage at RMS 0.075 | Mean uniform-action MSE | Proposal RHS256 |
|---|---:|---:|---:|
| Full RML | 3/10 | 0.00835286 | 21,504 |
| Shape-target RML | 3/10 | 0.00842637 | 21,672 |
| Bilateral-target RML | 2/10 | 0.00838969 | 22,152 |
| Fixed weak8 ray | 3/10 | 0.00873832 | 11,392 |
| Dynamic weak8 ray | 3/10 | 0.00859198 | 11,384 |
| Random66 control | 3/10 | 0.00838875 | 12,520 |

The shared anchor covers 2/10. All three new arms have 100% proposal qualification, yet the bilateral arm has no coverage gain over the anchor and one fewer covered parent than full RML. Its slightly larger mean diameter, 0.15257 versus 0.14620, does not establish more useful material possibilities. Uniform and exp(-Phi) actions are declared heuristics, not posterior expectations.

For bilateral minus full RML, the paired coverage difference is −0.1 with a 98.75% bootstrap interval [−0.4,0]; the uniform-action MSE difference is +0.0000368256 with interval [−0.000234365,+0.000361607]. Neither contrast supports added benefit. Full RML minus fixed weak8 has zero coverage difference and an empirical interval [0,0] because all ten paired coverage outcomes tie; this is not population equivalence. Its action-MSE difference is −0.000385460 with interval [−0.00148291,+0.000790055]. These 20,000-resample parent-level intervals are descriptive at N=10.

The three new arms share 12-forward/8-Jacobian caps and full-coordinate freedom. Historical controls retain 6-forward/4-Jacobian caps and their radius/progress rules. Thus the whole-six comparison has equal K, **not matched optimization costs or constraints**; only the three new arms isolate target changes under matched solver settings. Proposal charges in the table omit shared generation and anchor costs. The complete fresh campaign pays 118,760 RHS256 and 4,313 LU256, including 18,016/496 for anchors, 80/10 for generation, and 40/1 of interrupted work. Two software adapter failures are retained; failed physical solves are zero. Controlled wall time is 123.7383 seconds. Sources: [scores](results/candidates_fresh10/STATISTICS.json), [cost reconciliation](results/candidates_fresh10/COST_RECONCILIATION.json), and [independent bootstrap](audit/independent_qa/INDEPENDENT_PAIRED_BOOTSTRAP.json).

### C. Conditional references remain conditional

The earlier full-order transfer study reports Schur node medians of 0.552–0.781 ms versus LU medians of 1.963–2.282 ms, with setup amortization after 21–58 evaluations in that separate measurement window. It is not the same timing comparison as the paired heldout kernel above. The initial three conditional references failed qualification; a separately frozen revision obtained `NUMERICAL_FIXED_SHAPE_TASK_REFERENCE_PASS` for all three. These are numerical task references at fixed shapes, without full-support certificates. Total transfer charges, including the initial failed attempt, are 69,784 RHS256 and 270 LU256. Full66 posterior inference remains `UNQUALIFIED_NOT_ESTIMATED`. This result neither rescues the negative material-coverage comparison nor requires an accurate Full66 posterior before finite-set imaging can be studied.

A subsequent **executed probability-workload replay** uses the same exposed development shape and all 836 saved nodes and weights of its qualified revision-2 rule. Order 64 passes 836/836 guards without fallback. Rebuilding both shape operators and evaluating the responses and quadrature costs 0.117270 seconds for the reduced kernel versus 0.287601 seconds for exact Schur, a **59.2248% saving**. Evidence, amplitude means, and event probabilities meet their original task tolerances. This adds 6,688 RHS256 and 6,688 RHS64, with no LU256 or RHS576. Historical profile, Hessian, and event-grid construction are shared existing work and are not retimed; therefore the saving applies to this fixed-node workload, not the complete adaptive integration or imaging pipeline. This `NUMERICAL_FIXED_SHAPE_TASK_REPLAY` is neither a new conditional reference nor independent material evidence, and supplies no full-support probability certificate. Sources: [replay report](REDUCED_RESPONSE_PROBABILITY_WORKLOAD.md) and [saved results](results/reduced_response_probability_workload/RESULTS.json).

## 6. Limitations and interpretation

The model is a finite discrete scalar-TM setting with a declared material decoder, host, acquisition, and noise level. Conditional transfer reuse does not establish a low-dimensional all-material realization. Noncommutation rules out particular common diagonalizations, not every useful approximation. Finite response equivalence does not establish a material gauge. Energy constraints require correct physical channel normalization; an arbitrary finite measurement matrix cannot inherit a complete scattering operator's optical theorem. Material-amplitude variation at fixed frequency is not a causal frequency-dispersion model.

Accurate Full66 posterior sampling, full nonlinear amplitude profile optimization, transfer-aware candidate allocation, task-specific generalized-eigenvector proposals, integration of the reduced kernel into candidate inversion, and exact minimal-order certification are `NOT_RUN`. The implemented bilateral covariance is one heuristic, and its negative result does not rule out every bilateral or task-aware proposal. The positive result is a kernel benefit for repeated amplitudes at fixed shapes on the declared discrete problem. It does not establish end-to-end inverse-imaging speed, new theory, unseen-shape uniform accuracy, or population-level coverage.

## 7. Closest references

The central comparisons are [TSOM](https://doi.org/10.1088/0266-5611/25/8/085003), [FFT-TSOM](https://doi.org/10.1109/TAP.2010.2103027), [RTO](https://doi.org/10.1137/140964023), [rMAP](https://doi.org/10.1137/16M1060625), [LIS](https://doi.org/10.1088/0266-5611/30/11/114015), and [goal-oriented reduction](https://doi.org/10.1137/16M1082123). The theory also distinguishes [T-matrix SVD](https://doi.org/10.1103/PhysRevA.95.053834), [radiative corrections](https://doi.org/10.1103/PhysRevA.87.012504), and [frequency-dependent oscillator representations](https://doi.org/10.1038/s41467-023-43221-2). The detailed ledger states which claims were checked from full text, abstracts, or metadata and records retrieval failures.
