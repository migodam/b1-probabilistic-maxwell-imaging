# Research G: mathematical and interface handoff

Prospective design, 2026-10-15. The statements below are elementary derivations under their explicit assumptions, not results of a new Maxwell experiment or claims that a global feasible manifold has been characterized.

## 1. Physical feasible worlds

Let x in R^(2N) be the realification of a canonical complex material field, M its fixed quadrature mass matrix, eta the permitted SHARED nuisance variables, and F_h the fixed full-discrete multi-source forward map. Define

r_h(x,eta) = W(pack(F_h(x,eta)) - pack(y)).

V_y = {(x,eta): x in X_ad, eta in H_ad, ||r_h|| <= tau_data, solver defects <= tau_solver, available shared-constitutive defects <= tau_const}.

The material feasible set is the projection of V_y onto x. A material field is not a separate world merely because an inactive latent coordinate or a nuisance optimizer changes. Conversely, two fields with different internal material geometry remain different worlds even if their receiver responses are indistinguishable.

All thresholds refer to the declared numerical experiment. They are not an automatic simultaneous confidence region, continuum-Maxwell certificate, or posterior-mass statement. Analytic interfaces that do not expose internal currents/fields must report their constitutive metric as null/NOT_EXPOSED.

## 2. The corrector, rather than the raw generator, defines discovery

Let C_y(a;xi) be a possibly randomized and budget-limited corrector applied to proposal a, with output in V_y or failure. The distribution of accepted fields is the pushforward/restriction of q(a) under C_y, not q itself and not a Bayesian posterior.

For any measurable target region B of feasible material fields,

P(next accepted field in B) = integral P(C_y(a;xi) in B) q(a) da.

Thus high raw-proposal entropy cannot overcome a corrector for which the probability on the right is zero. Optimizing proposal diversity before correction can be unrelated to useful endpoint diversity.

### A defensive exploration bound

Assume every step allocates probability at least alpha > 0 to a fixed global proposal q0, and that the full pipeline (including screening, correction and validation) reaches region B with probability at least p_B > 0 under q0, uniformly over prior histories. Then

P(B never reached in n attempts) <= (1 - alpha p_B)^n.

Proof: condition on each previous history and multiply the conditional failure upper bounds. For finitely many audit regions, a union bound applies.

This does NOT require topology classes in the generator. It also does NOT certify p_B > 0, practical sample complexity, or coverage of infinitely many modes. Surrogate rejection can destroy the assumed lower bound; defensive proposals need nonzero access to full evaluation.

## 3. Continuous latent variables need not preserve object topology

A one-parameter family can already change topology. Inside a fixed disk, define material present where r^2 >= z. For z <= 0 the disk is solid; for z > 0 it has a central hole. The indicator-valued material fields vary continuously in finite Lp norms as z crosses zero because the removed volume tends to zero. The topology of the thresholded support changes at the critical parameter.

Smooth phase fields offer an analogous threshold crossing. A low-dimensional latent is therefore not inherently topology-preserving. The real restrictions may be a diffeomorphic decoder, fixed mesh connectivity, restricted basis support, insufficient spatial bandwidth, or the trained generator's image. This example does not imply that one latent parameter represents arbitrary shapes or preserves Maxwell data.

## 4. Distances and archive utility

Use fixed real/imaginary scales s_R, s_I from physical ranges, never candidate-dependent normalization:

d_M(x,x')^2 = sum_i m_i [((Re chi_i-Re chi'_i)/s_R)^2 + ((Im chi_i-Im chi'_i)/s_I)^2] / sum_i m_i.

A multiview distance may combine d_M, multiscale smoothed field distances, signed support-distance embeddings, boundary distances and bottleneck/Wasserstein distances between persistence diagrams. Only positive, fixed weights and versioned feature maps are used in a frozen comparison.

Maintain separate counters for material diversity and structural diversity. Amplitude persistence alone does not exclude a one-cell high-contrast hole. Add explicit spatial-scale/volume diagnostics. In 2D, H1 describes holes; in 3D, H1 and H2 respectively distinguish tunnels and enclosed voids. Do not quotient by receiver equivalence: exploration within such near-equivalence fibers is the objective.

One bounded, permutation-invariant utility uses fixed reference probes z_j in material-feature space:

U(A) = sum_j w_j max_{x in A} k(d_*(x,z_j)),  w_j >= 0,

where k is nonnegative and nonincreasing, such as exp(-d^2/(2 ell^2)). The empty maximum is zero. Probes are continuous material descriptors or TRAIN-only reference fields, not topology-class labels. Their distribution is still a declared preference/bias, not a posterior. Explicit unstructured novelty and defensive global proposals are needed to avoid complete dependence on a narrow reference bank.

For a fixed finite candidate pool and fixed probes, U is monotone submodular: adding a candidate can only increase each maximum, and its marginal increase decreases when the existing set grows. This supports greedy set selection. It does not prove that a dynamic search policy discovers the right candidate pool or that U improvement lowers truth MAE.

Keep an immutable success ledger for reward accounting. Pruning a working archive must not allow the same old field to earn discovery credit again.

## 5. Proposal-retentive restoration

A prospective constrained corrector is

min_(x,eta)  1/2 ||x-x_prop||_M^2 + lambda ||P_coarse(x-x_prop)||^2
subject to (x,eta) in V_y and declared trust/retention constraints.

This is a feasibility-restoration problem, not an instruction to drive every proposal to the same regularized MAP minimizer. A finite solver returns a best-found feasible point or an explicit failure.

For a nonzero proposed edit e = x_prop-x_parent, a diagnostic retention statistic is

retention = <x_out-x_parent,e>_M / ||e||_M^2.

This alone cannot certify geometry retention: supplement it with coarse support and local spatial-scale diagnostics. A constraint on retention must be reversible/relaxable under a frozen schedule; no hole-count target is imposed.

At a regular point of an exact equality constraint F(x)=y, the first-order M-orthogonal projection onto the tangent null space is

P_T = I - M^(-1) J^T (J M^(-1) J^T)^dagger J.

This is a local linear expression. Noise inequalities, active material bounds, rank changes and finite steps require the appropriate tangent cone/active-set treatment. Clipping a nominal null direction can destroy its weak-data property.

## 6. A17-compatible compensated-edit score

At a frozen base, let e be a structural edit, Bc an allowed compensation, and C a common compact coefficient set incorporating legal bounds, trust radius and any edit-retention constraints. Define

psi_J(e) = inf_(c in C) ||r + J(e+Bc)||.

Without constraints this residual is governed by projection away from Range(JB). In an implementation, retain the constraints rather than declaring all data residuals cancellable.

Suppose the frozen full and reduced Jacobians obey the correctly instantiated A17-style factorization

J_F - J_R = D = Y Gamma^(-1) N,

and both scores use the same r, B and C. Then

|psi_F(e)-psi_R(e)| <= sup_(c in C) ||D(e+Bc)||
                   <= ||Y Gamma^(-1)|| sup_(c in C) ||N(e+Bc)||.

Proof: reverse triangle inequality for each feasible c, followed by taking infima. An anchor residual discrepancy adds its norm to the bound. A finite nonlinear endpoint additionally needs a Taylor remainder bound, e.g. L||e+Bc||^2/2 when a valid local Lipschitz constant for the Jacobian is known.

Potential use: prioritize edits whose apparent compensability is sensitive to reduced-model defects, and request a full action when the ranking is unresolved. This is NOT a proof of nonlinear feasible transport, global topology discovery or a speed improvement. Compare the cost and ranking benefit against ordinary full JVP/VJP and Krylov-based correction. If Gamma is ill-conditioned or the finite-step remainder is unbounded, the score is a heuristic, not a certificate.

## 7. RL / contextual decision problem

State s = (y, legal base, permutation-invariant success/archive encoding, prefix-only failure/frontier encoding, physics summaries, remaining cost vector, frozen model/metric versions).

A set encoder can use rho(sum_i phi(x_i, diagnostics_i)); a parent-selection head should be permutation-equivariant. Time/order tags of failures can be meaningful and must not be discarded merely because candidate order is arbitrary. A finite pooled embedding is not guaranteed to be a sufficient Markov statistic; truncated-state implementations should acknowledge partial observability.

Action a = (parent or global restart, generic operator, continuous location/scale/amplitude, 8-32 edit/compensation coefficients, permitted solver allocation).

Transition = propose -> screen -> nonlinear correct -> full validate -> update ledger/archive/frontier -> deduct costs.

Primary finite-budget reward:

R_t = U(L_(t+1)) - U(L_t) - lambda_cost c_t.

Here U evaluates verified material endpoints in the monotone success ledger. Failures and duplicates incur cost without fictitious discovery credit. Hard physical constraints are not exchangeable for a larger diversity reward. For multiple resource limits, use a constrained finite-horizon formulation or declared scalarization, rather than adding incompatible state/RHS/GPU units.

A changing archive is not inherently a nonstationary environment if the full archive and fixed update rule are in the state. Updating the metric/encoder while training changes the reward definition; freeze per experiment or explicitly version and recompute rewards.

## 8. Group credit and the GRPO zero-signal trap

For a group of independent actions at the SAME pre-group state, let S be their validated corrected outputs. Define

R_group = U(A union S)-U(A)-lambda sum_i c_i.

Giving this identical scalar to every member and then subtracting the within-group mean produces zero advantages. Adding epsilon to the standard deviation only prevents division by zero.

A useful difference contribution is

Delta_i = U(A union S)-U(A union S_without_i)-lambda c_i.

For conditionally independent actions, a frozen environment/reward, and additive per-action costs, the leave-one-out baseline is independent of action i. Therefore

E[sum_i Delta_i grad log pi(a_i|s)] = grad E[R_group].

This is a score-function identity under those assumptions, not a general theorem for clipped/normalized GRPO. Random within-group normalization and PPO clipping change the estimator. If proposals are generated sequentially, later actions depend on earlier actions and the simple leave-one-out baseline need not be independent; use trajectory-level credit or a suitable critic.

Original DeepSeekMath GRPO does not by itself encode output diversity. GAPO's known-valid-set uniformity reward cannot be transplanted as a posterior over unknown continuous materials. A no-RL group CEM update is a mandatory comparator.

When all corrected endpoints fail, true discovery rewards contain no positive signal. Fit auxiliary feasibility/progress predictors from the failed trajectories or use a carefully specified curriculum, but do not report auxiliary loss reduction as new-world discovery. Persistent all-fail groups trigger a correction/representation NO_GO.

## 9. Amortization

At matched delivered quality, let Delta C_train be extra RL training/setup cost relative to the comparator and Delta C_query = C_baseline-C_RL its per-query saving. If Delta C_query > 0,

N_break_even = ceil(Delta C_train / Delta C_query).

If Delta C_query <= 0, there is no positive cost break-even. Costs must use the same measured time, money or declared full-solve-equivalent unit; publish raw factorization, RHS, GPU and wall counters separately.

Illustration only: 50,000 extra full-state equivalents and a reduction from 64 to 40 equivalents/query imply ceil(50000/24)=2,084 queries. At 5,000 extra equivalents the threshold is 209. These are assumptions, not measured performance forecasts for A30.

## 10. Proposed sidecar interface (not asserted to match existing schemas)

`GSearchEpisode/v1` references, rather than overwrites, native candidate package IDs. Required groups:

- target: backend/version, mesh, quadrature, frequencies, sources, receivers, units, material convention, noise packing/whitening, shared nuisance policy;
- candidate: canonical complex field, parent/action ID, raw field, corrected field, creation cost receipt, actual verified prediction;
- validation: status, full residual, material legality, solver defect, constitutive defect OR null plus exposure status;
- structure: metric version, multiscale descriptors, material/structural separation, retained-edit diagnostics;
- history: correction trajectory, raw/corrected feasibility, failure reason, proxy estimate, audit propensity, trust information;
- probability: raw proposal_log_density if genuinely evaluable, interpretation=PROPOSAL_ONLY; posterior_weight=null; unknown mass=null;
- cost: unique append-only receipt IDs, factorization/state attempts, RHS by type, wall/CPU/GPU, cache status;
- restrictions: no truth fields in online handoff, no fabricated digest, no repeated global-cost snapshot summation.

Research T should supply legal tangent cones, active-set conventions, local JVP/VJP operators, finite-step conditions, any certified low-rank defect factors, and a clear distinction between a path witness and a proven connected component. G supplies actual corrected paths, structural changes, failures and budgeted reachability statistics. Research B receives real fields and joint nuisance rules, not invented candidate posterior weights.
