# R12 mathematical scope and derivations

## Inherited physical contract

Time convention exp(-iωt); finite VIE contrast current J=mQ, Q=Einc+GJ. Material is real isotropic scalar shared across all sources and, in vector3D, all three components of each cell. With X=J/s the metric is ||ΔX||²+||ΔQ||². The Pro report derives exact box-constrained shared-slope TLS projection and fixed H=I+s²G*G background projection. These are inherited constructions, not new R12 theorems.

## Exact first-step control (DERIVED)

Let z0 be the full physical fields of the common legal m0. Then z0 belongs to C, hence PC z0=z0, regardless of a nonunique material label at dark cells. For ω=1,

DR(z0)=z0+PA(2PCz0-z0)-PCz0=PAz0=AP(z0).

Thus AP and DR have identical first updated state when supplied the same complete physical initialization, scale, data ball and exact projectors. Differences thereafter can be attributed to reflected memory, subject to identical implementation and stopping. A first-step discrepancy indicates different initialization/projection or numerical implementation, not a basin advantage. At a physical point inside the data ball both remain fixed. For ω≠1 the first DR step is a relaxed affine projection, so the equality is specific to ω=1.

This is an elementary control identity; no dedicated repeated precision experiment is justified. It guides checkpoint interpretation. It does not prove escape, convergence or material recovery.

Implementation metadata caveat: `run_projection` records iteration0 as a baseline physical audit before applying either projector. Its stored `gap=0` and A-state arrays are initial-state placeholders, not a measured A–C consensus gap or an evaluated A projection. If the initial material is outside the data ball, that placeholder A state is not in A. Exclude this baseline gap from feasibility/convergence claims and from the conditional shadow-gap bound below; use actual projected shadows from subsequent audits. This clarification changes documentation only, not any iteration or saved result.

## Conditional local projection regularity

For each cell set a=||x||², b=||q||², c=Re〈x,q〉. The stationary-slope equation is c α²+(b-a)α-c=0. Multiple closest admissible slopes can occur; near an eigenvalue degeneracy of the real 2×2 Gram, the selected slope can be unstable. Box-endpoint switches make a globally smooth material-shadow map unavailable. Consequently a local smooth tangent-space DR rate cannot automatically explain finite trajectories from boundary-saturated bad materials. Numerical diagnostics must distinguish shadow branch changes from physical-gradient stationary behavior.

## Selection boundary

The propagation/data projector enforces a data ball, not the original material prior. Candidates are therefore evaluated with the same original full Maxwell Phi and prior as controls. This protects comparison of search paths but does not prove that minimum Phi identifies ground truth. Noise-feasible incorrect materials, model mismatch and failed search are separate outcomes.

## Shared-material projection (inherited exact finite-dimensional construction)

Write the scaled cell state as complex vectors x,q containing all illuminations and, for vector Maxwell, all three components. Set alpha=m/s. For fixed real alpha the closest constitutive pair is

    qhat=(q+alpha*x)/(1+alpha^2), xhat=alpha*qhat.

The minimized squared distance is ||x-alpha*q||²/(1+alpha²). Its stationary equation is the quadratic above. Evaluate all real stationary roots in [lower/s,upper/s] and both endpoints, retaining the least distance. Degenerate cells can have nonunique minimizers. This is a small real shared-slope total-least-squares problem, not separate material fits per source or component. It applies to real isotropic contrast, not automatically to complex or anisotropic tensors. The scaling s changes the product-space metric and hence the algorithm; it is not a free physical identity that leaves trajectories invariant.

## Fixed propagation/data projection (inherited)

For input (x,q), put B=sG, H=I+B*B and xbar=H^{-1}[x+B*(q-Einc)]. Eliminating Q gives a strictly convex quadratic in X with Hessian H. For each source let T=WsS include its legitimate receiver mask and whitening, and y=WY. The data constraint across all sources is one global ball. The constrained minimizer is

    X_lambda = xbar - lambda H^{-1}T* (I+lambda T H^{-1}T*)^{-1}(T xbar-y).

Use lambda=0 if the unconstrained residual is inside the ball; otherwise determine the single nonnegative multiplier from the global residual norm. Different source port maps share lambda, not necessarily the same small eigensystem. Residual components orthogonal to range(T) cannot be removed: a ball below this floor is unreachable. The H factorization depends on fixed background and scale, not on the material iterate. Its setup, port decompositions, and every triangular RHS are nevertheless paid costs. Dense storage is quadratic and dense setup cubic; this implementation is not yet a large-scale FFT projection solver.

## What a small shadow gap does and does not certify (DERIVED)

Let c=(Jc/s,Qc) belong to C with material m, a=(Ja/s,Qa) belong to A_delta, and g=||c-a|| in the product Frobenius norm. Assume the finite Maxwell matrix I-Dm G is invertible. Define e=Qc-Einc-GJc and the full physical current Jphys=(I-Dm G)^{-1}Dm Einc. Then

    ||e|| <= sqrt(1+s²||G||²) g,
    Jc-Jphys = (I-Dm G)^{-1}Dm e.

The first inequality is the norm of the block map [-sG,I] applied to c-a, using Qa=Einc+GJa. The second follows by substituting Jc=Dm Qc, and is an exact finite-matrix relation. With a common port map and whitener,

    ||W(S Jphys-Y)|| <= delta + ||WS|| [s + ||(I-Dm G)^{-1}Dm|| sqrt(1+s²||G||²)] g.

For source-specific ports the same bound holds using the corresponding block-diagonal operator norm. It is a conditional bound, not a computable uniform certificate unless the resolvent factor is independently bounded. Near a poorly conditioned material response, a visibly small shadow gap can coexist with an unacceptable physical data error. No material-distance-to-GT bound follows. At g=0, uniqueness makes the constitutive shadow an exact finite-model physical solution; it can still be the wrong material or a poor approximation to the continuous experiment.

This explains why R12 audits material shadows with the original forward solver rather than accepting a small gap alone. It does not warrant repeatedly recomputing a norm bound or testing this identity to machine precision.

## Basin changes without convexification (scope statement)

The DR state z is not restricted to physical field/current pairs. The material shadow PCz and the propagation shadow PA(2PCz-z) therefore need not be a physical descent path. An ordinary material stationary point outside the data ball is not thereby a DR fixed point. This permits escape but proves neither that escape occurs nor that it finds GT. AP shares the exact projectors, so AP versus DR from identical physical initialization isolates the reflected iteration beyond the first step. The R11 endpoint experiments provide empirical evidence for this distinction in three exposed cases.

The scalar material objective Phi and its stationary points have not been convexified or removed. Nonconvex constitutive projections, boundary switches, and poor observability remain. If a wrong material already satisfies the same legal constitutive and data constraints, its physical pair is fixed under AP and DR. No reflection parameter can turn that statement into a GT-selection guarantee without additional information or a different estimator.

## Corrector versus constrained selection (DERIVED decision boundary)

For Phi=D+R, a successful descent may have Delta D>0 and Delta R< -Delta D. It can leave a previously satisfied noise ball while lowering Phi. Even if both D and R improve, material error need not improve. R12's exposed 3D and Austria diagnostics instantiate both phenomena; these are properties of the estimator and its data/model, not proof of an incorrect gradient.

The primary experiments keep inherited Phi selection. A separately registered confirmation rule selects the minimum same R among already paid, full-physical, legal, CN99-feasible candidates, breaking ties by Phi. This constrained candidate rule changes selection semantics and is reported separately. It supplies no missing candidates, no extra observations, and no GT oracle. If the candidate pool has no qualifying point the secondary estimate is unavailable. Its registration timing and exact cost eligibility are recorded in SECONDARY_SELECTION_FREEZE.json.

## Originality boundary

DR/RRR, total least squares, nonphysical auxiliary states, fixed quadratic elimination, and resolvent error propagation are established mathematical tools. R12 does not claim any of these identities as a new theorem or a global convergence result. Potential contribution lies in the physically shared-material implementation, the controlled evidence on actual wrong endpoints, and the identified separation of search success from physical material selection. Closest-method comparisons and literature-access limitations are recorded in PRIOR_ART_AND_NOVELTY.md.


## Finite-noise ambiguity diagnostic (elementary consequence, not a novelty claim)

After sealing, saved clean data allow checking whether both the generating legal material m_true and a returned legal material m_hat satisfy the same observed-data ball. If both do, the triangle inequality gives ||W(F(m_true)-F(m_hat))|| <= 2delta. This constructs two data-compatible materials for that observation and noise allowance. It does not assert equal noiseless responses, global nonuniqueness, or impossibility of distinguishing them using a different prior or new data. It is not an online GT selection rule. All24 saved GT observations in confirmation satisfy the frozen CN99 radius; the11 DR noise-feasible objects that miss the fixed IoU description therefore exhibit finite-noise compatibility with distinct materials. Their amplitudeNRMSE≤.15 and relativeL2.324–.600 must remain visible rather than calling them absent images.
