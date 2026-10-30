# Paper candidate: HOLD

**Working title:** *Finite-Step Resolvent Geometry and the Limits of Basin Navigation in High-Contrast Maxwell Imaging*

The scientific problem is to separate a locally accurate material derivative, a locally attracting correction region, and the ability to recover the true spatial material from a public initialization. A second question is whether data-conditioned nonlocal seeds improve the material image, rather than merely the measured-data fit.

## At most three candidate claims, with attribution

1. **Explicit qualification of finite material steps.** Exact receiving/illuminating field identities yield a sufficient invertibility/step-error ellipsoid and a source-adjoint exact Hessian action, with actual LU/RHS cost. These formulas are standard resolvent and second-order adjoint tools. The project-specific contribution is their consistent use and diagnosis in a 576-real-material VIE. No independent new Maxwell theorem is established.
2. **Deferred correction of data-conditioned rational minima changes finite-budget access to data-loss regions.** V2 historical7 and V3 fresh targets show that a higher raw uniform-response loss can lead to a better corrected data objective than the detected raw-score minimum. This belongs to the multistart/basin-hopping family and directly reuses R10 amplitude response structure. The fresh ordinary multistart beats the rational method on digit2 residual; universal search or cost superiority is unsupported.
3. **Lower residual can coexist with substantially worse material images under a fixed honest protocol.** All three new targets show worsened online material RMSE for the new bridge, and no best@K GT coverage gain. This is a useful diagnostic/negative result, not a successful new imaging method, a proof of nonuniqueness, or a statement about a posterior.

## Supporting evidence

- Parent T1–T8 and two independent theory appendices: precise finite-step order, conditional stability, rational ray, local contraction, and an injective rational-VIE box counterexample.
- Four exposed-history problems, six GT-path nodes each, twelve real endpoints; same initial Jacobian but different target-dependent directions and later trajectories.
- A fixed-face Hessian diagnosis of failed7 and true residual curvature at failed8, with failed eigensolver work retained.
- Two fixed7 shapes × three amplitudes, holding background/geometry/initialization fixed.
- Final frozen V3: 3 new anonymous targets × 6 methods, all failures visible, physical and material gates separate.

## Direct precedent risks

SMW/Schur transfer functions, Krylov projection, global line search, multistart, basin hopping, variable projection, and second-order adjoints are established mathematics. CSI has a direct material/current formulation; modern constraint-relaxation analysis already explains material-dependent Gram residual weights. R2/R10 are project-internal direct predecessors. The analytic counterexamples are mathematically elementary once the rational map is exposed; worldwide novelty is not established by the limited retrieval.

## Claims currently prohibited

No global convergence; no nine-digit recovery; no fresh GT-recovery improvement; no exact alternative data roots; no all-basin/all-root enumeration; no probability weights/posterior samples; no proven full576 local identifiability; no continuous/vector-Maxwell result; no measured cross-target-cache acceleration; no claim that a local contraction theorem reaches its own neighborhood.

## One most important missing step

Develop a non-GT spatial-structure proposal that survives full-data correction and improves **both** fresh best@K coverage and online top1 under a common, equally available structural prior and mature ordinary baseline. More uniform amplitudes or lower data loss do not fill this gap.

**Final judgment: not ready as a successful imaging-method paper.** The theory/diagnostic package is a real intermediate result. A methods paper requires the missing material-recovery evidence and closer prior-art comparison.
