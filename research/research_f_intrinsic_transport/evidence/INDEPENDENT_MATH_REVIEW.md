# Research F independent mathematics and implementation review

Reviewed 2026-10-17 against source commit `938be85331aa36fa70954337b7e307f2e3cdad59` and completed saved F1 artifacts. This review made **0 new full states, 0 LU factorizations/solves, 0 complex RHS, 0 Mie calls, and 0 GPU calls**. It imported only NumPy/SciPy for saved-array algebra; it did not import the project solver, execute `anchor`, or change source, configuration, gates, ledgers, or closed-budget records. Only `p`, `field`, and `prediction` were read from saved state archives; saved LU/pivots were not loaded or used. No truth/offline evaluator arrays were read.

## Review outcome and binding limits

1. **No formula, dimension, conjugation, or packing defect was found in the checked F0/F1 operations.** Saved F1 `C`, `A`, and a separately contracted packed Jv agree at double-precision rounding level. This validates the algebra using the saved transfer/field; it does not independently re-solve Maxwell or certify that a cached transfer is physically correct.
2. **Gate2 remains FAIL.** Both noise-ball tests pass, both frozen stationarity tests fail. Neither successful `ftol` termination nor a retrospective boundary interpretation changes the frozen 0.001 criterion. No F2 experiment or further paid calculation follows from this review.
3. **The large unconstrained weak-background Gauss–Newton material step is mainly infeasible at the weight-ratio support boundary.** It is not evidence of an available material descent direction or an exact fiber. Small constrained linear improvements are descriptive diagnostics, not an alternate anchor qualification gate.
4. **F0 establishes a sound numerical decomposition and conditionally positive total-cokernel lower bounds, but not the missing Maxwell Hessian vectors.** Parent Gate1's narrow decision to investigate finite-noise material retraction is compatible with this evidence. A statement that the full historical three-vector decomposition was reproduced would be unsupported.

## F0 formulas, units, and gauge

Let the real material coordinate be `m(chi)=[Re chi; Im chi]/24`, so its Euclidean norm equals complex-material RMS on the 576-cell grid. `D_m` is 1152 by 8. If its retained SVD is `D_m=Q S V^T`, then `Q` is 1152 by 7 and projects onto the local material tangent. The cached whitened material Jacobian `J` is 512 by 1152. `JQ` is 512 by 7. The implementation projects by `Q(Q^T x)` and `U_A(U_A^T y)`, without confusing parameter and material metrics.

The split `d2_parallel=QQ^T D2m`, `d2_perp=(I-QQ^T)D2m`, followed by `C_coordinate=J d2_parallel`, `C_prior=J d2_perp`, is correct. `Qout C_coordinate` should vanish because the coordinate acceleration lies in the tangent image. Recorded relative defects are 2.96e-15 and 1.53e-15. The raw-cokernel/prior-cokernel agreement is 1.18e-11 and 1.03e-14 respectively; the weak case is a small projected quantity, so this relative amplification is expected.

The two positive Gaussian weights share a genuine common-scale gauge because the rasterizer normalizes by its image maximum. In Gaussian coordinates the null direction is proportional to `weight_i/(d weight_i/du_i)`; in uniform coordinates it is proportional to `weight_i/(hi_i-lo_i)`. The SVD correctly removes this null direction for local geometry. It does **not** replace the original 8-dimensional probability law by a uniform 7-dimensional material law. The latter would require integrating the original density and Jacobian along gauge fibers, whose allowed lengths depend on the weight ratio and support boundaries.

The fixed-active-max analytic second derivative in `tools/research_f_geometry.py:60-84` agrees with direct differentiation: the log-Gaussian second derivative includes `-4(r dot c')s'/s^3`, and the normalized quotient has all four terms. The local coordinate transformation adds a tangent acceleration at unchanged initial material velocity, so its normal projection is unchanged. The reported quantity is the extrinsic second fundamental form of the prior image in a fixed Euclidean material space; it is not intrinsic Riemann curvature.

All three predeclared finite-difference steps were reported. Worst relative prior-cokernel discrepancy against the analytic derivative is 1.23e-4 for weak and 5.23e-7 for ring; worst reparameterized-cokernel discrepancy is 4.01e-4 and 1.15e-6. The weak sequence is roundoff-limited at small steps rather than monotonically convergent. Active-max margins and sampled side ranks support the local calculation, without proving constant rank or smoothness over an entire finite path.

Both old tangent responses have seven retained nonzero singular values (minimum 12.8487 / 297.076 in whitened-data per material-RMS units). Hence the caches supply strong numerical evidence of no nonzero local tangent kernel; this is not an exact-arithmetic rank theorem. An approximate weak direction with nonzero first derivative is not an exact-fiber obstruction instance. Its projected second derivative is a fixed-initial-velocity finite-noise diagnostic; allowing a different velocity changes the question.

## Missing Maxwell vector and the Gate1 scope

The original executed curvature code retained the vectors `q` and `bmax` only locally and persisted scalar Maxwell norms (`src/finite_feedback/diagnose.py:193-199`). The F0 inventory and source are consistent with the reported missing-vector status. No Hessian vector can be inferred from a filename, scalar norm, or finite endpoint alone.

For any orthogonal output projector, `||Qout C_Maxwell|| <= ||C_Maxwell||`. Therefore the reverse/forward triangle bounds used by F0 are valid, even without knowing its direction:

| Old calibration anchor | `||Qout C_prior||` | Historical `||C_Maxwell||` | Lower bound for `||Qout(C_Maxwell+C_prior)||` | Upper bound |
|---|---:|---:|---:|---:|
| 250033 weak | 0.000648559834 | 0.000095209015 | 0.000553350820 | 0.000743768849 |
| 250034 ring | 1.435993714 | 0.177677343 | 1.258316372 | 1.613671057 |

These use the original unit parameter-speed direction. Dividing by the squared material speeds gives intervals [2210.1901,2970.7565] and [5401.3328,6926.6955]. The lower bounds are positive by margins much larger than the reported decoder finite-difference/reparameterization discrepancies. Consequently, **conditional on the historical Maxwell norms being correct for the same cached direction and whitening**, no possible Maxwell-vector orientation can cancel the entire prior-cokernel term. The review checked the source's direction/norm provenance, not its old paid physical computation.

The missing vector still prevents recovery of Maxwell/prior/coordinate cosines, actual cancellation, exact total response, and the full three-vector reproduction. Gate1 must remain explicitly a conditional finite-noise experiment permission, not a complete historical-vector audit, exact-fiber theorem, or nonlinear finite-step guarantee. A local quadratic coefficient alone does not bound higher-order terms at a proposed displacement.

## F1 likelihood, dimensions, and saved-array checks

For circular complex iid noise with `E|epsilon|^2=sigma^2`, each real component has variance `sigma^2/2`. Packing source-major `[Re32,Im32]` and multiplying once by `sqrt(2)/sigma` is correct. F1 takes both observation and sigma directly from the online archive and does not change sigma or gain. The residual norm threshold is the declared 95% 512-real-channel noise-ball diagnostic; after fitting it is not a posterior credible region or coverage guarantee.

`T` is 32 by 576 and represents `S L^{-1}`. `field` is 576 by 8. The complex material Jacobian is `C[(q,s),p]=T[q,p]*field[p,s]`, with shape 256 by 576. The complex parameter Jacobian is `C D`; each column is reshaped to 32 by 8 before the source-major real packing. Thus `A` is 512 by 8. The response frame is `E=V_7 S_7^{-1}` (8 by 7), with `D_m E=Q`; consequently `A E` expresses data response per unit orthonormal material tangent coordinate.

Using only final saved arrays, this review reconstructed `C`, reconstructed all columns of `A`, and compared a random parameter Jv with the independent contraction `T@(field*(D@v)[:,None])`:

| Check | ring | weak |
|---|---:|---:|
| Relative `C` reconstruction error | 0 | 0 |
| Relative packed `A` reconstruction error | 7.55e-16 | 7.91e-16 |
| Relative direct-contraction Jv discrepancy | 1.03e-15 | 1.08e-15 |
| Material gauge defect, relative to `||D_m||` | 5.01e-17 | 4.65e-17 |
| Data gauge defect, relative to `||A||` | 4.23e-17 | 5.70e-17 |

This catches packing/indexing/conjugation mismatches within the saved first-order machinery. It does not constitute a fresh physical finite-difference check or solve residual audit.

## Uniform coordinates, independent initialization, and isolation

The optimizer variable is `p=Phi(u)`, not `u`. The physical Gaussian parameters are affine in p. `Prior.derivative` divides the Gaussian-coordinate derivative by `phi(u)` to produce `D chi / D p`; there is no Gaussian-CDF saturation suppression of the optimizer gradient. Likelihood-only optimization is appropriate to the original conditional uniform parameter density (constant density in p), whereas a Gaussian-coordinate MAP objective would require its Gaussian prior term. Neither parameter MAP is automatically a material-density MAP under the many-to-one pushforward.

Initialization in the sanitized online rows is fixed before fitting: centers (-0.052,0),(0.052,0), widths 0.0725, weights 0.65. It uses no unknown first-two-Gaussian values, GT material, N3 output, or target-specific fitting. Known global center, rotation, scale, axis ratio, family, background, and amplitude are revealed conditions and materially simplify the task. This is conditional reconstruction for two backgrounds of the **same exposed development parent**, not two independent objects or a blind evaluation. Redacted parent seed and preserved unused parameters do not feed the smooth two-Gaussian decoder.

The reviewed execution path reads online observation/sigma and sanitized metadata, then public background assets; it does not read `offline_only`. The module docstring's phrase “cannot access truth arrays” is stronger than the implementation guarantees: no OS/capability boundary prevents the same Python process from accessing sibling files. Evidence supports **no truth access in the reviewed path**, not an enforced filesystem isolation claim.

The numerical bounds [1e-10,1-1e-10] keep p interior. They exclude a tiny positive-measure boundary strip from optimization; this is a numerical search domain, not an exact claim that the full original probability law was sampled. No posterior samples or posterior mass are produced.

## Frozen stationarity failure and the weight-ratio boundary

For `r=whiten(F-y)`, `g=A^T r` is the exact gradient of `0.5||r||^2` with respect to p. SciPy TRF's reported first-order optimality is `||g*v||_inf`, with `v_i=upper_i-p_i` if `g_i<0`, and `v_i=p_i-lower_i` if `g_i>0`. This is a bound-scaled KKT residual. It is not `||g||`, material-gradient norm, or a projected Gauss–Newton material step. At an exact active bound, the feasible-descent sign matters: lower-bound `g_i>=0` and upper-bound `g_i<=0` satisfy the one-sided KKT condition even with a large raw component.

The `projected_gradient_inf` diagnostic in `anchor.py:48` uses an approximate active set p<=1e-7 / p>=1-1e-7. It is threshold-dependent and differs from SciPy's continuous bound scaling. At a strictly interior point, the literal tangent cone is still unconstrained; treating a near-bound coordinate as exactly active is a numerical convention. Neither heuristic substitutes for the frozen declared SciPy gate.

| Frozen F1 diagnostic | ring | weak |
|---|---:|---:|
| Final whitened residual | 22.2395549531 | 21.5121122563 |
| Noise threshold | 23.7854490041 | 23.7854490041 |
| SciPy bound-scaled optimality | 0.0452733958 | 0.00200098872 |
| Required maximum | 0.001 | 0.001 |
| Raw gradient infinity norm | 0.1013091923 | 0.0852956603 |
| Reported approximate projected-gradient infinity norm | 0.1013091923 | 0.00471301179 |
| Termination | max 12 evaluations, status 0 | ftol, status 2 |
| Qualified | **false** | **false** |

The largest binding scaled-gradient components arise in the six shape coordinates, so suppressing correctly signed weight-bound gradients does not rescue the gate. Ftol is a relative cost-change condition, not stationarity.

Actual weights are `(0.4000000563,0.8999979072)` and `(0.4000000003,0.8999999998)`, with ratios 0.4444455404 / 0.44444444486, near the original feasible minimum 0.4/0.9=0.44444444444. In both cases `g_weight0>0` and `g_weight1<0`: unconstrained descent wants to reduce that ratio beyond its support floor. A common-scale gauge cannot remove this ratio bound. At the corner, feasible weight perturbations satisfy `delta p6>=0, delta p7<=0`, and only the zero weight perturbation is a common-scale gauge direction.

## Material-normalized and constrained saved-array diagnostics

The saved `projected_GN_coordinate_step=||pinv(A) r||` is a minimum parameter-coordinate norm despite its name. It is gauge/convention/coordinate-metric sensitive and should not be interpreted as physical displacement or information. The material RMS of its linear displacement is `||D_m pinv(A) r||`. Independently forming the orthonormal material response `B=A E` gives the same vector norm `||pinv(B) r||` up to rounding. At the current arrays both use rank 7 at rcond=1e-7.

To inspect feasibility, this review solved only the **saved linear algebra** problem `min_delta ||r+A delta||`, subject to `1e-10-p <= delta <= 1-1e-10-p`. This is a bound-constrained linearized diagnostic on the same original parameter box; it evaluates no nonlinear endpoint. Its domain limits a full linear step and correctly respects the quotient weight-ratio floor. It is not a new frozen acceptance test.

| Saved-array linear diagnostic | ring | weak |
|---|---:|---:|
| Unconstrained minimum parameter step norm | 0.235612278 | 4.59835157 |
| Unconstrained material GN step RMS | 0.000172170894 | 0.0118470629 |
| Equivalent orthonormal material-frame GN norm | 0.000172170894 | 0.0118470629 |
| Minimum material response singular value | 265.425758 | 52.8186164 |
| Unconstrained predicted residual | 22.2395002885 | 21.5027770424 |
| Box-constrained linear step material RMS | 3.29631031e-6 | 4.95415614e-7 |
| Box-constrained predicted residual | 22.2395547341 | 21.5121122554 |
| Box-constrained residual norm reduction | 2.19e-7 | 8.56e-10 |

The unconstrained material-gradient norms `||B^T r||` are 21.4787 / 66.7088. Their size alone is misleading at this one-sided support boundary: they include an infeasible ratio direction. Conversely, the small constrained predicted gain does not prove nonlinear optimality or justify changing the frozen KKT tolerance. These endpoints remain bounded-pilot failures to meet the declared anchor qualification, while retaining valid noise-compatible P0 candidate status.

## Remaining implementation qualifications

- `Prior.frame` checks material rank 7 but F1 qualification does not separately test the response rank. The saved responses here do have seven clearly nonzero singular values, so this creates no observed false pass. Future claims should explicitly report response rank and scale, and not infer information from material rank alone.
- The fixed-max derivative is locally correct away from ties. F1 does not log active-max margin/rank for every trial; saved-array consistency cannot rule out all nonsmooth crossings. This limits a global convergence interpretation, not the correct final FAIL decision.
- The full first-order Jacobian is available at all 24 F1 states and uses 32 receiver-transfer RHS plus 8 state RHS per state, consistent with the parent ledger's 960 RHS. This review did not reserve or consume any of those RHS.
- `Book` prevents repeated stage attempts and reserves states/RHS before the paid actions. Its cumulative wall recovery relies on recorded `JOB_END` events; an unhandled process kill before that event would require explicit conservative recovery, not automatic resumption. The present completed F1 run has a normal end record. No recovery was attempted here.

**Required interpretation:** preserve Gate2 FAIL and no paid F2. Report valid noise-compatible original-prior conditional candidates separately from independently qualified anchors; retain missing-vector and historical-norm conditions on F0. Final novelty, physical conclusions, and any future change of experiment require parent review.
