# R7 A+ readout implementation note

Status: development candidates implemented; final scientific judgment belongs to parent Codex. No GT has been loaded by either candidate or development runner.

The callable contract is `readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost)`, unchanged from the original A material readout. The R7 alternating loop supplies `material_readout` only for `a_plus`; `a_material` remains the inherited exact-gradient L-BFGS-B control. `aplus_regularizer=model.reg` exposes the original prior geometry and coefficients. Material selection and the original physical corrector retain the original Maxwell objective.

## Unchanged objective and actual normalization

Both candidates minimize the same fixed-current loss

`f(m)=rho/2 * tr(R_A(m)^* R_A(m) Os) + reg(m)`.

The exact supplied Cholesky factor is used as `Jw=J@Ls,Qw=Q@Ls`, so the source-weighted state term becomes an ordinary Frobenius loss. The regularizer is the original `Regularizer` object. Crucially, its edge operator is `(m_j-m_i)/19.5`, and its smooth/TV terms are means over the full set of nearest-neighbor edges. The code extracts actual `D` entries, `smooth`, `tv`, `epsilon`, and `edges`; it does not substitute raw material differences or refit weights.

## Möbius four-color conditional MM

Set `b=beta-Gii`, `z=m/(1+b*m)`, `H=Qw+b*Jw`, and `V=Jw-z*H`. The retained residual is `R=beta*(V+z*(Goff@V))`. For one local coordinate `i`, every residual row is affine in `z_i`. Its coefficient is

- row i: `B_i=beta*(-H_i+(Goff@V)_i)`;
- row j != i: `B_j=-beta*z_j*Goff[j,i]*H_i`.

The code forms `A=R-z_i*B`, `q=rho*||B||²`, and `c=rho*sum(conj(B)*A)`. The nonconstant state loss is

`(n1*m+n2*m²)/(1+d1*m+d2*m²)`

with `n1=Re(c)`, `n2=q/2+Re(c*b)`, `d1=2Re(b)`, and `d2=|b|²`.

Its derivative numerator is `n1+2*n2*m+(n2*d1-n1*d2)*m²`. The prior tangent MM bound has coordinate derivative `h*m-s`, where each incident edge weight is

`(2*smooth+tv/sqrt((D_e*m)²+epsilon²))*D_ei²/edges`.

Thus the total stationary polynomial is the quadratic state numerator plus `(h*m-s)*(1+d1*m+d2*m²)²`, at most degree five. Batched companion eigenvalues provide all numerically real box candidates; the material is scaled by the largest box magnitude before root finding. Endpoints and the old coordinate are always included. Each color contains one coordinate per disjoint 2x2 block and shares no original prior edge. Residual updates affect every row of each block.

Candidate budget frozen before full trajectories: 2 four-color sweeps, early exit only for material max-change <=1e-8. The final exact original A loss is evaluated; any increase relative to the best of previous material and the inherited scalar seed returns that incumbent. Numerical polynomial roots are not an exact-arithmetic global minimization certificate. Observed fixed-current nonincrease does not imply original physical-objective descent, basin retention, uniqueness, or full-cycle convergence.

## Exact-adjoint projected control

The second candidate calls the inherited exact material gradient, including the derivative of the material-dependent residual transform. It takes at most 8 projected real-box BB steps with at most 16 Armijo reductions per step. The same original prior and scalar seed are used, with final exact-loss fallback. Its online criteria contain no truth metric.

## Development protocol and costs

The frozen development manifest contains the two old bars `obs_7b8d5df517f1497d` and `obs_b20d5a55677d4eca`, disclosed slit-shell `obs_40ac5b79e8a74ac3`, and disclosed easy control `obs_00f185ff5525406b`. For each, all three readouts execute the same native FFT current loop, mf3/6/10 (500 each), mf12 (1500), and original 3000-iteration Maxwell corrector. One child at a time, Apple CPU, one BLAS thread. This is disclosed development, not a new blind holdout.

`WIRING_CHECKS.json` records fixed-current tests using existing R6 stage states without any forward solve. Full runner outputs retain every auxiliary iteration loss/cost record, every stage joint current/material, stage physical states, the full corrector trajectory, original-Phi selected prediction, paid LU/RHS, and process CPU. Cost selection uses complete trajectories; local readout timings alone are not end-to-end speed evidence. Source snapshots reside in `source_freeze/`.

## One declared engineering refinement

After the frozen 12 full trajectories, parent Codex authorized one refinement: `aplus_inherited_gtol=True` (default remains False). Before solving conditional polynomials and after each completed sweep, call the original `onehop_material_value_gradient` including its full material-transform derivative and stop if the box-projected infinity norm is <= the inherited `material_inner_gtol=1e-8`. The maximum remains two sweeps. The exact evaluated final loss and old/scalar-seed fallback are unchanged. This avoids polishing well below the threshold the inherited L-BFGS-B uses.

This refinement has a new source snapshot in `source_freeze_guarded/`, a separate `GUARDED_DEVELOPMENT_MANIFEST.json`, and separate four-object trajectory output directories. It is not silently substituted into the original 12 results. All scientific claims still require the parent's independent scoring/review; original-Phi retention alone cannot establish material recovery.

## Final structural comparison: quadratic arc proposals

Parent Codex derived and authorized a final structural simplification after the guarded comparison. `aplus_arc_candidates=True` changes only the candidate-generation step: solve the exact conditional *state* stationary quadratic `n1+2*n2*m+(n2*d1-n1*d2)*m²=0` with a stable quadratic formula, then compare its feasible real roots, both endpoints, and old material under the **complete original-prior MM surrogate**. This is an approximate regularized conditional proposal; it does not globally minimize the regularized quintic surrogate. Including old material preserves surrogate nonincrease, and final exact A-loss fallback is retained. The complete-source-metric derivatives, prior, box, two-sweep maximum, inherited projected-gradient threshold, alternating current loop and physical corrector are unchanged.

`ARC_DEVELOPMENT_MANIFEST.json` declares all four predetermined disclosed objects, with the two old bars first. All four finish under the original budget. Parent Codex removed a worker-proposed arbitrary 1% relative-Phi gate before route selection: it has no physical qualification basis. Final interpretation uses full original Phi, public observation/noise residual contracts, order-of-magnitude changes, and complete paid CPU together. No sample is omitted due to a small relative objective difference. No further hyperparameter variations are authorized in this subtask. R7 primary results and truth remain unaccessed by this worker.
