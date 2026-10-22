# Saved-array decision bounds

Only `src/p3l/decision_bounds.py` and this audit directory were written. No physical arrays, truth, backend, source provider, old result, ledger, or authorization file was used. All checks are small synthetic complex matrices. The implementation does not perform or authorize saved-inverse setup; the parent owns any physical RHS and final assessment.

## API and interpretation

`stability_from_saved_inverse(L0,R0,Lnew)` computes:

- An induced2-norm upper estimate using `sqrt(||A||1 ||A||inf)`, including explicit entrywise residual-error inflation.
- `rho0_upper` for `||I-R0 L0||2` and an upper norm for `Delta=Lnew-L0`.
- When `rho0<1`, `B0=||R0||/(1-rho0)` bounds the old inverse norm under the stated numerical model.
- Only when `kappa=B0||Delta||<1`, an inverse bound `Bnew=B0/(1-kappa)` and positive lower bound `beta=1/Bnew` are returned. Otherwise `beta=None`, status `BOUND_UNAVAILABLE`. This does not say the new matrix is singular.

No inverse or solve is performed in this stability function. R0 may be an approximate saved inverse: its residual is explicitly assessed. It is not presumed exact because it is labeled an inverse.

`dual_corrected_prediction(L,rhs,S,V,Z0,beta,*,observation=None,sigma=None)` uses supplied finite matrices of dimensions `L(n,n)`, `rhs(n,m)`, `S(p,n)`, `V(n,r)`, `Z0(n,p)`. Beta may be a positive scalar or the stability-result dictionary. Unavailable beta returns no bound and performs no matrix products. A failed projected solve returns unavailable without asserting singularity of L.

It forms the numerical Galerkin state `j=V solve(V* L V,V* rhs)`, residuals `rp=rhs-Lj`, `Rd=S*-L*Z0`, and corrected output `fc=Sj+Z0*rp`. For exact residuals and exact correction arithmetic,

```
F-fc = Rd* L^-1 rp
||F-fc||F <= ||Rd||F ||rp||F / beta.
```

The implementation adds residual and correction rounding errors, including `||Z0*||2` times the primal-residual arithmetic error. Numerical projected solve accuracy requires no additional theorem: its returned j is simply an approximate state whose residual is bounded afterwards. The returned bound applies to that numerical state/correction, conditional on beta and floating-point assumptions.

Providing both observation and fixed sigma forms a likelihood interval using complex Frobenius norm and real whitening scale `sqrt(2)/sigma`. Residual-center arithmetic is added to the radius. With `t=||rwhite||`, `e=eps_white`, the endpoints are outward-inflated versions of `[max(0,t-e)^2/2,(t+e)^2/2]`. This agrees with the original complex noise likelihood `||F-y||²/sigma²`. No GT, exact candidate prediction, or offline reference error is accepted by the evaluator to decide fallback.

`same_uniform_decision(current_interval,candidate_interval,c,logu)` combines current/candidate likelihood intervals with the supplied log prior plus reverse/forward proposal ratio c. It returns ACCEPT, REJECT, or REFINE and echoes `logu`, requiring the caller to retain that same proposal/uniform on refinement. It generates no randomness. Unavailable bounds yield REFINE; nonfinite ratio arithmetic also yields REFINE. Boundary equality remains unresolved. It does not automatically call Full or spend unallocated RHS.

Every bound returns the label `CONDITIONAL_FLOATING_POINT_BOUND_NOT_INTERVAL_CERTIFIED`. An exact MH guarantee remains conditional on genuinely enclosing intervals and unchanged proposal/random number; observed toy agreement alone does not establish those properties for funded TM24 matrices.

## Floating-point inflation assumptions and limits

The code uses unit roundoff `u=eps_float64/2`, `gamma_k=ku/(1-ku)`, and the entrywise matrix product model

```
|fl(A B)-A B| <= gamma_(8k+4) |A| |B|
```

where k is inner dimension. The positive absolute-product calculation is inflated for its own multiplication/summation, and complex absolute values, additions/subtractions, norm sums and scalar divisions/square roots receive additional gamma factors and outward `nextafter` movement. Upper norms are computed with scaled Frobenius evaluation to reduce avoidable overflow. The factor8k+4 is a conservative operation-count model for complex real-arithmetic dot products, not a verified characterization of every possible BLAS kernel.

Assumptions are IEEE float64 round-to-nearest, ordinary stable matrix multiplication/dot-product evaluation, finite intermediates and no material underflow/overflow. Hardware, compiler/BLAS reassociation, subnormals/flush-to-zero, initial saved-array errors relative to a physical model, and the rounding mode are not independently certified here. An actual interval-arithmetic library or independently justified reproducible error analysis is required to upgrade this label. Nonfinite resulting bounds are unavailable, not usable infinity certificates.

Counts include ordinary matrix products, additional absolute matrix products for error estimation, matrix action labels/shapes/RHS-column counts, projected solve columns and projected linear solves, induced/Frobenius norm evaluations, and residual-center norm evaluations. Physical solve/RHS/factorization counts remain zero. Small projected array solves can cost wall time even though they are not new Maxwell RHS. The conservative absolute-product bound may be expensive or loose on TM24; no speed advantage is claimed.

## Controls and observed results

`check_decision_bounds.py`, `results.json`, `check.stdout`, and empty `check.stderr` preserve32 controls (seed2026102251):

- Analytic2x2 diagonal example: returned beta `1.8999999999999744`, below actual sigma_min2.1; kappa about.05.
- Stable `Lnew=3I` with saved `L0=I` has kappa>=1 and returns unavailable; it is not called singular.
- Deliberately bad saved inverse fails its residual sufficient condition.
-20 synthetic nonnormal complex6x6 cases: actual output error is below the bound and exact toy likelihood is inside its interval. Maximum observed actual-output-error/bound ratio about.244. Those checks are numerical evidence, not proof of all possible arrays.
-100 shared-uniform toy comparisons resolve consistently after tightening with analytic/small-matrix reference values. Initial wide-bound resolutions also match direct MH; unresolved cases retain the uniform. Test-only exact small solves are never passed into evaluator fallback selection.
- Exact-arithmetic dual cancellation still returns a positive roundoff allowance.
- Known accept/reject, uncertain/unavailable refinement, no-beta/no-projection-work, and zero physical counts pass.

From worktree root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ../../.venv_nn/bin/python \
  research/p3l_local_posterior_imaging/AUDIT/decision_bounds/check_decision_bounds.py
```

This is conditional saved-array algebra and a small control suite. It is not a TM24 certificate qualification, physical posterior result, finite-precision proof, performance measurement, or authorization for576-RHS setup.
