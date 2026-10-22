# Frozen conditional GN and importance module

`src/p3l/conditional.py` makes no calls at import time, instantiates no solver, owns no ledger, and imports no old optimization code. The parent supplies a live Target to explicitly requested operations and owns all RHS/wall accounting and execution authorization. The module's analytical controls instantiate only a synthetic linear target, with zero physical calls.

## API

```python
from p3l.conditional import make_frame, conditional_reference, importance_correct

# J and D are already obtained/counted by the parent; construction is array-only.
frame = make_frame(q0, J, weak_dimension=8, kind="prior_relative")
# Alternatives: material_generalized or posterior_material_risk, with D supplied.
# frame = make_frame(q0, J, 8, kind="material_generalized", material_derivative=D)

conditional = conditional_reference(target, y, frame, a, u,
    max_evaluations=6, initial_b=None, initial_point=None, callback=save_event)
reference = conditional["reference"]
corrected = importance_correct(target, y, frame, a, u, reference,
    samples=K, rng=rng, callback=save_event)
```

`make_frame` accepts the exact frozen native anchor `q0`, a real whitened measurement Jacobian with 64 or 66 columns, and optional complex material derivative. It selects only Gaussian64 spans. Prior-relative weak directions have smallest `Jz.T Jz` eigenvalues. Material-generalized directions have smallest generalized eigenvalues against a declared regularized material Gram matrix; its positive ridge is recorded and affects directions only. Posterior-material-risk directions are the leading eigenspace of `C^.5 M C^.5`, mapped by `C^.5`, where `C=(I+Jz.T Jz)^-1`. Material weights default to the full-cell mean squared norm; explicitly supplied weights must be nonnegative. These selection rules are declared frozen heuristics, not proved optimal or novel directions.

Every selected span is Euclidean-QR orthonormalized and completed in Gaussian64. The chart is exactly `z=q0z+Vw*a+Vs*b`, with `|det[Vw,Vs]|=1`; native `u` remains explicit. Original prior centers are `-Vw.T@q0z` and `-Vs.T@q0z`, recorded in the frame. No material Gram determinant is multiplied into the prior. Frames must be frozen outside production; no dynamic refresh is implemented.

`conditional_reference` performs bounded-count, unbounded-strong-coordinate TRF on `[white(F(q)-y), z64]`. The fixed outer coordinates are `(a,u)` and the original Gaussian prior is evaluated at the full shifted z, not at b alone. The requested 4–8 fresh evaluations are shared by this conditional optimizer; identical calls reuse points. Exact receiver Jacobians are requested from Target and counted separately. The final best point's Jacobian may require an additional derivative operation but no additional forward state. The parent must account for all actual receiver-J RHS. Callback events expose evaluation attempts, completed evaluations, accepted optimization iterations and final reference metadata.

The returned reference uses `H=Jb.T@Jb+I` and `C=H^-1` at the best finite mode candidate. **This is Gauss–Newton curvature, not exact `U_bb`.** It reports the strong gradient norm and finite-optimization convergence indicator; no small-evaluation cap implies mode convergence. The `U_mode` convention includes the original Gaussian normalization `32 log(2pi)` and omits only observation-only likelihood constants. Thus

`logZ_GN = -U_mode + (64-d)/2 log(2pi) - 1/2 logdet(H)`.

This is a conditional normal-volume approximation with respect to db, not a full posterior normalizer. The normalized Gaussian reference log-density is available through `reference_logpdf(reference,b)` and reference draws through `sample_reference(reference,K,rng)`. A failed Cholesky triggers explicitly reported eigenvalue flooring and consistently changes precision, covariance, log density and volume; it never silently claims the unmodified GN matrix. Extreme conditioning may still cause numerical failure, which propagates.

`conditional_material_moments` uses actual `prior.A`, `prior.config['field_mean']`, `prior.temperature`, fixed native property intervals, host and ROI mask. It obtains Gaussian scalar field means/SDs and calls the separately audited decoder moments module. These material moments are exact only under the declared Gaussian reference; they do not account for local truncation, GN approximation error, spatial covariance or property uncertainty beyond conditioning on u.

`evaluate_slice` returns an explicitly named `SLICE_PROFILE/SLICE` evaluation at fixed b. The conditional optimizer also returns a `SLICE_PROFILE/PROFILE` section value. Neither is labeled the original weak marginal, and a finite mode search is not an exact profile optimum.

## Explicit correction and support

`importance_correct` draws K independent Gaussian-reference b values and makes K visible Target evaluations. It computes full original log gamma (including Gaussian prior normalization) minus normalized Gaussian log r. `LOCAL_C` requires an explicit region predicate; a region cannot be silently attached to a `FULL` result. Native uniform properties remain fixed outer coordinates and contribute density one.

The ordinary average of unnormalized weights is a nonnegative unbiased estimator of the fixed-outer conditional normalizer under ordinary reference sampling and deterministic exact target evaluations. Its logarithm is generally biased. Self-normalized material means and variances are **finite-K biased** estimates and are labeled accordingly. Weighted reference particles are not claimed to be independent posterior draws. Importance ESS is a weight concentration diagnostic, with no inference about omitted support or teacher qualification. No MCSE or pseudo-marginal chain is implemented. All-zero local weights return `NO_POSITIVE_WEIGHTS`, zero ESS and no material estimate. No rescue clipping, prior change, CV numerator, sample expansion or free fallback is performed.

For `LOCAL_C`, the GN optimizer/reference remains based on the full conditional energy and is deliberately untruncated. Local-mode containment and the fact that GN volume ignores local retained mass are reported. Exact importance weights then enforce C. This provides a valid full-support importance reference when integrability holds, but does not make the uncorrected Gaussian moments or GN normal volume a local posterior calculation.

Results are dictionaries with target, reference-measure, approximation and cost-count fields; live `Point` objects are included for parent-owned cache persistence. Array-valued logweights may contain `-inf` for zero target weights and should be persisted in array artifacts. Callback events represent zero logweights/normalizers as `None` with explicit status/zero flags, allowing the parent's strict JSON serializer to record them without nonfinite values. Callback I/O and budget exceptions propagate; the module does not restart a failed physical operation.

## Control evidence

`check_conditional.py` uses a shifted Gaussian64 linear observation target and a synthetic five-cell saturating decoder with nondefault halfwidth `.37` and nonzero field mean `.13`. It checks all three orthonormal frame mechanisms, exact conditional mean/covariance, normalized Gaussian log density against SciPy, constant importance weights and normal-volume agreement, actual decoder field scaling, slice labels, and all-zero local weights. `results.json` records errors and five mode-state evaluations. These are probability implementation controls only, not Maxwell references or a claim that a physical conditional posterior is Gaussian.
