# Gaussian conditional decoder moments: zero-physics audit

Status: **PASS_GAUSSIAN_REFERENCE_MOMENTS_ONLY**. No Maxwell solve/RHS,
training, external service, or provider change was used. This audit validates
the supplied-reference integration module; it does not validate a posterior
approximation, a local-region normalizer, or an imaging-risk claim.

## Reading and source contract

The complete formal `Research_P3L_Report_ZH.md` (1,328 lines), complete
`technical_appendices/B_Probability_and_Algorithms.md` (796 lines), and
`src/p1_material/prior.py` were read. Source references below refer to the
unchanged package under `source_package/Research_P3L_Package/`.

- Report lines 363–405 and probability appendix lines 334–389 define a supplied
  Gaussian strong conditional reference, fixed outer properties, and the
  exact C2 quintic. Its general-halfwidth implementation equals the existing
  `saturation(field / halfwidth)`; the frozen P1 halfwidth is 0.5.
- Appendix lines 393–406 specify complex modulus variance and the limitation
  to pixel marginals/MSE trace. The module returns
  `E[abs(chi-Echi)**2]`, including the deterministic ROI mask squared in the
  variance, and `E[abs(chi)**2]` including the complex host cross term.
- Appendix lines 412–420 use property/field independence specifically for the
  original prior. This does not establish posterior independence.
- Report lines 419–445 and appendix lines 432–485 concern exact weighting and
  pseudo-marginal correction. The present module does not implement them and
  does not make their correctness claims.

## Implementation

Only `src/p3l/moments.py` was added outside this audit directory. Its public
interfaces are:

1. `gaussian_saturation_moments(mean, sd, *, halfwidth=.5)` returns fraction
   mean, second moment, and variance in a `SaturationMoments` dataclass.
   All numerical inputs broadcast. SD zero is deterministic. Well-conditioned
   cases use the truncated-Gaussian polynomial recurrence; difficult cases use
   positive, scaled one-dimensional quadrature with exact polynomial evaluation.
2. `gaussian_material_moments(mean, sd, inclusion_increment, *, host=.1+.04j,
   mask=1., halfwidth=.5)` returns complex mean, real nonnegative modulus
   variance, second absolute moment, and masked fraction moments.
3. `linear_gaussian_field_marginals(matrix, mean, covariance, *, offset=0.)`
   supports correlated and singular positive-semidefinite latent covariance.
   It returns field means/SDs without claiming spatial independence.

Fixed material properties must be finite. They may broadcast across **explicit
outer states**; the caller remains responsible for supplying the declared
bounded property values. Uncertain dependent properties must be conditioned on
and then combined using total expectation/variance. No independent-product
shortcut is provided.

The fallback computes centered polynomial differences and scales by the local
polynomial response. This retains a variance of `3.515625e-300` at field mean
zero and SD `1e-150`, although the raw first/second moments round to 0.5/0.25.
Reflection handles near-one means without subtracting raw moments to recover
variance. Underflow below float64 subnormal range is unavoidable. Gaussian
integration cuts at 38 standard deviations; this is numerical quadrature,
not an interval-arithmetic certificate. On this Apple runtime `longdouble`
has no extra mantissa precision, so stability does not rely on extended
precision being available.

## Independent checks and observed errors

`results.json` and `quadrature_cases.json` contain the complete final results.

- 183 independent SciPy quadrature cases using the original P1 decoder:
  maximum absolute first/second-moment error **3.3306690738754696e-16**.
- 400-digit independent polynomial integration validates three nearly
  deterministic interior variances, including the `3.515625e-300` case.
- 90-digit independent direct integration validates seam-adjacent tiny SD
  cases and a far-tail mean `6.629577372929671e-306`. One direct float64 SciPy
  reference issued roundoff warnings at mean `-0.50000001`, SD `1e-8`;
  those warnings are preserved in JSON. High-precision integration resolves
  the comparison to about machine relative precision at that seam.
- 600,000 draws of an eight-dimensional correlated, rank-five Gaussian,
  mapped to twelve pixels: maximum mean discrepancy **1.9612 standard
  errors**, maximum variance discrepancy **2.4043 standard errors**.
- Three additional 600,000-draw difficult-tail cases have maximum first/
  second-moment discrepancy **2.4604 standard errors**. Ordinary Monte Carlo
  cannot validate events with astronomically small probability; those cases
  are covered by scaled quadrature/high precision instead.
- Independent direct complex-variance quadrature agrees within `3e-12` in
  all three specified cases; actual differences are at rounding scale.
- Broadcasting, nondefault halfwidth, deterministic zero SD, extreme finite
  input scaling, nonfinite inputs, negative SD, asymmetric covariance and
  indefinite covariance are checked.

During development, the added high-precision seam check exposed a fallback
scaling error: at mean `-0.49999999`, SD `1e-8`, scaling only by SD allowed
relative errors around `4.50e-7` in mean and `1.20e-5` in second moment.
Scaling by the local polynomial response corrected it. Final high-precision
checks pass. The same review corrected multiplication order to preserve
subnormal second moments in the far tail.

## Mathematical cautions for final Codex review

The supplied-reference formulas themselves are correct under their stated
assumptions; no fatal mathematical issue was found. Their use requires care:

1. A Gaussian strong reference restricted to a nonlinear LOCAL-C region
   generally no longer has Gaussian scalar field marginals. Untruncated
   Gaussian moments must not silently become LOCAL-C posterior moments.
2. The module returns diagonal complex covariance only. Correlated linear
   latent fields do not make transformed pixels independent; simultaneous
   intervals, topology probabilities and full spatial covariance require
   additional integration.
3. Property/field dependence matters. The audit supplies the concrete bounded
   example `F=Z`, `U=Phi(Z)`, `A=.8+U`: correct `E[A*phi(F)]` is
   `0.7722348506749165`; the independence shortcut gives `0.65`.
4. A correct Gaussian decoder integral establishes neither conditional
   Gaussianity of the actual Maxwell posterior nor a net runtime advantage.
   Mode/Hessian/reference construction and actual posterior corrections still
   have to be charged and validated by their owners.

## Reproduction

From the worktree root, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ../../.venv_nn/bin/python \
  research/p3l_local_posterior_imaging/AUDIT/moments/check_moments.py
```

The audit requires NumPy, SciPy and mpmath (already available in that runtime).
The production module requires only NumPy/SciPy. Seed: **20261022**.
