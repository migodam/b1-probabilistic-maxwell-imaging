# G0 small-matrix stress metrics handoff

Status: implementation and small-matrix test cases written; **NOT_RUN**. The main
thread owns the metered test wrapper, physical work, stress-class assignment,
and scientific review. No physical solve, truth/reference read, or old-W edit
was performed by this worker.

## Interface

`bgga.stress_metrics.stress_metrics(a, b, k, *, prior_precision, operator_error)`
returns separate `data_only` and `proper_gaussian_prior` dictionaries.

- A is `(observations, wall coordinates)`, B is `(observations, nuisance
  coordinates)`, and K is `(task coordinates, wall+nuisance coordinates)`.
  Inputs must be real and finite. Caller owns realification, whitening, physical
  units, material coordinates, and task weights. Zero nuisance columns are allowed.
- `operator_error` is a mandatory finite nonnegative **absolute spectral-norm
  error**, common to A/B; caller supplies a conservative common bound/estimate.
  It is not divided by an operator norm or replaced by machine epsilon.
- `prior_precision` is a mandatory, exactly symmetric, positive-definite full
  joint precision in `[wall, nuisance]` order. A singular/improper prior or
  numerically noninvertible posterior raises `ValueError`; there is no jitter,
  pseudoinverse, fitted precision, or damping-as-prior.

Separate entry points: `svd_rank_audit`, `data_only_diagnostics`,
`proper_gaussian_diagnostics`, and `finite_pair_diagnostics`.

## Preserved diagnostics and limitations

1. SVD threshold rank is `count(s > operator_error)`. The lower/upper effective
   rank at that same threshold uses `max(0,s-e)` and `s+e`. A disagreement is
   `UNRESOLVED`; this is an error-resolution rank, not a true algebraic-rank
   certificate. Singular values, their error envelopes, threshold, support and
   excluded directions are retained. No machine-epsilon rank rule is used.
2. Data-only absolute wall information is `A.T@A`; Schur information is
   `R.T@R`, `R=(I-U_B U_B.T)A`, with B selected by that threshold. The full
   absolute/Schur spectra, generalized retention on A's effective support,
   principal angles, sum of angles, bare sin-squared values/sum, and retention
   trace are exposed. Retention outside A's support is undefined. Empty support
   has empty spectra and `None` sums, never an invented zero summary.
3. Threshold-rank-conditional numerical values remain available under
   `UNRESOLVED`. The input operator error is preserved; **no projector/Schur or
   covariance error envelope is claimed**. Raw principal cosines and any
   arccos-only domain clipping are retained. Floating-point roundoff is not
   certified.
4. The separate proper Gaussian label reports joint covariance, wall/nuisance
   cross covariance, wall Schur spectrum, task covariance, and the local
   wall-conditioned/omission task covariance. It is not bare sin-squared
   geometry, a finite-counting-prior posterior, or a physical certificate.
5. `finite_pair_diagnostics(current_difference, task_difference, *,
   current_error, task_error, task_norm, task_contract_id)` accepts only already
   computed caller differences. Errors refer to the pair difference in the
   declared norms. Task norm must explicitly be `euclidean` or `max_probe_l2`;
   the latter requires `(probes, real measurements)`. It retains vectors,
   magnitudes, errors, norm envelopes and task contract. Task/current ratios are
   withheld with a distinct status for a zero denominator, an error envelope
   reaching zero, or floating-point range overflow.
6. No class, gate, novelty, continuum, or true-posterior decision is made here.
   The tests only exercise finite supplied matrices and direct supplied finite
   differences. Their execution is reserved for the main thread's paid wrapper.

Files: `src/bgga/stress_metrics.py`, `tests/test_stress_metrics.py`, this summary.
