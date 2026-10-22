# Supplied-sample material estimators

Only `src/p3l/estimators.py` and this audit directory were changed. No physics, GT arrays, network, old results, ledger or authorization files were read or modified. Existing P1 public decoder configuration and geometry were read to construct a two-state supplied decoder control. Those states are synthetic controls, not Maxwell posterior samples.

Public interfaces:

- `posterior_mean_variance(samples, ddof=0)`: sample mean and pixelwise `E|chi-Echi|²`. Input shape is exactly `(n,576)`, finite and nonempty; rejection repeats retain their mass. Default variance describes the empirical distribution, not an iid variance estimate or MCSE.
- `complex_geometric_median(samples, tolerance=1e-10, max_iterations=2000)`: pixelwise complex-plane geometric median. Modified Weiszfeld handles coincident iterates and multiplicity, including all-identical samples and nonunique two-point medians. Returns per-pixel value, convergence, iterations, mean absolute distance, normalized subgradient residual. Convergence is numerical stationarity at the supplied tolerance (up to floating-point floor), not interval certification.
- `component_quantiles(samples, levels=(.025,.5,.975))`: separate empirical real/imaginary quantiles; these are neither joint complex credible regions nor simultaneous spatial bands.
- `material_diagnostics(samples)`: empirical variance trace per pixel, spread RMS, each sample's RMS distance from mean, real/imaginary spatial averages and maximum pixel variance. No truth risk, coverage, ESS, escape mass or posterior qualification is inferred.
- `project_legal_bayes_action(mean, prior, initial_q, posterior_latent_samples=None, maximum_starts=4, max_nfev=1000, tolerance=1e-9)`: constrained local decoder projection by SciPy least squares. Objective is `||g(q)-mean||²/576`, without prior penalty. Bounds are unrestricted on declared Gaussian coordinates, `[0,1]` on remaining coordinates. Starts include the supplied initial state and deterministic farthest-point choices from supplied latent draws in a fixed Euclidean metric; selection uses no GT or decoded fit. A real MaterialPrior or a compatible explicit supplied decoder is accepted. Analytic decoder derivatives are used when provided, otherwise numerical differences.

Projection reports q/material, distance squared/RMS, selected optimization success/status/gradient/projected-gradient, every start's objective and solver diagnostics. It preserves the original start if numerical optimization worsens it. In that fallback case selected-gradient/status fields are `None`, while the attempted optimization remains recorded. The output is always labeled `LOCAL_NUMERICAL_PROJECTION_NOT_GLOBAL_ARGMIN`. The projected gradient uses the gradient of the full squared-distance objective; SciPy's separately reported optimality uses its half-squared-residual convention. All latent posterior draws are checked for legal bounds. No p3l/core or physics backend is imported.

`check_estimators.py` runs 17 zero-physics controls; `results.json`, `check.stdout`, and empty `check.stderr` record the final run. Seed 20261022.

Observed results:

- Symmetric antithetic Gaussian ensemble: mean error `4.58e-16`, complex median error `0`.
- Triangle complex median differs from separate coordinate medians; analytic median error `1.59e-10`, normalized subgradient residual `9.73e-11`.
- Repeated majority, all identical and nonunique two-point controls pass.
- Linear decoder projection distance squared `2.77e-32`.
- Parabola decoder with supplied states ±1 has mean outside its image. Local projection reaches squared distance `0.3750000000200632`; the initial zero state is stationary but inferior, demonstrating why multiple starts and non-global status matter.
- Supplied P1 mixture mean has a spatial range of imaginary/real inclusion-increment ratio `0.02255751`; a single legal P1 sample has one shared ratio, so this particular mean is outside the single-sample law. The two-start, 40-evaluation projection returns distance squared `2.41857e-5` and `optimization_success=false`. That residual is a feasible local upper bound; it is not the exact minimum. No budget extension was attempted.
- An intentionally one-evaluation cap remains unsuccessful; shape/nonfinite/boundary violations fail explicitly.

Reproduce from the worktree root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ../../.venv_nn/bin/python \
  research/p3l_local_posterior_imaging/AUDIT/estimators/check_estimators.py
```

These controls establish estimator semantics on supplied inputs. They do not show posterior mean headroom, calibrated uncertainty, successful MAP optimization, scientific imaging improvement, or globally optimal legal Bayes action.
