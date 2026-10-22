# P3-L sampler and diagnostic audit

These are independent probability controls. They import only the new probability modules, NumPy and SciPy. They make zero Maxwell/PDE, adjoint, decoder, training, paid-compute, or operation-ledger calls. The target is an analytically specified Gaussian field times two bounded exponential property densities; it is not a P1 posterior or a replacement prior. Results do not establish posterior calibration, scientific acceptance, material-risk headroom, solver cost, or a qualified teacher.

## API contract

`src/p3l/kernels.py` accepts **full unnormalized log density with respect to the declared coordinates**, including the Gaussian prior. The caller supplies any region indicator by returning `-inf` outside it. No method calls a solver or performs numerical differentiation. A target callable may have its own cache and cost accounting; invocation counts here are not RHS costs or an independent proof of a callable's honesty. Reusing `ChainState` is valid only with the same target, transformation and gradient.

```python
import numpy as np
from p3l.kernels import initialize, PCN, PropertyReflection, EllipticalSlice, PreconditionedMALA

rng = np.random.default_rng(7)
state, init_counts = initialize(log_density, initial_q)  # exactly one target invocation
kernel = PCN(log_density, beta=.15, property_sd=.04)
state, transition = kernel.step(state, rng)

# Alternative fixed production sweeps, with explicit separate property block:
field = EllipticalSlice(log_density)
property_kernel = PropertyReflection(log_density, sd=[.03, .02])
field_mala = PreconditionedMALA(log_density, full_log_density_gradient,
                              step_size=.1, preconditioner=frozen_covariance)
for block in (field_mala, property_kernel):
    state, metadata = block.step(state, rng)
    # Persist state.q even when metadata['accepted'] is False.
```

`Transition` has `.state`, `.metadata` and supports tuple unpacking. State stores `.q`, `.log_density` and optional cached `.gradient`. Rejection returns unchanged coordinates and density; MALA may add a valid newly evaluated current gradient to the cache. Accepting any non-gradient block clears that gradient. Each transition reports `target_evaluations`, `gradient_evaluations`, `rejected_proposals`, `repeated_state` and exactness/status. Initialization evaluation is separate and chargeable. A MALA proposal outside hard support is rejected without a reverse gradient evaluation; current gradient evaluations are still counted. MALA performs no finite-difference calls.

`PCN` supports the product proposal or a field-only proposal (`property_sd=None`). Its MH log ratio adds the Gaussian-reference proposal term. `PropertyReflection` only supports fixed diagonal Gaussian scales: arbitrary correlated increments followed by coordinate reflection generally lack the symmetric box proposal required by its simple MH ratio. Reflection permits multiple boundary crossings and never clips to an endpoint.

`EllipticalSlice` uses the full target after removing the Gaussian log prior, the full initial angular bracket, and sign-based bracket shrinkage preserving angle zero. Default `max_evaluations=None` places no internal finite limit. An explicit cap returns the current repeated state with `SLICE_CAP_UNQUALIFIED` if hit. Do not silently resume or relabel a capped run as an exact numerical reference. A finite fallback/cap is an experimental decision requiring its own invariance treatment; this module does not claim finite-cap exactness.

MALA uses a fixed SPD preconditioner and optional orthonormal `(64,r)` basis. It updates Gaussian coordinates without clipping/reflection and evaluates the complete reverse drift. Default `dimension=64` leaves properties separate. `dimension=66` supports an explicitly declared **unbounded transformed** target, such as `eta=Phi_N^{-1}(u)` with the correct Jacobian included in the target and chain rule included in the gradient. This option implements no transform itself and must not be used to reinterpret native uniform coordinates as Gaussian. `ChainState` validates finite 66-coordinate states; canonical pCN, reflection, elliptical slice and default64 MALA additionally enforce native property bounds. An original-coordinate boundary `u=0` or `u=1` has no finite probit image and requires a declared initialization treatment.

`run_chain(target, initial, kernels, draws, rng)` is a convenience fixed-sweep recorder. It records **every block transition**, so shape is `(draws * number_of_blocks, 66)`, and includes all rejections; initial state/counts are separate. To diagnose a standard sweep chain, take a consistent sweep phase (e.g. every final block), not an adaptively selected or accepted-only subset. The caller owns persistence if an external budget exception interrupts a run. `run_chain` is not a budget manager.

## Diagnostics

Input shape is `(chain, draw, *observable_shape)`, after externally removed warmup. Complex arrays append a final real/imaginary component axis. At least two chains and four draws are required. Odd sample counts discard the single middle draw when splitting. All nonfinite input raises `ValueError` rather than being dropped.

- `rank_normalized_split_rhat(values)` returns the maximum of pooled rank-normalized and folded rank-normalized split Rhat, using average ranks for ties and Blom normal scores.
- `effective_sample_size(values, method='bulk'|'tail'|'mean'|'raw')` uses split multi-chain variance and FFT autocovariance with Geyer's initial-positive and initial-monotone paired sequence. Bulk uses rank scores; tail uses the minimum of the two quantile-indicator ESS values. Antithetic ESS may exceed stored draws, capped by `N log10(N)`. This is an autocorrelation estimator, not importance-weight ESS.
- `material_diagnostics(values)` returns material means/variances, Rhat, bulk/tail/raw mean ESS, mean MCSE, and variance-observable MCSE. Variance MCSE uses autocorrelation ESS for `(h - pooled_mean)**2`; it is an asymptotic plug-in approximation and not an exact finite-sample error bound.
- `crosscheck_chains(first, second, independent=True, tolerance=frozen_scale)` compares independently seeded mechanisms' mean estimates and combines their MCSE in quadrature. Its compatibility flag requires the approximate 95% mean-difference interval to lie inside the predeclared tolerance. It does not establish equivalence of targets or apply a multiple-comparison correction. Without independent runs it supplies no combined MCSE; without a tolerance it supplies no agreement flag.

Globally constant observables have undefined ESS/Rhat (`NaN`), zero algebraic MCSE and `CONSTANT_UNASSESSED`; such outputs must not auto-pass a reference gate. Chains stuck at different constant values have infinite Rhat, conservative ESS 1, undefined MCSE and `CHAINS_STUCK_AT_DIFFERENT_VALUES`. Discrete/tied tail indicators may be constant even for nonconstant input; that tail ESS remains undefined. A caller may exempt an analytically known deterministic material observable under an explicit rule; this module does not make that decision.

## Reproduction and caution

From repository root, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ../../.venv_nn/bin/python research/p3l_local_posterior_imaging/AUDIT/samplers/check_samplers.py
```

`results.json` records versions, seed, control-specific tolerances, per-chain invocation counts, moments, Rhat/ESS/MCSE and exact-proposal ratio comparisons. The tests check repeated rejection states, pCN prior cancellation, bounded reflection, independent normalized MALA proposal densities, transformed66 callable behavior, Gaussian/bounded moments, IID and AR(1) ESS, folded-Rhat sensitivity, constants, nonfinite rejection and complex-array handling.

The thresholds used for these low-stakes controls are written into the audit script and result; they are not the package's production acceptance thresholds. No tail accuracy, multimodal exploration, reference-effect scale, physics budget, material-risk reduction or prior/noise modification is inferred. Reference qualification still requires the parent-owned frozen protocol (including at least four overdispersed chains, material ESS/MCSE gates, independent mechanisms, boundary/saturation strata and all costs). Small Rhat and a finite ESS do not rule out shared missed modes. Fixed or history-dependent tuning during production is not analyzed here. Any property/target-dependent kernel selection, adaptive geometry, stopping rule or fallback requires separate statistical justification.
