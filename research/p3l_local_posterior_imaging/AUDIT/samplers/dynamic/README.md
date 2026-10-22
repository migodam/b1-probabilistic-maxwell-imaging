# Position-dependent Gaussian MH proposal audit

```python
from p3l.kernels import PositionDependentMALA
kernel = PositionDependentMALA(log_target, log_target_gradient,
    step_size=.1, precision=frozen_precision_function, dimension=66,
    callback=optional_geometry_event_writer)
state, metadata = kernel.step(state, rng)
```

The deterministic precision callable receives a copy of the full 66-coordinate state and must return a finite symmetric SPD matrix of shape `(dimension,dimension)`. Default `dimension=66` is for an explicitly declared unbounded transformed target (for example exact probit properties); all coordinate-Jacobian terms belong in the supplied target and gradient. `dimension=64` changes only the canonical Gaussian field block and keeps two native `[0,1]` properties fixed. The module performs no probit transform and imports no physics geometry. All parameters and the precision function must be frozen as functions of state during production.

At x, precision `H=LL.T`, drift `d=(h/2) H^-1 gradlogpi`, and proposal `x'=x+d+sqrt(h)L^-T xi` are used. The normalized proposal is evaluated at both endpoints:

`log Q(x,x') = 0.5 logdet H(x) - dimension/2 log(2*pi*h) - (x'-x-d(x)).T H(x) (x'-x-d(x))/(2*h)`.

The MH ratio is `logpi(x')-logpi(x)+logQ(x',x)-logQ(x,x')`. Varying log determinants and both quadratic forms are included. No divergence drift is added or needed to define this fully MH-corrected proposal. This implementation is not described as an ideal Riemannian diffusion, and it makes no scientific novelty or mixing claim. A deterministic approximate metric affects proposals; the supplied exact target/gradient still define the correction.

Precision factorizations and log determinants are cached in a two-entry LRU using the exact complete state bytes. Current geometry and current gradient survive rejection, preventing repeated current evaluations. Cache size is bounded even during long chains. `clear_geometry_cache()` explicitly discards this computation cache; it does not authorize changing the production metric. Sharing a kernel instance across targets or changing its precision function violates the cache contract.

Transition metadata includes separate `target_evaluations`, `gradient_evaluations` and `geometry_evaluations`, forward/reverse proposal log densities, endpoint log determinants, MH ratio, acceptance/repeated-state flags and cache size. Geometry invocation counts count callable executions, not new Jacobians, RHS or floating-point work. Callables own actual cost accounting. Optional geometry callbacks emit `GEOMETRY_EVALUATION_ATTEMPT` (geometry invocation count one) and `GEOMETRY_EVALUATED` (count zero, with completed log determinant); summing callback invocation counts therefore agrees with transition counts. Cache hits cause no geometry callback. Callback or geometry failures propagate rather than silently replacing the precision.

Hard-support proposals with `logpi=-inf` are rejected without evaluating an undefined reverse gradient or proposed precision. Their current gradient/geometry evaluations remain visible, and reverse proposal metadata is `None`. Every rejected transition retains repeated coordinates. Invalid/non-SPD precision raises rather than being silently regularized; any desired deterministic regularization belongs in the parent's declared precision function.

`check_dynamic.py` compares one variable-metric transition to independent SciPy normalized multivariate Gaussian densities, checks exact support rejection and caches, and runs four independently seeded overdispersed chains against an analytic standard Gaussian in 66 dimensions. The precision is an analytic bounded diagonal function plus a small fixed rank-one term. It deliberately changes determinant with state. Each chain retains 8,000 production states after 2,000 warmup transitions, with every rejection retained. Separate callable and callback invocation totals are checked.

Run from repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ../../.venv_nn/bin/python research/p3l_local_posterior_imaging/AUDIT/samplers/dynamic/check_dynamic.py
```

`results.json` records seed, versions, counts, errors and control-specific fixed tolerances. These checks make zero Maxwell/adjoint/Jacobian/decoder/training calls. Their tolerances are not the physical campaign's frozen qualification gates. They do not establish the parent's geometry quality, global exploration, physical posterior convergence, material-risk improvement, speedup or equal-cost benefit.
