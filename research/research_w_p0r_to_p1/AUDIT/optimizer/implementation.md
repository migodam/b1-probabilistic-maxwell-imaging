# Generic constrained optimizer implementation

Scope: `src/wr/optimization.py` only implements numerical optimizers. No Maxwell,
forward PDE, ground truth, physical experiment, or scientific acceptance claim is
included. Geometry legality, physical area constraints, authorization, accounting,
and validity of the prior remain the adapter's responsibility.

## Interface

- `trust_region(evaluate, x0, lower, upper, prior_lambda=0, prior_center=None,
  max_iter=30, radius=.3, callback=None, ...)`.
- `lbfgsb(...)`, also exported as `l_bfgs_b`, is the L-BFGS-B alternative.
- `evaluate(x, jac=True)` returns a real whitened residual vector and dense real
  Jacobian with shape `(number_of_residuals, number_of_parameters)`. The generic
  trust-region trial calls use `(r, None)` with `jac=False`; provisionally accepted
  candidates require a full-FJ refresh before becoming anchors. L-BFGS-B requests
  full FJ for every new point. Neither method accesses a forward model directly.
- Objective is `0.5*r.T@r + prior_lambda/2 * ||x-prior_center||^2`.
  Default center is zero. Bounds broadcast to `x0`; infeasible initial points,
  invalid bounds, complex/nonfinite outputs and incompatible shapes are rejected.
- Return value is a SciPy `OptimizeResult`, including `x`, `fun`, objective `jac`
  (gradient), `residual`, `projected_gradient`, `nfev`, `njev`, `nit`, `success`,
  `message`, `history`, and fully evaluated accepted `anchors`. `nfev/njev` count
  attempted calls / attempted Jacobian calls, including failed ones. No automatic restart is performed.

## Trust region and convergence

The model is the Full-J Gauss-Newton quadratic including the exact quadratic
prior. SciPy `lsq_linear` solves each linearized least-squares problem using the
intersection of parameter bounds and a coordinate trust-region box. Fixed
parameters are eliminated. By default a Euclidean radius also applies: repeated
bounded solves bracket and bisect a nonnegative Tikhonov damping multiplier until
the step is inside that ball. Damping is used only to compute the step, not in
actual or predicted objective reduction. No unconstrained-step clipping is used.

Actual/predicted reduction accepts a decreasing step at ratio >= 0.1 by default.
Ratio < 0.25 shrinks the radius by four; ratio > 0.75 near the boundary permits
radius doubling. Failed evaluations, failed/nonpositive models and rejected
trials backtrack by resolving at a smaller radius. A rejected/failed candidate
never overwrites the retained anchor. Every new trial goes through `evaluate`. Residual-only rejected trials avoid
Jacobian work. A provisionally accepted trial is refreshed through
`evaluate(..., jac=True)`; the adapter may reuse exact-point cached forward
fields. The refreshed full objective determines final acceptance. A failed or
call-capped refresh retains the previous anchor. This module does not estimate
PDE/RHS cost from Python call counts.

Convergence requires the true nonlinear objective gradient at the current anchor:

`||x - clip(x-gradient, lower, upper)||_inf <= gtol_abs + gtol*initial_mapping_norm`.

The criterion includes the prior and bound KKT conditions. Residual norm alone
cannot trigger success. Radius floor, objective floor, exhausted backtracking,
iteration limit, evaluation failure and call cap have distinct messages and are
not labelled gradient convergence. The final accepted iterate is checked even
when the iteration limit is reached. No objective floor is enabled by default.

## Evidence and failure recording

`callback(event_dict)` receives evaluation starts, completions, failures, accepted
and rejected trials, anchor updates, model failures and call-limit events. Every
row is also retained in `history`; rows identify calls and candidate coordinates.
A callback must provide its own persistent append-only write if durable recording
is needed; this module does not write campaign ledgers.

All adapter exceptions are logged. Ordinary exceptions are rejected trials in
trust-region optimization; L-BFGS-B terminates at its last accepted anchor because
its line-search interface does not expose a robust failed-evaluation recovery.
Pass budget/authorization exception classes through `fatal_exceptions`, e.g.
`fatal_exceptions=(BudgetStop,)`, so they propagate immediately after logging.
Callback exceptions always propagate. KeyboardInterrupt/SystemExit are not caught.

L-BFGS-B uses analytic full-FJ gradients and an optional exact, consecutive-point
cache (`np.array_equal`, no approximate match). Its accepted-iterate callback
labels all pending successful evaluations accepted/rejected; unfinished trials
are rejected on exit. Failed evaluations have their own failure event. SciPy's
function-reduction termination alone is not reported as gradient convergence;
the returned anchor must pass the same explicit projected-gradient criterion.

## Validation status

Source AST parsing passed. An explicitly synthetic identity illustration (zero
PDE calls) ran using the existing `Gaussian/.venv_nn/bin/python` runtime. The
three-parameter problem used target `[2, -0.5, 4]`, box `[-1,-1,1]` to `[1,1,1]`
(one fixed coordinate), initial `[0,0,1]`, and prior lambda 0.2 centered at zero.
Both methods reached the analytic box/prior optimum `[1,-0.4166666666666667,1]`,
objective `5.220833333333333`, projected-gradient norm `5.55e-17`.
Trust region used 7 evaluator calls / 4 Jacobian calls / 4 retained anchors;
L-BFGS-B used 3 / 3 / 3. Evaluator counters were checked against actual calls.
Both returned `projected_gradient_converged`. This is one elementary synthetic
illustration, not a physics validation or an optimizer robustness campaign.

The initial default/Homebrew runtime attempts failed at import because NumPy was
unavailable. No packages were installed. No physical correctness or performance
claim follows from this illustration.
