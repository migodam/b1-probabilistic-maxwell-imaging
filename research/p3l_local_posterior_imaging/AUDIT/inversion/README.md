# Independent inversion module audit

`src/p3l/inversion.py` owns only numerical optimization. It imports the public `Target`/`Point` API and never instantiates a physics backend or job/ledger. The caller owns observation and initial-cache provenance, RHS and wall budgets, independent starts and scientific interpretation.

```python
from p3l.inversion import reconstruct
best, report = reconstruct(target, y, parent,
    method="full_TRF", max_evaluations=96,
    initial_point=None, polish=False)
```

`full_TRF` uses SciPy bounded trust-region reflective least squares with exact callable receiver Jacobian. Augmented residual is `[target.white(prediction-y), q[:64]]`; augmented Jacobian is `[target.jacobian(point); I64, zeros]`. The squared least-squares cost is precisely `Phi + ||z64||²/2`, with native property bounds enforced by TRF. The default initial state is zeros64 plus property midpoint. Explicit online points are allowed without consulting truth; their existing cache belongs to the caller.

`coarse_TRF` first optimizes the original DCT modes `k,l<4` (16 Gaussian coordinates plus properties), then `k,l<6` (36 plus properties), allowing at most eight fresh states per coarse phase; it finally releases all 66 coordinates. Frozen higher modes retain the current base values, rather than changing the decoder or prior. Each phase starts from the best evaluated objective endpoint across prior phases. Accepted and rejected trial evaluations both count toward the shared fresh-state cap.

All identical native coordinate arrays reuse the point and residual cache within an invocation, including matching `fun`/`jac` calls. Receiver J is obtained through `target.jacobian`; source gradients are obtained only in optional polish through `target.phi_gradient`. Calls remain chargeable through Target's own ledger; the module's counters are invocation evidence, not replacement RHS accounting.

With `polish=True`, at most `min(polish_max_evaluations,32,max_evaluations//3)` slots are reserved from the same total fresh-state cap. L-BFGS-B is invoked only if scaled projected KKT exceeds `.001`. It uses source-adjoint likelihood gradients plus the exact Gaussian-prior gradient, native bounds, an independent hard objective/gradient-call cap (including the initial cached point), and preserves the best objective across all trials. A fresh terminal receiver J may still be needed to evaluate final KKT; its real cost is visible through Target.

KKT reports the active-box projected objective gradient, its raw infinity norm and `max(abs(projected_g)/sqrt(max(diag(J_aug.T@J_aug),1)))`. Lower-bound positive derivatives and upper-bound negative derivatives are suppressed at box tolerance `1e-8`. This is a declared diagnostic scaling, not a theorem about global optimality. The report retains SciPy status, success, message, iterations and evaluation counts for each phase, including status zero or state-cap exits. Failures preserve the best finite point when one exists; after an external budget failure the module does not request further physics for KKT.

Each invocation writes flushed JSONL evaluation-attempt/evaluation/Jacobian/iteration events and a best-point JSON under `book.jobdir/inversion/<parent>/<method>/`, followed by final `best_cache.npz` and `result.json`. Existing traces cause an explicit error to prevent overwriting earlier evidence. Target owns its own physical attempt/failure persistence. Initialization and final diagnostic Jacobians remain chargeable. A failed final cache write propagates rather than claiming persistence.

Noise-ball qualification is a separate descriptive check `Phi <= chi2.ppf(.95, real_measurement_dimension)/2`, with confidence and degrees of freedom written into the result. It is never a stopping rule, posterior hard support or optimization-convergence substitute. Changing that descriptive confidence does not change the objective.

`check_inversion.py` is a small zero-physics linear Gaussian control with two active box optima. It checks both method interfaces against the analytic bounded minimizer, invocation counts, small projected KKT, and preservation under a one-state cap. Run it with the project's Python environment. `results.json` contains its evidence. This control does not prove adequate convergence, robustness, independence of starts, or material-risk improvement on a real inverse problem.
