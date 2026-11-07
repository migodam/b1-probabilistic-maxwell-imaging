# Standard multiscale discrepancy control

This addition is a standard multiscale discrepancy control, not a new reconstruction algorithm. It tests a stronger ordinary noise-stopping baseline. Final-only stopping at24 may allow earlier scales to continue fitting noise, but that is a hypothesis; no observed failure cause or superiority is established by this code change.

## Entry points and retained results

The added configuration is `discrepancy_each_scale=False` by default. The existing primary strict TRF, final-only `--discrepancy`, and weakreg1e-5 branches retain their existing solver behavior, budgets, suffixes, and completed files. An old receipt missing the new key is treated asFalse. No historical source or result is rewritten.

The new explicit CLI is `experiments/run_anchors.py --discrepancy-all-scales`, or `experiments/run_baselines.py --kind discrepancy-all-scales`. Outputs use the separate suffix `_discrepancy_allscales`. The baseline runner's default remains the old two kinds, final-only discrepancy and weakreg1e-5. The all-scales flag selects an independent unregularized data-TRF control and cannot be combined with `--full-trf` or a nonzero `--weakreg`.

Example registration command, documented only and not executed by this implementation task:

```
python experiments/run_baselines.py --kind discrepancy-all-scales --split val --split test --workers 4 --core-mask 20560
```

## Stopping and endpoint semantics

The scales remain8→12→24. At every scale the physical state is a full576-cell nonlinear VIE, and the Jacobian is the existing counted analytic Full-J multiplied by the same material interpolation basis. The material box remains[.02,20], initial material.1, strict TRF tolerances1e-12, and the original cumulative phase fractions.20/.25/.55 remain unchanged. Each independent observation keeps45 seconds,500 LU and30000 complex RHS; features and shared reference setup retain their existing accounting.

At each accepted Jacobian point, after the actual forward and analytic Jacobian have completed, compare the relative data residual with1e-6 for clean data or.0011 for noisy data. A reached gate at8 or12 produces a `STAGE_DISCREPANCY_STOP` event, returns that accepted coefficient vector, and prolongs to the next scale. It does not set a global stop reason. At24 all576 material variables are released, its starting physical forward and analytic Jacobian are performed before any gate is considered, and a reached gate produces the final `DISCREPANCY_STOP`.

The new control records endpoints only after their analytic Jacobian succeeds. If TRF terminates without requesting its final accepted endpoint Jacobian, that endpoint is confirmed through the same counted forward/J API before use. Rejected residual trials remain trace entries and are never used as the endpoint. The extra endpoint confirmation belongs to the same scientific budget; there is no new solve budget or restart.

Budget and numerical stops remain valid negative outcomes. They are recorded honestly and do not bypass the requirement for a successful24-scale forward/J. In particular, the all-scales event field `accepted_forward_and_analytic_j_completed` and existing `all_unknowns_released` show whether that requirement was reached. Numerical termination without a reached discrepancy is not labeled a discrepancy success. Only the24 gate can set the global discrepancy stop; global budget/numerical stops retain their original semantics.

## Verification and authorization boundary

Only Python syntax compilation of the three edited Python files was performed. No numerical task, FD check, random validation, data regeneration, remote upload, queue interruption, or scoring was run. The parent owns registration under the method-freeze protocol and the later independent fullVAL/TEST execution. All previously incurred primary/final-only/weakreg costs remain chargeable.
