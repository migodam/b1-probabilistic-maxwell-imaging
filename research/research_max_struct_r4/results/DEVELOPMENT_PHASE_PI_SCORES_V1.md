# Development addendum: phase_pi

POST_SEAL_PUBLIC_DEVELOPMENT_NOT_BLIND

Separate versioned addendum; the previous13-row table remains unchanged. Original frozen RMSE/12, Otsu256, public99% CN radius, prior and material box. No new physics or GT selection; no confirmation truth access.

| Method | RMSE/12 | IoU | Relative residual | White residual | Joint | LU | Physical RHS | Wall s |
|---|---:|---:|---:|---:|---|---:|---:|---:|
|onehop_phase_pi|0.558703269|0.257812500|0.097692953412|2344.636820|FAIL|1019|24372|74.455652|

All displayed physical costs are original run top-level totals; do not add nested stage costs again. Solver wall/caller wall scopes remain explicit in JSON. Grid32 shared operator generation is a separate one-time receipt, not included in each solve. Corrector stop states are retained in JSON; no convergence claim. Scientific judgment stays with parent Codex.
