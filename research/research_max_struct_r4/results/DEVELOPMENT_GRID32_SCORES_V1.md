# Development addendum: grid32

DISCLOSED_DEVELOPMENT_GRID_DIAGNOSTIC_NOT_BLIND

Separate versioned addendum; the previous13-row table remains unchanged. Original frozen RMSE/12, Otsu256, public99% CN radius, prior and material box. No new physics or GT selection; no confirmation truth access.

| Method | RMSE/12 | IoU | Relative residual | White residual | Joint | LU | Physical RHS | Wall s |
|---|---:|---:|---:|---:|---|---:|---:|---:|
|onehop|0.053067749|0.935483871|0.000970516939|23.292542|PASS|1026|24540|128.757037|
|global|0.053102210|0.935483871|0.000970503354|23.292216|PASS|1027|24564|92.725427|

All displayed physical costs are original run top-level totals; do not add nested stage costs again. Solver wall/caller wall scopes remain explicit in JSON. Grid32 shared operator generation is a separate one-time receipt, not included in each solve. Corrector stop states are retained in JSON; no convergence claim. Scientific judgment stays with parent Codex.
