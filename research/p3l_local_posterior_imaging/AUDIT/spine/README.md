# Frozen nonlinear property-spine diagnostic

This is an optional Book-B algorithm diagnostic. The implementation has not
been run against Maxwell and supplies no material-accuracy, ESS, convergence,
normalizer, conditional-MAP, or reference certificate.

For physical probit coordinates `x=(z,eta)`, freeze a global quadratic
`beta: R2 -> R64`. Use `r=(v,eta)` and `T(r)=(v+beta(eta),eta)`.
The inverse subtracts the same beta. The derivative is block triangular
`[I, Dbeta; 0, I]`, so its determinant is one everywhere. All 64 residual
coordinates remain unrestricted. The transformed density is
`-Phi(from_probit(T(r))) - ||T(r)||^2/2`, and its gradient is
`DT(r).T @ (-gradient_Phi_x - x)`. A standard Gaussian prior on `r` would be
incorrect. The shear does not assert that the posterior lies on a manifold.

`tools/run_p3l_spine.py --build` is a separately charged Book-B job named
`build_property_spine_DEV030`. Before any target query it freezes protocol,
target configuration, and source snapshots. It reads the existing c0 y-only
probit mode endpoint, rebuilds its true state and full J, and obtains a linear
field compensation from the full probit GN blocks. Each of nine fixed eta
nodes uses at most 24 fresh field states and analytic full receiver-J actions.
Target writes every attempted and completed endpoint. Each node saves its best
cache, objective, field gradient, solver status, and endpoint history. A stopped
solver's endpoint remains approximate. Failed jobs preserve partial artifacts
and are never automatically overwritten/replayed. The full quadratic fit saves
all node arrays, coefficients, rank, residuals, and errors. Native transfer
caches are cleared between nodes; optimization retains only current and best
physical states.

`--chain c` uses Book B job `spine_hmc_DEV030_pilot_c{c}`. It starts from the
same paid y-only first source-chain row as ordinary HMC, with a freshly billed
full state. If floating-point inverse/forward cancellation changes represented
x, it bills that represented state separately and saves the round-trip delta.
No approximate cache alias is used. Both the source and represented initial
coordinates are recorded. Warmup is 512 and production 1024. GN mass refreshes
occur at 128/256/384 using true full J and `Hr=DT.T@(Jx.T@Jx+I)@DT`.
Step size starts/resets at .05, dual averaging targets .8, trajectory lengths
are drawn from 4/8/12, and the same ordinary-HMC seed `[2026102290,c]` is used.
All production parameters are fixed. Endpoint HMC uses the unchanged existing
full-dimensional kernel; every rejection repeat is retained.

`CHAIN.q` is native `(z,u)`; `CHAIN.x` is physical probit `(z,eta)`;
`CHAIN.r` is the separate transport coordinate. Saved `log_density` is the
full physical-probit target, equal to the r target because `dx=dr`. This
preserves the existing reference auditor's coordinate contract.

`check_spine.py` uses synthetic analytic functions only. It checks the
polynomial formula, global inverse, finite-difference transform derivative,
unit determinant, original-prior gradient, full target composition gradient,
leapfrog reversibility, and the existing HMC endpoint/energy rule against an
independent integrator. The intentionally wrong Gaussian-r prior differs on
the control point. These controls do not execute Book, Target, Maxwell, truth,
training, or network work.
