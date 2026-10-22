# Conventional surrogate-force HMC diagnostic

Only the C-owned proxy module, runner, and these controls are authored here.
Frozen `src/p3l` and locked estimator code are unchanged. No paid chain is
launched by authoring or by the pure controls.

The existing HMC kernel receives the unchanged Full probit target and a
deterministic frozen surrogate force. Each kick and drift preserves volume;
the symmetric sequence followed by a momentum flip is an involution. This
remains true for an inaccurate force. Acceptance uses the exact Full endpoint
Hamiltonian, with the complete original Gaussian66 prior and the same fixed
mass kinetic energy. Internal leapfrog steps never request Full gradients or
physics. This is conventional surrogate-force HMC, not a novelty claim or
posterior mixing certificate. Numerical failures remain explicit diagnostics.

The linear arm freezes the true white probit J and uses
`-Jx.T @ (white(F0-y)+Jx@(x-x0)) - x`.
The material arm freezes the paid receiver/source response
`Cchi=einsum('qp,ps->qsp',T,field).reshape(256,576)` and uses
`Fhat=F0+(Cchi@(chi(q)-chi0)).reshape(32,8)`. Its force is
`-2/sigma^2 Re(Dchi(q).conj().T@Cchi.conj().T@vec(Fhat-y))`
times the probit derivative, then minus x. Complex rows remain receiver-major,
source-minor, and whitening remains the original all-real/all-imag source pack.
The full decoder law is retained; no arbitrary current coordinates are sampled.

`tools/run_p3l_surrogate_hmc.py --arm {linear,material} --chain {0,1,2,3}`
runs one serial Book-C job with an immutable protocol and source snapshot.
All arms use the original mala-mode chain's first physical-probit state,
freshly rebuilding its Full state. Warmup1024/production4096 use lengths16/32/64,
seed2026102341+chain, initial step .05, and dual averaging target .8. At
256/512/768 the true current paid Full endpoint supplies J/T and GN inverse
mass; the surrogate, cached current force gradient, and adaptation are reset.
An exactly unchanged paid point may reuse its already paid J/T. Production
force, mass, and step size are frozen.

Every Full endpoint is saved by Target, including rejected proposals. Proposal
proxy-response/Phi discrepancies are computed only after the MH decision using
that paid endpoint. All rejected-state repeats are retained. Chain arrays,
force invocation/time counts, costs, current cache, and RNG checkpoints are
saved every128 transitions; per-transition metadata is also retained. Both
force CPU time and total controlled wall are disclosed. No Full gradient is
ever used for production dynamics or diagnostics.

The conservative bound is 164,384 RHS per four-chain arm, 328,768 for both:
each chain has at most5121 fresh Full forwards at8 RHS and four receiver-J
setups at32 RHS. Every setup, failed/rejected forward, and actual cache reuse
is owned by the shared Target/Book accounting, with a7200-second per-job wall
cap and the global one-paid-worker lock. No previous budget is transferred.

The pure controls deliberately use a wrong nonconservative fixed Gaussian
force and independently check endpoint Hamiltonian MH, detailed-balance joint
flux, involution, and the finite-difference phase-space volume determinant.
They then read previously paid DEV030 transfer/field arrays and online y,
with no Target/Book/backend, to check both physical proxy gradients against
finite differences, complex receiver/source packing against a direct einsum,
the white-J adjoint gradient identity, and the unchanged Gaussian66 prior.
An explicitly inward coordinate used in controls does not claim a new Full
center; its forward response is only arithmetic under the supplied frozen
matrix. Passing these controls does not establish actual sampler efficiency,
material accuracy, posterior convergence, or a qualified reference.
