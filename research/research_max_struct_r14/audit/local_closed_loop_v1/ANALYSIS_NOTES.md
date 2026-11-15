# Offline root mechanism audit

One declared8x8 finite-TM anchor; no ground-truth/test access and no new recovery claims. ROOT_FULL_STATE saves the Fourier/DPCS physical-root embedding. Every actual trial saves complete input/output state and PR/readout diagnostics, including zero-perturbation seed responses. IDEAL arrays use an independent orthonormal current space rather than singular redundant Fourier coordinates.

TABLES directional rows compare central responses with exact-current/exact-material ideal maps and split actual-minus-ideal into current-block and material-block error responses. The closure is an algebraic accounting check using defined residuals, not a certificate that finite epsilon is asymptotic or that the program has a Frechet derivative. ADDITIVITY arrays test the unnormalized h1+h2 response. Interpret plateau/nonadditivity together with readout guard, root selection, finite sweeps, root seed gaps and box events in TRIAL_LOG.

Counts separate one physical forward/setup, ideal normal factors/RHS, actual clone increments and actual-step increments. Dense diagnostic matrix products are offline work; inherited counters retain action/readout detail. Total CPU/wall include this program from main (plus process_CPU including imports), while individual trial times include state serialization. No search or hyperparameter tuning is performed. Parent retains final interpretation.

## Execution observations (no recovery interpretation)

The physical root embedding has 48 independent current dimensions among64native Fourier coefficients, relativeJgap7.24e-15. Zero-perturbation material steps move6.10e-14–8.16e-14. Ideal normal residuals are≤6.24e-16. The commonB00norm may exceed1even for11because its ownmaterialGram differs; do not call that a failedcommonpotential result.

At epsilon.001 current additivity relative errors are.01240257for00/01and.15447828for10/11, remaining close across the3declared epsilons. Arc readouts take0sweeps at the exactroot and1or2sweeps on perturbed trials; this guard/sweep dependence is retained. No single programJacobian is inferred from centraldifferences.

One failed initial attempt paidone64LU/sixRHS then hit the inherited arm_config explicitcontinuation validation(mf12only). Its directory/log remain siblingfailed_attempt01; the final source obtains inheriteddefaults then installs the declaredsinglemf4stage. Summary includes failedcosts.72diagnostic dense materialerrorLU solves were retrospectively instrumented exactly as36responses×2calls; no extra run was needed.
