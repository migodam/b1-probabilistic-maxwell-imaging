# Paper claim disposition

**HOLD for an improved GlobalDual imaging method.** This follow-up does not establish an original Maxwell theorem or an overall recovery/cost advantage.

Working descriptive title: **Contrast-Stratified Basin Accessibility in Nonlinear Maxwell Imaging: Strong Jacobian Baselines and Field-Gram Proposals**.

The useful result is an auditable contrast/noise-paired comparison: twelve full 24×24 material shapes, four frozen solvers, 288 i9 P-core runs and 72 paired RTX4060 runs. Mature native TRF reconstructs many intermediate-contrast targets quickly. GlobalDual preserves a few IoU/threshold gains, but has no practical NRMSE wins against the three comparators in the CPU collection under the predeclared tolerance. Larger intermediate-contrast errors are traceable to the shared data-L-BFGS coarse path; this is not an isolated causal indictment of the final Gram stage.

The fixed-background resolvent expansion, source-basis invariance, weighted root Hessian and linear unbiased covariance argument are standard algebra/least-squares consequences. They delimit the claim; they are not new inverse-scattering theorems. The GlobalDual construction retains the CSI/SOM/variable-projection precedent identified in frozen R3/R4 reviews. No new literature search or priority claim is made here.

The shortest missing experiment is a newly frozen matched TRF coarse initialization followed by data, fixed-reference Gram and dynamic Gram arms. This would separate coarse-path accessibility from source weighting and the material derivative of the weight. It was **NOT_RUN**, and the present disclosed objects cannot serve as its fresh holdout.

Not allowed: general GlobalDual superiority, high-contrast robustness, all-root enumeration, true alternative solution fibers, posterior sampling, end-to-end GPU speedup, continuous-PDE/3D vector Maxwell generalization, or a causal claim that only the Gram stage produced the intermediate-contrast failures.
