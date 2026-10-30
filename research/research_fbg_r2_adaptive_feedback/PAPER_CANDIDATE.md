# Paper candidate — evidence-qualified

**Title:** *Residual-Calibrated Feedback Metrics for High-Contrast Maxwell Basin Search*

**Scientific problem.** A data-generating material can exist while a mature local inverse solver cannot reach a useful reconstruction from public initial states. We ask whether changing a physically exact auxiliary residual geometry, without changing measurements or the legal material space, improves finite-budget accessibility.

## At most three candidate contributions

1. **Uniform-reference feedback gradient decomposition.** Under the stated fixed-source, full-rank, real-material and reference-normalization assumptions, the gradient shift is a nonnegative, residual-conditioned contrast drive. Its entries are column energies of the existing observable compatibility residual, are source-basis invariant, and activate quadratically with a small data perturbation. F09/F12 contain the proof. This is a specific algebraic mechanism, not a new Maxwell law; closest-prior novelty is not certified.
2. **A complete scalar-TM feedback-Gram imaging implementation and frozen comparison.** At 576 unknown real materials, fresh high-contrast clean success is 8/12 for fixed beta=.2, 6/12 for automatic scalar calibration, 3/12 for original GlobalDual and 0/12 for each of three ordinary controls. Fixed beta has six paired gains and one loss versus GlobalDual. Eight public starts and real 3/6/12/24-source acquisitions test conditional accessibility separately. All costs and failures are retained; the result is not a universal speed or recovery theorem.
3. **A limitation of tangent-energy feedback design.** The current-fit objective is exactly transformed-field tangent attenuation, yet the useful scalar feedback can increase both field tangent energy and projector movement. Diagonal/block regression, short physical pilots and a Schur-Galerkin approximation fail to outperform the scalar. The B=-G and scaling counterexamples explain why attenuation/rotation alone is not a design certificate. This is a mechanism/negative-result component, not an invention of ridge regression or Schur elimination.

## Claims not allowed

- The scalar fractional gauge is already CIE-I-related; standard Schur, ridge, variable projection and Grassmann derivatives are prior tools.
- A17 is not shown necessary for the positive images. The main successful arms use an exact reference adjoint; the Schur probe is one conditional Galerkin bank, not the full A17 exchange algorithm.
- Q-beta does not outperform fixed beta, and no feedback family is proved optimal.
- A prototype joint-current result is not a comparison against mature CIE-I/CSI/TSOM.
- Eight-initial-state frequencies are not a continuous basin volume. Saved path nodes are not independent exact roots or posterior samples. Fresh K1 best@K equals online top1.
- Same-grid, single-frequency scalar-TM simulations do not establish 3-D vector, continuous-PDE, cross-grid or experimental robustness.
- Nominal wall caps have measured finalization overruns. Physical counters and actual wall are separate; there is no general whole-solver speedup claim.

## Most important missing work

Compare against a mature CIE-I/CSI/TSOM implementation and strong fractional-material parameterizations on a new, prospectively frozen dataset with matched actual computational resources. This is the shortest test of whether the feedback-Gram mechanism adds something beyond existing CIE reconstruction.

## Disposition

**GO for further method development; HOLD for a TAP-ready superiority or novelty claim.** The frozen images and accessibility counts support a real research advance over this project's tested GlobalDual path. The mean paired NRMSE cluster interval crosses zero, the sample is small, failures remain, and the closest algorithmic prior is not yet fully controlled.
