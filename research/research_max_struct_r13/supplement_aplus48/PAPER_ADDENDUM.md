# Material Readout, Initialization and Observable Feasibility in Maxwell Inversion

## Abstract

We complete a previously partial comparison of one-hop material readout within a mature CIE/FFT-TSOM implementation. On the same48 exposed two-dimensional materials, unchanged R13 A+-MM achieves mean normalized material error0.15210 versus0.19454 for CIE and0.17985 for Dual-CIE. A controlled transplantation of the final R7 arc/guard readout preserves the mean error at0.15214 while reducing mean full-process CPU from119.37 to113.22 seconds. At common per-object paid CPU prefixes, its error is0.15282, compared with0.19468 and0.17978 for CIE and Dual-CIE. A fixed16-object, six-solver initialization study finds almost no final A+ benefit from a frozen neural warm start, whereas G19 has a concentrated rescue benefit and DR has a smaller mean benefit. Separately charged native-state diagnostics reveal large differences between auxiliary fitting and physical material feasibility, and material-readout-induced changes of the subsequent current gradient. These findings support a conditional recovery and efficiency contribution, rather than universal dominance, global convexification, or a new continuous-Maxwell identifiability theorem.

## 1. Scope and relation to the original R13 study

The original R13 mechanism arm included A+ on16 predetermined objects, while CIE, Dual-CIE, G19 and DR covered48. We complete the missing32 with the frozen runner, parameters, observations, prior and candidate-selection rule. These objects had already been exposed; the completion and subsequent controls are follow-up paired experiments, not a new blind cohort.

Eight material families cover separated objects, open and closed pores, graded random fields, eccentric layers, patchwork and thin bridges. The inverse material has576 independent real variables. Thirty-six cases use matched grids and12 use cross-grid data. Acquisition and noise variations are bundled rather than a fully crossed design. Historical three-dimensional and measured-data evidence retains its original limitations and is not enlarged by this supplement.

## 2. Methods and experimental controls

All methods retain the same legal box, original smooth/TV prior, public observations and original-Maxwell objective. Each run has a120-second whole-process CPU allowance, with at most72 seconds for proposals and48 for the common physical corrector. Original native iterations and every paid physical audit are charged. The primary estimate minimizes the original paid, legal Phi; no ground-truth metric selects candidates or stopping points.

A+-MM retains the global CIE current update but changes material readout through one-hop constitutive-residual propagation and a conditional Möbius/MM solve. The original R13 default omitted two final R7 flags. A separate controlled arm restores exactly the five frozen R7 readout options, including arc stationary candidates and inherited material-gradient stopping. It is a readout transplantation within the R13 contract, not a reproduction of every historical R7 condition.

The frozen24,737-parameter R8 initializer is reused on the same16 materials for CIE, Dual-CIE, G19, DR and both A+ variants. Initialization maps differ: cold auxiliary current, projected physical current, and complete physical field–current state are distinct interfaces. The experiment therefore measures complete deployed pipelines, not a pure initialization-map causal effect. Existing supervised training information remains an additional prior.

## 3. Completed reconstruction results

|Method|Mean material NRMSE|Mean IoU|Mean actual process CPU, s|
|---|---:|---:|---:|
|CIE|.19454|.69456|106.15|
|Dual-CIE|.17985|.71592|106.26|
|G19|.16758|.68577|119.08|
|DR|.17672|.68726|84.50|
|A+-MM/R13|.15210|.75063|119.37|
|A+-Arc/R7-readout|.15214|.75153|113.22|

For each object, common-CPU evaluation uses only candidates whose complete physical audits had finished before the minimum actual process consumption across the compared methods. It does not interpolate checkpoints or use future candidates. The resulting mean errors are.19468/.17978/.18222/.17672/.15746/.15282 in the table order. This is a paid-prefix comparison, not a newly measured interrupted deployment time.

Benefits concentrate in difficult cases. On0017,0026 and0045, A+ recovers main material structure missed by CIE and Dual-CIE. DR instead succeeds on0002,0038 and0042 where A+ remains incorrect. A+ and ordinary CIE preserve0005 and0032 while DR fails. All failures and all48 images remain available. Existing global-Dual16 still fails0045, providing a material-readout control on that object; the full48 effect has not been isolated from the global-versus-sourcewise current distinction.

The arc/guard and MM final errors differ by only0.000041 on average. Their version difference is real but does not explain the final ranking in this cohort. The useful arc contribution is lower cost and better early prefixes.

## 4. Neural initialization and the common corrector

On the same16 objects, uniform-to-neural final errors are.14920→.14980 for CIE,.12786→.12818 for Dual,.15008→.12656 for G19,.12313→.11936 for DR,.10966→.10943 for MM and.10970→.10996 for arc. G19's main rescue occurs on0024; neural DR still fails that object. Neural initialization improves the initial material error from.39078 to.27702, but does not automatically improve the final auxiliary-optimization outcome. Arc has an early neural benefit at15 CPU seconds that largely disappears at the endpoint.

The corrector improves mean error for Dual, G19, DR and A+, with damaging cases in every group. CIE improves33 of48 individual errors while its mean worsens slightly. Physical-objective descent, material improvement and noise-ball feasibility must therefore be reported separately. Selecting whichever phase looks closer to ground truth would define an invalid estimator and is not done here.

The historical same16 direct reference improves from.40966 to.27223 with neural initialization, but actual CPU rises from55.92 to101.39 seconds. This is not evidence of equal-cost neural acceleration. Training plus training/validation data generation previously cost402.558889 Mac CPU seconds; the reused R13 feature/inference preparation cost1.285002 Mac CPU seconds. These are separate from i9 solver times.

## 5. Mechanism evidence and mathematical limits

Let u=J−Dm Q, v=Q−Einc−GJ and Am=I−Dm G. For invertible Am, the physical-current discrepancy is Am⁻¹(u+Dm v), and its receiver-visible component is SAm⁻¹(u+Dm v). This known finite-model relation organizes the observed failure: an A+ state on0002 fits auxiliary data within0.891 noise radii while the same material has physical residual49.037 radii. A small auxiliary fit or defect energy alone does not establish physical feasibility.

At fixed postcurrent fields, material readout changes the next current coefficient gradient. Its data term cancels in the before/after difference, leaving the material-weighted constitutive residual and Green adjoint. Large rotations occur on both successful and failed paths; rotation magnitude is not a recovery predictor. Among475 retained material readouts,139 increase the physical objective, including111 with decreasing counterfactual one-hop loss. Retained checkpoints are not independent samples or all iterations.

These observations motivate a falsifiable next question: when does readout-induced current-gradient change reduce the port-visible physical defect? They do not establish a new global-basin theorem. A smooth invertible Möbius coordinate change alone preserves strict interior local minima of the same physical objective. Reflected memory also does not resolve an already feasible wrong material without selection information, prior preference or additional measurements.

## 6. Reproducibility, cost and limitations

The supplement adds184 runs,20657.390625 process-CPU seconds,156439 physical LU factorizations and3040452 complex solve right-hand sides. Native diagnostic replays, offline physical gradients, fixed-background projection work, failed operational attempts and archival work are separately charged. Three-dimensional/experimental-data evidence, true runtime bundles, legal-current initialization, observer overhead and exposed-cohort status are explicitly documented.

The results justify further development of a material-readout algorithm/mechanism paper with complete paired failure evidence. They do not yet justify a broad strong-TAP claim across measured and three-dimensional physics, an original global projection theory claim, or continuous-Maxwell stability/uniqueness claims for Inverse Problems. Prior CSI/CIE, WRI, variable-projection and DR/RRR mechanisms remain prior art, as documented in the parent novelty ledger. No further algorithm loop is started in this supplement.
