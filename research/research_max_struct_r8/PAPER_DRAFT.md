# Transport and Selection of Learned Warm Starts in Nonlinear Electromagnetic Inversion

**Research manuscript, not an accepted-paper or universal-superiority claim.** All 750 primary/ablation predictions and 126 affine-offset predictions on 42 independent materials were completed and sealed before offline test scoring. Classical, measured and vector-three-dimensional experiments are reported with their distinct contracts and limitations.

## Abstract

Learned initializations can supply useful material information for nonlinear electromagnetic inverse scattering, but that information must survive the constraints and candidate selection of the subsequent solver. We study this interaction in contrast-inversion equations (CIE), FFT-TSOM, and a nonlocal material readout with analytic Möbius updates, denoted A+. A 24,737-parameter convolutional predictor is followed by full-wave physical optimization. Across 24 in-distribution, 12 held-out-family and six cross-grid materials, neural direct optimization attains mean amplitude-normalized material errors of 0.0965, 0.1624 and 0.0843; neural global A+ attains 0.2961, 0.5939 and 0.3293. A projection decomposition separates hard deterministic-current mismatch from Fourier truncation. Retaining the physical ambiguity-current component as a fixed affine offset repairs some trajectories but also produces counterexamples. Selecting the lower physical objective after paying for two complete branches improves the in-distribution cohort yet worsens the other two, exposing a further limitation of residual-based material selection. Neural direct optimization has lower mean deployment CPU and error than the tested auxiliary pipelines in each hardware/cohort stratum, but not universal equal-budget superiority; it costs more than a uniform direct optimizer that often stops at poor materials. Classical Austria, measured Fresnel and vector three-dimensional experiments do not establish a stable A+ recovery advantage. We therefore position the work as a reproducible study of learned-state transport, selection and model adequacy, rather than a generic invention of neural warm starts or a universally superior Maxwell solver.

## 1. Introduction

High-contrast inverse scattering combines nonlinear propagation with incomplete and imperfect observations. Auxiliary-current formulations can make alternating updates tractable, but their geometry is not identical to that of the original material-to-data map. Replacing one material update, changing current coordinates, or supplying a better-looking initialization can therefore change a reconstruction trajectory without adding physical information or guaranteeing material accuracy.

The present work develops an evidence-led investigation of this issue. Earlier experiments found that A-material and its cheaper analytic A+ readout recovered two exposed difficult curved targets that conventional CIE variants did not recover. The same advantage did not persist uniformly across other materials. That evidence motivates a stronger question than whether the new auxiliary loss decreases: which information supplied by a material initializer is retained by the current representation, and which candidate reaches the physical corrector?

We investigate three distinct mechanisms. First, a learned material is mapped to a physical current and then subjected to deterministic-current and Fourier restrictions. Second, sourcewise and global Polak–Ribiere updates interact with the same ordinary or A+ material readout. Third, selecting a candidate before nonlinear physical correction can yield a different result from retaining and correcting multiple candidates. These mechanisms are compared with a direct physical optimizer; an auxiliary algorithm receives no credit when its final result simply comes from the unchanged learned initial material.

The evidence is separated into exposed development targets, new independent synthetic materials, original-paper acquisition contracts, measured data, and genuinely vector three-dimensional problems. This separation is necessary because a geometry-inspired toy problem is not a reproduction of a published experiment, and a small dyadic model is not a demonstration at the original electrical size.

## 2. Related Work and Scope of Novelty

The nonlinear integral equation (NIE) formulation of Zhong et al. [1] and the contraction integral equation CIE-I of Zhong and Xu [2] motivate the inherited physical and auxiliary formulations. Fourier/TSOM and multi-resolution CIE are established methods. The present work retains their original coefficient coordinates in its qualified comparisons; exact current elimination and SVD-whitened coordinates are separately labelled exploratory controls.

Learning an initialization for a physical inverse solver is also established. CS-Net [3] learns radiation/noise current components and feeds SOM/TSOM/CSI; it is a particularly close precedent for learning-assisted current initialization. Learned warm starts for inverse obstacle scattering [4], initialization followed by distorted Born iteration [5], NeuralBIM [6], and CIE-SOM-NET [7] further preclude a generic claim of inventing neural initialization or learned CIE. A change from learned current to a material network followed by a physical forward solve does not by itself establish a new physical principle.

Möbius transformations, conditional minimization, majorization-minimization, orthogonal projection, Maxwell resolvents, and effective-medium polarizability are likewise not new here. Potential contributions are narrower: an explicit diagnostic of how a learned material state is changed by the inherited current constraints; a simple affine repair and its counterexamples; sourcewise/global material-readout interactions; and reproducible recovery-versus-cost evidence across physically qualified contracts. Final novelty and publication claims must be limited to what those comparisons establish.

## 3. Forward Model and Experimental Objectives

For the two-dimensional scalar-TM finite model, let X=diag(m), with real contrast m=εr−1 and fixed μ=1 under the exp(−iωt) convention. The full-wave current and data are

\[
J=X(E^{inc}+GJ),\qquad F(m)=S(I-XG)^{-1}XE^{inc}.
\]

The primary synthetic comparison uses 576 independent real material parameters, common bounds [0.02, 20], and

\[
\Phi(m)=\tfrac12\|F(m)-Y\|_F^2/\|Y\|_F^2+R(m),
\]

where the same smooth/TV-type explicit prior R is used within each comparison. This normalized objective is not the Gaussian likelihood whitened by the declared noise standard deviation. Noise consistency is evaluated separately. Sources and receivers lie outside the unknown region.

Complex noise is CN(0,σ²): the standard deviation of each real or imaginary component is σ/√2. Synthetic material error is RMSE divided by a predeclared contrast amplitude, stored in private scoring metadata. Predicted morphological support uses its own Otsu threshold. Because the inherited ROI contains a nonzero host, this morphology is not the complete contrast support relative to vacuum.

The original-paper Austria contract, measured-data contract, and vector-DDA contract differ and are not pooled. Austria includes an air initialization and no additional prior in the paper-faithful setting. Measured data use source-specific observation operators and masks. Vector inversions retain three-component currents and dyadic interactions rather than stacking scalar-TM slices.

## 4. Methods

### 4.1 Lightweight material initializer

The primary network has 24,737 parameters and two convolutional scales of width 16. Its ten input channels contain complex Born diagonal least-squares and deterministic-current material features, current/field energy summaries, aperture, source count, and spatial coordinates. Feature calculation uses public operators and observations only. It requires receiver factorization and Green/receiver actions, but no nonlinear Maxwell LU solve.

Training uses 1,536 independent materials from eight families; 192 separate validation materials select the epoch. Contrasts, noise, aperture and source count vary. A same-architecture Born-only control removes four deterministic-current/field channels. Both networks are trained for 120 epochs; the physics-feature network selects epoch 86 and uses 387.69 CPU seconds. Training plus its own training/validation data generation costs 402.56 CPU seconds, before other research development costs. A separate acquisition-domain repair network is trained for the classical paper experiments and never substitutes for the primary confirmation weights.

The network only supplies m0. It does not replace the full-wave solver, reduce the 576-dimensional material class, or provide test labels to the solver. Its training distribution is additional prior information, even when the explicit R is unchanged.

### 4.2 Current constraints and affine transport

Let the retained receiver singular vectors define Pd=VdVd*, and let Jd be the deterministic current obtained from the retained data singular components. The inherited truncated family is

\[
J=J_d+UA,\qquad U=(I-P_d)\mathcal F_k.
\]

Material-only initialization leaves A0=0. A joint initialization projects the paid physical current jp=Jphys(m0) into this same family:

\[
A_0=U^\dagger(j_p-J_d).
\]

The proposed diagnostic repair instead fixes

\[
J_a=J_d+(I-P_d)j_p,\qquad J=J_a+UA,\quad A_0=0.
\]

The offset is never recomputed during material updates. It changes the early truncated affine family while preserving its optimized dimension. At the final full Fourier space, the feasible current family is the same full APCS family. The initial physical current is reused from an already charged objective evaluation; one complement projection and full fixed-current storage remain chargeable costs. This is an initialization/representation modification, not a new observation or a claim of removing Maxwell nonlinearity.

### 4.3 A+ and sourcewise/global current updates

A+ retains ordinary CIE current updates and changes the fixed-current material readout. In the inherited finite model, a real material coordinate is mapped to z_i=m_i/[1+(β−G_ii)m_i]. The zero diagonal of the off-self propagation makes the conditional residual affine in z_i. Analytic real-arc stationary candidates, endpoints and the old value are checked against the inherited prior majorizer and original fixed-current material loss.

This acceptance concerns the auxiliary material block; it does not guarantee descent of Φ after a full physical solve or descent of material error. Sourcewise A+ uses a separate current PR coefficient for each illumination; global A+ uses one coefficient for the joint current block. The ordinary/A+ by sourcewise/global factorial retains the same weights, material class, Fourier schedule and physical corrector.

For three-dimensional isotropic real materials, a scalar self block sI3 permits the same real coordinate structure after stacking three field components and all illuminations. General complex material occupies a two-real-dimensional feasible region and is not treated by the real circular-arc routine. The complex-material experiments use a legal constrained alternative.

### 4.4 Selection and direct physical controls

Each stage endpoint and the initialization are evaluated by the original full-wave Φ. The inherited pipeline corrects its smallest-Φ candidate. We also retain a direct physical trajectory from the neural initializer and compare a paid two-branch portfolio selected by final Φ. Both complete branches are charged; this is a conventional retained-candidate control, not a free extra corrector or a new global convergence theorem.

The primary auxiliary schedule has at most 3,000 iterations and a common physical corrector of at most 3,000 accepted steps. Direct physical optimization has a 6,000-step maximum. We therefore report both final outcomes and comparisons at shared actual CPU/LU/RHS allowances from paid incumbents. Sparse checkpoints are not interpolated into unexecuted trajectories.

## 5. Mechanism Analysis

**Projection decomposition.** With Pu=UU†, the change caused by the original joint initialization is

\[
\delta j=\underbrace{J_d-P_dj_p}_{\delta j_D}
+\underbrace{P_uj_p-(I-P_d)j_p}_{\delta j_N}.
\]

The two components are Euclidean-orthogonal. The first equals VdΣd⁻¹Ud*(Y−Sjp), so hard deterministic-current mismatch remains even with a full ambiguity subspace. The second is Fourier representation loss. The constitutive defect is the exact finite quantity (I−X0G)δj, not merely the current L2 error.

**Affine projection property.** Ja is the closest current to jp in the full affine set PdJ=Jd. Every feasible displacement has fixed Pd component δjD; its orthogonal complement is minimized by zero. The fixed-offset initialization therefore removes δjN at the starting state without adding optimized degrees of freedom.

This property does not imply smaller constitutive defect. For a general invertible B, BδjD and BδjN need not remain orthogonal and may cancel. For example, B has diagonal 1 and off-diagonal −3/4, while δjD=(1,0) and δjN=(0,1). Removing δjN changes the transformed norm from sqrt(2)/4 to 5/4. This is an algebraic counterexample to a norm inference, not an asserted Maxwell material counterexample. Actual physical effects require the executed full-wave experiments.

**Candidate-ordering limitation.** A nonlinear descent map C need not preserve initial objective ordering. It is possible that Φ(a)<Φ(b) but Φ(C(a))>Φ(C(b)). Consequently, C(argmin Φ) and argmin Φ(C) need not agree. The latter requires paying for multiple corrections. No material truth is needed to define either choice, and neither choice certifies material uniqueness.

**Dimensional limitation in the original Austria class.** The 64² inverse material has 4,096 real parameters, while 20 illuminations × 40 complex data provide at most 1,600 real observations. Its real Jacobian nullity is at least 2,496. Under local constant-rank and interior well-posedness assumptions, the constant-rank theorem yields a local equal-data manifold of dimension 4,096−r. This standard argument rules out a uniqueness claim over the full open material class; it does not prove a particular constrained Austria truth has a same-box material fiber. Air-boundary activity and forward/inverse discretization differences must be retained. The argument also cannot be transferred unchanged to the other data contracts.

## 6. Experiments

### 6.1 Independent synthetic confirmation

The completed frozen study contains 42 new independent materials: 24 within the training families, 12 from four held-out families, and six generated on 36² and inverted on 24². The original 750 tasks contain 630 method/factorial comparisons and 120 initialization ablations on eight predefined objects. A subsequent 126-task fixed-offset extension was frozen from development evidence before test labels were read. All 876 predictions were completed and both sets were sealed before offline scoring. They are 42 independent material objects, not 876 statistical samples.

| Method | ID, n=24 | Held-out families, n=12 | Cross-grid, n=6 |
|---|---:|---:|---:|
| Uniform direct physics | .24966 (2/24) | .40474 (1/12) | .19593 (2/6) |
| Neural direct physics | .09647 (9/24) | .16245 (6/12) | .08431 (2/6) |
| Uniform Dual-CIE | .36843 (7/24) | .63825 (1/12) | .32821 (2/6) |
| Neural Dual-CIE | .32822 (7/24) | .56779 (2/12) | .32835 (2/6) |
| Neural Global-Dual | .28919 (8/24) | .60481 (1/12) | .32834 (2/6) |
| Uniform global A+ | .34062 (8/24) | .63748 (1/12) | .32833 (2/6) |
| Neural global A+ | .29609 (7/24) | .59393 (2/12) | .32933 (2/6) |
| Neural sourcewise A+ | .32182 (7/24) | .59164 (2/12) | .32811 (2/6) |
| Neural affine-offset Dual | .29007 (7/24) | .52299 (4/12) | .25988 (2/6) |
| Neural affine-offset global A+ | .31395 (5/24) | .51831 (4/12) | .22800 (2/6) |
| Neural affine-offset sourcewise A+ | .31293 (6/24) | .50160 (4/12) | .22833 (2/6) |

Entries are mean amplitude NRMSE with the inherited composite material-success count in parentheses. That label requires NRMSE≤.15, prediction-only Otsu IoU≥.7, and legality. It is not an experimental stopping gate. Smooth and nested materials make morphology an imperfect summary: 12 of the 15 ID neural-direct label failures fail only the IoU component. Full continuous rows, all primary arms and all 42 images are retained.

The neural initializer alone has mean NRMSE .14690/.21898/.12052 across the three cohorts. Direct physical refinement improves these values, whereas a simple attachment to an auxiliary solver often loses the useful trajectory. On predetermined object 0027, the neural starting objective is .30965 and a global-A+ candidate reduces it to .01565 before correction. Nevertheless, neural direct optimization terminates by projected-gradient tolerance after 1,814 accepted steps at Φ=.00042866 and NRMSE=.12082. Global A+ reaches its 3,000-step corrector cap at Φ=.00046350 and NRMSE=1.01370. Affine-offset Dual stops by the numerical projected-gradient tolerance after 2,789 accepted steps but still has NRMSE=1.06415. This is a concrete trajectory counterexample; it does not prove an inescapable basin.

Neural direct optimization costs 1.06–2.74 times the mean final CPU of uniform direct optimization across the five hardware/cohort strata. It is therefore a quality improvement, not a measured final-time acceleration in that comparison. At common per-object CPU/LU/RHS allowances, its mean material errors are lower in these strata. Relative to the five uniform and five neural auxiliary methods, its final mean deployment CPU and error are lower in all five strata; these ten detailed contrasts were added after inspecting the primary results and are explicitly descriptive. For example, ID/Mac mean neural-direct CPU is 32.03 seconds versus 62.12 for uniform Dual; the corresponding held-out-family/Mac values are 39.12 and 69.89. Training is additional. Equal-budget gains are not universal: ID/i9 neural direct versus neural Global-Dual has mean NRMSE difference +.00281 with descriptive 95% interval [−.06069,+.06401]. Sparse paid checkpoints, unused allowances and unadjusted small-stratum bootstrap intervals are reported.

The primary physics-feature model and its training/validation generation cost 402.56 Mac CPU seconds. Against uniform direct optimization there is no finite mean-CPU break-even because neural direct costs more per deployment. Against uniform Dual, the observed Mac average savings imply 14/14/13 warm deployments to amortize that training cost across the three cohorts. This is conditional CPU bookkeeping, excluding cold loading, not an economic or equal-quality guarantee; other research-only network costs remain in the full ledger.

Paying for neural direct and neural global A+ and selecting by final Φ gives mean NRMSE .07754/.34498/.22882, with 13/24, 6/12 and 2/6 success labels. It improves ID but worsens held-out-family and cross-grid mean error relative to neural direct alone, at 2.38–3.17 times its CPU cost. This negative result is not replaced by an offline best-material candidate. Global-Dual portfolios behave similarly, so the ID benefit is not specific to A+.

The fixed-eight additional initialization study is reported separately as ID4/held-out2/cross2. A Born-only network obtains 4/4 ID global-A+ labels versus 3/4 for the physics-feature network, but both obtain 0/2 on that held-out subset. This observation does not justify retrospectively replacing the frozen primary network or assigning the eight-object controls a 42-object denominator.

### 6.2 Development mechanisms and exposed historical regression

The eight fixed validation materials show why a single favorable initialization experiment is insufficient. On validation object 0003, ordinary neural global A+ and its truncated physical-current projection yield amplitude NRMSE 0.4249 and 0.7420. The fixed-offset Dual/global A+/sourcewise A+ variants yield 0.0315/0.0326/0.0326, while a uniform offset remains 0.4249. However, direct neural physical optimization is about 0.0297 and cheaper. On objects 0002 and 0005, the offset pipelines recover only because the selected auxiliary stage is −1: they retain the original neural material and reproduce direct correction. These outcomes cannot be credited to a superior auxiliary-current solution.

Object 0006 is a decisive negative control. Original joint initialization retains the initial material and reaches 0.0527, whereas the fixed-offset Dual/global/sourcewise variants yield 0.6846/0.5844/0.6036. All eight initial constitutive defects decrease when the truncated component is removed, yet final material recovery does not uniformly improve. This refutes that initial diagnostic as a standalone success predictor.

The four exposed historical targets contain two earlier successful curved objects and two failures. With fixed material-only neural initialization, global and sourcewise A+ each retain 2/4 historical recoveries, Dual recovers 0/4, and direct physical optimization recovers 1/4. The same two failures remain. These are regression/mechanism observations, not new confirmation samples; old and new software versions are labelled in the provenance table.

### 6.3 Classical Austria experiments

The study reconstructs the published 2 m domain acquisition rather than only drawing an Austria-shaped target. The original 150² generating grid and 64² inverse grid are actually used, alongside clearly labelled 32² adaptations and longer corrector follow-ups. For the εr = 4, 400 MHz original-grid comparison, amplitude NRMSE is approximately 0.441 for CIE, 0.427 for Dual, 0.501 for global A+, 0.440 for sourcewise A+, and 0.441 for the independent real-NIE implementation. Residuals are near 0.10, but reconstructions contain inner-ring artifacts and contrast peaks close to the upper bound 20 although the target contrast is 3.

The independent NIE code is a formula-based reproduction, not author-certified mature software. A separate domain-augmented initializer improves some early global-A+ endpoints, but longer physical correction can lower residual while worsening material error across methods. The result does not establish a stable A+ advantage on the original-paper class.

### 6.4 Measured Fresnel data

The 2001 single/twin-cylinder data and the 2005 FoamTwinDiel data are real acquired files with recorded provenance. Empty-target source-free multipole fitting and source-specific masks replace the earlier single-point cylindrical illumination approximation. Empty-field held-out error improves substantially, but cannot separately identify the transmit and receive horn responses.

For FoamTwinDiel at 8 GHz on 60², the completed CIE/Dual/Global/A/A+ held-out residuals are 0.26594/0.26589/0.26609/0.26593/0.26587. Their similar ring/spike material artifacts do not justify a recovery claim. A+ reduces A's material-inner evaluation count from 11,135 to 2,710 and its recorded auxiliary CPU from 76.93 to 57.74 seconds, but total cost also includes different physical corrector trajectories and shared-host timing. Dual alone reaches its relative-function stopping condition; the other four reach the 1,000-step cap. No broad speed or accuracy victory follows.

Frozen neural initialization of the 2001 cylinders does not establish held-out prediction improvement. The independent real-NIE Foam run uses a separate air (zero-contrast)/no-prior/no-corrector contract, reaches a smoother lower-amplitude estimate, and has worse held-out residual around 0.348; without voxel truth or matched contracts it is not a material-accuracy winner. A separately labelled 6 + 8 GHz follow-up compares added information with an 8 GHz continuation from the same paid anchor. Both runs reach the cap of 300 accepted steps and 306 objective evaluations. Adding 6 GHz data reduces the 6 GHz held-out residual from 0.33652118 to 0.18644632 and the two-frequency RMS residual from 0.30330496 to 0.23243573, while slightly increasing the 8 GHz held-out residual from 0.26597236 to 0.27072220. Recorded CPU time rises from 285.20 to 538.92 seconds. Peak contrast falls from 7.6393 to 6.81016 but remains high. These results establish a prediction tradeoff across frequencies, not material restoration or an A+ advantage.

A subsequent one-change model diagnostic removes the inherited positive contrast floor, which excludes air in the measured ROI. From the same paid anchor, replacing the lower bound 0.02 by zero reduces the 8 GHz held-out residual from 0.26597 to 0.13613, while peak contrast increases from 7.639 to 7.964. The improvement therefore concerns prediction and background feasibility, not restored material accuracy. Fresh-start Dual-CIE and A+ controls then use the corrected air-zero bound consistently in every auxiliary and physical update, with the original initialization 0.1, prior, calibration and 3,000/1,000 iteration limits unchanged. They obtain held-out residuals 0.137237 and 0.137128, CPU times 741.18 and 750.06 seconds, and peaks 6.929 and 6.912; both correctors reach the 1,000-step limit. Their estimates differ by 2.35% in relative material norm and by 0.1215% in acquired-channel predictions. This does not establish an A+ imaging advantage. All 25 measured-data endpoints and the earlier failed attempts are retained.

### 6.5 Genuine vector three-dimensional inversion

TensorDDA/dyadic Green experiments retain three-component currents, 60 external electric-dipole illuminations and actual three-dimensional materials. The final original projected-IDFT PR comparison for a small four-cube object gives relative-L2 material errors about 0.4632/0.4620/0.4637/0.4630/0.4624/0.4608 for LSIE/CIE/Dual/Global/A/A+. A+ costs 40.84 CPU seconds versus approximately 17.47–22.09 for the other arms. The small error differences do not establish a cost-competitive advantage.

An 8³ inverse-grid refinement gives approximately 0.249/0.251/0.249/0.249 for LSIE/Dual/A+/direct physics, with appreciable cross-grid discrepancy remaining. A literal hollow-cube case remains difficult. A separate Mie test improves spherical discretization error using polarizability-preserving partial-cell treatment, but does not certify all imaging examples. A real 3λ-domain, 16³-forward/12³-inverse attempt is retained as a 20-auxiliary/3-corrector-step truncation with 379.59 CPU seconds, not as an original 60³/30³-scale reproduction or converged reconstruction.

## 7. Discussion

The central distinction is between useful material information, compatibility with an auxiliary-current representation, and the endpoint chosen after nonlinear correction. A good initializer can be lost when its current is constrained; retaining more of that current can help, do nothing because the initializer is selected anyway, or steer the auxiliary solver toward a lower initial objective that corrects poorly. The observed counterexamples require all three effects to be reported.

Similarly, a more efficient conditional A+ update has a legitimate engineering benefit even when it does not improve material recovery. Such a benefit should be compared with mature ordinary readouts and a direct physical control at actual cost, rather than presented as a new information source. Conversely, two historical objects show that a material network does not automatically replace the useful trajectory effect of A+.

The classical and measured experiments highlight regularization and model adequacy. When the material class is underdetermined or calibration/discretization error is non-negligible, an increasingly accurate data fit can produce a less credible material. Training priors, auxiliary restrictions and stopping rules all participate in material selection. Their influence must not be concealed behind an unchanged explicit regularizer.

## 8. Limitations and Reproducibility

The independent cohorts are modest and heterogeneous; repeated algorithm endpoints are not independent statistical samples. Hardware-stratified CPU comparisons and object-paired resampling are reported without converting a small cohort into a universal superiority claim. Training and feature-generation costs are separated from deployment, and any amortization threshold uses mean paired savings and an explicit hardware-cost convention.

The primary scalar-TM model inherits a finite disk-pulse discretization with radius 0.45h. Cross-grid tests reduce a same-grid inverse crime but remain within that discretization family. Measured horn calibration is incomplete and exact voxel truth is unavailable. The initial five-arm measured comparison used a positive lower bound that excludes air; the corrected two-arm comparison removes this defect but retains uncertain illumination/receive response and spiky material estimates. Three-dimensional meshes remain small, the original large electrical-size attempt is truncated, complex A+ circular-arc updates are not claimed, and no 3D neural initializer or 3D measured reconstruction is reported.

All new code, configurations, observations, frozen network weights, raw predictions, paid trajectories, failures and cost receipts are archived additively under MAX-STRUCT-R8 for the final public release. Third-party measurement files and papers are not repackaged; acquisition scripts and provenance are provided. Test truth is read only after public predictions and candidate selection are sealed. Failed and interrupted runs are retained with unknown timing explicitly marked rather than replaced by zero.

## 9. Conclusion

The completed evidence supports a study of learned-state transport, candidate ordering and conditional-readout cost in Maxwell inversion. A lightweight learned material followed by direct physical optimization provides the clearest synthetic quality/cost result against the tested auxiliary pipelines. Simply attaching that material to sourcewise or global A+ often loses its benefit. A fixed physical-complement offset supplies a simple partial repair, together with confirmed counterexamples, and post-correction lowest-objective selection introduces another material-selection failure. Historical A+ gains remain conditional; classical, measured and three-dimensional tests do not establish broad superiority. Neural initialization, projection and candidate retention are established ideas, so this manuscript's publication case rests on the specific transport/selection evidence and its reproducibility. A reliable new data-only predictor or mechanism that preserves neural-direct gains while adding recovery beyond it remains missing. The evidence does not justify presenting the work as a universally superior A+ algorithm.

## References

1. Y. Zhong, M. Lambert, D. Lesselier, and X. Chen, “A New Integral Equation Method to Solve Highly Nonlinear Inverse Scattering Problems,” *IEEE Transactions on Antennas and Propagation*, 64(5), 1788–1799, 2016. [DOI: 10.1109/TAP.2016.2535492](https://doi.org/10.1109/TAP.2016.2535492).
2. Y. Zhong and K. Xu, “Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems,” *Journal of Imaging*, 5(2), article 27, 2019. [DOI: 10.3390/jimaging5020027](https://doi.org/10.3390/jimaging5020027).
3. Y. Sanghvi, Y. Kalepu, and U. K. Khankhoje, “Embedding Deep Learning in Inverse Scattering Problems,” *IEEE Transactions on Computational Imaging*, 6, 46–56, 2020 (early online 2019). [DOI: 10.1109/TCI.2019.2915580](https://doi.org/10.1109/TCI.2019.2915580). Full author manuscript inspected.
4. M. Zhou, J. Han, M. Rachh, and C. Borges, “A neural network warm-start approach for the inverse acoustic obstacle scattering problem,” *Journal of Computational Physics*, 490, article 112341, 2023. [DOI: 10.1016/j.jcp.2023.112341](https://doi.org/10.1016/j.jcp.2023.112341).
5. A. D. Magdum, H. Shimoga Beerappa, and M. Erramshetty, “Deep learning based distorted Born iterative method for improving microwave imaging,” *Frequenz*, 78(1–2), 1–8, 2024 (online 2023). [DOI: 10.1515/freq-2023-0074](https://doi.org/10.1515/freq-2023-0074). Abstract-level evidence; no unverified implementation details claimed.
6. T. Shan, Z. Lin, X. Song, M. Li, F. Yang, and S. Xu, “Neural Born Iterative Method for Solving Inverse Scattering Problems: 2D Cases,” *IEEE Transactions on Antennas and Propagation*, 71(1), 818–829, 2023. [DOI: 10.1109/TAP.2022.3217333](https://doi.org/10.1109/TAP.2022.3217333). Abstract-level evidence.
7. K. Xu, Z. Qian, R. Song, X. Ye, N. Xu, X.-M. Pan, P. Zhao, S. Chen, G. Wang, and W. Li, “Physically Unrolling Network Under Contraction Integral Equation for Limited-Aperture Inverse Scattering Problem,” *IEEE Transactions on Antennas and Propagation*, 71(11), 9130–9135, 2023. [DOI: 10.1109/TAP.2023.3305834](https://doi.org/10.1109/TAP.2023.3305834). Abstract-level evidence.
8. A. Dubois, K. Belkebir, and M. Saillard, “Retrieval of inhomogeneous targets from experimental frequency diversity data,” *Inverse Problems*, 21(6), S65–S79, 2005. [DOI: 10.1088/0266-5611/21/6/S06](https://doi.org/10.1088/0266-5611/21/6/S06). [Primary manuscript](https://www.fresnel.fr/perso/belkebir/Articles/IP_Dubois05.pdf). Phase convention and nominal Foam geometry inspected.

Full bibliographic and claim-level evidence, including IMSA-CIE, FBE-CIE and LUCIE, is retained in `PRIOR_ART_AND_NOVELTY.md`; title/abstract evidence is not represented as full-paper reproduction.
