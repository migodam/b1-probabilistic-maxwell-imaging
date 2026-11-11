# Reflected Shared-Material Feasibility for Maxwell Inversion: Basin Intervention and the Limits of Physical Selection

**Research manuscript draft, MAX-STRUCT R12; original study 11 November 2026, paired follow-up 12 November 2026.** This is an evidence-constrained draft, not a claim of acceptance or a camera-ready manuscript. The independent cohort has been executed, sealed and scored; detailed tables and figures accompany the confirmation report. The implementation, exposed diagnostic cases and confirmation cohort must remain distinguished.

## Abstract

Nonlinear electromagnetic inversion can fail because an optimization trajectory remains in a poor region, because the available data admit inaccurate material estimates, or because a discrete model does not explain the experiment. We study these mechanisms using reflected feasibility in the space of total fields and contrast currents. One set enforces background propagation and a measured-data ball; another enforces a real scalar material shared across all illuminations and vector components at each cell. The constitutive projection reduces to a constrained shared-slope total-least-squares problem, while the propagation projection reuses one fixed background factorization. From three exposed poor endpoints of inherited CIE/Green-feedback solvers, reflection returns more accurate publicly selected materials than direct continuation, alternating projections and a specified smooth-restart control at paid CPU prefixes. This conditional benefit does not establish a generally superior default solver: post-seal scoring of 24 independent objects gives mean amplitude-normalized RMSE 0.1846 for reflection, 0.1855 for Dual-CIE and 0.1439 for the inherited Green-feedback method. Vector experiments reveal an additional failure mode: a common data-plus-prior corrector can degrade an accurate, noise-feasible reflected material while reducing its objective. Measured-data and cross-grid tests expose model and selection limitations. The resulting contribution is a controlled study of search, physical feasibility and material selection, rather than a new global convergence theorem or universal replacement for CIE.

A frozen follow-up on these now-exposed objects adds 44 same-start runs. Triggering rescue only when the original Green-feedback endpoint lies outside its declared observation ball improves three of four triggered objects; the other twenty are retained without rescue. The complete 24-object mean error decreases from 0.143891 to 0.106261 at 7.85% additional actual CPU cost. A matched legal G19 restart does not reproduce these three improvements, although it outperforms reflection on the older ring endpoint. Four frozen-neural-start experiments further show comparable G19 and reflected primary means (0.086768 and 0.086094), with about 31% lower measured CPU for reflection. These follow-up results support a conditional rescue role, not a new blind confirmation or attribution of every neural-interface gain to reflection.

## 1. Introduction

High contrast changes both the forward response and the behavior of inverse algorithms. A low auxiliary residual need not imply a physical material, and a low full-wave data residual need not identify the correct material. A method that escapes one bad endpoint may still be inferior when applied from a uniform initialization across a new material population. These distinctions motivate an intervention study rather than a sequence of favorable reconstructions.

The R4–R11 implementation lineage supplies mature CIE/FFT-TSOM, Dual-CIE, global-current controls, A-material/A+ and Green feedback, together with saved failure endpoints. The present study preserves those records. We implement a field/current feasibility scheme suggested by the supplied Pro research package, reproduce its small experiments, add the missing mature baselines, intervene at actual saved endpoints, and test the method on an independent cohort and additional physical models. Ground truth is excluded from online selection in the new solver cohorts. Original exposed toy replay retains its original diagnostic truth logging and is explicitly not a blind experiment.

Our questions are whether reflected memory changes accessible physical regions, whether that helps beyond additional computation, and whether the material selected after search retains the benefit. We do not require victory on every object, but do require a stated domain of benefit and explicit counterexamples.

## 2. Related work and novelty scope

Contrast-source inversion already introduces auxiliary currents [1]. CIE changes the integral-equation representation and its nonlinear material coupling; FFT-TSOM controls the current representation [2]. WRI enlarges the search space through a wave-equation penalty [3], and DCTMC permits nonlocal intermediate interaction potentials while seeking a data-compatible local material [4]. Electromagnetic ADMM and PADMM-CIE establish relevant splitting precedents [5,7]. Douglas–Rachford and RRR are established feasibility tools; projection-based phase retrieval already distinguishes iteration fixed points from recovered objects [6].

Consequently, neither auxiliary variables, nonphysical intermediate states, reflection, nor quadratic elimination is claimed as a new principle. The specific real-material shared projection and finite-channel implementation are an application-level construction. The experimentally testable contribution is the controlled distinction between endpoint escape and final material selection. A full equation-level comparison and numerical benchmark against PADMM-CIE and WRI remain missing; no superiority over those methods is asserted. Primary sources and access limitations are documented in [the prior-art audit](PRIOR_ART_AND_NOVELTY.md).

## 3. Method

### 3.1 Finite Maxwell contract

For the adopted time convention exp(-iωt), the finite volume-integral model is

\[
Q=E^{inc}+GJ,\qquad J=D_mQ,\qquad F(m)=SJ.
\]

The unknown m is real isotropic contrast in a declared box. In scalar TM each cell shares m across sources; in vector DDA it additionally shares m across three components. These are distinct physical discretizations. Comparisons within a cohort use its identical operators, data, legal box and physical regularizer. Initialization is declared by arm: NN arms use a learned material initializer while ordinary arms use the specified uniform material. In particular, the saved-endpoint AP/DR/direct intervention controls share the same complete physical initialization. The frozen NN adds learned prior information even though its downstream physical regularizer is unchanged. Continuous Maxwell accuracy is not inferred from finite-matrix agreement.

Set X=J/s with fixed s=6. In the product norm, define the constitutive set C by X_i=(m_i/s)Q_i and the propagation/data set A_delta by Q=Einc+sGX and the global whitened or explicitly raw data-ball constraint. Source-dependent masks are represented by their actual port operators. One common multiplier enforces the total data ball, rather than separate balls per illumination.

### 3.2 Projections and reflection

For a cell with scaled vectors x,q and real alpha=m/s, elimination of its projected field gives

\[
\widehat q=(q+\alpha x)/(1+\alpha^2),\qquad
\varphi(\alpha)=\|x-\alpha q\|^2/(1+\alpha^2).
\]

All admissible roots of c alpha²+(b-a)alpha-c=0 and both box endpoints are evaluated, where a=||x||², b=||q||² and c=Re(x*q). This is an exact finite-dimensional minimization, subject to floating-point ties and tolerances; complex or tensor-valued constitutive classes require different projections.

The propagation projection eliminates Q and factors H=I+s²G*G once. A small port spectral system and a scalar multiplier impose the data constraint. Unreachable residual components are retained. Grouping identical source-port maps into multiple RHS preserves the mathematical problem and improves the tested small projection workload, without a claim of bitwise identical nonlinear trajectories.

DR uses c=PCz, a=PA(2c-z), z+=a-c. AP instead projects c directly onto A. Every 50 steps a material shadow is evaluated by the original full-wave solver. Primary candidates minimize the inherited physical objective Phi=D+R. The common corrector optimizes the same objective. All evaluations and background projection operations are charged.

### 3.3 Secondary constrained candidate selection

Before any independent confirmation inversion or scoring, a secondary rule was registered after exposed 3D diagnostics: among already evaluated, legal, cost-eligible materials inside the raw CN99 ball, select the smallest same R, breaking ties by Phi. It performs no extra forward solves and does not alter trajectories. It is a separate estimator, not a replacement for primary results. No qualifying candidate means unavailable. Means over available secondary cases cannot be compared directly with primary means over all 24 objects.

## 4. Mechanism and conditional analysis

A complete physical initialization belongs to C. For relaxation one, AP and DR therefore have the same first state PAz0. Their subsequent difference is reflected memory, not a hidden change in the initial material. A physical point already in the data ball is fixed for both methods, even when its material is inaccurate. Thus reflection cannot solve an information or selection problem merely by continuing iteration.

For constitutive and propagation shadows c,a, let g=||c-a||, e=Qc-Einc-GJc and suppose I-DmG is invertible. Then

\[
\|e\|\le\sqrt{1+s^2\|G\|^2}\,g,
\qquad J_c-J_{phys}=(I-D_mG)^{-1}D_m e.
\]

This gives a conditional bound on physical data error involving both g and a material-dependent resolvent norm. Small shadow gap alone is insufficient near poor conditioning, and no bound on error to ground truth follows. The iteration-zero archive gap is only an initialization placeholder and is excluded from this argument.

Reflection does not change the stationary points of the material objective Phi. It searches in another state space and may temporarily increase Phi. In the exposed ring replay the audited objective initially increases before the trajectory reaches a low-residual region, whereas AP retains the bad incumbent. This is evidence of a different search route, not proof of global convexification.

## 5. Experimental design and results

### 5.1 Toy reproduction and mature baselines

The original four-object code is replayed without changing its physical model. Forty mature-baseline runs cover two explicitly distinguished mode schedules. The toy ring remains difficult for the added CIE-family controls, but the shell is reconstructed better by those controls than by the original DR replay. The new controls allow longer correction and adapted current schedules, so this is a recoverability control rather than an equal-compute attribution. This prevents interpreting the supplied toy result as universal superiority over CIE. Small BLAS-dependent differences in the bridge trajectory and corrector degradation are preserved.

### 5.2 Intervention at saved R11 endpoints

Three real saved bad materials initialize every arm with freshly solved full physical fields. At paid 240 body-CPU prefixes, amplitude NRMSE for direct continuation / smooth restart / AP / DR is respectively:

|Exposed object|Direct continuation|Smooth restart|AP|DR|
|---|---:|---:|---:|---:|
|Ring|0.42288|0.74888|0.42288|0.16585|
|Split shell|0.51514|0.45633|0.51514|0.11377|
|Bridge|0.43673|0.43673|0.35144|0.10617|

Ring and shell endpoints have near constrained-local-minimum evidence; the bridge has not converged equivalently and is reported as an endpoint/path intervention. Some old solvers recover these objects from other initializations. The claim is therefore rescue of specified bad endpoints, not discovery of previously unrecoverable objects. Early archive overhead caused some full runs to exceed nominal windows; the table uses archived paid prefixes, and full costs remain separate.

Post-seal oracle checks of the actually saved, cost-eligible physical candidate pools give best control NRMSE0.37454/0.45202/0.33477, still above publicly selected DR. This diagnostic does not enter online selection and does not cover unaudited intermediate shadows or other restart distributions. It strengthens the distinction between a search difference and merely choosing poorly from an equivalent recorded pool.

### 5.3 Independent matched 2D confirmation

Twenty-four independently generated materials cover level sets, branching structures, islands and continuous random fields, with high8/high18 and mixed regimes. The 576 real unknowns, 12 sources, 48 receivers, 1% circular complex noise, legal box and original prior are fixed. Eight methods receive 240 whole-child CPU caps, with prescribed proposal/corrector windows and 30/60/120/240 prefixes. All 192 predictions were sealed before scoring. Six additional development objects are not pooled with confirmation.

The frozen primary scores are:

|Method|Mean amplitude NRMSE|Descriptive success /24|
|---|---:|---:|
|CIE|0.202932|9|
|Dual-CIE|0.185456|9|
|Global-Dual|0.184007|9|
|A+|0.171692|9|
|G19|0.143891|9|
|DR|0.184631|7|
|Physical smooth multistart|0.619361|4|
|Frozen NN→direct|0.149514|9|

These data do not support replacing mature G19 by DR as the default solver. Pairwise uncertainty, per-family effects, feasibility classifications and full galleries are supplied by the post-seal scoring report; the old success threshold is descriptive, not an experimental stopping gate. All cases, including shared failures, remain in the denominator.

The paired mean difference DR minus Dual is−0.000825 (95% material-bootstrap interval[−0.040982,0.053017]); DR minus G19 is0.040740 ([0.001934,0.095095]). These unadjusted exploratory intervals neither establish equivalence with Dual nor universal inferiority. Eleven of the17 DR threshold failures are already noise-feasible and have amplitudeNRMSE≤0.15; their remaining failure is the fixed IoU threshold, with relative material L2 errors0.324–0.600. They must not be described as wholly absent reconstructions. The other six lie outside the ball. At60CPU seconds meanNRMSE is0.14988 for G19 and0.18396 for DR, so the final average difference cannot simply be attributed to G19 spending longer.

### 5.4 Neural initialization

Four exposed R8 materials compare uniform→direct, NN→direct, NN→CIE and NN→DR using the frozen 24,737-parameter model. Mean amplitude NRMSE is approximately 0.4470, 0.1419, 0.3089 and 0.0862, respectively. NN→DR costs more than NN→direct and is not better on every object. This supports a conditional interface benefit, not a new learned algorithm or independent generalization claim. Historical training is inherited and disclosed. The independent cohort tests NN→direct, not every neural combination.

### 5.5 Vector Maxwell and measured data

Small vector-DDA experiments use 216 real material cells, 648 field rows, 60 dipole illuminations and 180 receive components. With a longer common corrector, moderate-contrast method differences largely disappear. In two high-contrast matched cases, DR raw relative material errors near 0.0616 and 0.0783 become approximately 0.1298 and 0.2272 after correction. The raw candidates are noise-feasible; the corrected candidates leave that ball while decreasing the smoothing prior enough to reduce Phi. A cross-grid response mismatch of about26% precludes a noise-only continuous-accuracy claim.

Fresnel single/twin measurements use calibrated incident fields and actual masks, with held-out prediction. Full native Dual continuation and DR end near the same corrector-dominated fits; no measured-material superiority is demonstrated. Austria150→32 exposes significant model mismatch: the coarse GT itself lies outside the chosen data ball. A matched32 follow-up makes GT feasible but DR still fails, so mismatch is not a complete explanation. All these negative controls remain part of the study.

### Post-hoc structural interpretation requested during review

The historical conjunction (amplitude NRMSE, thresholded IoU, and legality) is not a definition of all useful imaging. After the user highlighted case0009, the lead reviewer inspected all 24-by8 primary panels using the same three descriptive labels: S for recognizable principal overall structure, P for recovered dominant components with important weak-structure limitations, and F for an incorrect principal structure. This unblinded, post-hoc visual interpretation does not replace frozen quantitative scores or create a new confirmatory success rate. DR has10 S and8 P cases, versus19 combined S/P for Dual and20 for G19 or NN-direct. Case0009 recovers the dominant high-contrast branch with IoU near0.5, while weak branches and material values remain biased. G19 in case0007 and several methods in case0013 show genuine structural increments hidden by the binary cutoff; DR in0019 and all but NN in0022 supply corresponding negative controls. The full per-object labels and neutral plots are provided in STRUCTURAL_REVIEW_ZH.md. No algorithm, candidate, or online choice was changed. Independent multi-reader agreement has not been established.

## 6. Computational costs

Dense projection setup has cubic work and quadratic storage; iteration uses fixed-H multiple RHS and Green/port actions. These must not be replaced by a count of material-dependent LU factorizations alone. The 192 confirmation children consume28,147.796875 process-CPU seconds. Failed startup attempts, generation, scoring and transport are separate paid entries. Recorded under-load i9 telemetry reported low active-core frequency; baseline code explicitly sets BelowNormal priority whereas search inherits its controller, with only a single effective-priority snapshot. CPU caps remain recorded, but strict identical-scheduling wall-throughput claims are unavailable. No CUDA scientific execution is claimed for this campaign.

### Use-specific positive evidence

Default-solver ranking does not exhaust algorithmic value. Saved-endpoint intervention is the strongest positive DR evidence; the frozen neural interface adds a second, exposed-development use case. Across four NN objects, NN–DR mean NRMSE is0.0862 versus0.1419 for NN–direct and0.3089 for the inherited NN–CIE interface, with recorded body CPU totals78.43,43.44 and104.46 seconds. The gain is concentrated in direct-stagnation repair on0022 and preservation of useful structure on0027; representation and reflection are not separately identified. Conversely, G19 supplies independent-cohort evidence in routine inversion: mean NRMSE0.1439 versus Dual0.1855, with positive control-minus-G19 paired mean differences and an advantage already present at the60-second paid prefix. Its higher setup/full cost and failures remain. These observations support complementary, conditional roles rather than dismissing both methods for lacking universal dominance. COMPREHENSIVE_EVALUATION_ZH.md provides the full evidence matrix.

## 7. Discussion

The decisive distinction is between reaching a physically useful candidate and choosing the final estimate. Reflected memory can escape particular bad trajectories without globally reshaping Phi into a convex landscape. A common corrector can remove an apparent short-budget advantage, or can degrade a better noise-feasible material through the data/prior tradeoff. The latter is not automatically a derivative error. A constrained selector cannot repair a pool that never contains a suitable candidate.

This suggests a conditional role for reflection as an intervention after demonstrated physical search failure, rather than an unconditional front end. Such a cascade still requires a new frozen selection rule and independent evaluation; the present paper does not claim that it has been established by retrospectively choosing the best arm per material.

## 8. Limitations

The matched confirmation cohort fixes acquisition and noise; development studies on other apertures and models are not an independent factorial confirmation. Vector grids are coarse at high index. The material class is real isotropic. A uniform resolvent bound and nonconvex convergence theorem are absent. WRI/EM-ADMM/PADMM-CIE are not qualified numerical baselines here. The current dense implementation is not a scalable matrix-free CUDA solver. CPU clocks, priorities and incomplete preparation timers limit global speed claims. Real measurements lack pixelwise truth and therefore have no fabricated material NRMSE.

## 9. Conclusion

Reflected shared-material feasibility has demonstrated conditional endpoint-escape capability and has exposed a consequential mismatch between physical-objective descent and material selection. It has not established general superiority over mature CIE-family solvers, stable high-index continuous3D accuracy, or measured-data imaging gains. A credible manuscript should focus on these controlled mechanisms and their boundaries. Further algorithm claims require new, explicitly frozen tests rather than reclassification of exposed examples.

## Paired follow-up: shared endpoints and accounting

The [supplement](supplement_paired/RESEARCH_REPORT_ZH.md) contains three old bad materials, all 24 original G19 endpoints with an observation-only out-of-ball trigger, and four identical frozen neural materials. The four restarted methods are G19, reflection, alternating projections and direct original-physics refinement. Their whole-child CPU caps are 240 seconds for old/NN starts and 120 seconds for incremental rescue. Original parent CPU and newly spent CPU are reported separately. The G19 restart uses the already implemented least-squares current initialization in the unchanged affine DPCS/Fourier family; PR history is reset. It is not a continuation of unavailable optimizer history or a full-current affine-space extension.

|Triggered R12 object|Original G19 NRMSE|G19 restart|Reflection|AP/direct approximately|
|---|---:|---:|---:|---:|
|0001|0.4318|0.4318|0.1072|0.4318|
|0003|0.4297|0.4297|0.1191|0.4297|
|0010|0.3658|0.3656|0.3656|0.3656|
|0022|0.2959|0.2959|0.0283|0.2959|

Three reflected candidates become noise feasible and improve materially; only two satisfy the inherited joint IoU/NRMSE descriptive criterion. The0001 image recovers the principal curved structure with IoU0.659. The failure0010 and all twenty bypassed objects remain in the denominator. Actual added CPU is336.875 seconds for reflection versus355.000 for G19 restart,239.328 for AP and62.500 for direct; equal caps are not equal expenditures. Total G19-parent CPU is4290.1875 seconds. Parent plus reflection therefore costs4627.0625 seconds. Candidate selection uses public original Phi, with a separately reported noise-feasible minimum-prior rule and explicit fallback.

The old ring can also be recovered by G19 restart (0.114232 versus reflected0.140651); reflection remains better on the shell and bridge (0.112699/0.097269 versus0.515137/0.371690). In the four neural-start cases, AP also repairs the direct failure0022, while AP alone produces a poor structure on0001 at essentially the same actual cost as reflection. G19 no longer reproduces the earlier native single-beta CIE interface failure on0027. Readout, PR and initialization differ, so this comparison does not isolate current projection as the sole cause.

The common corrector improves the reflected mean from0.200448 to0.116873 on the old starts, from0.110176 to0.106261 across the24-object rescue strategy, and from0.100916 to0.086094 on neural starts. Neural cases0004/0022 nevertheless worsen materially. The previous vector counterexamples remain. The method contribution is consequently a controlled conditional intervention and a separation of search from estimation, not a monotone relation between objective and material accuracy.

The follow-up is posthoc on exposed objects, even though its policies and budgets were frozen before its own runs. Fresh confirmation of the complete conditional pipeline and qualified continuum/vector generalization remain open. See [paired figures](supplement_paired/GALLERY.html), [raw/primary/secondary tables](supplement_paired/RESULTS/offline_v2/) and [cost ledger](supplement_paired/COST_LEDGER.md). No new algorithm or neural training was introduced in this supplement.

## References and reproducibility

See [PRIOR_ART_AND_NOVELTY.md](PRIOR_ART_AND_NOVELTY.md), [REPRODUCE.md](REPRODUCE.md), [BENCHMARKS.md](BENCHMARKS.md), [COST_LEDGER.md](COST_LEDGER.md), original task and source package under `sources/`, and the final artifact manifest. No omitted experiment is implied by this draft.


### Bibliography

1. P. M. van den Berg and R. E. Kleinman, “A contrast source inversion method,” *Inverse Problems*13,1607–1620(1997). [DOI](https://doi.org/10.1088/0266-5611/13/6/013). Metadata verified; the new full-text reader was blocked.
2. Y. Zhong and K. Xu, “Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems,” *Journal of Imaging*5(2),27(2019). [Primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/), DOI10.3390/jimaging5020027. The paper itself describes distinct current/material contraction parameters; this is not a new R12 idea.
3. T. van Leeuwen and F. J. Herrmann, “Mitigating local minima in full-waveform inversion by expanding the search space,” *Geophysical Journal International*195,661–667(2013). [DOI](https://doi.org/10.1093/gji/ggt258); [author formulation](https://slim.gatech.edu/research/full-waveform-inversion).
4. Levinson and Markel, “Solution of the nonlinear inverse scattering problem by T-matrix completion. I. Theory,” *Physical Review E*94,043317(2016). [Accepted primary manuscript](https://link.aps.org/accepted/10.1103/PhysRevE.94.043317). The inspected abstract/introduction directly establish the nonlocal intermediate-potential precedent.
5. J. Liu, H. Zhou, L. Chen and T. Ouyang, “Alternating direction method of multiplier for solving electromagnetic inverse scattering problems,” *International Journal of Microwave and Wireless Technologies*12,790–796(2020). [DOI](https://doi.org/10.1017/S175907872000015X). Publisher abstract inspected; full update equations not obtained.
6. V. Elser, “Phase retrieval by iterated projections,” *JOSA A*20,40–55(2003), [author preprint](https://arxiv.org/abs/math/0111080), DOI10.1364/JOSAA.20.000040. Primary abstract inspected; not used to assert nonconvex global convergence for Maxwell.
7. “Proximal Alternating Direction Method of Multipliers for Solving3-D Electromagnetic Inverse Scattering Problems,” *IEEE TMTT*72(2),981–995(2024), [DOI](https://doi.org/10.1109/TMTT.2023.3298203). Metadata inherited from the supplied research package and discovery record; closest full-text comparison remains unresolved.
8. Y. Zhong, M. Lambert, D. Lesselier and X. Chen, “A New Integral Equation Method to Solve Highly Nonlinear Inverse Scattering Problems,” *IEEE TAP*64(5),1788–1799(2016). [DOI](https://doi.org/10.1109/TAP.2016.2535492). The Austria setting is inherited through R8; our adapted experiment is not a literal reproduction of every original-paper choice.

The unpublished supplied Pro report and the frozen R4–R11 repository records are implementation provenance under `sources/pro/` and the corresponding repository study directories. They are not substituted for external peer-reviewed novelty evidence.
