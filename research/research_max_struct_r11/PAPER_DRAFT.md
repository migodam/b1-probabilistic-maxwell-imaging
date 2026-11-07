# Green residual metrics and block interactions in high-contrast Maxwell inversion

**Working manuscript — partial results, not a submission-ready claim.** The independent sixteen- and forty-eight-object studies, eight-object paired-grid study and exposed mechanism experiments below are complete. The new vector-Maxwell study is running; its results are not presumed here. Per the owner's instruction, no additional algorithm or experiment will be launched after these frozen computations.

## Abstract

A material-dependent Green denominator alters contrast-source inversion trajectories without adding measurement information or enlarging the Fourier state space. We separate constitutive correction from variation of the residual metric, and show why freezing a convex material block does not convexify the Maxwell inverse problem. Actual saved states do not support attributing all gradient cancellation to the new Green term. Applying the Green residual to both current and material repairs an exposed bridge but loses a new annulus in a separately frozen sixteen-object study. In a larger forty-eight-object comparison, CIE, Dual-CIE, the matched global-current control, material-only Green, and a partially frozen feedback rule meet the material criterion on 21, 22, 23, 23, and 23 objects. The latter three succeed on exactly the same objects. Partially freezing feedback repairs two substantial Green-induced error increases and reduces its valid-run CPU by about 5.5 percent, but does not establish superior recovery to Dual-CIE at common measured CPU prefixes. Conditional positive evidence remains substantial: in eight additional paired-grid materials, Green recovers an annulus and a slotted shell missed by all three controls, attaining six versus four material successes under both matched and cross-grid observations. Complete-run Green CPU is about21 and35 percent higher than Dual-CIE in these two conditions. The common-prefix advantage persists in count but changes one contributing object, showing why endpoints, costs and stopping must be reported together. All cohorts retain failures, independent material units and the original-Maxwell corrector. The evidence supports conditional recovery gains and explicit failure mechanisms rather than a universal denominator improvement or global convexification.

## 1. Introduction

High-contrast electromagnetic inverse scattering couples material contrast and total fields through repeated scattering. CIE and Fourier-TSOM introduce auxiliary currents and contraction variables to make this interaction more tractable. Changing the residual used to read out the material can create useful trajectories, but it can also stabilize a poor interpretation of an auxiliary current.

The issue considered here is narrower than finding a universally better optimizer. We ask which block is responsible for observed recovery changes, whether a simpler fixed-weight material update repairs detrimental behavior, and whether improvement survives independent objects and explicit computational accounting. A lower auxiliary objective, a lower Maxwell data residual, and a more accurate material are distinct outcomes.

Earlier project results motivate the study but are treated as exposed development evidence. In particular, a low-rank Green denominator recovered annuli and damaged a shell and a bridge. Those successes do not define a confirmation set. All new cohorts retain the full material parameterization and use public-data selection only.

## 2. Related work and scope of novelty

The comparison builds on Zhong and Xu's three-dimensional CIE formulation [1], established Fourier/current subspace methods, and earlier material preconditioning in extended CSI [2]. Matrix inverses, Woodbury reduction, lagged weights, quadratic material fitting, source Krylov spaces, and inertia under low-rank changes are established tools.

The vector-Maxwell volume operator has a substantially different spectral structure from a scalar compact-kernel model. Essential-spectrum results [3,4] motivate caution about transferring a small scalar Green truncation to a three-component dyadic discretization. The exact hypotheses of those continuous results are not replaced by a finite DDA experiment.

Our prospective contribution is an auditable mechanism and failure analysis, matched baseline comparisons, and a physically explicit assessment of the material-only and two-block residual modifications. We do not claim a new global convergence theorem, a first continuous Maxwell ambiguity construction, or a new spectral theory. Full publication positioning must also account for the ongoing confirmation study and the remaining closest-prior-art limitations documented in PRIOR_ART_AND_NOVELTY.md.

## 3. Model and methods

For a finite scalar-TM model, write

\[
J=D_m(E^{\rm inc}+GJ), \qquad Y=SJ+n.
\]

All inversions retain 576 independent real contrasts, the same box, observations, smooth/TV prior, initialization, Fourier/APCS stages, and original-Maxwell corrector. The online material selection objective is

\[
\Phi(m)=\frac{\|F(m)-Y\|_F^2}{2\|Y\|_F^2}+R(m).
\]

This data-normalized objective is distinguished from the declared iid circular-complex noise statistic. For \(E|n_j|^2=\sigma^2\), the reported noise radius is
\(\epsilon_{.99}=\sqrt{\chi^2_{2n,.99}/2}\) in units of \(\|F(m)-Y\|/\sigma\).

At fixed auxiliary current, let \(Q=E^{\rm inc}+GJ\), \(e=J-D_mQ\), and

\[
A_K(m,\beta)=I+\beta D_m-D_mK,\qquad w=A_K^{-1}e.
\]

Here \(K\) is a fixed public low-rank approximation. The designation “full” refers to the complete derivative of this material-dependent residual; it does not mean \(K=G\).

The compared material rules are the inherited ordinary CIE readout, the full derivative Green rule, a fixed-\(A\) rule, and a contracted-coordinate rule which freezes only
\(T(z)=I-D_zK\), where \(z=m/(1+\beta m)\). The latter preserves the prior as \(R(m(z))\). All extra operations are charged.

The two-block control additionally replaces the current residual derivative by

\[
H_K=\beta_j A_K(m,\beta_j)^{-1}(I-D_mG)\mathcal U,
\]

where \(\mathcal U\) is the same Fourier/APCS map. Its exact adjoint and quadratic direction line search are used. The inherited PR aggregation, representation capacity, and material space are unchanged. Since \(\beta_j\) and \(\beta_m\) differ, using the same residual family does not imply a common alternating objective.

## 4. Mechanism and limitations

### 4.1 The material derivative changes the effective illumination

For a real material direction \(h\),

\[
dw[h]=-A_K^{-1}D_h[Q+(\beta I-K)w].
\]

With \(a=A_K^{-*}(w\Omega)\), the data-part material gradient is the sum of a constitutive term proportional to
\(-\operatorname{Re}\sum_s\bar a_{is}Q_{is}\), and a metric-variation term proportional to
\(-\operatorname{Re}\sum_s\bar a_{is}[(\beta I-K)w]_{is}\).
The second term vanishes to second order in the constitutive residual near an exact root. It need not be harmless far from that root.

In ten actual exposed states, the ordinary \(\beta\) contribution is 2.288–5.390 times the Green contribution in norm. Early cancellation predominantly involves the ordinary contraction term. Thus the initial explanation that the new Green derivative alone causes the observed cancellation is not supported.

### 4.2 A convex material block does not convexify inversion

Freezing \(A_{\rm ref}\) gives

\[
q_{\rm ref}(m)=\frac{\rho\beta^2}{2}
 \|A_{\rm ref}^{-1}(J-D_mQ)\Omega^{1/2}\|_F^2+R(m).
\]

Its data Hessian is positive semidefinite, and the entire block is convex when \(R\) and the material set are convex. This is not an MM majorizer of the moving-denominator objective: equality of values at the reference point does not give matching gradients or a global upper bound. Refreshing the metric also does not guarantee monotonicity of the physical objective.

The exposed fixed-state experiment makes this distinction concrete. Fixed-\(A\) steps reduce the raw constitutive defect in all ten states but reduce the original physical objective in only seven. A coated-object step lowers \(\Phi\) while worsening material error. Full inverse replays preserve annulus recovery but fail to repair the shell and bridge.

### 4.3 The current effect is not merely a larger state space

In the restricted uniform-material, full-\(G\), eigenchannel model, the scalar and Green current factors are

\[
h_0(\lambda)=\frac{\beta(1-m\lambda)}{1+\beta m},
\qquad
h_G(\lambda)=\frac{\beta(1-m\lambda)}{1+\beta m-m\lambda}.
\]

The second factor saturates for large \(|m\lambda|\), whereas the first grows. This motivates a current-curvature mechanism. It is not a bound for nonnormal \(G\), nonuniform material, or finite \(K\).

Actual same-state diagnostics show only moderate conditioning improvement but large changes in reset direction. The shell's conditioning also improves without recovery. Furthermore,

\[
H_K\delta J=\beta\delta J-\beta^2A_K^{-1}D_m\delta J
-\beta A_K^{-1}D_m(G-K)\delta J.
\]

On saved reset directions, raw Green omissions can be large while their weighted contribution to the current derivative is only 1.2–5.85 percent. Static Green energy coverage is therefore insufficient to diagnose recovery.

## 5. Completed experiments

### 5.1 Exposed four-cell block ablation

The following amplitude-normalized errors are from three public development objects. They are not confirmation successes.

|Object|ordinary current/ordinary material|ordinary current/Green material|Green current/ordinary material|Green current/Green material|
|---|---:|---:|---:|---:|
|Annulus|.42256|.09890|.42203|.10517|
|Slotted shell|.11032|.51514|.51514|.51514|
|Bridge|.09629|.33700|.55462|.10088|

The bridge repair requires an interaction of the two changed blocks in this experiment. Neither individual change explains the successful two-block trajectory. The annulus instead benefits primarily from the material modification. Every real-shift Green variant in this four-cell ablation fails the exposed shell. These are trajectory interventions, not a unique causal identification of all intermediate nonlinear effects.

### 5.2 Independent sixteen-object study

Eight geometric families, each with a high and a mixed contrast object, were generated after freezing all five algorithms, rules, and budgets. All eighty outcomes and public selections were sealed before truth access.

|Method|Material criterion|CN99 compatible|Mean amplitude NRMSE|Mean actual child CPU, s|
|---|---:|---:|---:|---:|
|CIE|4/16|12/16|.23134|47.96|
|Dual-CIE|5/16|13/16|.21068|48.57|
|Global-Dual control|6/16|13/16|.21344|47.39|
|Material-only Green|7/16|14/16|.14485|60.85|
|Two-block Green|5/16|13/16|.17357|62.84|

The two-block method loses a new high-contrast annulus that material-only Green recovers. It also crosses the fixed IoU criterion on a shell with only a small continuous-error difference; these two regressions should not be described as equally large failures. At a common achieved physical-LU budget the success counts are unchanged. Auxiliary costs remain separate.

The paired mean difference between two-block and material-only Green is \(+.02873\), with an object-level high/mixed-stratified bootstrap interval \([+.000176,+.084961]\). This small-cohort comparison does not establish a universal ordering. It does reject promoting the exposed bridge result into a generally superior replacement.

The matched global-current control already recovers six objects. The additional success unique to the Green material modification is the new high-contrast annulus. The other two net gains over sourcewise CIE are also attained by ordinary-material Global-Dual. For this annulus, actual child CPU is 54.48 s for Green, 73.70 s for CIE, 64.79 s for Dual-CIE, and 47.40 s for Global-Dual. Thus the substantial annulus gain is not purchased with more CPU than CIE or Dual-CIE, although it remains more expensive than the matched global-current control.

A posthoc budget control continued eligible CIE/Dual endpoints up to the corresponding Green actual CPU expenditure. Seven publicly selected restarts consumed 62.86 additional CPU seconds and 7370 additional 576-dimensional factorizations; material success counts remained four and five. Selection used only original objective values attained within the remaining budget. Normal gradient stopping and all atomic overruns were retained, rather than padding runtime. L-BFGS history was reset, so this is a limited restart control and not a claim about every possible longer run.

### 5.3 Why data compatibility is not a material guarantee

Every ground-truth prediction lies inside its realized CN99 observation ball. Sixty-five of eighty candidates also lie inside, including thirty-eight material failures across eight independent objects. All sixty-five have both a lower data term and a lower prior than the ground truth.

These facts do not prove exact nonuniqueness or Gaussian statistical indistinguishability. In particular, the standardized prediction gaps for the material-failing compatible pairs are 4.30–6.09, not nearly zero, and candidates were selected from the same noisy data. The audit establishes that the adopted discrepancy criterion and lower \(\Phi\) do not certify material recovery. The remaining fifteen candidates are outside the ball although the ground truth is inside; their failure cannot be excused by this noise-compatibility argument.

### 5.4 Physical barriers between existing public endpoints

We evaluated nine preselected points on each legal material segment joining the higher-objective CIE/Green endpoint to the lower-objective endpoint for the three exposed objects. No truth material defined the direction. The annulus and shell had 467 and 562 of 576 variables at the lower bound, respectively, and projected-gradient infinity norms of approximately \(10^{-8}\). Their initial directional derivatives toward the better endpoint were positive: 0.01514 and 1.33849. Objectives rose from 0.00939 and 0.17393 to sampled maxima of 0.15283 and 1.59608, before reaching approximately \(5\times10^{-5}\). The bridge also showed a sampled barrier, but its projected gradient was larger.

Almost all of each sampled peak came from the physical data term. These measurements distinguish a physical path barrier from a purely auxiliary-loss issue. They establish neither local minimality nor barriers along every possible curved path. They show why a local descent correction at the wrong bound-active endpoint and an auxiliary trajectory that reaches a different endpoint are distinct capabilities. The low-objective endpoint was itself obtained by a paid alternative algorithm; it is not a free proposal generated by this diagnostic.

The exact physical Hessian restricted to the free box face further differentiates these endpoints. The annulus/CIE and shell/Green endpoints have positive smallest face eigenvalues (6.69e-9 and 3.59e-7), correctly signed active multipliers, and small projected gradients. This is numerical evidence consistent with constrained local minima, not an interval certificate. The bridge/Green endpoint has 24 negative face eigenvalues and eight wrong-signed lower-bound multipliers; it must not be classified as equally converged. The residual-weighted nonlinear Hessian and the exact smoothed-TV curvature are retained, rather than substituting a Gauss-Newton or IRLS matrix.

### 5.5 A bounded-current control does not ensure material recovery

For a passive outgoing finite Green operator and strictly positive real contrast, write Z=D_m^{-1}-G. With the fixed auxiliary shift c=-i beta and complete K=G, the current derivative is beta(Z-i beta I)^{-1}Z. Dissipativity gives its Euclidean norm an upper bound beta without assuming normality or uniform material. This is an application of established complex-shift/Cayley structure, not a new general contraction theorem; complex-beta CIE itself has prior art [5]. A finite Galerkin K introduces the separately bounded G-K remainder. Near physical resonance the current-error gain may be small, so an upper bound alone does not guarantee useful material excitation.

We fixed the phase and retained the original beta magnitudes, priors, physical corrector and initialization. Four exposed objects were run with K19 and K0, with all eight predictions sealed before offline scoring. Absorptive Green recovered only the shell; K0 recovered only the bridge. Both lost the annulus benefits, yielding one material recovery out of four each. All runs completed without reaching the CPU cap. This control rejects the fixed absorptive shift as a default repair; no phase scan was used to manufacture a positive result.

A separate zero-forward diagnostic of ten saved states did not support a universal loss of material excitation in the shell. The exact material-Jacobian column energy factors into the inverse-denominator column norm and effective-illumination row energy. Relative to ordinary CIE, their product has median ratios 1.0602/1.0560 for the shell's early/late states, rather than widespread attenuation. These diagonal energies do not resolve column correlation or full physical observability.

### 5.6 Completed three-dimensional development inversions

Twelve genuine vector-Maxwell inversions were completed on two exposed designs with six methods. For v04, Global and full Green had amplitude NRMSE 0.02120 and 0.02038, but IoU 0.8551 and 0.8169. For v09 all methods retained a blurred, attenuated high-contrast core, with relative contrast L2 error about 0.49--0.50 despite noise-compatible predictions. Reporting amplitude NRMSE alone would hide this sparse-volume normalization effect. All methods exhausted the 320-accepted-step corrector allowance; these are finite-budget outcomes, not convergence certificates. They show actual integration, but no clear three-dimensional material advantage.

### 5.7 Completed forty-eight-object comparison

The primary cohort contains eight geometry families, high/mixed material regimes, and high contrasts 8, 12 and 18. Every combination has a separately generated material, rather than reusing one geometry at three contrasts. All six algorithms were frozen before generation. All 288 public predictions, 144 conditional-deployment choices and 384 common-actual-CPU pair selections were sealed before truth access. The secondary selection rules were registered before this cohort was scored, but after the separate sixteen-object analysis, and are labelled accordingly.

|Method|Material criterion /48|CN99 /48|Mean amplitude NRMSE|Mean valid child CPU, s|
|---|---:|---:|---:|---:|
|CIE|21|42|.159082|159.681|
|Dual-CIE|22|44|.138267|155.749|
|Global-current control|23|44|.138219|158.381|
|Material-only G19|23|43|.146752|207.786|
|Frozen complete denominator|16|31|.225933|163.026|
|Frozen nonlocal feedback T|23|45|.135053|196.396|

The global-current control, G19 and contracted T succeed on exactly the same 23 objects, including the same five of six annuli. Hence the additional successes over sourcewise CIE cannot be attributed to the denominator. A new high-contrast annulus fails for every method, so the annulus label alone cannot predict the earlier benefit.

The strongest pure-high-contrast subgroup, chi=18, has zero material successes out of eight for every method. Its mean amplitude errors are .23639, .29022 and .21926 for Dual-CIE, G19 and contracted T, respectively. The mixed chi=18 subgroup gives four successes out of eight for each. Different contrast levels use independently generated geometries, so these descriptive strata are not a controlled contrast-only causal experiment. Among all 24 mixed objects, the median low-phase errors are 1.3210, 1.3474 and 1.3331 in absolute contrast units, against a low-phase truth of 2.5. Small globally normalized error differences therefore do not establish recovery of low-phase detail.

There are nevertheless useful continuous improvements. On a high-contrast bridge, contracted T reduces G19's error from .36949 to .08307, close to Dual-CIE's .08388. IoU improves from .1959 to .6762, still below the frozen .7 criterion. On a comb, error falls from .48253 to .22211, again near ordinary baselines but without recovering the thin structure. These are partial repairs, not additional material successes or superiority to Dual-CIE. Complete-denominator freezing is substantially worse.

For each method pair, a public rule takes B as the smaller measured valid child runtime and selects only paid minimum-Phi states within B and the 480-second cap. At these common CPU prefixes, Dual-CIE versus G19 recovers 22 versus 21 materials, while Dual-CIE versus contracted T recovers 22 versus 20. Mean paired error differences, Dual minus comparator, are -.011780 (exploratory 95 percent object-bootstrap interval [-.029557,.002132]) and +.001285 ([-.003681,.007889]). All objects and continuous outcomes remain included; these are not exact equal-spend claims or multiplicity-adjusted universal orderings. Excluding the one engineering-interrupted object gives the same qualitative conclusion.

The three conditional deployments, CIE followed if necessary by G19, Dual or Global, each recover 23/48 after six discrepancy-triggered fallbacks. G19 provides no extra recovery over these paid controls. Valid scientific children total 49968.922 CPU seconds, while the complete campaign total remains partly unknown because of the interrupted attempts described below. All 48 primary objects and a 47-object cost sensitivity are reported separately.

### 5.8 Completed paired-grid study: a conditional recovery gain

Eight new materials cover annuli, slotted shells, bridges and multiple islands, each in high and mixed contrast regimes. The high contrast is12, with low contrast2.5 and host.1. Each material produces matched24 observations and fine48-to-coarse24 observations, sharing standardized noise and a common2-by-2 volume-averaged scoring target. The target is a specified material representation, not an exact scattering-equivalent homogenization. Four methods yield64 inversions but only eight independent materials. All576 public selections and128 common-CPU pair definitions were sealed before the separate truth authorization.

|Condition|CIE criterion /8|Dual /8|Global /8|Green /8|Mean NRMSE: Dual / Green|Valid CPU sum: Dual / Green, s|
|---|---:|---:|---:|---:|---|---|
|matched24|4|4|4|6|.226128 / .065615|1308.422 /1580.828|
|cross48-to24|4|4|4|6|.231978 / .066069|1194.203 /1608.313|

The complete-run additional successes are the pure-high annulus and slotted shell. The Green annulus errors are.07357/.07731 versus Dual-CIE's.93779/.99572; the shell errors are.07027/.07218 versus.48558/.47850. Both gains are absent from the ordinary-material Global control. All methods preserve four ordinary successes, and all final results miss the frozen IoU criterion on two multiple-island objects. The latter still show recognizable locations and shapes; their moderate boundary errors must not be equated with the controls' catastrophic annulus errors.

The common-B Dual comparison also gives four versus six successes in each grid condition, but not the same two incremental objects. In matched data they are the annulus and shell. In cross-grid data the shell has insufficient IoU at the shorter Dual endpoint budget, while an earlier mixed-island Green state meets the criterion and later loses it as physical Phi decreases. CIE- and Global-matched prefixes give seven Green successes, but their different B values are reported separately and do not replace the six-success primary endpoint. No truth-dependent budget or stopping rule is adopted.

The cached fine/coarse clean-data discrepancy is1.45--5.30 percent of the fine signal norm. This is a two-discrete-model discrepancy, not a rigorous continuous-PDE error bound. The conditional gains therefore survive this specific inverse-discretization change, without establishing arbitrary-grid or universal high-contrast superiority. A separate diagnostic search inadvertently printed one private generation-metadata field; the process deviation is retained. Online solvers and already frozen selections did not use it, no materials were redrawn, and the whole research workflow is not described as perfectly blinded.

### 5.9 Ongoing vector-Maxwell study and engineering accounting

The true vector-Maxwell study uses 1728 shared real material variables, 5184 complex vector-current components per source, sixty external dipoles, and a finite dyadic CM+RR operator. Five methods are scheduled on eight designs, with two separate cross-grid observations. These are real ongoing computations, not completed numerical evidence.

After 132 completed two-dimensional runs, a controller/Windows receipt failure interrupted four attempts on one object; three subsequent preparation attempts failed before any forward solve. Completed science was retained and the remaining work resumed with unchanged algorithms. The original interrupted CPU totals are partly unknown and are not replaced by the successful retry times. The completed 48-object results retain this object, with a separate 47-object cost sensitivity analysis and explicit lower-bound or missing end-to-end costs.

The paired-grid study above was frozen and generated before the primary48 truth was opened, and executed after the primary288 completed. Its results were scored only after separate complete public sealing; it does not inherit three-dimensional truth permission.

Before any fresh three-dimensional truth was read, complete object blocks were assigned to the existing CPU and CUDA backends to use available resources. Every object's five methods and its matched/cross-grid variants remain on the same host. Timing is reported by platform; raw Mac and Windows seconds are not pooled to assert an algorithmic speedup. Scientific precision, acquisition, iteration budgets, and final physical validation are unchanged.

## 6. Discussion

Material-only Green has meaningful conditional positive results. Its increase from four to seven material successes over CIE in the sixteen-object cohort comes at approximately27 percent additional actual CPU, with one success unique relative to Global. The paired-grid cohort adds an annulus and a shell missed by Global, and preserves the gains after discretization mismatch. The larger forty-eight-object cohort does not reproduce a Green-specific success increment. Universal victory is not required, but no cohort may be omitted to favor either a positive or negative claim. Partially freezing the feedback repairs some detrimental trajectories while remaining more expensive than Dual-CIE and unproven as a better common-cost solver in the primary48 comparison.

The proposed two-block repair is not supported as a replacement. Its development success helps isolate block interactions, but the independent annulus counterexample is decisive against a broad matching-block explanation. Likewise, freezing a material-dependent metric can make an inner problem simpler without improving the material implied by an incorrect auxiliary current.

The use of a long common physical corrector deserves separate scrutiny. When materially inaccurate candidates have lower data and prior terms than the truth, further objective decrease is not a sufficient scientific outcome. In a posthoc diagnostic, all 65 first-discrepancy-hit states had larger material NRMSE than their final counterparts. Using this fixed stopping rule yielded no success gains and nineteen regressions across method–object pairs. A separate ground-truth oracle audit found no saved physical-objective incumbent that met the frozen material success criterion and subsequently lost it at the final state. These observations concern retained incumbents only; they do not rule out other stopping rules or unsaved iterates, and no ground-truth oracle is used for deployment.

## 7. Limitations and provisional conclusion

The primary48 and sixteen-object confirmations use matched finite forward/inverse operators. The supplementary eight-object study adds a paired cross-grid condition, but remains synthetic rather than measured data and does not certify continuous-Maxwell accuracy. The current three-dimensional Green space is incomplete, and its prospective results must not be conflated with the two-dimensional two-block modification. The inherited real-data and large-electrical-size failures remain part of the evidence.

The present manuscript supports a mechanism study with explicit counterexamples and a conditional material-only Green benefit. It does not yet support a general faster solver or a broadly superior three-dimensional imaging algorithm. The full publication claim will be chosen after the remaining frozen computations, without changing criteria to rescue the proposed repair.

## References

1. Zhong and Xu (2019), *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*. https://doi.org/10.3390/jimaging5020027
2. van den Berg, van Broekhoven, and Abubakar (1999), *Extended contrast source inversion*. https://doi.org/10.1088/0266-5611/15/5/315
3. Budko and Samokhin (2005), *Spectrum of the Volume Integral Operator of Electromagnetic Scattering*. https://doi.org/10.1137/050630660
4. Costabel, Darrigrand, and Sakly (2012), *The essential spectrum of the volume integral operator in electromagnetic scattering by a homogeneous body*. https://doi.org/10.1016/j.crma.2012.01.017
5. *Wavelet-Based Subspace Regularization for Solving Highly Nonlinear Inverse Scattering Problems with Contraction Integral Equation* (2020). https://doi.org/10.3390/electronics9111760
