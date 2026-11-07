# Green residual metrics and block interactions in high-contrast Maxwell inversion

**Working manuscript — partial results, not a submission-ready claim.** The sixteen-object study and exposed mechanism experiments below are complete. The independently frozen forty-eight-object comparison and vector-Maxwell study are running. Their results are not presumed here.

## Abstract

A material-dependent Green denominator can substantially alter contrast-source inversion trajectories without adding measurement information or enlarging the Fourier state space. We examine why such a modification can recover a high-contrast annulus while degrading other objects. An exact derivative split separates constitutive correction from variation of the residual metric. Actual saved states show that strong cancellation cannot be attributed solely to the Green contribution: the ordinary CIE contraction term is larger. Two lagged-denominator controls preserve annulus recovery but do not repair the exposed shell and bridge failures. Applying the Green residual to both current and material updates repairs an exposed bridge, whereas changing the current alone does not. This interaction nevertheless fails to generalize: in sixteen independently generated high-contrast and mixed-contrast objects, ordinary CIE, Dual-CIE, the inherited material-only Green modification, and the two-block modification satisfy the fixed material criterion in 4, 5, 7, and 5 cases, respectively. Material-only Green requires approximately 1.27 times the measured CIE CPU cost. All methods retain the same observations, material degrees of freedom, prior, and full-Maxwell corrector. A postseal audit further finds materially inaccurate candidates whose predictions and the ground truth are both compatible with the realized noise ball; this does not establish statistical indistinguishability or exact nonuniqueness. These results support a conditional trajectory mechanism, not universal convexification or unconditional superiority.

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

The bridge repair requires an interaction of the two changed blocks in this experiment. Neither individual change explains the successful two-block trajectory. The annulus instead benefits primarily from the material modification. Every Green variant fails the exposed shell. These are trajectory interventions, not a unique causal identification of all intermediate nonlinear effects.

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

### 5.3 Why data compatibility is not a material guarantee

Every ground-truth prediction lies inside its realized CN99 observation ball. Sixty-five of eighty candidates also lie inside, including thirty-eight material failures across eight independent objects. All sixty-five have both a lower data term and a lower prior than the ground truth.

These facts do not prove exact nonuniqueness or Gaussian statistical indistinguishability. In particular, the standardized prediction gaps for the material-failing compatible pairs are 4.30–6.09, not nearly zero, and candidates were selected from the same noisy data. The audit establishes that the adopted discrepancy criterion and lower \(\Phi\) do not certify material recovery. The remaining fifteen candidates are outside the ball although the ground truth is inside; their failure cannot be excused by this noise-compatibility argument.

### 5.4 Ongoing, independently frozen studies

The forty-eight-object study compares six methods on independent high/mixed geometries and three high-contrast levels. Its complete scoring is pending. The true vector-Maxwell study uses 1728 shared real material variables, 5184 complex vector-current components per source, sixty external dipoles, and a finite dyadic CM+RR operator. Five methods are scheduled on eight designs, with two separate cross-grid observations. These are real ongoing computations, not completed numerical evidence.

## 6. Discussion

Material-only Green has a meaningful conditional positive result in the completed study. Its increase from four to seven material successes over CIE comes at approximately 27 percent additional actual CPU cost and includes both a substantial annulus recovery and threshold-sensitive cases. Universal victory is neither required nor demonstrated.

The proposed two-block repair is not supported as a replacement. Its development success helps isolate block interactions, but the independent annulus counterexample is decisive against a broad matching-block explanation. Likewise, freezing a material-dependent metric can make an inner problem simpler without improving the material implied by an incorrect auxiliary current.

The use of a long common physical corrector deserves separate scrutiny. When materially inaccurate candidates have lower data and prior terms than the truth, further objective decrease is not a sufficient scientific outcome. In a posthoc diagnostic, all 65 first-discrepancy-hit states had larger material NRMSE than their final counterparts. Using this fixed stopping rule yielded no success gains and nineteen regressions across method–object pairs. A separate ground-truth oracle audit found no saved physical-objective incumbent that met the frozen material success criterion and subsequently lost it at the final state. These observations concern retained incumbents only; they do not rule out other stopping rules or unsaved iterates, and no ground-truth oracle is used for deployment.

## 7. Limitations and provisional conclusion

The completed two-dimensional confirmation uses the same finite operator for generation and inversion. It is not a measured-data result or a continuous-Maxwell accuracy certificate. The current three-dimensional Green space is incomplete, and its prospective results must not be conflated with the two-dimensional two-block modification. The inherited real-data and large-electrical-size failures remain part of the evidence.

The present manuscript supports a mechanism study with explicit counterexamples and a conditional material-only Green benefit. It does not yet support a general faster solver or a broadly superior three-dimensional imaging algorithm. The full publication claim will be chosen after the remaining frozen computations, without changing criteria to rescue the proposed repair.

## References

1. Zhong and Xu (2019), *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*. https://doi.org/10.3390/jimaging5020027
2. van den Berg, van Broekhoven, and Abubakar (1999), *Extended contrast source inversion*. https://doi.org/10.1088/0266-5611/15/5/315
3. Budko and Samokhin (2005), *Spectrum of the Volume Integral Operator of Electromagnetic Scattering*. https://doi.org/10.1137/050630660
4. Costabel, Darrigrand, and Sakly (2012), *The essential spectrum of the volume integral operator in electromagnetic scattering by a homogeneous body*. https://doi.org/10.1016/j.crma.2012.01.017
