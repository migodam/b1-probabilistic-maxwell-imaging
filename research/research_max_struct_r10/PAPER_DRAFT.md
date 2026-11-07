# Projected Green Matrix Denominators for Contrast-Source Inversion: Mechanisms, Costs, and Limits

Evidence-complete mechanism manuscript, MAX-STRUCT R10, 7 November 2026. All numerical statements refer to sealed object-level tables. This is a research draft, not an accepted or peer-reviewed paper. It distinguishes finite-model algorithm evidence from continuous-Maxwell validation and does not claim general convexification of inverse scattering.

## Abstract

Contrast-source inversion separates auxiliary-current updates from material readout. We investigate whether the material residual should transport spatial coupling through a matrix denominator rather than a local, one-hop correction. With the current fixed, the proposed residual is `R_K = beta (I + beta D_m - D_m K)^(-1) [J - D_m(Einc + GJ)]`, where `K = PGP` is a low-rank compression of the public background Green operator. All material degrees of freedom and the original Maxwell forward model remain unchanged. A Woodbury implementation and its exact material derivative retain the denominator dependence on the material. For the exact two-dimensional disk-pulse Green matrix used here, a radiative-Gram argument establishes nonsingularity of the projected denominator for nonnegative contrast and beta. The result is qualitative: it does not bound floating-point conditioning or establish inverse uniqueness. We separate current-representation error, auxiliary/physical-current mismatch, material recovery, and actual cost. In a frozen24-material2D study, the matrix denominator recovers one high-contrast annulus missed by CIE but worsens other paths; mean NRMSE is0.17669 versus0.16287 for CIE, with19% more child CPU. A separately frozen four-material discrepancy-triggered fallback improves recovery from2/4 to3/4 at1.46 times CIE CPU, whereas a similarly priced Dual-CIE fallback remains2/4. True vector3D and cross-grid experiments identify budget, finite-noise ambiguity and model mismatch, without establishing an A+ recovery advantage. The evidence supports conditional complementarity and a limited mechanism claim, not general dominance.

## 1. Introduction

A small auxiliary residual need not imply an accurate material. This distinction is particularly consequential for high-contrast inverse scattering: an auxiliary current may fit the measured channels while remaining inconsistent with the material's true full-wave current. Increasing the current representation or reducing its auxiliary objective can therefore improve a surrogate without reaching a better material solution.

The present study asks a narrower question than global convexification. Holding the mature CIE/FFT-TSOM current machinery fixed, can a physically motivated material-residual transport change the search trajectory usefully, and is its cost justified? We compare the inherited local A+ readout with a low-rank, nonlocal matrix denominator. The low rank is a computational device for the denominator, not a reduction in material dimension or additional measurement information.

The evidence has three distinct parts: exposed historical mechanism controls; a frozen confirmation cohort of 24 independent two-dimensional materials and seven methods; and four three-dimensional vector materials with two additional cross-grid acquisition conditions. The three-dimensional study evaluates the inherited isotropic A+ method, not a purported unimplemented vector extension of the new matrix denominator.

## 2. Related work and scope of novelty

CIE and its FFT-TSOM implementation already exploit contrast-source reformulations, Fourier representations, and distinct current/material contraction parameters [1,2]. Twofold subspace optimization already uses Green-operator spectral structure in current-space continuation [3]. Material preconditioning, contrast updates, and regularized contrast-source inversion also have longstanding precedents [4]. Woodbury inversion, Möbius reparameterization, radiative positivity, and Rellich uniqueness are established tools [5,6].

Our candidate contribution is the specific material-readout residual, its implementation within an otherwise matched current solver, and the evidence identifying its useful and harmful regimes. We do not claim invention of these underlying identities. The available publisher abstract for [3] was verified, but its full formula-level comparison remains incomplete because the full text was unavailable through the accessible sources. The novelty claim consequently remains narrower than a first-of-kind assertion.

## 3. Methods

We use the exp(-i omega t) convention. In the finite two-dimensional scalar-TM model, `Q = Einc + GJ`, `J = D_m Q`, and the measured scattered field is `Y = SJ + noise`. Sources and receivers lie outside the unknown region. Positive real contrast is optimized over the same box and with the same smooth/TV prior for every method. The online selection objective is the original full-wave `Phi = ||F(m)-Y||_F^2/(2||Y||_F^2) + R(m)`, not the auxiliary objective and not a sigma-whitened chi-square statistic.

The inherited A+ material update uses `a_i = 1 + (beta-G_ii)m_i`, `z_i=m_i/a_i`, `u=beta e/a`, and `R_A=(I+D_z G_off,block)u`, where `e=J-D_mQ`. Its off-diagonal action is confined to disjoint 2 by 2 spatial blocks. The Möbius/MM solver reduces the cost of this material subproblem; it does not insert the entire Green operator into its scalar denominator.

We instead define `A_K=I+beta D_m-D_mK` and `R_K=beta A_K^(-1)e`. A fixed orthonormal basis U is computed only from the public operators. The primary candidate uses the strong right singular subspace of G (rank19 at the frozen threshold0.1). A second candidate adds the receiver-based twofold construction (rank30 in confirmation; rank28 on the historical operator). Both use `K=U(U*GU)U*`. An equal-rank Green control and a self-term-retaining variant are restricted to exposed development experiments.

The main confirmation retains all576 real material unknowns. The current updates, actual projected Fourier spaces, acquisition, prior, initialization, and common original-physics corrector are inherited unchanged. The two matrix-denominator candidates use at most eight bounded L-BFGS-B material iterations. Old material, ordinary-CIE seed, and the optimized candidate are compared by the same new auxiliary objective. The final trajectory output is selected only by paid original-Phi evaluations.

## 4. Theory and mechanism

Write `K=BU*`, `a=1+beta m`, `z=m/a`, `F=D_z B`, and `T=I-U*F`. Then

`A_K^(-1)b = D_a^(-1)b + F T^(-1) U* D_a^(-1)b`.

For fixed J and Q, `w=A_K^(-1)e` obeys

`dw[h] = -A_K^(-1)D_h[Q+(beta I-K)w]`.

For a fixed Hermitian source metric Omega, the real material gradient of `rho beta^2 Re<w,w Omega>/2 + R(m)` is obtained with one adjoint small solve: `v=A_K^(-*) (w Omega)` and `g_i=-rho beta^2 Re sum_s conj(v_is)[Q+(beta I-K)w]_is + partial_i R`. Freezing the material-dependent denominator would omit a required derivative. The present implementation does not do so.

If both inverses exist, the exact discrepancy to the full shifted denominator is

`R_K-R_G = -A_K^(-1)D_m(G-K)R_G`.

Thus a small static Green tail alone does not guarantee a small material-update error. Contrast, omitted coupling, and the residual resolvent all matter. Since material multiplication and projection generally do not commute, replacing complex Green action by scalar singular values is not an equivalent material formula.

For H=U*GU and M=U*D_zU, the small denominator is T=I-MH. If Im_H H is strictly positive, T cannot have a nonzero null vector: `x=MHx`, with `y=Hx`, implies `x*Hx=y*My` is real, contradicting strict dissipation. The exact disk-pulse model satisfies `Im_H G_ij = c J0(k|x_i-x_j|)`, `c=(pi k a/2)J1(ka)>0`. The Bessel integral makes this a full-angular radiation Gram. A zero quadratic form would yield a finite point-source field with identically zero far field; an outgoing Fourier-Hankel expansion and unique continuation force every logarithmic singularity coefficient to vanish. Hence the exact Gram is strictly positive for distinct centers, and any full-column compression inherits this property. The actual ka is below1, where the stated sign of J1 follows directly from its alternating series.

This is an exact-model nonsingularity result, not a useful uniform condition-number bound or a floating-point interval certificate. It does not transfer automatically to complex contrast, negative materials, arbitrary self terms, PML, or the vector DDA. Nor does full-angular radiative positivity make finite receiver data injective in the material.

The new residual does not guarantee convexity. Already `f(m)=m^2/(1+m)^2` has `f''(m)=2(1-2m)/(1+m)^4<0` for m>1/2. The unchanged original full-wave objective retains its own landscape. The proposed mechanism changes the lifted auxiliary search path; it does not prove removal of original-objective local minima.

## 5. Experiments

### 5.1 Frozen acquisition and comparison contracts

The two-dimensional confirmation uses eight independently sampled geometry families—annulus, slotted shell, islands, curved bars, comb, coated inclusion, perforated star, and bridge—with one independently generated low, high, and mixed object per family. These are24 materials, not168 independent samples. The24 by24 inverse grid has576 real unknowns,12 exterior sources,48 exterior receivers, and circular complex Gaussian noise with expected relative amplitude1%. Contrast levels are2.5 and12/18 with a0.1 host; all methods optimize the full box[0.02,20]. The same finite discretization generates and inverts the data.

Seven methods are evaluated per object: CIE-single, Dual-CIE, Global-Dual, A+global, A+sourcewise, Green19, and Two-fold30. Global-Dual is an internal matched-current ablation, not a separate literature method. All methods allow the same Fourier schedule and original-physics correction, without a wall-time cutoff. We report fixed-schedule endpoints and saved paid-cost prefixes separately. Full child-process CPU and wall are measured externally; nested component timers are not added again.

The vector study uses12^3 real material cells,5184 complex field components per illumination,60 exterior electric dipoles, and180 receiver components in a0.6-wavelength region. It contains an epsilon9 cube, an epsilon9/epsilon2 dimer, an epsilon5 hollow cube, and an epsilon2-shell/epsilon9-core cube. The solid and coated cubes additionally use independent18^3 generation and12^3 inversion. Each condition compares sourcewise CIE, Global-Dual, and A+global. The CM+RR dyadic DDA, three field components, same-LU adjoint, and padded nonperiodic FFT are retained. The common corrector permits80 accepted steps; reaching that limit is not a stationarity certificate. Fine/coarse differences and the inherited Mie limitations remain separate from algorithm performance.

### 5.2 Exposed mechanism controls

On historical b20d/40ac, the inherited A+ final amplitude-NRMSE is0.0406/0.6260. Green19 obtains0.0409/0.0548, and Two-fold28 obtains0.0400/0.0557. Retaining every original diagonal Green term also succeeds, limiting an explanation based on an incorrectly implemented or simply discarded self term. A control that only changes the old one-hop material solver to the same eight-step L-BFGS-B gives0.3846/0.6235; merely switching the inner optimizer is insufficient on these two exposed cases. The failed first dispatch of this control and its costs are retained separately.

The improvement on40ac includes a substantial contribution from the common physical corrector: matrix-denominator endpoints before correction have NRMSE about0.15–0.16. An earlier full-current low-penalty reference in R9 had also recovered this object. We therefore claim neither the first recovery of40ac nor dominance over all optimization routes.

At24 saved states, the local one-hop feedback spectral radius is below0.168, whereas the projected global feedback exceeds1 at18 later states. The small denominator nevertheless has minimum singular value0.312–0.688. Failure of a Neumann series is not the same as a near-singular direct inverse. Green19 is closer to the full shifted residual than one-hop at every checked state; this is auxiliary-residual evidence, not a material recovery theorem.

### 5.3 Independent confirmation results

All168 primary predictions were sealed before scoring. There were no execution failures. The predeclared material descriptor combines amplitude-NRMSE at most0.15, common-support IoU at least0.7, and legal contrast; continuous errors remain the main evidence.

| Method | Material successes/24 | Mean NRMSE | Median NRMSE | Mean full-child CPU (s) |
|---|---:|---:|---:|---:|
| CIE |11|0.16287|0.13801|157.71|
| Dual-CIE |10|0.16346|0.13829|153.98|
| A+global |8|0.18044|0.14531|169.49|
| A+sourcewise |10|0.18045|0.14479|174.98|
| Green19 |8|0.17669|0.14638|187.66|
| Two-fold30 |8|0.18099|0.14544|192.53|
| Global-Dual, internal ablation |10|0.16369|0.13894|150.13|

The matrix denominator is not an overall improvement. On the high-contrast stratum, mean NRMSE is0.21316 for CIE,0.26514 for A+global, and0.25363 for Green19. Low/mixed continuous differences are small; some success-count changes reflect threshold crossings near0.15 rather than large image differences.

Three high-contrast objects are more informative than the aggregate. On annulus0001, CIE/A+global have NRMSE about0.423 and relative residual13.7%, while Green19 reaches0.09889 and0.960% residual. On slotted shell0004, CIE reaches0.11109 and0.968% residual, whereas all A-family readouts reach approximately0.515 and58.98% residual. On bridge0022, CIE/A+global reach approximately0.097/0.096; Green19 and Two-fold30 instead reach0.337/0.437 with substantial residuals. This is conditional path complementarity, not universal dominance.

All24 objects, including all failures, appear in the three regime-specific image plates. Object-paired intervals and paid CPU/LU prefixes are retained in the source tables.140 correctors stop by projected-gradient tolerance,4 by relative objective change, and24 reach the3000-step cap. A termination flag is not itself a certificate of global or even exact local optimality.

CIE produces22 noise-compatible endpoints but only11 material successes. The true material and the estimated material both lie inside the actual observation's per-case CN99 noise ball in all11 data-compatible/material-failed cases, and both obey the box. These are finite-noise pairs, not exact multiple solutions or an information-theoretic lower bound. The common soft prior may assign them different preferences.

A posthoc two-candidate portfolio selecting CIE or Green19 by the lower original Phi reduces mean NRMSE to0.14960, but yields only10/24 material successes and costs2.19 times the CIE CPU. It gains the annulus but loses two threshold-near CIE successes. The unselected candidate remains fully charged. A separate, newly frozen conditional fallback study tests whether invoking the candidate only after CIE fails the public discrepancy criterion can retain ordinary-case behavior at lower cost; it is not retroactively treated as a predeclared rule on these24 objects.

### 5.4 A single posthoc structural correction

The strong right singular space is an input space. The two-sided approximation drops the additional output `(I-P)GP`, whereas `GP` retains the complete Green response to that input space. Its discarded output amplitude is16.83% in the measured Frobenius norm, or2.83% of the corresponding squared norm. We tested the already implemented one-sided form on the three exposed positive/negative objects, changing no other solver setting, rank, prior, or budget.

One-sided GP preserves the annulus result (NRMSE0.10167), leaves the shell failure unchanged (0.51514), and worsens the bridge to0.43608 versus0.33700 for PGP. The bridge Phi decreases from0.01421 to0.01040; GP reaches the3000-step cap while PGP previously stopped by relative objective change. All701.64 child-CPU seconds are charged. Thus a better standard low-rank approximation of G is not sufficient to improve material recovery. These three exposed diagnostics do not replace the original confirmation or supply an additional independent cohort.

### 5.5 Vector and cross-grid results

All18 primary predictions were sealed before the vector truth was opened. The following errors are relative L2 errors in contrast, not the amplitude-normalized2D metric. Every physical corrector reaches its80-accepted-step budget.

| Object/condition | CIE | Global-Dual | A+global | CIE/A+ relative data residual |
|---|---:|---:|---:|---|
| epsilon9 cube, matched12 |0.06082|0.05486|0.05758|2.979% /2.978%|
| epsilon2 shell/epsilon9 core, matched12 |0.49222|0.48836|0.49656|2.993% /2.993%|
| epsilon9/epsilon2 dimer, matched12 |0.39200|0.39929|0.38306|3.082% /3.069%|
| epsilon5 hollow cube, matched12 |0.31816|0.32140|0.32124|3.004% /3.005%|
| epsilon9 cube,18-to12 |0.30701|0.30372|0.30005|3.173% /3.173%|
| coated cube,18-to12 |0.50825|0.50199|0.50806|2.998% /2.999%|

The matched solid cube is reconstructed accurately by all three methods. The coated and hollow objects retain substantial contrast error despite near-noise data residuals. A+ differences are small and mixed; they do not establish a distinct vector recovery advantage. These are four independent materials, not18 independent samples. The relative clean-data discretization gaps are24.80% for the solid cube and6.41% for the coated cube, exceeding the3% measurement-noise scale. Algorithm comparisons on the matched finite model therefore do not certify continuous-Maxwell accuracy.

The inherited vector selection objective is Phi3=mean_s[||F_s-Y_s||^2/||Y_s||^2], without one-half and with R=0. It is common to all three methods, but is neither the2D global-normalized objective nor the ML objective for the declared iid scalar-sigma noise. True-sigma whitened discrepancies are reported separately. This implementation-contract limitation is retained rather than silently changing the frozen runs.

Saved fields distinguish noise ambiguity from model mismatch. For the matched coated cube, the true and CIE materials have whitened residual squared10806.71 and10754.59, both below the per-case CN99 bound11043.23, despite49.2% contrast error. For the two cross-grid cases, the unchanged coarse-reference material instead gives749142 and60516, far outside that same bound. Returning the coarse-reference truth and fitting the fine-generated observation within the stated measurement noise are therefore incompatible requirements at this particular reference. This does not prove that the entire coarse material class lacks a fitting solution, nor is it an exact nonuniqueness theorem. No additional forward solve was used for this check.

A single additional, explicitly exposed budget diagnostic continued the dimer's three public80-step endpoints with the same objective and reset-history rule. It was selected for its large public projected gradients (about1.3e-5 to1.9e-5) and residual above the noise level, not by selecting a favorable truth score. After80 additional accepted steps each, contrast errors decrease from0.39200/0.39929/0.38306 to0.35572/0.36594/0.35225 for CIE/Global/A+. All three enter the actual CN99 ball, while all again reach the iteration budget with projected gradients5.4e-6 to9.1e-6. Thus budget matters, but the A+-versus-CIE difference shrinks from0.00894 to0.00347. The extra243 vector LU and29160 complex RHS are charged separately; the original18 endpoints remain unchanged. Historical probes and abandoned slow deployments are charged and are not extra independent materials.

### 5.6 Independently frozen conditional fallback

After analyzing the24-object study, we froze a new four-material cohort with independent geometries: annulus, slotted shell, bridge and coated inclusion. Each deployment first runs CIE. Only a whitened residual norm above the predeclared CN99 radius25.16542 triggers one cold-start fallback; the endpoint is selected by the lower original Phi, with ties retaining CIE. Green19 and ordinary Dual-CIE are evaluated as separate fallback deployments using the same trigger, candidate count, prior and per-run budgets. No truth enters the trigger or selection.

Two objects trigger. On the annulus, Green19 reduces NRMSE from1.05000 to0.10958, improves IoU from0.04899 to0.71098, and reduces the whitened residual norm from135.70 to23.10. Dual-CIE does not recover it. The nontriggering shell and bridge retain exactly the successful CIE predictions (NRMSE0.12797 and0.09071). The triggering coated object remains a shared failure, with NRMSE about0.46018 and whitened residual about702.

CIE, CIE-to-Dual and CIE-to-Green19 produce2/4,2/4 and3/4 material successes and mean NRMSE0.43222,0.43222 and0.19711. Their deployment child-CPU totals are692.84,1015.55 and1013.44 seconds: the new fallback costs1.46 times standalone CIE and approximately the same as the Dual fallback control. Physical576-dimensional LU counts are8054,11560 and10225; small19-dimensional LU are additional. These deployment views overlap the CIE prefix and must not be added. All eight actual scientific inversions and the0.125-second zero-physics startup failure are retained.

This is a positive independent example of conditional path complementarity, not a statistically established population gain: only four objects and two triggers were observed. By construction the rule preserves CIE when it already meets the noise criterion, but it cannot repair a noise-compatible yet inaccurate material. The new cohort is not pooled retroactively with the original24 or with exposed structural diagnostics. No Global-Dual fallback was run on these four objects, so this deployment comparison alone does not uniquely attribute its gain to the denominator rather than the entire alternative path. The original24-object annulus does have a matched-global-current control, whose NRMSE remains0.42256.

## 6. Discussion

Representation and optimization must be distinguished. On exposed b20d and40ac, the actual FFT/APCS affine current floors at the final Fourier capacity are about0.232 and0.205 in the declared squared norm, despite very different A+ material outcomes. The current representation is identical across the matched methods. This rules out an explanation based solely on granting A+ extra Fourier modes; it does not establish a general relationship between a static Green projection and nonlinear material observability.

Auxiliary and physical currents also differ. For a shared material, `e_aux=e_phys+g` yields `||e_aux||^2=||e_phys||^2+||g||^2+2 Re<e_phys,g>`. One historical state has200.610=781.243+624.985-1205.618. Omitting the signed cross term would conceal cancellation and misidentify auxiliary fitting as physical recovery. Physical-corrector outputs are not constrained to the auxiliary affine current family and may legitimately beat its oracle floor.

The scientific comparison is regime-specific, not an all-cases dominance requirement. A reproducible improvement on a defined hard class can be valuable, provided ordinary-case regressions, failed objects, extra cost, and the matched baselines are also reported. Positive historical controls alone cannot supply that evidence.

## 7. Limitations

The24-object two-dimensional study is a finite-model algorithm confirmation with a fixed geometry family design, not a population-wide or real-measurement guarantee. Bootstrap intervals are descriptive object-level summaries and are not corrected for multiple exploratory comparisons. The vector study has only four independent materials and uses a coarse CM+RR discretization; its cross-grid cases do not certify continuous-Maxwell accuracy. The new matrix denominator is not implemented in3D, and no new Fresnel or neural initialization result is claimed in R10. Inherited R1–R9 positive and negative evidence is not overwritten.

The public-operator dense SVD used to form the low-rank basis is O(N^3) and is charged. The small readout solve costs approximately O(Nr^2+r^3+Nrs+r^2s) per value/gradient, in addition to the unchanged current/prior work. A scalable vector extension would require its own compatible basis construction and physical qualification. Qualitative dissipativity does not remove numerical-conditioning risks. Access limitations prevented a complete full-text comparison to every closest prior method.

## 8. Conclusion

The primary confirmation supports a precise distinction between local residual transport, projected nonlocal feedback, current representation, and material recovery. It contains a genuine independent high-contrast recovery gained by the matrix denominator and substantial opposing failures. The method therefore cannot yet be recommended as a general replacement for mature CIE. A publication claim must be restricted to the identified mechanism and conditional recovery evidence, with the independently frozen fallback and vector results reported separately.

## References and reproducibility

1. Zhong and Xu, *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*, Journal of Imaging5(2),27,2019. https://doi.org/10.3390/jimaging5020027
2. Xu, Zhang and Wei, *Fourier Bases-Expansion Contraction Integral Equation for Inversion Highly Nonlinear Inverse Scattering Problem*,2020. https://doi.org/10.1109/TMTT.2020.2978165
3. Zhong and Chen, *Twofold subspace-based optimization method for solving inverse scattering problems*,2009. https://doi.org/10.1088/0266-5611/25/8/085003
4. van den Berg, van Broekhoven and Abubakar, *Extended contrast source inversion*,1999. https://doi.org/10.1088/0266-5611/15/5/315
5. NIST DLMF10.9.1, Bessel integral representation. https://dlmf.nist.gov/10.9.E1
6. Colton, *Inverse Acoustic and Electromagnetic Scattering Theory*, MSRI Publications47,2003; the cited Rellich theorem is stated in3D, while the text above supplies the2D outgoing Fourier-Hankel argument. https://library.slmath.org/books/Book47/files/colton.pdf

The exact contracts, source, unselected per-object outcomes, failures, and accounting are provided in BENCHMARKS.md, METHODS_AND_ABLATIONS.md, THEORY_DERIVATIONS.md, REPRODUCE.md, COST_LEDGER.md and the complete same-repository evidence release. Model and result scope are part of the reproducibility contract.
