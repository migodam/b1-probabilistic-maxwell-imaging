# Metric Placement and Finite Recovery Dynamics in Maxwell Inverse Scattering

**Research manuscript draft.** The evidence supports the claims below, rather than a claim of universal superiority or global convergence. The companion Chinese report, contracts and raw release archives are authoritative for implementation details and qualifications.

## Abstract

Auxiliary field–current formulations of inverse scattering use residual transformations to improve numerical search. A transformation may preserve the physical zeros while changing which material solutions a finite alternating solver can access. We investigate this distinction in a retained Fourier/DPCS Maxwell solver through a two-by-two factorial experiment: the current and material blocks use either the original scalar transformation or a local Green one-hop transformation. Material dependence and its derivatives are retained. On 48 previously exposed objects, the asymmetric A+-type configuration recovers 35 objects under an inherited joint descriptive criterion, versus 31 for a consistent nonlocal configuration. A common material inner solver preserves a decisive asymmetric advantage on one difficult object. On a separately frozen set of 24 new materials and acquisition conditions, A+ recovers 12, CIE and Dual-CIE each recover 9, the consistent nonlocal method recovers 11, and field–current reflection recovers 14. Importantly, the consistent nonlocal method alone recovers a high-contrast object missed by all older comparators. A simple late consistency switch preserves the A+ recovery set but adds no recoveries. Exact-best-response root spectra do not predict the measured finite program: incomplete current solves, PR memory, readout guards and non-root metric-derivative forces must be distinguished. A small genuine vector-Maxwell analogue shows limited additional imaging benefit. Full states, counterexamples, paid candidate histories and computational costs are released to support a more faithful theory of material accessibility.

## 1. Introduction

At high dielectric contrast, recovering a material from finite electromagnetic measurements is difficult for several distinct reasons. A material can be outside the effective auxiliary-current representation; the finite optimizer can follow a poor trajectory; different materials can satisfy the noise tolerance; or the forward model can be inaccurate. A lower auxiliary objective does not distinguish these explanations. Neither does a lower physical data residual alone establish a more accurate material.

This work concerns the placement of a known local Green residual transformation inside a CIE-derived solver. The existing A-material/A+ construction retains an ordinary global current update while modifying the material readout. It has recovered difficult materials missed by ordinary CIE, but has also failed on objects recovered by Green feedback or field–current reflection. The question is whether its gains arise from a generally better nonlocal residual norm, or from the particular interaction of two blocks with different norms.

We do not identify these explanations by inspecting one final image or an idealized spectral radius. We implement both residual families in both blocks, retain the actual material derivatives, reproduce the original A+ updates, and intervene from identical complete states. Separate controls address the material inner solver, the off-diagonal Green terms, scalar-beta mismatch and the order in which metrics are applied. A new independently generated cohort tests whether the observations extend beyond the exposed cases. The resulting evidence supports selective, state-dependent accessibility rather than a universal prescription to enforce or remove metric consistency.

## 2. Related work and scope of novelty

NIE and CIE already exploit electromagnetic integral-equation structure, auxiliary currents and normalized constitutive residuals [1,2]. CSI and CC-CSI connect data fitting and constitutive fitting, including normalization and iterative-search issues [3]. Relaxed PDE-constrained inversion and variable projection provide a broader context for auxiliary-state elimination and Gram weighting [4]. Convergence results for alternating or inexact descent require assumptions that do not automatically hold when two blocks optimize different potentials [5]. Nonlinear conjugate-gradient methods also have a finite algorithmic state; a one-step response is not synonymous with an exact quadratic minimizer.

The local Green transformation and the Möbius/Arc readout are inherited project results. We do not claim their first invention here, nor do we claim that weighted block Gauss–Seidel contraction, the chain rule, or reflection methods are new mathematics. The contribution is a controlled Maxwell study connecting actual metric functions, finite inner programs, complete-state interventions and recovery on new objects. The most useful negative result is that a favorable ideal root statement can coexist with substantially different finite recovery behavior. The positive numerical result is that both asymmetric and consistent nonlocal configurations possess material recoveries absent from relevant controls.

## 3. Forward model and residual transformations

For the finite scalar-TM model, write

\[
Q=E^{inc}+GJ,\qquad e=J-D_mQ,\qquad Y\simeq SJ.
\]

Here m is real susceptibility contrast and J is the normalized contrast-polarization current, not the impressed electric current. Time dependence is exp(−iωt); the permeability is fixed. The two residual transformations are

\[
P_0(m)=\operatorname{diag}\frac{b_0}{1+b_0m},\quad
P_1(m)=(I+D_aG_o)D_p,
\]
\[
d_i=1+(b_1-G_{ii})m_i,\quad a_i=m_i/d_i,\quad p_i=b_1/d_i.
\]

G_o contains exactly the inherited disjoint four-cell off-diagonal Green blocks. The metric is M_j=P_j*P_j. The diagonal self interaction remains in d. This construction is neither a full dense Green inverse nor a higher-order scattering expansion. Actual whole-box bounds ensure invertibility of the residual transformations for the frozen finite kernels; these bounds do not certify conditioning of the physical resolvent or accuracy of a continuous PDE discretization.

The current block minimizes its data loss plus a transformed constitutive loss. The material block minimizes its transformed constitutive loss plus the unchanged prior. We use 00, 01, 10 and 11 for the two choices. The same index means the same material-dependent function, including beta: in two dimensions b0=1 and b1 follows 1/.2/.2/.2. Consequently 00 is not a relabeling of the original Dual-CIE, whose two scalar beta values differ. Both derivatives of P1 and P0 with respect to m are retained.

The physical selection objective remains Φ=D+R, with D=.5||F(m)−Y||F²/||Y||F² and the original smooth/TV prior. This physical normalization is separate from the source-normalized auxiliary objective. Every selected material is evaluated through the original complete forward model. The primary selector uses the minimum eligible, completed, paid Φ; no truth participates in selection.

## 4. Ideal and finite block dynamics

### 4.1 Exact-best-response root model

At an exact root, in realified independent current coordinates, let K=L_*U, Ch=D_hQ_* and D=WSU. Define Hc=KᵀMcK+DᵀD and Bm=CᵀMmC. The ideal material loop is

\[
T=B_m^{-1}C^TM_mK H_c^{-1}K^TM_cC.
\]

When Mc=Mm, a positive-semidefinite block Hessian implies that T is self-adjoint in the corresponding material Gram and has spectrum in [0,1]. This is a standard weighted block-Gauss–Seidel property. It requires an exact physical root, no prior, a fixed regular active face, independent coordinates and exact block minimization. Inconsistent metrics need not possess the same structure.

### 4.2 Inexact-current and material errors

For actual local updates define rc=Hc δa−KᵀMcCh and rm=Bm δm−CᵀMmKδa. Then

\[
\delta m=Th+B_m^{-1}C^TM_mKH_c^{-1}r_c+B_m^{-1}r_m.
\]

This accounting identifies how an incomplete current solve is transmitted to material space. It is not a new recovery bound. One exact line-search steepest-descent step already depends nonlinearly on its direction through its step length. Finite PR adds history and stage transitions. Redundant native Fourier coordinates must not be substituted for an invertible independent-current Gram without qualification.

A declared 8×8 Maxwell anchor gives an ideal asymmetric eigenvalue 1.002712449. Material-only perturbations followed for 32 actual updates contract at the tested scale, whereas the ideal eigenvalue power is 1.09055. This alone does not falsify asymptotic expansion: an incomplete-current model can have a contracting material transient before an expanding full-state mode dominates. We therefore also initialize the compatible material/current direction. At the declared ±.001 scale the actual A+ readout performs zero material iterations in all 128 observed updates; its projected gradient is below the inherited guard. This records a finite-scale deadband, not asymptotic stability.

### 4.3 A non-root term absent from the ideal root model

For fixed actual J and Q,

\[
\nabla f_m=g_{const}+g_{metric}+\nabla R,
\]
\[
g_{const}=-\rho\Re C^*P^*(r\Omega_s),\quad
\langle g_{metric},h\rangle=\rho\Re\langle r\Omega_s,(DP[h])e\rangle_F.
\]

The second term is quadratic in the constitutive defect under a bounded derivative of P. It disappears in the ideal first-order root model but cannot be dropped at a non-root state. For scalar P0 its cellwise contribution is −ρb³||e_i||Ω²/(1+bm_i)³. This familiar normalization effect is not new measurement information. In 56 fixed first-step interventions, the metric/constitutive gradient-norm ratio ranges from .1003 to 2.0673. Constitutive force alone has positive work on the actual material step in 38 cases; the full gradient has no positive-work cases. These observations qualify interpretations based on a frozen-metric material projection.

## 5. Experimental design

We perform 576 full inverse runs: 192 runs on the 48 exposed materials, 48 self/beta controls, 32 stage-order controls, 16 common-material-solver controls, 264 runs on 24 new materials, and 24 vector-Maxwell runs. Runs are not independent material samples. The new two-dimensional cohort uses new seeds within eight inherited geometry families, three contrast regimes, varying illumination/noise conditions and four cross-grid cases. It is new-object confirmation, not an unseen-family claim.

All two-dimensional material spaces retain 576 real unknowns. No true support, font, shape label or test material is used online. The new cohort is scored only after all 264 predictions and public selectors are sealed. Complex noise uses a declared complex RMS convention; cross-grid model error is not included in the raw CN99 noise ball.

All new-cohort methods share a maximum 120 process-CPU-second allocation. Auxiliary proposals have an absolute 72-second deadline; unused time can be used by the common physical corrector. Direct physical L-BFGS receives the entire 120 seconds. Whole-child OS CPU is measured separately. A documented initial contract discrepancy concerned redistribution of unused proposal time; older strict-48-corrector results remain separate historical references. Paid-prefix and matched-physics-pool views do not refund computation or masquerade as stopped reruns.

Windows execution uses eight P-core and four E-core lanes, one BLAS thread per child, with all methods for an object executed sequentially on the same core. The small vector and common-inner-solver studies use a separate Apple CPU environment. We do not infer a cross-host speedup. No new CUDA backend or neural training is introduced.

## 6. Results

### 6.1 Exposed-object mechanism evidence

The 48-object four-pair experiment gives 29/35/30/31 joint descriptive recoveries for 00/01/10/11. A+ has four gains and no losses relative to 11: objects0017,0024,0026,0045. This does not prove that inconsistency as an abstract property is beneficial. Within the fixed material column it shows that imposing the nonlocal current metric does not automatically preserve material gains.

On the 16-object self/beta control subset, full01 recovers eight objects, self-only01 and global-Dual recover six, and full11/self11 recover five. The off-diagonal Green terms are decisive for0017 and0045. Object0024 is also recovered by global-Dual, showing why beta and ordinary readout effects must remain controls. With a common material L-BFGS60 program,0017 is still recovered by01 and missed by the other three pairs. The gain is therefore not exclusively an Arc inner-program effect.

Applying01 in the first three original stages and11 only at full Fourier capacity preserves the eight01 recoveries. The reverse order retains the five11 recoveries. This implicates early finite trajectories, while leaving a state-independent rule for choosing the better trajectory unresolved.

### 6.2 New-object recovery and counterexamples

| Method | Joint criterion /24 | Mean amplitude NRMSE | Mean IoU |
|---|---:|---:|---:|
| CIE | 9 | .243742 | .647744 |
| Dual-CIE | 9 | .237757 | .649437 |
| global-Dual | 10 | .221394 | .669424 |
| A+ /01 | 12 | .188710 | .697226 |
| consistent nonlocal11 | 11 | .205952 | .680138 |
| 01→11 | 12 | .188162 | .697643 |
| G19 | 11 | .197078 | .645968 |
| DR | 14 | .151093 | .715065 |
| direct L-BFGS | 2 | .862482 | .233693 |

A+ adds three recoveries and loses none relative to either CIE reference under the inherited joint descriptive criterion. Its paired mean NRMSE difference from CIE is −.0550, with a descriptive, unadjusted object-bootstrap interval [−.1086,−.0040]. This is not a population or novelty certificate. Low-contrast cases are mostly comparable; the stronger differences occur in high and mixed contrast.

Separate two-arm paid-pool comparisons retain12 versus9 for A+ against both CIE references at common process-CPU endpoints and at common componentwise physical-LU/RHS endpoints. The corresponding CPU-matched mean NRMSE differences are−.05570 and−.04937. These are retrospective budget views with truth-free selection, not measured stopped reruns or cost refunds. Auxiliary actions and DR background projection costs remain charged separately. An all-eleven-method common endpoint is additionally reported; an early-terminating direct method makes that a substantially different budget question.

The consistent nonlocal method alone recovers new high-contrast patchwork0021 (NRMSE .087, IoU .925, residual/δ .928), while A+, CIE, Dual, global-Dual, G19, DR, direct and the late switch miss it. Conversely01 recovers0012 and0015 while11 misses them by large material margins. Thus neither consistency nor inconsistency is a universal accessibility rule. The late switch has exactly the same joint recovery set as A+; its small mean error change is not evidence of added recovery capability.

DR provides strong complementary recovery on0009 and0020. Some membership differences are threshold-sensitive rather than basin changes: on0002, A+ and DR both have NRMSE near.053 although their support IoUs lie on opposite sides of.70. All continuous metrics and complete galleries accompany the binary descriptive table.

The pre-frozen public01→G19→DR saved-result replay recovers14/24, spending3846.796875 child-CPU seconds, compared with2797.609375 seconds for DR alone. It is not a120-second method or a GT-selected union. It can preserve a good A+ candidate and rescue others, but equal recovery counts do not imply equal cases or equal cost.

### 6.3 Corrector, difficult failures and vector Maxwell

Raw/corrected results show that physical correction often improves recovery, but not invariably. Across the exposed four-pair cohort, raw joint recoveries19/20/22/22 become29/35/30/31. This aggregate benefit does not justify discarding raw candidates. On the common-inner-solver subset, material error increases after correction in14/16 comparisons, with the size and structural relevance of each change retained.

The five inherited severe failures remain unresolved. Their true materials have lower physical Φ in the matched model, so this is not simply a prior that globally prefers the observed endpoints. Nevertheless19/20 tested four-pair endpoints admit a small declared physical descent. We do not label them certified strict local minima or prove uniqueness of the true inverse solution.

The vector analogue uses genuine dyadic Green coupling, CM+RR self terms,216 material variables,648 current components and60 external dipoles. Four public conditions represent three geometries, one repeated across grids. The main material errors are largely comparable across metric placements. Some support differences are real array changes, but ordinary CIE also performs well on the bracket case; no general new three-dimensional recovery advantage is established. The inherited Mie/discretization limits and the frame cross-grid mismatch remain explicit. No new experimental Fresnel inversion is claimed.

## 7. Discussion and limitations

The evidence rejects a simple prescription that a nonlocal material norm should also be used in the current block. It also rejects the opposite universal statement that a shared nonlocal norm cannot add recovery capability. Small invertible transformations can change finite paths substantially in a nonlinear, representation-limited problem. Our root experiment explains why ideal block spectra, although valid under their assumptions, are insufficient diagnostics of the actual PR/readout program. Non-root metric derivatives provide a second omitted mechanism.

The present controls do not isolate an abstract benefit of every nonpotential vector field. They concern two concrete metric families, specific finite programs, a fixed prior and declared data distributions. The new24 sample is small, subgroup analyses are descriptive, and initializations are uniform. Existing NN warm-start evidence is inherited rather than rerun. The vector study is a finite-discretization analogue with limited independent physical accuracy, and no new public measured-data result is included. These limitations prevent a claim of a broadly superior or globally convergent Maxwell algorithm.

The most useful next theoretical question is a finite-time, state-dependent accessibility statement that accounts for incomplete current solves, the nonlinear material force and representation transitions. The released complete states and counterexamples allow that statement to be falsified. A rule proposed after inspecting0021 would require a separate future confirmation cohort; we do not retrofit it into the present blind result.

## 8. Conclusion

Local Green feedback and the placement of residual metrics change material recovery beyond a scalar reparameterization. The asymmetric A+ configuration has reproducible two-dimensional gains relative to the tested CIE implementations, while a consistent nonlocal configuration unlocks a different independently generated material. Neither ideal local expansion nor metric consistency alone predicts this behavior. A late consistency switch preserves but does not enlarge the measured A+ recovery set. The resulting contribution is a state-resolved mechanism study and a new qualified nonlocal-current solver variant, with full costs and negative evidence; stronger general theory and external physical validation remain open.

## References

[1] Zhong et al., NIE, IEEE TAP2016, DOI10.1109/TAP.2016.2535492.

[2] Zhong and Xu, CIE-I, Journal of Imaging2019, DOI10.3390/jimaging5020027.

[3] Sun et al., CC-CSI, IEEE TAP2017, DOI10.1109/TAP.2017.2673758; public manuscript arXiv1906.10864.

[4] van Leeuwen and Yang, An analysis of constraint-relaxation in PDE-based inverse problems, Inverse Problems2025, DOI10.1088/1361-6420/adab89; preprint arXiv2403.15292.

[5] Attouch, Bolte and Svaiter, Mathematical Programming2013, DOI10.1007/s10107-011-0484-9.

Exact source verification, claim scope and retrieval failures are recorded in `PRIOR_ART_AND_NOVELTY.md` and `audit/prior_art/`. Historical project formulas are inherited and cited as such; no claim of exhaustive prior-art exclusion is made.
