# Conditional Möbius Updates and Finite-Action Rational Search in Electromagnetic Inverse Scattering: Algorithms, Costs, and Limits

**Research draft, 2026-11-04. Submission status: not submitted.** This manuscript reports only executed experiments. Broader algorithmic superiority and continuum-qualified three-dimensional imaging are not assumed. Numerical tables are generated from the final sealed-score artifacts described below.

## Abstract

We study two uses of rational material dependence in electromagnetic inverse scattering. The first preserves a contrast-source current update and replaces the material readout by conditional Möbius proposals under the original regularizer. A one-hop constitutive-residual metric is multi-affine in single-site Möbius coordinates, but the regularized joint problem remains nonconvex. We derive a four-color update whose candidate set consists of the unregularized state's quadratic stationary points, box endpoints, and the current material. Evaluating these candidates with a tangent majorizer of the unchanged smoothed-TV prior gives fixed-current descent without solving a regularized quintic. The second uses exact background-conditioned finite material responses and precontracted port Grams to evaluate finite, nonlocal actions without repeated complete forward solves. We distinguish evaluating an action accurately from selecting an action capable of leaving a poor reconstruction. Experiments include a thirty-object synthetic confirmation benchmark, a separately reported twenty-object selected-family confirmation, historical regression and initialization controls, genuine vector-DDA material inversions, and two measured Fresnel targets with withheld-channel prediction. Analytic structure reduces specific material-readout work, while the evidence places substantial limits on universal recovery, finite-action basin escape, and physical transfer. No learning model or test-truth-based online selection is used.

## 1. Introduction

Auxiliary-current formulations can follow material paths that differ from direct physical optimization. An inherited one-hop material readout (A-material) recovers two disclosed difficult curved objects where established CIE variants fail. That result does not establish robustness: previous independent materials, grid refinement and other geometries remove or reverse the advantage. Moreover, its material inner loop adds substantial work. The relevant questions are whether its useful behavior can be retained cheaply, whether it improves new objects at matched actual computation, and whether exact finite rational actions provide a genuinely different search capability.

Our contribution is a bounded algorithm-and-mechanism study. We implement a prior-preserving conditional proposal method; compare it with original A and strong established controls; and test finite-action search at actual poor full-wave anchors. Exact Schur/Woodbury identities, fractional-linear geometry, MM and least-squares contraction are established tools. The value of their application here is assessed by complete physical inversion, counterexamples, cost, and transfer, rather than by restating identities or increasing Taylor order.

## 2. Related work

NIE and CIE introduce established nonlinear integral-equation representations and FFT-TSOM implementations [1,2]. A method that estimates currents before material readout must also be distinguished from LUCIE [3] and classical CSI. Eigenpermittivity Green expansions already study material-spectrum poles [4]. DCTMC uses nonlocal intermediate interaction potentials while imposing data compatibility and locality [5]. Rational approximation and Lippmann–Schwinger–Lanczos ROMs compress parameter-dependent transfer information or estimate internal fields [6–8]. LNO learns a Laplace-domain operator representation [9]; its temporal poles are not automatically material-direction poles. Concavity-based majorization and MM descent are established [10]. Our evidence ledger states which comparisons use full text, abstracts, metadata or inherited repository evidence; it does not certify exhaustive novelty.

## 3. Physical model and objective

For two-dimensional scalar TM fields with exp(−iωt), real contrast m, self-inclusive outgoing disk-pulse Green G and external sources/receivers,

J=diag(m)(Einc+GJ), F(m)=SJ, Phi(m)=||F(m)−Y||F²/(2||Y||F²)+R(m).

The full finite model is solved with one retained complex LU for all illuminations. R is the shared normalized edge-smoothness plus smoothed-TV penalty. All inverse candidates obey the same box and are selected by original Phi. Source-wise data/state metrics in auxiliary CIE are distinct from the final observation normalization and the declared complex-noise discrepancy.

A-material keeps the ordinary global-current PR update and changes only its material block. Thus A versus Global-Dual is the cleaner material-update ablation; A versus sourcewise Dual-CIE also changes the current PR geometry. With b=β−diag(G), z=m/(1+bm), H=Q+bJ, and actual zero-diagonal Goff inside fixed disjoint cell blocks,

RA=β(I+diag(z)Goff)(J−diag(z)H).

This changes off-shell material updates. It does not replace physical Maxwell propagation or add measurement information.

## 4. Conditional material algorithm

Because Goff has zero diagonal, RA=A_i+z_iB_i when other coordinates are fixed. For real material, z lies on a legal Möbius circle arc. The conditional state loss is a quadratic-over-quadratic rational function of m. Its stationary numerator is quadratic, whereas adding a quadratic prior majorizer generally gives a quintic.

The proposed A+ avoids the quintic. For each color, it evaluates legal quadratic stationary points, box endpoints and the old coordinate using the exact conditional state loss plus the original prior-MM surrogate. A complete fixed-current loss check retains the old/seed value if numerical error causes degradation. Two sweeps and the inherited material-gradient tolerance are fixed in development. Color independence follows from disjoint four-cell propagation blocks and nearest-neighbor regularization. No material class/support reduction is used.

**Proposition (conditional descent).** Under independent same-color coordinates, a prior majorizer Qprior that touches the original prior at the old point, pole-free legal coordinates, and a candidate set containing the old point, selecting the smallest complete surrogate value yields nonincrease of the original fixed-current material-block loss. The proof is f(new)≤M(new|old)≤M(old|old)=f(old), where M=fstate+Qprior. This does not prove physical-Phi descent of the entire alternating method, retention of a reconstruction basin, or convergence to the true material. Floating-point fallback and observed loss histories are reported separately.

The full derivation, actual prior normalization, exact dP term, numerical restrictions and pseudocode-equivalent source are in THEORY_DERIVATIONS and src/aplus_readout.py. A projected-gradient readout is an important competing simplification: a lower block loss alone cannot justify it as an equivalent replacement.

## 5. Exact finite-action rational search

At physical anchor X0, let A0=I−X0G, O=SA0⁻¹, C=GA0⁻¹ and E0=Einc+GJ0. For ΔX=PDPᵀ, the exact finite response is

ΔF=OP(I−DKP)⁻¹DEP, KP=PᵀCP.

Two scalar cells yield a bilinear numerator divided by a bilinear determinant. Four port arrays give an exact small Gram representation of the original data loss. For general q, let T=(I−DK)⁻¹D and precompute Aq=OPᴴOP/||Y||², Bq=EP EPᴴ, Zq=OPᴴ(F0−Y)EPᴴ/||Y||². The data loss is Phi_data(0)+Re tr(ZqᴴT)+.5tr(TᴴAqTBq). This contraction eliminates repeated port loops at an unchanged anchor; setup, selected columns, Gram construction and all small solves remain paid.

We test distant residual-curvature cell pairs, adjacent and random controls, two four-cell patches (rank eight), and one structural refinement to rank-thirty-two actions. The latter uses weak-Hessian/translation/gradient support and signs with box-slack amplitudes; it is not an exact Hessian-eigenvector trajectory. Each selected endpoint is verified with full Maxwell and the common corrector. Continuing the untouched anchor and local exact search in the same action family are retained as controls. A free-face Hessian does not certify a box-constrained basin. No global search theorem is claimed.

## 6. Experiments

### 6.1 Frozen synthetic comparison

Thirty independently generated materials span six analytic families and five prespecified acquisition/contrast/grid conditions. Twenty-four are same-grid finite-model data; six are generated on a36-grid and inverted on a24-grid with the same disk-pulse backend. This breaks same-discretization inverse crime in that stratum but is not an independent continuum solver. All inversions retain576 real material coordinates. Conditions jointly vary contrast, complex noise, source count, aperture and generation grid; their descriptive differences do not isolate single-factor effects.

CIE, Dual-CIE, Global-Dual, original A and frozen A+ use3000 inherited auxiliary updates and a common3000-step original physical corrector. Direct physical L-BFGS-B has a declared6000-step cap. Equal nominal iterations are not equal work: final and paired overlapping CPU/LU/RHS frontiers are both reported using paid physical-Phi incumbents. Snapshots every ten observation events and float32 offline storage impose an explicit matching resolution. Bootstrap resamples material objects, not illuminations or iterations.

|方法|材料成功/30|同网格材料+CN99 /24|中位NRMSE|总CPU秒|LU|复RHS|
|---|---:|---:|---:|---:|---:|---:|
|成熟CIE|21/30|15/24|0.08722|1543.60|54873|1154196|
|成熟Dual-CIE|20/30|14/24|0.08584|1539.15|53408|1137060|
|Global-Dual|20/30|14/24|0.08653|1512.97|53923|1137684|
|原A-material|21/30|14/24|0.08177|2495.14|54384|1148820|
|A+|21/30|14/24|0.08233|1711.34|52156|1103040|
|原物理L-BFGS-B|10/30|9/24|0.31489|838.18|47754|1009764|

Material NRMSE, IoU, physical residual and Phi are reported separately. Same-grid CN99 discrepancy has a declared probabilistic interpretation. Cross-grid raw noise thresholds omit model mismatch and cannot certify total error. Public development choices used only disclosed objects; all main predictions were sealed before independent truth scoring. Repeated initialization and historical nine-digit regression do not increase the thirty-object denominator.

### 6.2 Selected-family secondary confirmation

After the primary benchmark, we freeze the same algorithms and draw twenty new objects: ten curved bars at contrast12/six sources and ten slit shells at contrast18/half aperture. A single-CIE control is added before target-truth access; all80 predictions seal before scoring. This is conditional, post-selection family confirmation, not an expanded pre-registered50-object population benchmark.

A+ succeeds on12/20, versus9/20 each for Dual/global and8/20 for single CIE; data-qualified counts are9/8/8/6. Curved-bar success counts are9/8/8/5; shell counts are3/1/1/3, with identical A+/single shell success sets. A+ mean-NRMSE minus Dual is−.06831 (paired20-object95% bootstrap interval[−.14213,−.01242]), while the success-rate difference interval[0,.30] touches zero. CPU totals are1089.69/1061.57/1027.61/855.19s. Under equal CPU versus single, material counts are10 versus8 and mean-error CI includes zero. One curved bar is recovered only by A+ among four controls, but seven shells still fail. All eight A+ failures have a better feasible offline objective witness; two report numerical convergence and six hit their allowance. No target is removed or used for algorithm reselection.

### 6.3 Readout mechanism development

Twenty complete disclosed runs compare original A, quintic MM, projected updates, inherited-tolerance MM and the final quadratic-candidate A+. Projected8 loses the old b20d outcome: Phi rises from about6.47e−6 to.005832. Final A+ retains the public-objective behavior of both old curved bars and the failed/easy controls. Four-object CPU falls from284.21s to197.78s, while one case remains4.8% slower. Actual quintic roots drop to zero. All60,000 retained material readouts across development show no retained block-loss increases under the recorded numerical check. This is not a universal speed or material-recovery theorem.

### 6.4 Finite-search anchors

The disclosed slit shell has no substantive improvement in tested rank2/8/32 families; its smallest computed free-face Hessian eigenvalue is about4.34e−9, a diagnostic only. Thin ring and42cb have budget/slow-convergence explanations: continued physical correction reaches approximately.00492946 and.00218249, comparable to or better than finite-action branches. Accurate small responses therefore do not establish new basin access. Twelve rank2/rank8/local/anchor-control endpoints and three separately sealed rank32 endpoints all fail the material-image criterion; these fifteen endpoints concern only three disclosed anchors. All setup, rejected candidates and continuation costs are retained.

### 6.5 Genuine vector-DDA material inversion

The vector benchmark has its own quadratic-plus-l2 prior, data normalization and material-relative-L2 NRMSE; equality of priors/objectives is within each benchmark, not between2D and3D. We perform actual125- and216-coordinate isotropic material inversions, with375/648-component electromagnetic systems, six transverse plane waves and external vector receivers. Vector-adapted ordinary, local-resolvent, A-material and direct controls are explicitly distinguished from a mature published3D CIE implementation. Longer common physical correction largely removes initial differences. On125 coordinates, final A NRMSE is about.17844 for the sphere, .19965 for two spheres and.13300 for a bracket, comparable to other adequately continued methods. No durable A recovery advantage appears.

Independent actual-grid Mie checks expose geometry aliasing: binary sphere error increases from6.58% to27.59% under the tested refinement. One fixed volume-fraction correction gives8.40%→8.20%. These finite-DDA inversions do not establish continuum-qualified3D superiority or high-contrast vector robustness.

### 6.6 Measured Fresnel data

Four original dielectric TM files are freshly obtained from the canonical IOP supplement. Two targets at2GHz are inverted with Dual/global/A using actual missing channels, explicit exp(+iωt)-to-Hankel1 conjugation, target-absent incident calibration and held-out scattered-field prediction. The horn interior field remains approximated by a cylindrical source. Single/twin training residuals are.14193/.10062 and withheld residuals.26755/.23444; all three methods yield identical saved material endpoints per target. Both hit the contrast cap; the twin reaches its corrector limit. No pixel GT errors or calibrated absolute permittivity recovery are asserted. The first common-corrector objective equals the host objective in all six measured runs, consistent with rejecting auxiliary endpoints. Because individual auxiliary candidate arrays were not retained, the exact auxiliary argmin cannot be certified from these records. A adds cost without measured benefit in this pilot. Raw data/fulltext redistribution is excluded; official retrieval scripts and derived results are supplied.

## 7. Discussion

The useful structure is conditional rational material dependence and its interaction with an off-shell constitutive metric. The final readout is simpler than exact regularized quintic minimization and can preserve a disclosed useful path. However, neither fixed-current descent nor exact finite propagation supplies missing information or a useful global material action by itself. Direct physical continuation is essential as a control: otherwise extra solve budget can be misreported as a basin mechanism.

Learning is unnecessary for the exact small evaluator and analytic coordinate proposal tested here. A learned proposal generator might choose better coherent actions, but this study does not test that hypothesis or claim that learning is generally useless. There is no independent train/test learning dataset, training cost or neural superiority result.

## 8. Limitations

Six material families and thirty new realizations do not establish population-wide robustness. Acquisition factors are bundled. Most synthetic observations have inverse crime; different-grid cases share the same backend. Equal-budget comparisons use stored checkpoints rather than continuous-time synchronized stopping. CPU is measured on a shared local machine; full operation counts accompany timings. The measured pilot lacks a fully calibrated horn model and exact material truth. Coarse3D geometry error remains substantial and vector controls are adapted rather than a mature3D CIE reproduction. NIE/LUCIE, LSL field lifting, Fresnel3D and learned proposals are not independently reproduced.

## 9. Conclusion

Conditional Möbius state proposals with unchanged prior-MM acceptance provide a concrete lower-cost A-material implementation on the tested synthetic benchmark. The observed success set is retained at approximately31.4% lower CPU than original A, while average material error versus Global-Dual improves with a limited, heterogeneous-sample interpretation. A separate selected-family confirmation repeats limited curved-bar gains, while overall success-rate superiority, measured-data benefit and continuum-qualified vector3D advantage remain unestablished. Exact rational finite-action evaluation is useful machinery but did not provide a distinct basin escape in the tested anchors. A paper should state these bounded algorithmic and mechanistic conclusions rather than claim a universal rational-learning inversion framework.

## References

1. Zhong et al., NIE, IEEE TAP2016, DOI10.1109/TAP.2016.2535492.
2. Zhong and Xu, CIE-I, J.Imaging2019, DOI10.3390/jimaging5020027.
3. Bevacqua et al., LUCIE, IEEE LAWP2026, DOI10.1109/LAWP.2026.3660123 (primary abstract evidence).
4. Chen et al., GENOME, Phys.Rev.Applied2019, DOI10.1103/PhysRevApplied.11.044018.
5. Levinson and Markel, DCTMC Part I, Phys.Rev.E2016, DOI10.1103/PhysRevE.94.043317.
6. Druskin et al., LSL2021, DOI10.1088/1361-6420/abfca4.
7. Baker et al., regularized LSL2025, DOI10.1016/j.jcp.2025.113725.
8. Gosea and Güttel, matrix-valued rational approximation2021, DOI10.1137/20M1324727.
9. Cao et al., LNO2024, DOI10.1038/s42256-024-00844-4.
10. Hunter and Lange, MM tutorial2004, DOI10.1198/0003130042836.
11. Belkebir and Saillard, Fresnel2001, DOI10.1088/0266-5611/17/6/301.

Full source identity, evidence tiers, failed retrievals and bounded novelty interpretation: PRIOR_ART_AND_NOVELTY.md. This draft is supported by the delivered experiments, not by a claim of journal acceptance.
