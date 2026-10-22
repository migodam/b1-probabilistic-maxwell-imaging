# Research W source-package evidence map

Date: 2026-10-22. Role: mechanical source/claim extraction and focused bibliographic verification. No new physical computation, no theory/novelty verdict, no hash checks. Paths below are relative to research-w-p0. The source report's scientific judgments are attributed, not adopted.

## Coverage and deliverable status

Read all 14 `source_package/report_parts/*.tex` files, the PDF through `pdftotext -layout`, `README_ZH.txt`, `P0_P1_PROTOCOL_ZH.txt`, `SOURCE_MANIFEST.json`, `REPRODUCIBILITY.json`, `build_report.sh`, font README, and saved analytic results. Extracted PDF is `docs/SOURCE_REPORT_EXTRACTED.txt`. PDF table of contents and source inputs agree on 12 numbered sections, mathematical appendix, and bibliography. Physical PDF has 52 pages (metadata); report's printed body numbering ends at 51 because title is unnumbered. This extraction is not a rendered visual audit or a recompilation.

The delivered package reports **executed exact analytic 1D Maxwell verification only**. It explicitly says no 2D/3D inverse imaging, full posterior reference, network training, P0-B, P1, or speed benchmark. The P0/P1 source describes a **proposed protocol**, requiring freeze before fresh data/seeds; it is not by itself a complete, frozen executable specification.

## File / section / claim map

| Source | PDF printed pages | Content and claim boundary |
|---|---:|---|
| `main.tex` | title, 1–3 contents | Article input order, formatting and declared evidence cutoff. |
| `overview.tex` (§1) | 4–5 | Prioritize propagation identifiability; distinguish known tools, conditional provable corollaries, analytic witness, untested mechanism, unresolved scope. Explicitly disclaims broad novelty. |
| `model.tex` (§2) | 6–8 | e^(-iωt); passive Im εc≥0; total material/nuisance model; real/imag stacking; wall/target label gauge; fixed query T with X=m restricted to T; level-set/Debye examples; joint prior and task types. |
| `identifiability.tex` (§3) | 9–14 | Set-theoretic fiber criterion; finite-dimensional kernel inclusion and KA† constant; infinite-dimensional instability counterexample; constant-rank extension; nuisance compensation; nonlinear task modulus; Gaussian two-point lower bound; exact planar slab/PEC witness and bandwidth failure; local localization error. |
| `probability.tex` (§4) | 15–19 | Four blocks d,A,B,C and exact finite-feedback response; bounded contrast and resonance-margin error bound; task pseudometric; same-joint-prior likelihood/TV bound; tails; wall/target conditional dependence; task Bayes risk; Laplace volume; repeated-likelihood prohibition. |
| `feedback.tex` (§5) | 20–23 | Resolvent/Woodbury; local-update four-block identity; A17 task I/O extension and same-space Galerkin identity; changed-wall GN terms; primal/dual residual bound; four error ledgers; total cost/strong controls. |
| `frontier_algorithm.tex` (§6–7) | 24–27 | Schur effective source/receiver/offset; uncertain Schur tuple; reversible empirical task frontier; candidate vs probability modes; prior-supported proposals; lazy exact MH / delayed acceptance / correct SMC stage targets; budget stop; acquisition separately identified. |
| `evidence.tex` (§8) | 28–32 | Five fixed repositories; original reading scope/missing drafts; A17 closure positive but slower results; A26/A30/P2 failures; internal A17-X/F overlap. Historical numbers are source-record summaries, not reruns. |
| `literature.tex` (§9) | 33–36 | Ten requested works and precise source read-depth; Marchenko and ambiguity prior art; Bayesian goal reduction; pathway constraints. No exhaustive novelty certification. |
| `experiments.tex` (§10–11) | 37–41 | P0-A accomplished; P0-B/P1 proposed; independent-grid caveats; finite 2D TM fixture; 2×2 inference ablation; all-cost comparison; realistic measurement/3D limitations. |
| `judgment.tex` (§12) | 42–44 | Source author recommends bounded I→II branch; six candidate contributions and stop/redirect conditions. Final scientific judgment remains parent Codex's responsibility. |
| `theory_addendum.tex` | 45–47 | Same-prior TV proof, branch+numerical cross term, conditional mixture, deterministic fixed-wall TV=1 caveat, dimensionless dangerous directions, analytic-statistics conventions. |
| `repo_bib.tex`, `literature_bib.tex` | 48–51 | 15 project-reference keys and 33 literature-reference keys, 48 total. |

## Assumptions that must accompany extracted formulas

1. **Physical forward model:** all inverse/resolvent uses assume a uniquely solvable field problem. Passivity alone is not a uniform cheap inverse-norm certificate. Shared material and nuisance across frequency/source are required; per-block free complex gains can absorb material information.
2. **Semantic decomposition:** arbitrary χW+χX relabeling is a parameter gauge, not a new Maxwell counterexample. Operational T is fixed and contains all local material, including wall if present. Joint target q must be a total-material functional unless an object prior is separately stated.
3. **Factorization:** H=h∘F is a set-level fiber statement only; continuity/measurability/stability do not follow automatically. Local derivative inclusion alone does not produce a finite nonlinear witness. Constant rank and legal interior parameters are substantive assumptions.
4. **Nuisance:** unpenalized local compensation permits a whitened orthogonal projection; finite prior penalty produces a Schur metric, not generally an orthogonal projector or full marginalization.
5. **Four blocks:** fixed discretization/reference material; V includes proper frequency/metric factors; physical complex operators precede real whitening. Strong-target bound requires ||V||≤v, ||A||≤a, ||B||≤b and ||(I−CV)^−1||≤β^−1 over the declared task class. PEC requires a separate boundary treatment, not v→∞.
6. **Probability compression:** same underlying joint prior and support; preserve X,γ conditionals within wall clusters; uniform likelihood error only controls stated region, with both posterior tail masses otherwise needed. Candidate weights/short chains do not establish posterior coverage. A fixed wall cannot approximate a non-atomic full material posterior in TV, although a task pushforward can be close.
7. **Error and cost:** numerical approximation, wrong fixed model, posterior uncertainty, and unmodeled structure are separate. Full-rank material changes can arise from one scalar latent. All setup/proposal/rejection/reference/check/fallback cost must be charged; equal current rank is not equal cost.
8. **Statistics:** data-distance uses two noiseless means with declared shared covariance. d≤1 is a chosen ambiguity gate (equal-prior Gaussian test error ≥0.309), not a 95% residual noise ball.

## Saved 1D evidence (not re-executed here)

`analytic_results/results.json` status is `executed_analytic_1d_Maxwell_check_only`. Walls: n=2,d=1/4 and n=3,d=1/6, common front face, μ=μ0, εr=4,9; λ0=1. Fixed comparison plane D=1. PEC planes at L1=1.5 and L2=1.4166666667. Nine frequencies per bandwidth, independent quadrature σ=0.01.

| Saved check | Value |
|---|---:|
| Center maximum reflection modulus | 1.6328623989e−16 |
| Same-position task response difference | 1.0000000000 |
| Round-trip phase difference | π/3 |
| Shifted-target joint center difference | 2.4492935983e−16 |
| Target shift | λ0/12 = 0.0833333 λ0 |
| Maximum lossless unitarity defect | 8.8820603167e−16 |
| Independent Fresnel discrepancy | 4.5182803599e−16 |

At 0.5% bandwidth saved CAL distance is 0.88713775 but shifted-target JOINT distance is 3.03824626; at 1% it is 6.05507116. This is a documented negative boundary: the selected broadband joint pair is distinguishable even while CAL remains near. Infinite planar normal-incidence geometry has no finite lateral support or aperture. Neither this pair nor its 1/12 wavelength target shift meets proposed 2D finite-wall/0.1 wavelength criteria. Different bandwidth configurations are nonnested nine-point grids.

## P0-B specifics available to freeze

These are **source proposals**, not a statement that the run is authorized or frozen.

- 2D TM Ez, z invariant, μ0, e^(-iωt), full multiple scattering; domain [−2,2]×[0,3] λ0, antennas y=−0.3, aperture 3 λ0, same-side 4 Tx×12 Rx.
- Wide frequencies f/f0={0.8,0.9,1,1.1,1.2}; second narrow configuration frozen before scientific output.
- Full-domain/level-set material representation, no true wall mask or hole count; Re εr∈[1,9], explicit small-loss interval; proposed inverse grid 32×48. Forward grid chosen by physical accuracy, independent refinement and PML extension; proposed whitened CAL/task consistency tolerance ≤0.1 each.
- Two fixed backgrounds of tilted variable-thickness dielectric with irregular air; task T fixed. CAL fixture makes T actual reference material; JOINT uses actual total material and nuisance. Distinct observation identities, no truth subtraction.
- Predeclared 12 task locations × two fixed dielectric contrasts with fixed probe scale/amplitude; full feedback, separate construction/validation probes.
- CAL candidate task separation dT≥3 while mean distance d≤1. JOINT additionally d≤1 and fixed total-material task threshold (example location ≥0.1 λ0 or declared occupation change).
- Local dangerous-direction+finite corrector versus Full-J/adjoint, random material direction+same corrector, multistart joint FWI. Same warm start, caches, prior, solver, cost rights. CAL first; passed pairs must face JOINT and independent refinement.
- Suggested whole-P0-B cap: 8,000 complex RHS; 2,000 material-frequency factorizations; 2 CPU hours, first limit binds. 40% search / 30% full endpoint+task validation / 20% independent grid / 10% failure recheck. A single published solver benchmark can adjust before scientific results only.
- With rigorous errors require dnum+εF1+εF2≤1 and dTnum−εH1−εH2≥3. Without PDE error certificates only finest declared discrete-model conclusion. Mesh/PML differences are consistency evidence.
- Budget exhaustion with no witness means search unresolved, not global uniqueness. Only CAL or ultra-narrow success requires scope redirect; coarse-grid/unphysical-gain-only success is mechanism NO-GO. Disclose robustness neighborhood; isolated tuned pair supports existence only.

### Parameters still requiring a concrete freeze before P0-B

Exact Tx/Rx x coordinates; precise narrow frequencies; absolute/source amplitude and calibration rule; noise covariance/σ and shared nuisance ranges; domain/PML/air padding and discretization; loss/dispersion interval; prior correlation length/occupancy/minimum feature; actual T and probe centers/radii/contrast values; two numeric seed identities/background generator; total-material task q/threshold; optimization/step acceptance/stopping; exact allocation across two acquisitions/seeds/controls; strict definition and ledger of RHS/factorization and CPU budget; numerical validation acceptance margins and robustness neighborhood. The prose examples cannot silently become an already frozen protocol.

## P1 deferred conditions

Only after actual JOINT ambiguity: low-dimensional independently integrable branch reference with explicit branch prior/evidence weights. Single vs multiple hypotheses crossed with joint vs progressive computation; same prior/data/proposal library. Four development backgrounds then frozen eight × two noise samples; aggregate by background. Proposed each-method cap 256 multifrequency/multisource full-forward equivalents plus wall time; reference cost disclosed. Reference event/CDF error≤0.02; same-accuracy total cost reduction≥20% or explicit same-cost task-loss improvement. Coverage remains pilot. Truth-selected Best-of-K, oracle wall and full-reference fallback remain diagnostics. No general unknown-topology posterior claim.

## Focused scholarly verification, 2026-10-22

Applied the `scholarqa-research` skill's Evidence QA/citation verification rules. A single `scholarqa-cli verify` batch (5 records) ran in approximately 0.67 seconds tool wall time, within the requested 60-second CPU scope. All five IDs resolved, unresolved_ids=[]; raw output: `SCHOLAR_VERIFY.json`. No broad search or citation-graph expansion. These records verify identity; abstract support below is Tier B, missing abstracts Tier C. Detailed full-text assertions are **the source report's prior-read claims**, not newly verified by this extraction.

| Boundary | Original source | What current evidence supports / does not |
|---|---|---|
| Closest through-wall application | [Watson 2024](https://doi.org/10.1088/1361-6420/ad5b83), [author manuscript](https://arxiv.org/abs/2403.10354) | Retrieved abstract: full-wave obstacle effects with approximately known wall and unknown nuisance fitted during imaging. Source §9.1 adds BEM/ROM/variable projection, scalar and target Born limitations, based on its stated full-text read. Abstract alone does not verify those detailed limitations. |
| Reflection/transmission ambiguity | [Waagaard–Skaar 2007](https://doi.org/10.1137/060672224), [arXiv](https://arxiv.org/abs/0707.3252) | Retrieved abstract explicitly studies 1D interacting multimode coupling from reflection, matrix factorization and ambiguities resolved by priors. Source §9.3 gives §2–3 mechanism. This is not itself a finite irregular Maxwell-wall theorem. |
| T-matrix and dynamic background | [Jakobsen–Tveit 2018](https://doi.org/10.1093/gji/ggy197) | Title/authors/year verified, abstract unavailable (publisher elision). Source report says publisher text search read. No new substantive scope verification here. Related original [2015 direct iterative T-matrix](https://doi.org/10.1088/1742-2132/12/3/400). |
| Marchenko full wavefield | [Diekmann–Vasconcelos 2021](https://doi.org/10.1103/PhysRevResearch.3.013206) | Retrieved short abstract supports a 3D full-wavefield Marchenko integral framework only; source claims first-arrival requirements and decompositional distinctions from its earlier inspection. |
| Partial data uniqueness | [Brown–Marletta–Reyes](https://doi.org/10.1016/j.jde.2016.01.002), [preprint](https://arxiv.org/abs/1410.8480) | Abstract supports local-boundary recovery, C1,1 domain and coefficients known near boundary. Does not equate continuous boundary response to finite antennas. Semantic Scholar year returned 2014; bibliography cites journal 2016. Resolve preprint/publication metadata before final bibliography. |

Task factorization is explicitly labeled a known basic fact in source §3.2; no dedicated mathematical original citation is supplied nearby. Do not claim a newly verified original provenance for it. Nearby goal-oriented works are [Lieberman–Willcox 2014](https://doi.org/10.1137/130928315), [Spantini et al. 2017](https://doi.org/10.1137/16M1082123), [Cui et al. 2016](https://doi.org/10.1016/j.jcp.2016.03.055), and [Zahm et al. 2022](https://doi.org/10.1090/mcom/3737); these four were extracted from source bibliography, not live verified in this limited batch.

Source Marchenko comparisons that parent should retain as nearest boundary: [2014 multicomponent single-sided](https://doi.org/10.1093/gji/ggu313); [2022 forward-scattered waveforms](https://doi.org/10.1109/TUFFC.2022.3163906); [2017 electromagnetic lossless representation](https://doi.org/10.1109/IWAGPR.2017.7996033); [2016 dissipative EM](https://doi.org/10.1109/ICGPR.2016.7572702); [2021 1D dissipative GPR](https://doi.org/10.2528/PIERM21020901). Source reports differing complete-reflection/focusing-field/transmission/two-sided requirements; the 2022 wave example is acoustic. No global novelty judgment is made here.

## Potential inconsistencies / parent-review flags

1. **Execution status vs freeze:** source gives suggested values, example thresholds and missing numbers; cannot call P0-B/P1 executed or fully frozen from this package. Metadata is a Linux run record, not this Mac's current runtime.
2. **Strong-target signature / whitening:** §4.1 warns whitening may be real-linear; §4.1 last paragraph also says whitening can be absorbed in R. Parent should maintain the earlier realified-operator qualification whenever presenting the bound.
3. **Gate definitions:** dT denotes both task metric and full target-response class; JOINT q is a separate total-material functional. Freeze these separately to avoid treating CAL separation as JOINT location risk. P0-A target shift < proposed 0.1λ0 threshold is no contradiction because it is not P0-B PASS.
4. **Grid checks:** ≤0.1 mesh differences do not certify PDE error; source expressly restricts no-certificate claims to declared discrete model. Parent should preserve margins in both d and dT gates.
5. **Counterexample scope:** no finite-support lateral structure/aperture and joint compensation breaks at modest bandwidth; any broad 'finite single-sided walls are nonidentifiable' conclusion needs separate math review.
6. **Measurability/stability/constant rank:** source expressly limits basic set factorization and local kernel checks; these conditions cannot be omitted when transforming them into scientific claims.
7. **Probability certificate:** uniform log likelihood error assumes original joint prior; data-dependent clustering can still require measurable representative selection. Tail unknown means conditional certificate, not safe global posterior. Source proposition is conservative and makes no sharp-constant claim.
8. **Internal overlap:** A17-X has finite local T-matrix, four port objects, bilateral feedback; F has primal/dual and likelihood interfaces. Source identifies these as reuse. Parent owns any assessment of independence.
9. **Bibliographic metadata:** five-source batch recovered Brown preprint year vs journal year mismatch; T-matrix abstract missing; Marchenko retrieved abstract too short to freshly support detailed requirements. Source's 48-source count is not 48 full-paper fresh verification.
10. **Build reproducibility:** font binaries intentionally omitted; `build_report.sh` needs two Noto Serif CJK OTFs plus XeLaTeX/DejaVu. No package versions lock file; stored versions and no-warning metadata are provenance, not a rebuild performed here.
11. **Reading claims:** §8 describes original report author's prior reading. This extraction checked those statements and their listed URLs, not all remote historical project files. L/Z/X ZIPs, P3-L, and several A-series originals explicitly missing in the source.

## Reproducibility metadata extracted

Report date 2026-10-22; declared 52 PDF pages; Python 3.12.14; Linux x86_64; NumPy 2.3.5; SciPy 1.17.0; Matplotlib 3.10.8; XeTeX 0.999995 / TeX Live 2023 Debian. Metadata reports zero Overfull, missing character, undefined, multiply defined or warning entries. Source list 48. These are stored records, not newly rerun measurements. Analytic reproduction command is `python analytic_witness.py`; PDF build makes three XeLaTeX passes. Neither was run by this worker.

## Fixed internal evidence links and complete reference-key index

The following tables are generated directly from SOURCE_MANIFEST, preserving fixed commits. Their presence is provenance extraction, not remote readback verification.

- **R_A17** (repo_bib.tex): [1](https://github.com/migodam/a17-maxwell-certified-exchange/blob/df61a1b533682a883cdaf65259b4528c354bd832/closure_final/reports/A17_CLOSURE_REPORT_ZH.md) ; [2](https://github.com/migodam/a17-maxwell-certified-exchange/blob/df61a1b533682a883cdaf65259b4528c354bd832/closure_final/START_HERE.md)

- **R_A17core** (repo_bib.tex): [1](https://github.com/migodam/a17-maxwell-certified-exchange/blob/df61a1b533682a883cdaf65259b4528c354bd832/closure_final/manuscript/CORE_THEORY.tex) ; [2](https://github.com/migodam/a17-maxwell-certified-exchange/blob/df61a1b533682a883cdaf65259b4528c354bd832/closure_final/manuscript/CONDITIONAL_EXCHANGE_THEORY.tex)

- **R_A17eff** (repo_bib.tex): [1](https://github.com/migodam/a17-maxwell-certified-exchange/blob/df61a1b533682a883cdaf65259b4528c354bd832/efficiency_followup/A17_EFFICIENT_VERIFICATION_FINAL.md)

- **R_A26** (repo_bib.tex): [1](https://github.com/migodam/a26-om-feedback-prior-completion/blob/2603f3bd16c48e48587fedc4363975779def17f5/A26_FINAL_DECISION.md)

- **R_A30** (repo_bib.tex): [1](https://github.com/migodam/a30-maxwell-observable-material-manifold/blob/717d3d26a290248bd2b9aeb55f3d03694b407b84/A30_E0_E1_DECISION.md)

- **R_B1** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/START_HERE.md) ; [2](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/joint_p2/GATE_DECISION.md)

- **R_C1** (repo_bib.tex): [1](https://github.com/migodam/c1-physics-state-generative-imaging/blob/83056881e89329212b16f5be47a4d0dbc4465943/p2/C1_P2_DECISION.md)

- **R_P2** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/p2_material_fiber_20261019/RESEARCH_REPORT_ZH.md)

- **R_F3** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_f3_quotient_exploration/RESEARCH_REPORT_ZH.md)

- **R_A17X** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/a17x_f3_competition_audit/pro_source/report/RESEARCH_A17X_FULL_REPORT_ZH.md) ; [2](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/a17x_f3_competition_audit/RESEARCH_REPORT_ZH.md)

- **R_F** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_f_intrinsic_transport/pro_source/report/RESEARCH_REPORT_ZH.md) ; [2](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_f_intrinsic_transport/RESEARCH_REPORT_ZH.md)

- **R_B** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_b/RESEARCH_REPORT_ZH.md)

- **R_A** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_a_20261015/RESEARCH_REPORT_ZH.md)

- **R_T** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_t_20261015/RESEARCH_REPORT_ZH.md)

- **R_G** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_g_20261015/RESEARCH_G_REPORT_ZH.md)

- **R_LXZ** (repo_bib.tex): [1](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/064d0aaed9772cadd573866eec84f159895b317b/research/research_lxz_20261016/RESEARCH_REPORT_ZH.md)

- **L_Somersalo** (literature_bib.tex): [1](https://doi.org/10.1088/0266-5611/10/2/017)

- **L_DBIM** (literature_bib.tex): [1](https://doi.org/10.1121/1.3458856)

- **L_Marks** (literature_bib.tex): [1](https://doi.org/10.1364/OPTICA.4.001514)

- **L_Tmatrix** (literature_bib.tex): [1](https://doi.org/10.1093/gji/ggy197)

- **L_Tmatrix15** (literature_bib.tex): [1](https://doi.org/10.1088/1742-2132/12/3/400)

- **L_Watson** (literature_bib.tex): [1](https://doi.org/10.1088/1361-6420/ad5b83) ; [2](https://arxiv.org/abs/2403.10354)

- **L_Lidar** (literature_bib.tex): [1](https://doi.org/10.1109/JSEN.2025.3576879)

- **L_Building** (literature_bib.tex): [1](https://doi.org/10.1109/JIOT.2025.3545739)

- **L_LayerUQ** (literature_bib.tex): [1](https://doi.org/10.1016/j.amc.2025.129453) ; [2](https://arxiv.org/abs/2504.05776)

- **L_Minsley** (literature_bib.tex): [1](https://doi.org/10.1111/j.1365-246X.2011.05165.x)

- **L_Topology** (literature_bib.tex): [1](https://doi.org/10.1137/17M1113850)

- **L_Marchenko14** (literature_bib.tex): [1](https://doi.org/10.1093/gji/ggu313)

- **L_Marchenko21** (literature_bib.tex): [1](https://doi.org/10.1103/PhysRevResearch.3.013206)

- **L_Marchenko22** (literature_bib.tex): [1](https://doi.org/10.1109/TUFFC.2022.3163906)

- **L_EM2017** (literature_bib.tex): [1](https://doi.org/10.1109/IWAGPR.2017.7996033)

- **L_EM2016** (literature_bib.tex): [1](https://doi.org/10.1109/ICGPR.2016.7572702)

- **L_EM2021** (literature_bib.tex): [1](https://doi.org/10.2528/PIERM21020901)

- **L_Dissipative18** (literature_bib.tex): [1](https://doi.org/10.1103/PhysRevApplied.10.044022)

- **L_Waagaard** (literature_bib.tex): [1](https://doi.org/10.1137/060672224) ; [2](https://arxiv.org/abs/0707.3252)

- **L_Golden** (literature_bib.tex): [1](https://www.math.utah.edu/~golden/docs/publications/GOLDEN_etal_inverse_IEEE_1998.pdf)

- **L_Chaderjian** (literature_bib.tex): [1](https://doi.org/10.1137/S0036139992241157)

- **L_NonlinearGoal** (literature_bib.tex): [1](https://doi.org/10.1137/130928315)

- **L_Goal** (literature_bib.tex): [1](https://doi.org/10.1137/16M1082123)

- **L_Cui** (literature_bib.tex): [1](https://doi.org/10.1016/j.jcp.2016.03.055) ; [2](https://arxiv.org/abs/1510.06053)

- **L_Zahm** (literature_bib.tex): [1](https://doi.org/10.1090/mcom/3737)

- **L_Cotter** (literature_bib.tex): [1](https://arxiv.org/abs/0909.2126)

- **L_Lie** (literature_bib.tex): [1](https://arxiv.org/abs/1712.05717)

- **L_BayesStrip** (literature_bib.tex): [1](https://doi.org/10.1088/1361-6420/ab6f9e)

- **L_GreenBayes** (literature_bib.tex): [1](https://doi.org/10.1029/2018JB016499)

- **L_Zhai** (literature_bib.tex): [1](https://doi.org/10.3934/ipi.2026060)

- **L_Partial** (literature_bib.tex): [1](https://doi.org/10.1016/j.jde.2016.01.002) ; [2](https://arxiv.org/abs/1410.8480)

- **L_SOM** (literature_bib.tex): [1](https://doi.org/10.2528/PIER10020903)

## Runtime record

Extraction/verification work completed 2026-10-22T08:28:16.611723+00:00. Five-record scholarly command tool wall time 0.67 seconds. Local reading/extraction calls completed in under 1 second each; artifact writing experienced a longer tool/session wait. No CPU-heavy work launched. Overall worker elapsed approximately 4 minutes; exact start was not instrumented, so no precise total CPU claim is made.
