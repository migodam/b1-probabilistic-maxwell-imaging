# Independent R7 implementation review

Review date: 2026-11-04. Reviewer: independent native Codex workstream. Scope: ordinary code reading and analytic derivation of `src/rational_search.py`, `experiments/run_rational.py`, `experiments/run_rational_coherent.py`, `experiments/generate_benchmark.py`, `experiments/run_case.py`, `src/benchmark_solver.py`, `src/cie_fft_r7.py`, together with inherited R6 interfaces/FFT operators and the original scalar regularizer. I inspected only public physical receipts, public configs and these sources. **No sealed truth, new material metrics or GT arrays were opened; no forward solve, experiment, hash or numerical test was executed for this review.**

## Assessment

The reviewed finite-response algebra, signs, prior update and coherent true Hessian are consistent with the finite scalar-TM model. R7's ordinary arms retain the R6 FFT loop, beta schedule and readouts; sourcewise versus global PR dispatch is explicit. Public receipts support accurate evaluation of the tested reduced finite endpoints against the original physical objective. This is implementation evidence, not a continuum theorem, global search certificate, statistically independent benefit or resource-matched superiority.

There is no blocking core mathematical error identified in the scoped code. Several reporting/cost limitations must remain explicit. Findings requiring parent attention were messaged immediately; no core code was silently changed.

| Finding | Severity / required interpretation |
|---|---|
| Serial rational comparator cost offsets include earlier comparators' correction work. | **Material cost limitation:** do not use these per-endpoint CPU/LU traces as standalone method-matched budgets without subtracting the unrelated earlier paths and preserving common setup. Study-total cost remains valid. |
| Positive free-coordinate Hessian is only a face restriction. | **Material scope limitation:** no full box-constrained local-minimum or basin certificate; active-bound feasible/critical-cone directions are not tested. |
| rank2/rank8/rank32 are selected internal-state dimensions. | Correct labels for maximum selected-cell resolvent dimension; not two poles, not two large regions having rank two, not measured minimal realization rank. |
| `setup_cpu_seconds` omits the anchor physical solve and per-action port-Gram/prior construction. | It is a narrowly scoped selected-resolvent setup timer; use whole-process costs or separately charge the omitted work. |
| New30 uses 24 same-grid and six different-grid same-backend objects. | Do not describe all 30 as inverse-crime-free, independent-physics, or 30 unrelated shape families. |
| Public-loader path guard omits literal `sealed_truth`. | Nonblocking hardening gap: actual manifests point to public files and current truth payload keys would be rejected, but the guard is not a complete private-path policy. |
| Rank32 searches directions derived from curvature/phase but changes their amplitudes to bound-distance weights. | Do not call the searched templates exact Hessian eigenvector paths or exact translation directions. |
| No executable release exclusion manifest was found in the scoped files. | Publication must explicitly exclude the third-party/raw artifacts listed below. This review does not assert that release filtering has already occurred. |

## 1. Finite response: derive first, then compare code

Let `A0=I-X0G`, anchor current `J0=A0^-1 X0 Einc`, total field `E0=Einc+GJ0`, observation residual `R=SJ0-Y`. Set `O=S A0^-1`, `C=G A0^-1`. For selected-cell matrix `P` and diagonal increments `D`, subtraction of the two discrete Maxwell equations gives

`A0 ΔJ = P D Pᵀ(E0 + G ΔJ)`.

Define `K=Pᵀ C P`, `OP=O P`, `EP=Pᵀ E0`. Then, on the nonsingular branch,

`ΔF = OP (I-DK)^-1 D EP`.

`RationalAnchor` obtains O using a transpose solve (`lu_solve(..., S.T, trans=1).T`), which is exactly `S A0^-1`; using an adjoint solve there would be wrong. The selected columns obtain `G A0^-1 P`, not `A0^-1 G P`. `FiniteAction` slices its selected rows to form K in the same order. These orders agree with the derivation.

For two distinct scalar cells and `D=diag(s,t)`, the determinant is

`δ=(1-sk11)(1-tk22)-st k12k21`.

The numerator is `sU+tV+stH`, where

`H=-k22 U-k11 V+k12 O1 E2+k21 O2 E1`.

The actual residual is `[δR+sU+tV+stH]/δ`. Therefore four numerator features are `R`, `U-k11R`, `V-k22R`, `H+(k11k22-k12k21)R`, with real coefficient vector `(1,s,t,st)`. `rational_search.py:53–60,70–73` implements precisely these signs, normalization by `||Y||`, and `0.5 pᵀGram p/|δ|²`. For real s,t, taking the real part of the Hermitian feature Gram is correct. This is an exact finite-discrete identity under invertibility; it is not a Born truncation or a frequency-global finite rational theorem.

For general selected state dimension q, let `T=(I-DK)^-1D`, `A=OPᴴOP/||Y||²`, `B=EP EPᴴ`, `Z=OPᴴR EPᴴ/||Y||²`. Expansion gives

`Phi_data = Phi_data0 + Re tr(ZᴴT) + 0.5 Re tr(Tᴴ A T B)`.

`np.vdot(Z,T).real` has the correct conjugation and sign, as does the trace term (`rational_search.py:40–42,75–77`). The inverse is `(I-DK)^-1D`; reversing D and K without the corresponding push-through change would be invalid. The code uses the correct side.

The general-q online evaluation costs a q×q solve plus q-scale matrix contractions, typically O(q³), independent of original source/receiver counts **after paid setup**. Rank2 evaluation is constant-sized after its four-feature Gram. Neither statement makes setup independent of n or ports. The coherent run explicitly obtains all 576 resolvent columns, constructs full port/Hessian contractions and pays those in whole CPU/RHS counts. A small q candidate is not a free n-dimensional anchor representation.

The current constructor/search assumes exactly two nonempty disjoint groups, selected cells present in the pool, nonzero legal weights, and the inherited [.02,20] box. Actual rank2 calls use two singleton groups, rank8 uses disjoint four-cell blocks, rank32 uses disjoint positive/negative sets. The reusable API does not enforce all these invariants; overlapping groups or a rank2 two-cell single group would not be valid general inputs. This does not invalidate the scoped calls but prevents an unrestricted API qualification.

## 2. Original prior and true Hessian

The original regularizer is `smooth*mean(d²)+tv*mean(sqrt(d²+eps²)-eps)`, `d=Dgrid m`. `FiniteAction.prior` identifies all edges touched by selected cells, updates those edges exactly, and adds their original-prior difference to the anchor prior. Cross-group edges are included once. It neither freezes a TV weight nor replaces the prior with a quadratic surrogate.

For the actual Hessian the edge weights are

`[2*smooth + tv*eps²/(d²+eps²)^(3/2)] / number_of_edges`.

`run_rational_coherent.py:29–31` uses exactly this expression and `Dgridᵀ diag(w) Dgrid`. This is the **true** prior Hessian. In contrast, the third return of inherited `Regularizer.evaluate` is `[2*smooth+tv/sqrt(d²+eps²)]/edges`, a useful surrogate/IRLS weight rather than the true TV Hessian. The coherent code correctly does not reuse that surrogate as an exact Hessian.

For the data Hessian, the Jacobian column is `Oi Ei`, so GN has entries `Re[(OᴴO)ij (conj(E) Eᵀ)ij]/||Y||²`. The residual-weighted second derivative is `Re[Tij+Tji]`, `T=C ⊙ conj(Z)`. Both `run_rational.py:31–33` and `run_rational_coherent.py:27–31` symmetrize correctly. The initial pair selector uses data curvature without the prior only as a heuristic ranking; candidate scoring itself retains the original prior. It must not be labeled the full physical Hessian including prior. The coherent variant does include both components.

Free coordinates are selected with `.02+1e-5 < m < 20-1e-5`. The four computed eigenpairs are the lowest of that principal submatrix, not a KKT/critical-cone spectrum. At the shell public receipt, dimension=50 and lowest reported eigenvalue≈4.34e−9 is positive; this is evidence only for the chosen finite face and numerical matrix. The negative values for ring and old curved bar show negative curvature on their selected faces; their anchors need not be stationary, and curvature alone does not establish a traversed barrier. These eigenvalues are not global Maxwell frequency/material poles.

## 3. Search scope, rank labels and available public receipts

`run_rational` evaluates 48 two-single-cell actions and 24 two-four-cell actions per disclosed historical anchor: selected state dimensions 2 and 8, each with two amplitudes. `run_rational_coherent` selects 16 positive plus 16 negative cells, giving state dimension 32, still two amplitudes. A zero amplitude can lower actual perturbation rank; “rank” here labels maximum selected-state dimension, not a certified minimal rank everywhere.

The coherent template selects cell positions using Hessian/gradient/interface/random directions, but weights them by the distance to the material box boundary (`wa=20-m`, `wb=-(m-.02)`). It therefore preserves a sign/selected-support idea, **not** the original eigenvector component ratios. “Curvature-informed weighted templates” is accurate; “exact negative-eigenvector line search” is not.

`FiniteAction.search` uses a fixed finite grid followed by L-BFGS-B from origin and the best three grid points. This is a multistart search on a two-dimensional legal domain. It does not certify a global minimum, exhaust all finite actions or construct a constrained-basin theorem. The 101-point straight path is sampled only; a sampled positive height would not be a certified continuous barrier, and absence of a sampled barrier does not exclude another path/basin geometry. The explicit code comment limits it appropriately.

Public `RESULTS/rational/*/RECEIPT.json` and `RESULTS/rational_coherent/*/RECEIPT.json` show finite predicted/original Phi endpoint differences on the order of 1e−17 for the executed candidates. This independently supports the stated **endpoint wiring** validation but not exactness of every unverified grid point, continuum discretization accuracy or scientific recovery.

Public outcomes inspected: the shell rank2/rank8/rank32 actions change original Phi only at tiny scales; ring/old-bar common corrector continuation accounts for most later improvement. Their records retain iteration-capped stops. No positive resource-adjusted escape claim follows merely from these reductions. No material-GT result was inspected in making this assessment.

## 4. Cost accounting and fairness limitations

`RationalAnchor.setup_cpu_seconds` starts **after** `model.state(m)` and base-data/prior work, and ends before every `FiniteAction` port-Gram/prior-edge setup. It measures O/C column setup only. Whole `run_rational`/coherent process CPU includes public load, anchor construction, proposal construction/ranking/maps, physical validation and correctors; that is the proper study expenditure. Physical RHS counters correctly record 48 receiver-transpose RHS and the selected-column pool size; coherent charges 576 selected columns. Existing LU is reused for these solves, while new physical endpoints use original full LU.

However, `run_rational.py:74–85` runs anchor continuation, rank2, rank8 and local comparator **serially using one model/process**. Each observer's `budget_cpu_offset` equals elapsed CPU since study start. Later endpoint budget curves therefore include earlier unrelated comparator correctors; LU/RHS counters are also cumulative. They are valid cumulative-study accounting but cannot be plotted as isolated per-method matched CPU/LU curves. Per-endpoint physical deltas in its receipt remain useful, although candidates may inherit cached factorization when identical. Any resource-matched comparison must define whether setup is common/reused or method-specific and correct these offsets without discarding paid work. No new experiment is required to avoid the unsupported claim; label the existing result as same-corrector-cap mechanism evidence.

The six full-forward checks are on already-ranked proposals and are not the same workload as all 14k reduced objective queries. They confirm endpoint evaluation and a narrow observed full-forward cost. They are not alone an equal-candidate-count throughput benchmark or an amortized inverse-algorithm advantage. Reduced solve counts/ranks and total CPU should accompany any engineering statement.

Primary direct physical L-BFGS is allowed 6000 accepted iterations while CIE-family correctors are allowed 3000 after 3000 auxiliary iterations. Those are declared caps, **not equal compute**. `BUDGET_STATES.npz` samples actual physical-Phi incumbents every ten paid physical observations with a process-prefix CPU offset, so a separate overlapping-cost comparison is possible for the primary pipelines. Its snapshots are float32 for offline material grading while actual optimization and exported final state remain float64; this numerical distinction must be disclosed. Do not extrapolate a capped trajectory as converged or assign a later incumbent to an earlier budget.

## 5. FFT dispatch and inherited algorithms

A direct source diff against `research/research_max_struct_r6/src/cie_fft.py` shows only the supported-arm dispatch expansion and material-readout callback hook in `cie_fft_r7.py`. The main R6 loop is retained. It uses the same `NativeFourierAPCS` and zero-padded `DiskPulseToeplitzGreen` classes from R6. Those retain the disk-pulse self term, finite nonperiodic convolution, forward/adjoint action and native Fourier coefficient metric. S and final physical LU stay dense. Native action omits only the historical roundoff/rank SVD; it does not imply bitwise long-trajectory equivalence to explicit dense APCS.

`benchmark_solver.arm_config` delegates the frozen settings to R6: common real box [.02,20], host .1, source metrics frozen from Y/E, material inner maxiter60, ftol1e−12/gtol1e−8, beta/material stages (3,1), (6,.2), (10,.2), optional mf12. Single-CIE uses equal current/material beta; Dual/global/A retain current beta1 while material beta becomes .2. `sourcewise_cg=True` for `cie_single`/`cie_dual`, False for global/A/A-plus; receipts inspected for object00 Dual and A confirm True/False respectively. A vs Global-Dual isolates the material change better than A vs sourcewise Dual.

A-plus substitutes only its declared material readout hook; original A retains the inherited exact one-hop derivative and L-BFGS readout. This review did not re-audit `aplus_readout.py` or certify that A-plus retained all its scientific benefits; that is outside the requested scope. `cfg.aplus_options` is merged for every arm, so frozen per-run configs should continue being inspected for non-A-plus override keys. The reviewed primary config has no such override.

The primary runner validates host and all physical stage endpoints using the original Phi and exports the lowest paid incumbent, never the substitute auxiliary objective. Its retained-candidate consistency check compares exact m arrays and objective values. The exact physical corrector is the original R6-observed/R4 implementation. Observers copy already computed states; no extra physics is hidden in ordinary checkpointing.

## 6. New30 generation and GT isolation

Generator source declares 30 unique geometry draws: six named analytic families across five acquisition/contrast/noise strata. Each object has its own public Y with independently drawn noise; repeats/noise views are not counted as separate materials. Public files contain observations/operators/noise scale/acquisition metadata, not true material. The solver imports no generator or private mapping and opens only task-declared public files. The public manifest points to 24-grid inversion operators and observation NPZs; the separate generator writes material/noise/clean-data/mapping into `sealed_truth`. I did not open those private files.

The loader rejects known truth payload keys and common private path components. Its component list omits the actual new directory name `sealed_truth` and does not resolve symlink aliases; this is a defense-in-depth gap, not observed GT leakage. Existing truth files include forbidden alpha/mask/Y_clean keys and actual tasks do not reference them. A strict claim that arbitrary private paths are impossible would overstate the implementation.

The first four strata (24 objects) generate and invert on the same 24-grid disk-pulse backend, explicitly `inverse_crime=true`. Only the last six use generation36/inversion24. Eightfold subpixel shape averaging reduces binary pixel aliasing but does not change the forward solver family. Those six qualify as different-grid same-backend robustness, not an independent Maxwell solver or continuum-reference experiment. “New independent objects” can refer to independently drawn geometries within these families; “new independent shape families” or 30 external physical validations would be incorrect.

Primary config status explicitly separates baseline/A freeze from A-plus development. The data are fresh relative to earlier R6 shapes, but an A-plus method changed after observing these public inverse outcomes is development on this cohort unless a separate method-freeze/confirmation boundary is documented. No-GT selection alone does not establish an untouched independent method-confirmation protocol.

## 7. Required publication exclusions

This review found source acquisition scripts and written rights notes, but no scoped executable publication allowlist/exclusion manifest. Before copying R7 to the public site/repository, preserve these exclusions:

- `data_public/fresnel/*.exp` and original archive/supplementary downloads: scientific use allowed, unrestricted redistribution license not established.
- `RESULTS/fresnel_*/endpoint.npz`: includes `measured_Y`, hence raw measured data despite residing under RESULTS. Use `RESULTS/fresnel_public_derived.npz`, metrics/configs and generated plots instead.
- Third-party `delegated/public_literature/*.pdf`, `*.html`, and full-text extracted `*.txt` (e.g. CIE/LSL); local access does not grant publication rights.
- Any copied full Pro handoff archive, full third-party paper source or external copyrighted figure. The review inventory/prose, DOI/URLs and short evidence summaries are the publishable handoff.
- Sealed new truth and private mapping until the parent's authorized scoring/release stage; do not copy them merely because their directory is within R7.

Original R7-generated figure PDFs are the project's own artifacts and need not be confused with third-party paper PDFs. Acquiring public papers for local review is distinct from republishing them. Preserve the official downloader, stable DOI/data links and access/error receipts as reproducibility evidence.

## Review conclusion for parent use

Implementation-level finite algebra: **consistent in the scoped scalar discrete model**. Original prior and coherent Hessian: **consistent**. Inherited FFT/beta/sourcewise dispatch: **consistent, with disclosed roundoff and readout-hook differences**. Public endpoint validation: **supported for executed tested endpoints**. Full-space basin/global novelty/equal-cost superiority/continuum truth: **not established by this review**. Outstanding items are reporting boundaries, cost-curve interpretation, optional loader hardening and explicit release filtering; no new scientific computation is requested by this review.
