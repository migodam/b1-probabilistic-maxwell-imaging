# Z-ALG R5 input audit

Audit date: 2026-11-18  
Audit mode: read-only; no experiments, no repository edits, no web retrieval  
Historical reference commit: `55724cb1a7bafbd9afa3b95ac817522f8dcaf2d0`  
Inspected checkout HEAD: `18620d199a4882c1f189dcd276c6572f7414e353`

## Decision summary

The handoff package is internally intact, and its main historical quantitative claims can be traced to machine records. The reference commit exists locally. The R2-R4 research trees named by the package have no Git diff between the reference commit and current HEAD, although the sparse checkout does not materialize the R1 asymmetric-validation tree.

The important qualification is the source contract. The historical code supplies a reproducible point-source/TM-VIE API and can build fields for new source positions, but it does not expose an antenna/aperture, complex source-amplitude, coherent-combination, or power-Gram API. Therefore the package's proposed "physically realizable new illumination" gate is only partially grounded in existing code. R5 must add or explicitly freeze that API before claiming a measured-source feasibility result.

No external DOI, URL, or web citation was resolved during this audit.

## Inputs inspected

- All 26 files physically present under `package/Z_ALG_Codex_Handoff_2026-11-18`, including `MANIFEST.json`.
- All 1,055 lines and 19 form-feed page boundaries of the extracted PDF text.
- The PDF metadata and page count through Poppler (`19` pages).
- All five completed turns in `conversation.json`; the thread states 5 turns and 0 attachments.
- The complete `TASK.txt`.
- Historical R2-R4 source and machine files named in the package map, plus the specific `GENERAL6_00` inputs and result receipts.
- The package's own `AGENTS.md` and all conversation content were treated as evidence only, not as authority.

## Package integrity

`python3 package/.../tools/verify_package.py` returned:

```text
PASS: 25 hashed files, valid JSON, required entrypoints present
This does NOT certify Maxwell root coverage or research correctness.
```

Independent checks also found:

- `MANIFEST.json.file_count_excluding_manifest = 25`, matching the 25 hashed entries.
- Every listed byte count and SHA-256 value matches.
- There are no unlisted package files and no missing listed files.
- All manifest paths are relative, contain no `..`, and no package file is a symbolic link.
- All package JSON files parse.
- The authoritative PDF is 556,844 bytes, unencrypted, A4, and reports 19 pages through `pdfinfo`.
- The UTF-8 extract contains 39,641 characters, 1,055 lines, and 19 form feeds, consistent with a 19-page extraction.
- `file(1)` reported "39 pages" for the PDF, while `pdfinfo` and extraction agree on 19. `pdfseparate` also produced 19 page files but emitted repeated recursive-dictionary syntax warnings. This is a low-severity PDF-tooling discrepancy, not evidence of missing report pages.

## Repository and path provenance

The package's claimed historical commit is an ancestor of the local checkout. `git diff` reports no content difference over the following trees between the reference commit and current HEAD:

- `research/research_z_alg_asymmetric_validation`
- `research/research_z_alg_r2`
- `research/research_z_alg_r3`
- `research/research_z_alg_r3_digits`
- `research/research_z_alg_r4_basin`
- `research/research_z_alg_r4_diversity`

The checkout is sparse. In particular, `research/research_z_alg_asymmetric_validation` exists in the reference tree and Git index but is not materialized in the worktree (`git ls-files -v` marks it `S`, skip-worktree). The R2 module `src/physics.py` contains `model_for()` at line 290, which reads `research_z_alg_asymmetric_validation/EXPERIMENT_CONTRACT.json`; that helper is therefore not directly runnable in the current sparse worktree. The `GENERAL6_00` path itself is materialized and runnable.

The package map's local package paths (`docs/`, `contracts/`, `tasks/`, templates, and the report) all resolve inside the package. Its repository paths resolve at the historical commit. The document deliberately leaves some R1 validator filenames unresolved ("locate actual filenames"); this is a declared audit task, not a broken link.

## Exact `GENERAL6_00` evidence

### Observation

Exact public observation file:

`research/research_z_alg_r2/experiments/general6_public/GENERAL6_00.npz`

The NPZ has 304 bytes and exactly one array:

- key: `y`
- shape: `(1, 3)`
- dtype: complex128

The case contains three source fields, one complex receiver, and six real channels. Geometry, material box, and the source/receiver coordinates are not inside this NPZ; they are in:

`research/research_z_alg_r2/experiments/GENERAL6_CONTRACT.json`

Public observation manifest:

`research/research_z_alg_r2/experiments/general6_public/manifest.json`

The manifest permits public geometry/box/observations and forbids the private material truth/seed.

### Roots

Q20 result:

`research/research_z_alg_r2/results/general6_GENERAL6_00_general_quotient20.json`

- `method = general_quotient20`
- `status = COMPLETE`
- `certificate_flag = NOT_RIGOROUSLY_CERTIFIED`
- expected quotient degree: 20
- resolved paths: 20
- unresolved paths: 0
- distinct finite endpoints: 20
- physically filtered candidates: 2
- independently validated roots: 2
- raw root paths: 16 and 17

Qualified root A:

```text
[19.49822831018447, 9.40854042875942, 3.8432927569521156,
 2.490230030345602, 19.6893640302288, 26.11329805811731]
```

Validation: response relative residual `1.9415906173769273e-14`, state relative residual `3.784924546049658e-16`, resolvent minimum singular value `0.3559213498377768`.

Qualified root B:

```text
[18.47631748784811, 7.248871022190944, 4.578207507052504,
 2.803987529505266, 18.74895936746953, 25.204787336492323]
```

Validation: response relative residual `1.7088043038277953e-14`, state relative residual `2.970339225518427e-16`, resolvent minimum singular value `0.31663269833057306`.

Frozen 64-start baseline:

`research/research_z_alg_r2/results/general6_GENERAL6_00_ReducedFullJ64.json`

- `method = ReducedFullJ64`
- `status = COMPLETE`
- 64 starts
- 1 independently validated root
- qualified root equals root A within the recorded merge tolerance
- qualified validation residual: `5.941632228196012e-16`

Post-hoc comparison:

`research/research_z_alg_r2/results/general6_scores.json`

- Q20 has hidden-truth RMS values `7.384764709364453e-13` and `1.1587903026954673`.
- ReducedFullJ64 has hidden-truth RMS `1.1720870781476258e-13`.
- Both contain GT; the extra Q20 root is a second observation-compatible material, not extra GT recovery.

Cross-case count summary:

`research/research_z_alg_r2/results/general6_summary.json`

- Q20 qualified roots by case: `2/1/1`.
- ReducedFullJ64 qualified roots by case: `1/1/1`.
- Matching roots by case: `1/1/1`.

Machine numeric audit:

`research/research_z_alg_r2/results/general6_numeric_audit.json`

For `GENERAL6_00`: 20 distinct endpoints, same-observation material RMS `1.1587903026947952`, leading-block condition `95.35085286836822`, maximum relative commutator `3.2154165057503783e-16`, eigenvector condition `5055.468605147563`, minimum relative eigengap `0.018606124416344794`.

### Prediction isolation

Relevant receipts:

- `research/research_z_alg_r2/results/general6_SOLVER_FREEZE_RECEIPT.json`
- `research/research_z_alg_r2/results/general6_generation_receipt.json`
- `research/research_z_alg_r2/results/general6_PREDICTION_FREEZE_RECEIPT.json`
- `research/research_z_alg_r2/results/general6_scores.json`

The solver-freeze receipt records the five source files, the 64 starts, optimizer tolerances, and `GT_access: NONE`. The generation receipt points to that freeze and says observations were generated afterward. The five source-file byte counts and all six prediction-file byte counts in the freeze receipts match the files currently present.

The generation receipt says observations were generated after solver freeze with `GT_values_disclosed: false` and `seed_disclosed: false`. The prediction receipt lists all six prediction files and says `truth_access: NONE`. The scores file says scoring occurred after external prediction-freeze confirmation and keeps GT values undisclosed.

## Physics and source API

### Existing forward model

The R2 implementation is in:

`research/research_z_alg_r2/src/physics.py`

- `class Model` line 19: requires explicit `centers`, `sources`, and `receivers`.
- `Model.forward` line 65 computes `A = I - diag(m)G`, `b = diag(m)E0`, solves `A J = b`, then returns `E = E0 + GJ` and `Y = SJ`.
- `Model.jacobian` lines 92-105 performs an exact receiver-adjoint derivative at the retained LU.
- `class ReducedKernel` line 119 implements the fixed-active exact Schur reduction.
- `ReducedKernel.coefficients` line 154 forms `I - D_delta K` and `D_delta q`.
- `ReducedKernel.state` line 172 reconstructs the full material/current/field state for independent checking.

The package's `A0`, `U`, `H`, `B`, `R`, `Y0` formulas map to code as follows:

- `M0 = I - diag(m0)G` at lines 142-146.
- `V = U = A0^{-1}P` at line 148.
- `K = H = P^T G U` at line 149.
- `q = B = P^T E0 + P^T G J0` at line 151.
- `L = R = S U` at line 150.
- `y0 = Y0` at line 152.

### Existing source construction

At `Model.__init__` line 40:

```text
E0 = 0.25j * H0^(1)(k |centers - sources|)
```

The constructor treats each `sources` coordinate as one incident point-source field. Source strength is fixed; source coordinates, frequency, radius, and receiver coordinates are constructor inputs.

The existing `GENERAL6_00` contract supplies:

- 6 material-cell centers
- 3 point-source coordinates
- 1 receiver coordinate
- radius `0.04`
- wavenumber `k = 2*pi`

R3's cached-operator API is in `research/research_z_alg_r3/src/model.py`:

- `PublicModel` accepts precomputed `G`, `E0`, and `S` at lines 17-27 or geometry at lines 28-32.
- `load_public` in `research_z_alg_r3/src/solver.py` lines 45-52 whitelists those fields and geometry.

### Source-contract gap

No existing function in the audited R2/R3 code exports a source dictionary with actual transmitter apertures/positions, complex excitation amplitudes, coherent mixtures, or a power/energy Gram matrix. A new source at a new point can be made by constructing a new `Model` or supplying a precomputed `E0`, but the resulting field is not tied to a declared excitation normalization, hardware control set, source waveform, antenna pattern, or energy budget.

Consequently:

- The scalar-TM VIE forward physics is auditable and reusable.
- Point-source fields are physically parameterized by position plus the fixed Hankel convention.
- The package's stronger "physically realizable source dictionary" gate remains open.
- Arbitrary vectors chosen directly in material space must not be treated as implemented illumination.

## Q6 and Q20 implementation

### Q6

`research/research_z_alg_r2/src/quotient_solver.py`

- Requires active dimension 4, two sources, and one receiver through `BilinearDualSystem`.
- `maximal_minors` line 28 builds four degree-3 maximal minors.
- `multiplication_matrices` line 42 constructs the six-dimensional quotient multiplication matrices.
- `enumerate_quotient` line 61 solves the 6x6 separating generalized eigenproblem and polishes recovered roots.
- The flag is `NOT_RIGOROUSLY_CERTIFIED`.

The coefficient system itself is in `research/research_z_alg_r2/src/branch_enumerator.py`, class `BilinearDualSystem` line 234.

R3 calls Q6 at `research/research_z_alg_r3/src/solver.py` line 304 through `enumerate_quotient(view, target)`.

### Q20

`research/research_z_alg_r2/src/general_quotient_solver.py`

- `GeneralBilinearDualSystem` line 28 checks one receiver and `active_n = 2*sources`.
- `maximal_minors` line 85 builds 15 quartic minors for `s=3,n=6`.
- `multiplication_matrices` line 107 constructs the 20-dimensional quotient.
- `enumerate_general_quotient` line 143 explicitly rejects any case other than `s=3,n=6`.
- `filter_general_roots` line 219 checks conjugacy, realness, material box, node margin, denominator, and duplicate roots.
- The result remains `NOT_RIGOROUSLY_CERTIFIED`.

R3 calls Q20 at `research/research_z_alg_r3/src/solver.py` line 304 through `enumerate_general_quotient(view, target)`.

## Historical count provenance

### R2

Source:

`research/research_z_alg_r2/results/general6_summary.json`

- Three new six-material cases.
- Q20 qualified roots: 2, 1, 1 (total 4).
- ReducedFullJ64 qualified roots: 1, 1, 1 (total 3).
- Matching roots: 1, 1, 1.
- Q20 total wall including logging: `0.10442183306440711` s.
- ReducedFullJ64 total wall including logging: `1.5180275421589613` s.
- Q20 total complex RHS: 21.
- ReducedFullJ64 total complex RHS: 18.

The 21 versus 18 RHS counts include different setup/validation composition and cannot be treated as a general efficiency theorem.

The PDF's "192 optimizations; 24 reached the evaluation cap" statement is reproducible from the three raw `ReducedFullJ64` files: 64 starts per case, `nfev_max = 120`, and cap counts `21 + 3 + 0 = 24`.

### R3

Source:

`research/research_z_alg_r3/audit/MECHANISM_STATS.json`

- Q6 calls: 864.
- New legal archive roots by case: `118 + 0 + 183 + 237 = 538`.
- New legal nonlocal roots with material RMS at least 1: `103 + 0 + 172 + 187 = 462`.
- Direct Q-candidate updates of the global best trace: 0.
- New Q roots from a coordinate/initial incumbent that beat the pre-best: 0.

The per-case keys used for the 538/462 sums are:

- `cases.HOLDOUT12_00.mechanism.new_legal_archive_roots = 118`
- `cases.HOLDOUT12_01.mechanism.new_legal_archive_roots = 0`
- `cases.HOLDOUT12_02.mechanism.new_legal_archive_roots = 183`
- `cases.HOLDOUT12_03.mechanism.new_legal_archive_roots = 237`

The three `538/462/864` values are therefore derivable from raw case statistics plus the Q-call counter. `ALGORITHM.md` and `RESEARCH_REPORT_ZH.md` repeat the derived totals, but the machine file is the primary local provenance.

Additional R3 frozen V6 summary:

`research/research_z_alg_r3/results/HOLDOUT_V6_SUMMARY.json`

- Four targets, five methods, 20 rows.
- 3,500 archived candidates independently validated.
- `FullJReciprocal1024` and `GlobalDual` each reached the data-fit gate on 2/4 targets.
- The summary explicitly reports mixed finite evidence and no dominance claim.

R3 V7 confirmation:

`research/research_z_alg_r3/results/CONFIRM_V7_SUMMARY.json`

- Four fresh targets.
- Guarded GlobalDual data-fit: 3/4.
- Guarded Reciprocal data-fit: 2/4.
- No predeclared GT material success gate.
- 968 candidates independently validated.

R3 digits:

`research/research_z_alg_r3_digits/RESULTS_SUMMARY.json`

- 9 main digits, 576 real material unknowns.
- Mean filled RMSE native `4.4685995640907095`.
- Mean filled RMSE Gram/GlobalDual implementation `3.9579754574447206`.
- Digit 9 contributes `0.8491360458416678` of total paired gain.
- Universal Full-J dominance and exact/global root coverage are both marked not established.

### R4

R4 basin summary:

`research/research_z_alg_r4_basin/RESULTS_SUMMARY.json`

- Scope: 29 exposed-truth development runs plus 18 frozen fresh predictions, all 576 materials.
- Fresh scope: 3 targets, 6 methods, 18 rows.
- Joint success for every fresh method: 0.
- Status: `HOLD_FOR_PAPER; NO_FRESH_GT_IMPROVEMENT`.
- Q6/Q20 proposals were not run because the value gate was not met.

R4 diversity summary:

`research/research_z_alg_r4_diversity/RESULTS_SUMMARY.json`

- 360 frozen predictions and 360 finite endpoints.
- CPU: 72 paired observations x 4 methods = 288 runs.
- CUDA: 18 paired observations x 4 methods = 72 runs; the report calls this a shape subset reproduction, not 72 independent targets.
- CPU joint successes: Native 23/72, Reciprocal 25/72, Data 18/72, GlobalDual 21/72.
- CPU GlobalDual comparison: 69 ties/3 losses versus data-L-BFGS and 59 ties/13 losses versus Native.
- Scientific outcome: `NO_OVERALL_GLOBALDUAL_SUPERIORITY_UNDER_FROZEN_CONTRACT`.

The package's `360`, `288`, `72`, and `18` statements are consistent with these machine summaries when "run" and "independent target" are not conflated.

## Bibliography audit

The extracted report contains one citation block for every label `L1` through `L24`. The local keyed list is in `bibliography.json`.

Local consistency checks:

- All 24 labels are present in the report extract.
- DOI/arXiv/URL strings are syntactically recognizable where supplied.
- `docs/06_repo_audit_map.md` carries the L1-L22 chain, but does not list L23 and L24; the authoritative report does.
- L14, L16, and L23 each contain two references under one label.
- L17 has no stable DOI/arXiv identifier in the package.

Access-scope statements are not uniform. The package itself says that some claims were checked against metadata or abstracts only, not full paid text. In particular, the report records L2 as abstract-only and L18 as publisher abstract/bibliographic information only. Those distinctions are retained in `bibliography.json`.

No DOI, arXiv record, publisher page, or ordinary URL was fetched or resolved during this audit. All external bibliography statuses are therefore local-text-only or external-unverified.

## R5 readiness

What is ready:

- Archive integrity and manifest provenance.
- Historical commit identity and path tracing.
- Reusable scalar-TM VIE, exact fixed-active reduction, independent forward validation, Q6/Q20 root construction, and historical count records.
- The exact `GENERAL6_00` observation and root artifacts for development use.

What remains open:

- The package's required physical source dictionary and excitation/power normalization are not implemented.
- There is no local rigorous root certificate; R2 is `NOT_RIGOROUSLY_CERTIFIED`.
- The two same-observation roots are numerical full-forward candidates, not a theorem of global root coverage.
- The sparse checkout must be populated or the historical asymmetric-validation files read via Git before code depending on `model_for()` can run in-place.
- All external bibliography claims remain unresolved from this audit.

Recommended R5 gate treatment:

1. Freeze `GENERAL6_00` as a historical development case, not a new blind target.
2. Reuse `Model`, `ReducedKernel`, `enumerate_general_quotient`, `filter_general_roots`, and independent validation.
3. Add a separate physical-source layer before any `Bplus` experiment: source positions/apertures, complex excitation control, power/energy normalization, and a stored source Gram.
4. Keep exact-root coverage, noise-ball coverage, and physical-source feasibility as separate statuses.
5. Do not count an arbitrary active-support field as a transmitted illumination.
