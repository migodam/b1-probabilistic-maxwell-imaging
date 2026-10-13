# B1 Source Audit - Reusable Physics and Prior Negative Evidence

Scope: `Gaussian/B1/b1-probabilistic-maxwell-imaging`.
Audit date: 2026-10-13. Author: B1 source-audit worker (Codex).
Status: READ-ONLY INVENTORY. No new physics, no solver calls, no dataset
arrays, no hashes, no installs. Codex parent owns new theory, algorithm choice
and final science.

This document catalogues what homogeneous / concentric-layered vector-Mie
formulas and implementations already exist in the tree, their exact
conventions, the tests that back them, license/provenance status, and the
prior negative results that bound any B1 probabilistic multi-hypothesis
imaging claim. It is a reuse map, not a validation of B1.

## 0. Constraints honored

- Only `research/delegated/source_audit/` was written.
- No script was executed that invokes physics; inspection used `cat`, `sed`,
  `rg`, `ls`, `find`, and `json`/`numpy` reads of already-frozen JSON only.
- No historical/locked dataset arrays were opened; no truth arrays read.
- No SHA256 checks were run (existing hash fields are historical provenance).
- No sub-worker or Agentic-AI-Scientist was spawned.

## 1. Sources inspected

Provided A32 package (input, read-only):

- `inputs/source_package/A32_Research_Package/a32_checks.py`
- `.../A32_CHECKS.json`, `.../A32_CHECKS_INITIAL.json`, `.../A31_RECHECK.json`,
  `.../A32_RESOURCE_LEDGER.json`, `.../MANIFEST.json`
- `.../A32_THEORY_REPORT_ZH.md`, `.../A32_NEXT_DECISION_ZH.md`,
  `.../README_ZH.md`
- `.../upstream_a31/a31_matrix_checks.py`

Executed / historic project sources:

- `Gaussian/A13/REVISION_R2/src/mie_validation.py`
- `Gaussian/A13/REVISION_R2/src/sphere_mechanism.py`
- `Gaussian/A13/REVISION_R2/results/mie_validation/mie_validation.json`
- `Gaussian/A9/code/sphere_reversal.py`
- `Gaussian/A12/deps/a9/sphere_reversal.py` (independent homogeneous reference)
- `Gaussian/A12/REVISION_R1/code/r1_sphere_mie.py`
- `Gaussian/A9/code/experiment_support.py`
- `Gaussian/A15/research/delegated/layered_reference/run_reference_checks.py`
- `Gaussian/A15/research/delegated/layered_reference/SUMMARY.md`
- `Gaussian/A15/research/delegated/layered_reference/results_20260930/reference_checks.json`
- `Gaussian/A15/research/delegated/layered_reference/results_20260930/reference_arrays.npz`
- `Gaussian/A15/src/validate_layered_dda.py`
- `Gaussian/A15/results/layered_dda/results.json`
- `Gaussian/A15/FAILURES_AND_NEGATIVES.md`
- `Gaussian/A6/vendor/treams` (0.4.7) and `treams-0.4.7.dist-info/LICENSE`
- `Gaussian/A17/CLOSURE_R1/manuscript/CORE_THEORY.tex`, `.../MIE_ZH.tex`
- `Gaussian/A32/a32-material-risk-directed-feedback/A32_E0_E1_DECISION.md`
- `Gaussian/A32/a32-material-risk-directed-feedback/src/a32/a17_factors.py`
- `Gaussian/A32/a32-material-risk-directed-feedback/src/a32/frozen_physics.py`
- `Gaussian/A30/a30-maxwell-material-manifold/A30_E0_E1_DECISION.md`
- `Gaussian/A30/a30-maxwell-material-manifold/src/a30/physics.py`
- `Gaussian/A30/.../vendor/frozen_a29/PROVENANCE.json`, `src/a24/physics.py`
- `Gaussian/A28/a28-e1-material-information/E1_DECISION.md`,
  `Gaussian/A29/a29-material-hypothesis-coverage/E1_DECISION.md`

## 2. Reusable homogeneous vector-Mie (the strongest reuse asset)

Primary kernel: `Gaussian/A13/REVISION_R2/src/mie_validation.py`.

Reusable functions:

- `_riccati(n, z)` - Riccati-Bessel `psi_n, psi'_n, xi_n, xi'_n`.
- `mie_coefficients(nmax, m, x)` - Bohren-Huffman `a_n, b_n`.
- `mie_s1_s2(mu, a, b)` - angular `S1, S2` via `pi`/`tau` recurrences.
- `mie_farfield_vector(...)` / `mie_farfield_angles(...)` /
  `mie_farfield_bank(...)` - physical vector far field and per-receiver
  algebraic projections; arbitrary incidence and transverse polarization.
- `mie_efficiencies(a, b, x)` - `Qext, Qsca, Qabs` (passivity signs).
- `_rayleigh_reference_checks`, `_mie_protocol_checks`, `field_metrics`.

Frozen configuration used and reported: domain edge `1.5`, sphere radius
`0.35`, `k = 2`, contrast `chi = 1 + 0.05j`, 6 illuminations, 64 receivers at
radius `5`, 2 transverse polarizations, `nmax` in `(20, 30, 40)`.

Exact conventions (declare these in B1 before reuse):

- time factor `exp(-i omega t)`;
- incidence `E_i = e_i exp(+i k u_i . x)`;
- outgoing far field `E_s = exp(+i k R)/R * A(n_hat)`;
- Mie amplitude matrix in the scattering plane
  `[E_perp, E_par_out] = exp(+i k R)/(-i k R) [S1 E_perp, S2 E_par_in]`;
- DDA receiver rows are algebraic projections `b_q . E_s` (no conjugation);
- `1/(-i k) = i/k` prefactor; forward (backward) caustic handled by the
  deterministic `S0` (`S2`) limit rather than a scattering-plane basis.

Independent cross-checks already recorded (floating-point, not a proof):

- local vs existing `sphere_reversal.mie_coefficients` abs diff `0.0`;
- forward oracle vs existing `S0` formula relative diff `9.24e-17`;
- optical theorem `Cext` relative diff `1.28e-16`;
- numeric 64-direction `Csca / Csca_efficiency = 0.99970`;
- Rayleigh limit `a1` rel err `7.67e-9`, `S0` rel err `2.29e-7`;
- `nmax = 20 / 30 / 40` bank-identical (`0.0` aggregate rel).

Related homogeneous implementations that reuse the same formulas:

- `Gaussian/A9/code/sphere_reversal.py` and
  `Gaussian/A12/deps/a9/sphere_reversal.py` - `mie_coefficients`,
  `forward_scalar`, `efficiencies` (this is the "existing" reference the
  A13 kernel is checked against; treat as a second writer, same convention).
- `Gaussian/A12/REVISION_R1/code/r1_sphere_mie.py` - independent co-polarized
  forward-amplitude check with a `nmax=30` vs `40` truncation indicator.

## 3. Reusable concentric-layered vector-Mie (multi-domain)

Primary kernel: `Gaussian/A15/research/delegated/layered_reference/run_reference_checks.py`.

- `coefficients(lmax, k, radii, epsilon)` assembles achiral multilayer `a, b`
  by calling `treams.coeffs.mie(l, k*radii, epsilon, mu=ones, kappa=zeros)`
  per `l = 1..lmax`, then converting helicity to the inherited `e^{-i w t}`
  outgoing Mie convention with the explicit plain-array map
  `a = -(h00 + h01)`, `b = -(h00 - h01)`.
- Input contract: `radii` strictly increasing positive; exactly one `epsilon`
  per region including the vacuum embedding; achiral finite/symmetric guard
  (`rel(h, h.T) <= 1e-12`, `|h00 - h11| <= 1e-12`).
- Angular projection is imported from `mie_validation.mie_farfield_bank`
  (Section 2); the layered far field is therefore NOT a second independent
  angular implementation.

Frozen configuration / checks reported (all 20/20 PASS, wall 0.512 s):

- radii `[0.18, 0.35]`, `eps = [2.4 + 0.08j, 1.7 + 0.04j]`, embedding `1`,
  `mu = 1`, `kappa = 0`, `k = 2`, `lmax` in `(20, 30, 40)`.
- collapsed homogeneous layer vs independent Mie coefficients `3.21e-15`
  rel; projected far field `3.23e-15` rel.
- vacuum outer shell vs smaller homogeneous core `4.30e-15` rel;
  all-vacuum layered scatter exactly `0.0`.
- homogeneous uniform-real derivative, multilayer engine vs independent Mie
  central FD: `2.8e-12 .. 2.9e-11` rel.
- layered complex-derivative linearity `<= 8.7e-10` rel; layered FD step
  convergence `<= 3.8e-10` rel; layered `lmax 20/30` vs `40` `0.0`.
- `treams.TMatrix.sphere(poltype='helicity')` vs direct block assembly `0.0`.

The kernel supports arbitrary strictly-increasing concentric radii, but only a
two-region (core + shell) target was exercised.

## 4. Reusable acquisition / observation conventions

`Gaussian/A9/code/experiment_support.py`:

- `acquisition(tx, rx, rotation)` -> Fibonacci illumination directions, unit
  transverse polarizations, receivers at `5 * unit`, two transverse receiver
  basis vectors. Used unchanged by both the A13 Mie harness and the A15
  layered harness (6 tx / 64 rx / R = 5).
- `grid(n, edge)`, `contrast(...)`, `haar(...)`, `Tangent` helpers.
- DDA far-field observation is the strict leading `sqrt(volume) * G_recv`
  receiver with algebraic (non-conjugated) projection.

These give B1 a drop-in geometry/observation contract for a sphere problem
without inventing a new convention.

## 5. Reusable 2-D scalar TM forward/cache backend (not 3-D vector)

- `Gaussian/A30/a30-maxwell-material-manifold/src/a30/physics.py` and the A32
  copy `.../a32/src/a32/frozen_physics.py` are adapters over an unmodified
  vendored A29 backend (`vendor/frozen_a29/src/a24/physics.py`, `a25`, `a27`,
  `vendor/a2/physics.py`).
- Contract: 2-D scalar TM, `c = sqrt(v) j`, `e = sqrt(v) E + G c`,
  material law `a(chi) = chi`, complex128, single SciPy LU per RHS batch.
- Frozen geometry: `24 x 24`, cell volume `1/576`, `adapter = "a24_2d_tm_a2"`,
  8 sources / 32 receivers, descriptor-validated. `PROVENANCE.json`: exact
  byte copy from A29 commit `5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4`,
  "no file digests computed or recorded".
- `Gaussian/A32/.../src/a32/a17_factors.py`: A17 bounded dense algebra
  (`gn_step`, `exact_step_change`, `paired_support`, `realify`); reusable for
  the paired likelihood-ratio signature theory, not a Maxwell solver.

This backend defines B1's "cheap/affine response" and budget-contract idioms,
but it is 2-D and cannot be presented as 3-D vector Maxwell.

## 6. License / provenance summary

- `treams` 0.4.7 (vendored at `Gaussian/A6/vendor/treams`): MIT License,
  Copyright (c) 2025 Dominik Beutel, https://github.com/tfp-photonics/treams.
  Distribution Beta; ships cp313 Darwin compiled `coeffs`/`config` extensions.
- In-house Gaussian modules (A9/A12/A13/A15/A17/A24/A30/A32): the user's own
  research artifacts, no separate LICENSE file inside those trees; published
  sub-repos under this program (A20/A21/A22) carry MIT License, Copyright (c)
  2026 migodam. Treat all in-tree research code as internal provenance.
- External analytic Mie formulas (Bohren-Huffman) are textbook/standard and
  carry no license obligation; only the `treams` binary/library dependency
  needs attribution and an MIT notice if redistributed.

## 7. Negative-evidence boundaries (do not overwrite these)

| # | Source | Negative fact | Consequence for B1 |
|---|--------|---------------|--------------------|
| N1 | A13 `mie_validation.json` | DDA-vs-Mie PARTIAL: field gate 3% fails at grids 16 (0.0417) and 28 (0.0644); derivative gate 5% fails at grid 28 (0.0616); 20/24/32 pass. Non-monotone. | 28^3 is a known bad grid; a single grid pass is not DDA convergence. |
| N2 | A17 `MIE_ZH.tex` | 28^3 sphere mask holds only 93.14% of true sphere volume; 32^3 = 99.56%, field 0.252% / derivative 0.731%. | Grid/mask volume bias is a first-order error source; report discrete volume ratio. |
| N3 | A15 `results/layered_dda/results.json` | Layered DDA qualifies only grids 20/24/32; 12/16/28 fail; 28 fails both gates. Scope: concentric radial reference, shared angular projection, uniform real derivative only. | Use 20/24/32 for layered DDA; no full voxel Jacobian certified. |
| N4 | A15 `SUMMARY.md` / `FAILURES_AND_NEGATIVES.md` | `treams` 0.4.7 ABI-incompatible with current SciPy (`__pyx_fuse_2sph_harm`); `poltype='parity'` path emits warnings; only `helicity` path clean; must run in a fresh process with vendor prepended. | Any B1 layered reference must isolate the treams process and lock the helicity path. |
| N5 | A30 `A30_E0_E1_DECISION.md` | E0 E0 PASS for offline capacity only (hull +13.12%, sparse64 +35.69%, PCA64 +21.91%); E1 NO-GO: min compensated full-F whitened distance 8.84 vs target <=1, 0/12 complex parents; 161/192 full-wave calls. | Capacity != measurement-selected recovery; do not reuse A30 numbers as online performance. |
| N6 | A32 `A32_E0_E1_DECISION.md` | E0 FAIL: main method worsens complex-material MSE by 39.46% vs N3 (target >= +5%); offline same-chart capacity exists (PCA64 +38.68%, Sparse64 +49.85%); E1 and 3-D NOT_RUN; 0 new full-wave states. | Global Gaussian/chart MAP overwrote N3 structure; a B1 estimator must preserve N3-class structure. |
| N7 | A28/A29 `E1_DECISION.md` | Material information exists but 6-step F/L TV improves complex MSE only +1.01% / +0.40% (bootstrap crosses zero) while residual drops 27.33%; budget compliance FAIL (>= 807/768, 8x8 backend forwards bypassed reservation). | Separate data-residual improvement from material recovery; enforce budget reservation. |
| N8 | A14/A15 reference program | Reference posterior not established (A15 `Gate_A = UNRESOLVED_REFERENCE_INSUFFICIENT`, `Gate_B = FAIL`, `Flow NOT RUN`, overall NO-GO; A14 reference/Flow not run; earlier layered integral 94/384 queries). | Do not use one sampler's output as reference truth; build a trusted E0 reference first. |
| N9 | A32 theory `sec-3.3/3.4` | No uniform topology detection without a minimum hole scale (`P_err -> 1/2` as `a -> 0`); topology event needs an explicit physical definition and `P(K=0) > 0` atomic prior mass. | Declare minimum cavity scale and a proper discrete/atomic topology prior, or the posterior question is ill-posed. |
| N10 | A32 `a32_checks.py` (checks 5,7,8,10,13) | Constant field offset changes material odds with zero Jacobian defect; equal best residual != equal evidence; support exclusion forces `TV >= excluded mass`; missing Ito divergence changes the invariant law; non-Gaussian prior breaks Gaussian mean-support. | Do not equate residual/field accuracy with posterior probability; keep the counterexamples as B1 falsification tests. |
| N11 | A32 package `README_ZH.md` | A22 independent archive not obtained; DiffusionEMIS and several papers only abstract/metadata; no first-ever priority claim. | Reading boundary for the B1 literature map; do not claim non-existence or priority. |
| N12 | A31 `A31_RECHECK.json` scope | 22/22 PASS but "floating-point verification and explicit counterexamples; not interval-arithmetic proof or imaging evidence"; only 1 tiny vector-dipole material state, 2 RHS. | Keep A31/A32 checks as algebra/implementation evidence, never as imaging or 3-D recovery. |

## 8. Reuse recommendation for B1 E0 (reference posterior)

The B1 plan (`A32_NEXT_DECISION_ZH.md`, E0) is a 3-D concentric layered sphere
vector-Mie reference with two topology branches and two continuous material
variables. The reuse stack that already exists:

1. Homogeneous branch/limit and angular projection:
   `A13/REVISION_R2/src/mie_validation.py` (Section 2).
2. Concentric multi-layer radial coefficients:
   `A15/.../run_reference_checks.py` (Section 3), cross-checked against (1) in
   the collapsed and vacuum-shell limits.
3. Acquisition / observation geometry: `A9/code/experiment_support.py`
   (Section 4).
4. Budget/state-reservation idioms: A30/A32 `frozen_physics.py` and
   `PROVENANCE.json` (Section 5) - reuse the reservation contract, not the 2-D
   physics.

Recommended first reuse actions (parent's call, not done here): copy the two
kernels into a B1 module, freeze the `exp(-i w t)` / `exp(+i k u.x)` /
`exp(+i k R)/R` / algebraic-projection conventions in code, and keep the
collapsed / vacuum-shell / vacuum-zero / multilpole-convergence checks as
regression tests before any B1 probability is computed.

## 9. Uncertainties / gaps (for the parent)

- U1: No second independent *layered* solver is present. The layered checks are
  step convergence + complex linearity of `treams`, not a two-implementation
  cross-check; the angular projection is inherited from the homogeneous
  harness. A genuinely independent layered angular implementation does not
  exist yet.
- U2: `treams` is an external MIT binary dependency with a current-SciPy ABI
  trap; B1 must either pin the vendor process or own a pure-NumPy multilayer
  implementation. Neither was executed here.
- U3: Only a two-region (core + shell) layered example was exercised; behavior
  for 3+ shells and for strongly absorbing complex contrast is unverified in
  this tree.
- U4: All prior *inverse* physics is 2-D scalar TM (`24 x 24`); there is no
  in-repo 3-D vector inverse solver. A 3-D reference is analytic-forward only.
- U5: A32 (E0) tolerances (P(H) within 0.02, log-BF within 0.1) and the 384
  material-frequency cap are research tolerances/budgets, not proven error
  bounds; E0 may legitimately return `REFERENCE_UNRESOLVED`.
- U6: The provided A32 package and the executed `Gaussian/A32` experiment share
  a label but are different objects (theory package vs. E0-FAIL experiment);
  do not merge their statuses.
- U7: This audit did not open any locked/development truth array, so B1 data
  hygiene rules remain exactly as the parent states them; prior exposed
  development results are descriptive only.
