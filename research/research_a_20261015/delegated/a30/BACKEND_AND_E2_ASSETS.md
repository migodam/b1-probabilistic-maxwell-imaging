# A30 backend and E2 implementation handoff

2026-10-15. Read-only audit of A30; only the Research A worktree received new artifacts. No Maxwell constructor, new state, new LU, or new RHS was executed by this worker. No truth values or locked payload were read. The only imported code was definitions, plus existing online case metadata/paid prediction arrays. One NumPy private-header-reader error was repaired with the public API; its cost was zero physics. Historical tests below were read, not rerun.

## Source and physical identity

A30 checkout: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A30/a30-maxwell-material-manifold`.

Observed clean HEAD `717d3d26a290248bd2b9aeb55f3d03694b407b84`. Its self-contained frozen A29 vendor commit is `5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4`; A28 source is `c7819fa142bb70f9920835e8e91ba145f15ed8e3`. These are Git identities, not file digests. No hashes were checked. Original decisions, state ledger, results and truth were not changed.

The backend is **2D scalar TM outgoing Helmholtz VIE**, not three-dimensional vector Maxwell. Frozen geometry:24×24 cells on a unit square,8 sources on radius.75,32 receivers on radius.8, full aperture. Time convention `exp(-i omega t)`; outgoing Green `i H0^(1)/4`; saved Green includes `k²` and the cell quadrature. Original frequency is `C0`Hz, so free-space `k=2π` in unit coordinates. Use the actual saved G/S/incident rather than rebuilding a nominally similar kernel.

The material is complex susceptibility χ with **linear constitutive map a(χ)=χ**, scalar shared across all8 illuminations. It is not the Clausius-Mossotti polarizability used by some 3D dipole backends. Actual box:Reχ∈[−.5,4], Imχ∈[0,.5], numerical tolerance1e−10. An object-background increment can have negative imaginary part while the full material remains passive. Current and exciting-field normalization is `c=sqrt(v)j`, `e=sqrt(v)E`, with v=1/576. The discrete equations are

```
Lχ = I - diag(χ)G
Lχ c = diag(χ)incident
e = incident + Gc
y = Sc
c = diag(χ)e
```

Material mass coordinates are `sqrt(v)*[Reχ576,Imχ576]`. **Data packing is per-source [Re32,Im32], source0 through source7**; it is512 real coordinates, unlike B1's all-real-then-all-imaginary packing. See `vendor/frozen_a29/src/a24/contracts.py:104` and `src/a30/compensation.py:21`.

Saved sigma is the **complex** noise RMS:independent Re/Im variance σ²/2. Proper complex log likelihood uses `-||F-y||²/σ²`; whiten a real packed residual with **sqrt(2)/σ**. Reuse each case's saved sigma and observation. Do not candidate-rescale, inflate noise, redraw measurement noise, or import B1's .12 real-channel convention. Original A30 has **known calibration**, so E2 keeps gain1. Introducing free gain is a new target, not a free corrector.

## Assets and data separation

The extracted A29 public inputs exist at `A30_ROOT/inputs/a29_release`. No download is necessary. Historical release names/size receipts are preserved in `research/delegated/a30_inputs/ACQUISITION.json`; no new network access was used.

|Asset relative to extracted release|Fields usable online|Separation|
|---|---|---|
|`inputs/a25/background/weak.npz`, `ring.npz`|χb, eb, yb, Db, G, S, incident, volume/grid/geometry descriptors|Saved known background; no new solve on import|
|`results/a29/pools/<case>.npz`|observation(32,8), sigma, n3(576), background, existing candidate fields|No truth field; metadata/method artifacts are exposed historical|
|`inputs/a28/results/a28/exact_states/<case>.npz`|paid N3 chi/mass/current/field/data/LU/pivots|Whitelist exact same material and background; missing/stale cache refuses|
|`results/a29/OFFLINE_TARGETS.npz`|only case/parent/family/background metadata online|Truth(32,576) **offline scoring only**|
|`inputs/a25/train.npz`|TRAIN labels only to freeze representations if needed|64 parents250001–250064; not development truth|
|`inputs/a25/validation.npz`|contains observation metadata **and truth/current/clean**|Avoid broad dictionary load in a proposal; use sanitized pool fields|
|`inputs/a28/results/a28/anatomy/<case>.npz`|contains both CNN and truth-centered oracles|Avoid wholesale online loading|

All16 parents250065–250080 have weak/ring paired cases at a1/SNR20. These32 conditions are historically exposed development, not fresh validation. Family blocks are shell250065–68, asymmetric_piecewise250069–72, smooth250073–76, mixed_scale250077–80. E2 runner freezes the **first parent of each family**, weak then ring:250065/250069/250073/250077. This is chosen by family/index before truth values are read, not by reported success. All8 required same-state caches exist. All32 pool counts58–91; native E2 can operate without trusting their old hypothesis probabilities.

`ONLINE_ASSET_AUDIT.json` records32 online metadata entries, sigma, shapes, cache existence and successful definition-only imports under the existing `Gaussian/.venv_nn/bin/python`. No package install is necessary. A30's unused helper `data.source_root` references an absent protocol path; avoid that helper. `release_root`, explicit pool/background/cache paths work.

## Reusable runtime interface and exact cost

Import `a30.physics.FrozenMaxwellAdapter` with `A30_ROOT/src` added to PYTHONPATH. It injects only its own frozen vendor `a20/a24/a25/a27` modules and refuses preloaded foreign modules of those names. Start a clean subprocess if another study has imported different A24 modules. Importing definitions does not initialize geometry.

**Do not instantiate `a30.budget.Book` for Research A.** It writes the historical ledger and hardcodes the old192-state budget and old E0 gate. The adapter accepts any new isolated book with `counts`, `span`, `reserve_state`, `reserve_rhs`. Research A's `src/research_a/e2.py` provides one, fixed to384states/9216RHS/3600executionseconds and append-only accounting. Old A30's161 calls remain closed.

|Method|New state attempts|Actual complex LU RHS|Other paid/timed work|
|---|---:|---:|---|
|`from_background_npz`|0|0|Disk read; saved equation/descriptor checks; dense matrix products; no geometry/LU|
|`cached_anchor(chi,path)`|0|8 perload|Saved packed-LU import and same-material equations; LU source audit|
|`full_state(chi)`|1 perattempt|8|One576×576complex LU, matrix construction, G/S/source8 applications and backward residual|
|`anchor_jvp(U,state)`, m complex directions|0|8m|Existing-state LU solves and S products|
|`receiver_transfer(state)` firstuse|0|32adjoint|Computes `S Lχ^-1`; same-state reuse0RHS|
|`cached_material_dot_identity`|0|0 if transfer ready, otherwise32|Raw/material-mass real-adjoint identity|
|`solve_anchor(rhs,state)`|0|number of rhs columns|Existing-state forward/conjugate-adjoint LU and residual audit|
|`L_response(chi)`|0|0|Saved affine `yb+Db((χ−χb)eb)`; not a full-Maxwell state|
|Compact J from transfer/field/basis|0|0 after transfer|`einsum(T,U,e)`; count material directions and wall separately|

Torch/LAPACK cache pivots are one-based swaps; the adapter converts to SciPy by subtracting1 and audits the actual saved LU. It does not create a replacement if a cache is missing. Historical backend tests002 completed with0actual Maxwell states; that historical arithmetic result is not an E2 real-state check.

## What A30 actually failed

The decisive records are `A30_E0_E1_DECISION.md`, `results/a30/FINAL_DECISION.json`, `results/a30/compensation_corrected/{E1_SUMMARY.json,COMPENSATION_PAIRS.csv,PAIR_CHOICES.json}`, and `results/a30/FAILURE_LEDGER.jsonl`.

E0 offline capacity passed:full convex hull improvement13.12%; sparse64 achieved35.69%; TRAIN-PCA64 achieved21.91%. These are **truth-fitted offline capacities**, not online reconstruction. Real measurement L-weighted mean worsened3.32%, ESSmedian1. Atlas/levelset nonconvex failures mix representation and optimization limitations.

E1 used a **single N3-anchor tangent**:four online magnitude regions, independent complex contrast coordinates8real, boundary-gradient coordinates4real, both endpoints12real each. Contrast radius.30 and boundary displacement.35pixel; the actual hole interior was protected. SLSQP optimized the frozen linearized response, after which endpoints were independently full-solved. This was **not iterative nonlinear GN/LM**, so sequential fresh linearizations are a meaningful falsification change, not a rerun of the same optimizer.

29/32 legal pairs completed;22/24 complex conditions. Missing old edits:p250069_weak,p250072_ring,p250074_ring. Median compensated full distances:shell19.18, asymmetric20.92, mixed20.21, smooth11.79. Global best8.84, complex best9.18. **No parent crossed nominal≥3 to compensated≤1**;0/6 required complex-parent witnesses;0 useful uniform bounds; old E2 NOT_RUN. Local optima are not global impossibility results. F/L solvers each hit100iterations6times, jointly9conditions. The q≥1 Neumann sufficient bound failing does not mean the Maxwell system is physically singular.

First implementation protected every phase-projection changed cell, leaving zero compensation columns.15cases×5calls=75charged;15nominal remained valid,60were redundant. Never delete those failures. Corrected followup reused valid nominal values; extra L endpoints were full-checked only at four preregistered first parents. Total161new states,3048complexRHS,146.97s charged-at-report-generation (not all interactive elapsed time). Original A29 NO-GO and A30 NO-GO remain unchanged.

Original nominal edit included a global phase projection, so it was not a pure causal cavity intervention. The new E2 runner uses direct object-removal-to-known-background edits with protected actual cells and no such phase projection.

## Implemented Research A E2 runner; not executed by worker

File:`src/research_a/e2.py`. Entry:

```python
from research_a.e2 import run
run(B1_RESEARCH_A_ROOT, e1_gate_path=E1_DECISION_JSON, a30_root=A30_ROOT)
```

Parent's JSON must explicitly contain `"E2_allowed": true`; the worker has not made that scientific entry decision. No force flag exists. Completed/failed RESULTS prevents accidental repeat; original ledger is never used. A clean new subprocess is recommended.

The runner freezes8conditions, no free gain, σ unchanged,95% χ²512 residual radius. Cached N3 residuals are28.0347/32.9870(shell),23.2129/31.7447(asymmetric),23.0519/37.8959(smooth),24.6407/31.5446(mixed), weak/ring order. An unfit N3/base is **not counted as a valid archive member**. A shared single constrained material-GN update improves the coarse base when it helps; no truth determines this choice.

Three corrected arms each get4trajectories×3fullstates=12states/condition. Each trajectory evaluates its actual finite proposal then performs two fresh nonlinear relinearization/correction attempts. Round1 uses16real DCT coefficients; round2 expands to32real multiscale Haar coefficients. Material field support is global; A1 protects only its edited hole, so neither support nor compensation is fixed to truth. Box-constrained material QPs preserve passive full fields; rejected nonlinear corrections remain charged. SLSQP statuses/iteration caps are logged.

1. `finite_history`:holes of1.5/2.2pixelradius at online material-based sites, parent selected from previously full-valid archive if available; actual edit cells are held at χb during compensation. This is finite-history search, **not a learned policy**.
2. `multistart`:independent global multiscale material starts under the same materialbox, followed by the identical corrector/budget.
3. `weak_null`:SVD directions from paid N3 material tangent (declared anchor, not falsely current-state exact after base changes), followed by the identical corrector. Pixel-box projection defects are recorded; after projection this is not claimed an exact null direction.
4. `same_edit_no_compensation`:the exact4A1initial proposals, independently full-evaluated with no correction. It costs4calls/condition, not12gratuitous repeats; comparison uses matched proposals and actual cost.

All accepted trial checkpoints may enter archive only after full residual and constitutive checks. Material-field RMS dedup is.03 relative to1.2; structural masks are separately reported with fixedonline |χ−χb|>.18 and≥6changedcells. Different ε amplitudes do not automatically count as different structures. Online TopK is farthest-material coverage starting near the base; data-selected candidate is minimum full residual. Best-of-K uses truth only after every search is finished. Full Re/Im/complexMAE, NMAE, ROI/cavityMAE, IoU/boundary/holes and Coverage@K diagnostics are then computed offline. Probability is null/candidate-only.

One newlycharged3-cell actual-backend state checks exact finite current/data secant, nonlinear remainder,3×3Schur, dual correction, real mass adjoint, and `current=chi*field`. This costs1state/51complexRHS and is distinct from the package's synthetic8dipole algebra. It is invoked only after conditional E2 entry. All failures stop while retaining reservations and partial arrays.

Planned worst case is **329states/9139RHS**, below384/9216. This includes8N3LU audits, all32receiver adjoints per needed anchor,8coarse endpoints,288corrected-arm states,32no-compensation states and1audit. Existing paid caches are shared setup, not newly billed material states; their fresh LU audits and timed loads are billed. Wall guard is cumulative over recorded executions; a previously interrupted execution refuses automatic restart until accounting repair. Runtime if nearing limit stops and preserves partial results.55states and77RHS intentionally remain; no forced budget spend.

## Parent review limits

The worker did not run E1, decide its gate, run E2, find a new solution, establish method superiority, or claim a paper contribution. `e2.py` was syntax-parsed; actual physics/interface identities await its one charged backend audit. It is a deliberately small falsification panel, not all32A30conditions and not a new holdout. Its two correction steps per trajectory may remain optimization-limited. A null result must retain that limitation. Existing A30 local≤1 pair objective and new y-compatible archive objective are distinct metrics; neither can silently replace the other.

The runner emits `results/E2/{PROTOCOL.json,ACTUAL_BACKEND_AUDIT.json,SEARCH_HISTORY.json,MATERIAL_AND_PREDICTION_ARRAYS.npz,METRICS.json,COST_LEDGER.jsonl,RESULTS.json}`. Parent remains responsible for final scientific decision, whole-iteration comparisons, uncertainty/statistics, gallery, native A package and versioned A→B adapter.
