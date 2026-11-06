# Frozen R4 GN-family endpoint asset map

Read-only discovery; no inverse solves, SVDs, source mutations, private unscored truth reads, SHA calculations, or network. Paths below are relative to `/Volumes/migodam's-external-brain/Research/Inv_SLAM`. Let `P = research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r4_diversity`.

## Primary cohort and exact loader contract

Recommended cohort: `method == native_trf`, `backend == cpu`, 72 endpoints = 12 shapes × 3 amplitudes × 2 noise levels. These are **GN-family bounded trust-region reflective nonlinear least-squares endpoints, not pure unconstrained GN iterates**. The implementation uses analytic Jacobians, common regularization, and 8→12→24 material stages. Exported endpoint is the online top1 selected by normalized full data objective plus common generic regularizer. `reciprocal_trf` is a separate coordinate-system control; its exported `m` is ordinary material, not reciprocal coefficients.

```python
from pathlib import Path
import json
import numpy as np
P = Path('research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r4_diversity')
rows = {r['task_id']: r for r in json.loads((P/'RESULTS_SUMMARY.json').read_text())['rows']}
with np.load(P/'predictions/endpoints.npz', allow_pickle=False) as z:
    ii = np.flatnonzero(z['finite'] & (z['method']=='native_trf') & (z['backend']=='cpu'))
    cases = [(int(i), z['m'][i].copy(), rows[str(z['task_id'][i])]) for i in ii]
with np.load(P/'public_data/operators.npz', allow_pickle=False) as z:
    G_internal, S_receiver, E_incident = z['G'].copy(), z['S'].copy(), z['E0'].copy()
for endpoint_index, mhat, row in cases:
    with np.load(P/'truth_after_freeze'/(row['observation_id']+'.npz'), allow_pickle=False) as z:
        m = z['m'].copy()
    material_error = mhat - m
```

Join metadata by `task_id`. `RESULTS_SUMMARY.json['rows']` supplies shape/kind/font, amplitude, background, noise/noise_relative, backend/method, residual, solver disposition, and old offline material metrics. All 360 task IDs, observation IDs, methods, backends and finite flags match the endpoint archive. Every released truth reference exists. `ASSET_INDEX.json` lists every endpoint index and metadata row.

## Keys, dimensions and conventions

| File | Keys and shapes |
|---|---|
| `P/predictions/endpoints.npz` | `m float64 (360,576)`; `finite bool (360,)`; `task_id <U40 (360,)`; `observation_id <U20 (360,)`; `method <U14 (360,)`; `backend <U4 (360,)` |
| `P/public_data/operators.npz` | `G complex128 (576,576)`; `S/H complex128 (48,576)`; `E0/U0 complex128 (576,12)`; centers float64 (576,2), sources (12,2), receivers (48,2); shape/grid_shape int64 (2,); scalar grid_spacing, radius, k; material_bounds float64 (2,) |
| `P/truth_after_freeze/<observation_id>.npz` | `m float64 (576,)`; `alpha float64 (24,24)`; `mask bool (24,24)`; `Y_clean complex128 (48,12)` |
| `P/public_data/<observation_id>.npz` | `Y complex128 (48,12)`; scalar Unicode observation_id |

All 360 endpoints are finite. CPU 288 and CUDA 72: each of native_trf, reciprocal_trf, data_lbfgs, globaldual has 72 CPU and 18 CUDA results. Do not count backend repeats as independent targets. Shape labels: digit1–digit9, ring, two_islands, curved_bar. Amplitudes: 0.4, 2.5, 12.0; backgrounds: 0.08–0.15; noise: clean 0, noisy relative 0.001. Source: `P/configs/CONTRACT.json`.

**Truth alpha is a normalized shape amplitude map**, not ordinary material and not reciprocal material. Across all 72 released truth files, verified `m = background + amplitude * alpha.ravel(order='C')` (maximum absolute difference recorded in JSON). Load `m` directly. Both endpoint/truth m and Green columns use the same flattened cell order.

## Physics and Green identification

Executable authority: `P/src/physics.py:1-29,82-101`:

`j = diag(m) (E0 + G j)`, `E = E0 + G j`, `y = S j`.

Material m is real scalar contrast; induced current j and fields are complex. The retained solve uses `A = I - diag(m) G`, RHS `diag(m) E0`. Twelve source illuminations are columns. G maps 576-cell current to internal field; S maps the same current to 48 receivers. Hence **G_internal = operators['G']; S_receiver = operators['S']**. No incident-field/material dressing is built into these Green assets. `H` aliases S, `U0` aliases E0; prefer canonical keys. The stored matrices already include disk-pulse integration/self term; do not multiply by an extra area/k² factor. Exact assembly is `physics.py:16-28`.

`P/src/solver.py:135-160` gives the analytic residual Jacobian and SciPy `least_squares(method='trf', tr_solver='exact', bounds=..., x_scale='jac')`. Reciprocal variables have the chain `-m²`; exported m remains ordinary material. `P/export_replay.py:18-31` exports `result['m']` and metadata; `P/score_results.py:102-126` independently scores these same fields. `P/PROVENANCE.md` identifies fixed R3 operator ancestry and retained physical model.

## Released truth provenance

`P/truth_after_freeze/README.md` explicitly discloses these synthetic materials after planned prediction closure. The historical basename `mapping_PRIVATE.json` **inside this public released-truth directory** has status `PUBLIC_SYNTHETIC_TRUTH_DISCLOSED_AFTER_PREDICTION_CLOSURE`; export strips generator entropy. `P/audit/ALL_PREDICTIONS_SEALED.json` is `PLANNED_RESULTS_CLOSED`; `RESULTS_SUMMARY.json` states truth_disclosed_only_after_seal. No `private_truth` data were read.

## Established local Twofold definitions and evidence limits

No literal mathematical `Gin/Gob` definitions were found in targeted project theory/code searches. Local explicit notation uses `G_D` (current → internal field) and `G_S`/sensing S (current → receiver). The R4 G/S mapping above follows domain/codomain; it is not a claim R4 implements classical TSOM.

- `Theory/Questions/A2.md:2332-2359`: existing local account of author's monograph printed pp.161–163, Eq.(6.70)–(6.71). V_D+ consists of leading right singular vectors of G_D. Implementable approximation is `P_S- V_D+ beta`, distinguished from exact intersection.
- `Theory/Questions/A3.md:733-747`: `j- ≈ (I - V_S+ V_S+*) V_D+ beta`; Eq.(6.71) belongs to the 2018 monograph, not the unverified 2009 original-paper PDF.
- `research/delegated/a3_rom/basis.py:1-20,104-171`: executable sensing chart uses leading right singular vectors of stacked sensing rows. twofold_chart constructs sequential `[V_S, (I-P_S)V_D]` with V_D from right singular vectors of domain_operator. This is a different geometry/model study; reuse its declared definition, not its matrices for R4.
- `research/delegated/a3_prior/EVIDENCE.md:153-179` retains original-paper identity Tier C and monograph equation attribution as prior-pass evidence. No fresh primary-literature verification was done here.

Projecting **material error** into Green right-singular directions is a parent-chosen diagnostic. Do not silently equate material error with current error or label the diagnostic a native SOM/TSOM implementation. Parent owns ranks, thresholds, mathematics, and interpretation. R3 fallback discovery is unnecessary because complete R4 released assets exist.
