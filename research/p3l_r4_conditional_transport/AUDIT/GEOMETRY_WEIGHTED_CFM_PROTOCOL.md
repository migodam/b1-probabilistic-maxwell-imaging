# Bounded geometry-weighted CFM: prepared, unreleased

New module `src/p3r4/train_geometry_cfm.py`; six-fit controller `tools/queue_r4_geometry_cfm.py`. Registry `configs/GEOMETRY_WEIGHTED_CFM_V1.json` remains `released:false`. All VI-frozen source identities were verified unchanged. No optimizer update, MPS allocation, child fit or new Maxwell query occurred in preparation. This is a separate training-loss experiment, not another GN coordinate transform.

## Pre-results coordinate clarification

The original whitened native Jacobian is converted to Gaussianized coordinates by `Jx=Jq*probitDerivative(anchor_x)` on its two property columns. Let `rawx=center+T*r`, where r is the actual TRAIN-standardized residual learned by baseline CFM and `T=frame_inverse.T @diag(TRAIN_x_scale)`. Build **Jr=Jx T**, not a raw-unit Jacobian mislabeled as a standardized metric. SVD `Jr=U diag(s) V.T` gives

`W_r=I+3 V diag(s²/(1+s²)) V.T`.

Stable SVD avoids negative Gram-roundoff eigenvalues; the implementation uses `(s/hypot(s,1))²`. All66 network velocity-error components remain present. Metric eigenvalues lie in[1,4], so no velocity direction is removed. The corresponding raw metric is `T^{-T}W_r T^{-1}`, bounded **relative to** the baseline raw metric `T^{-T}T^{-1}`. No raw-Euclidean[1,4] eigenvalue claim is made. This clarification was supplied by the parent before preprocessing/results.

The random control is `R W_r R.T` in the same network coordinates, with one fixed TRAIN-independent orthogonal R, seed20261025191. Every parent's metric spectrum is preserved. The conditioner, computational base, data normalization and coordinate map do not change. Only velocity ERROR is weighted: `E[error.T W_r(y) error]/66`; velocity itself is neither projected nor penalized toward zero. An SPD y-only metric preserves the ideal conditional mean-velocity optimum. That statement does not prove finite-model calibration or improvement.

## Frozen cache and cost

`inputs/geometry_weighted_cfm_v1/` contains TRAIN/VAL geometry+random float32 matrices, shrinkage spectra, original baseline normalization and fixed rotation. `CACHE_COMPLETE.json` records byte/mtime identities, original input sources, basis/formula, spectrum checks, and preprocessing cost. The preprocessing reads only parent/context_x/context_J plus the existing baseline's global TRAIN statistics; it never reads endpoint x/material labels, decoder features or GT. TRAIN and VAL metrics use their own observation-derived anchors/J, with no VAL tuning. The cache is shared by all seeds.

1024 TRAIN and256 VAL metrics were actually prepared. Double-precision geometry eigenspectra lie approximately[1,3.999999]; random-spectrum discrepancy<=1.12e-14. Float32 cached eigenbounds are checked with numerical tolerance at fit load. Metric preprocessing cost:2.384s wall,2.251s CPU, zero physical RHS. The full analysis-only preparation/audit Book cost:2.481s wall,2.312s CPU, peakRSS495927296B. Existing cached anchor/J acquisition costs remain owned by the prior strong-context pipeline; they are not claimed free or newly acquired here.

Readiness: `AUDIT/geometry_weighted_cfm_prepare_1792879193757865000/READINESS.json`, PASS_PREPARED_INERT. Identity-weight loss matches unweighted baseline loss exactly in the mock; weighted-error gradient discrepancy1.39e-17; fixed-y conditional-mean gradient residual1.31e-15. Same-seed fresh model initialization is bitwise replayed; all three original baseline normalizers match bitwise. All VI-frozen source identities remained unchanged. Initial preparation invocation had a corrected duplicate-key syntax error before any Book/physics; the successful audit retains the actual bounded preparation cost. A scalar-reporting Torch warning in the audit does not affect its gradient or stored numerical results.

## Six matched fresh fits, not launched

Two arms (geometry, random alignment) × seeds2026102501/2/3. Same original anchorCFM_none architecture and initialization algorithm:691022 parameters, full66, raw512+anchor66 conditioner, width128, velocity hidden547, four blocks. Fresh per-seed initialization, base/noise/minibatch streams, full-rank raw anchor Gaussian computational base, residual translation and TRAIN normalizers match the original baseline. Original completed checkpoint weights are not used to give the new arm an extra3000-step continuation.

TRAIN1024 / VAL256 only.3000 updates, LR5e-5, batch128, AdamW decay1e-4, gradient clip5. Same validation as baseline: raw66 latent energy, K4/NFE16/every100, validation batch128, fixed generator seed per checkpoint, update0 eligible. Weighted FM loss is a diagnostic, never the selection score. TRAIN statistics are replay-checked; validation never fits statistics. Logs, normalization, BEST/LAST and random states use the existing anchor full66 CFM inference format. No new J or W is needed at deployment; the conditioner stays rawy+anchor, without J summaries.

Parent must review and set the new registry `released:true`. Both controller and trainer additionally require the paired V2 VI controller to be CLOSED COMPLETE; they do not block waiting. Every fit uses the existing StageI training-MPS Book wrapper,600s outer quota,540s inner guard,2GiB sampled RSS policy, zero physical actions permitted. Controller stops on partial/failed fit and never retries or expands. Maximum6 fits/3600s outer reservation; no launch occurred here.

```text
<research-python> tools/queue_r4_geometry_cfm.py
<research-python> tools/queue_r4_geometry_cfm.py --execute
```

Outputs will be under `models/geometry_weighted_cfm_v1`. Existing raw/anchor/GN/VI models and prior reports remain untouched. CFM density remains numerical and unqualified; this weighting does not fix the earlier failed trace-agreement gate, establish an exact MH proposal, or qualify a Full posterior. Future DEV/Full comparisons need a separate parent-frozen evaluation release.
