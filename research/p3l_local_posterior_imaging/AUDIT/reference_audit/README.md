# Four-chain reference audit

Implemented only `src/p3l/reference_audit.py` outside this audit directory. Callable:

```python
from p3l.reference_audit import audit_reference
stats = audit_reference(four_chain_paths, original_prior, warmup=1024,
                        output_dir=new_audit_directory)
```

Each input must be a distinct explicit `CHAIN.npz` containing `q(draw,66)` and `Phi(draw)`; no other input files are opened. Saved warmup must agree when present. Exactly four chains, equal retained lengths and at least four production draws are required. Rejection repeats remain; the audit never thins, deduplicates or truncates to the shortest chain. Explicit FULL metadata is required for the within-mechanism gate. Target families cannot be mixed. Original prior public configuration/geometry are supplied by the caller; no backend, book, ledger or truth module is imported.

Material decoding is vectorized in batches of1024, using `field_mean+qz@prior.A.T`, the original C2 saturation, original mask, host and property intervals. It computes all576 complex pixel samples, not independent pixel approximations. Current original P1's256 ROI requirement is checked. A check against17 supplied states decoded through the original per-state method has maximum difference `3.5580223032530555e-15`; the observed BLAS summation difference slightly exceeded an initial `3e-15` tolerance, so final comparison uses `1e-14`. No physical calculation was involved.

The fixed material-function manifest contains52 real functions:

- ROI real/imaginary means and real/imaginary mean squares (4).
-16 fixed sorted ROI pixel indices chosen by rounded evenly spaced index positions, both components (32).
-8 unit-Euclidean-normalized Gaussian ROI projection weights from seed2026102250, both components (16).

Indices and weights are stored in NPZ. Material functions, all64 Gaussian coordinates, both bounded properties and Phi receive existing `diagnostics.py` rank/folded split Rhat, bulk/tail/mean ESS, mean MCSE, variance ESS and variance MCSE. Pixel diagnostics are also preserved in blocks of32 pixels to bound diagnostic working memory. Full-grid complex mean/empirical modulus variance and separate real/imaginary 2.5/50/97.5percent quantiles are saved. Centered accumulation preserves deterministic sample means, including fixed host pixels. Quantiles are marginal component intervals, not joint complex or simultaneous spatial credible regions.

Frozen within-mechanism thresholds:

1. All named latent/property/Phi/material-function Rhat values are finite and `<1.01`, with assessed status.
2. All ROI pixel component Rhat values are finite and `<1.01`; all ROI components are assessed. Constant ROI components stay `CONSTANT_UNASSESSED` and block the gate.
3. All52 material-function bulk AND tail ESS values are finite and `>=1000`.
4. Full-grid RMS mean MCSE, `sqrt(sum_ROI(mcse_Re²+mcse_Im²)/576)`, is finite and `<=.001`.
5. All chain targets explicitly say FULL.

Outside ROI,320pixels/640components are known deterministic host and explicitly counted. They are excluded only from mixing gates; their constant diagnostic statuses remain visible. Zero MCSE in an unassessed ROI never makes the overall audit pass. The extra ROI pixel Rhat gate is conservative relative to a material-function-only gate; no lower threshold or silent constant exemption is applied.

Even if all gates pass, `reference_qualified` remains false and status is `WITHIN_MECHANISM_GATES_PASS_CROSSCHECK_REQUIRED`. Separate independently seeded mechanism agreement remains required. Otherwise status is `UNQUALIFIED_WITHIN_MECHANISM_DIAGNOSTICS`. The audit does not infer seed independence, target implementation correctness, global mixing, escape mass, GT risk or stationarity from a file.

Outputs:

- `REFERENCE_STATS.json`: all gate decisions, exact thresholds, input paths/counts, named and full-pixel diagnostics, constant statuses. Nonfinite numbers become JSON null, with status fields retained.
- `REFERENCE_ARRAYS.npz`: material mean/variance/quantiles, production q/Phi, all52 material-function traces/names, fixed projection weights/indices/seed and every numeric/status diagnostic array. Raw full material trace is reconstructible from retained q and the original prior, avoiding another large archive.

The function holds decoded material draws in memory: for4x4096draws this is about151MB, with vectorized batch and diagnostic working arrays in addition. NPZ excludes LU/current/field and does not touch funded chain files. MCSE is the existing asymptotic plug-in estimator; it is not a rigorous bound.

`check_reference_audit.py` and `results.json` record11 supplied-chain controls. Separate synthetic IID and constant-chain inputs are retained inside this audit directory; none is a funded Maxwell posterior. Checks cover vectorized decoding, full-array moments/quantiles, fixed manifest, warmup, cross-mechanism hold, host constants, strictJSON, ROI constant failure, null diagnostics, chain count and warmup mismatch.

Reproduction from the worktree root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ../../.venv_nn/bin/python \
  research/p3l_local_posterior_imaging/AUDIT/reference_audit/check_reference_audit.py
```

No new physics, RHS, network access, GT arrays or old ledgers were used. Parent owns funded execution and final qualification.
