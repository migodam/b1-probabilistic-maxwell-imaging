# Separate explicit anchor-conditioned full66 family

Implemented `src/p3r4/anchor_density.py` and `src/p3r4/train_anchor_density.py`. No campaign training or physics was launched. `ALGEBRA_REPLAY.json` and `CHECKPOINT_REPLAY.json` record zero-update algebra/checkpoint audits. The new format is `p3r4_anchor_full66_density_v1`; raw-y first-family formats, CLI defaults and queued semantics are preserved.

## Conditioning and reference measure

NPZ contains original supervised `x[N,66]`, raw source-major `y[N,512]`, explicit y-derived `context_x[N,66]`, parent IDs, optional `base_mean`, `base_chol`, and cached native whitened `context_J[N,512,66]`. Only these fields are read. Exported material/GT helper fields are unused by this trainer. Normalization is fitted only on TRAIN and saved. All512 raw observation channels and all66 anchor coordinates are retained: default conditioner dimension578. Optional `--j-summary log_column_norm` adds66 columns, yielding644, from `log1p(sum_rows Jx^2)` with `Jx=Jq*diag(1,...,1,pdf(anchor_eta1),pdf(anchor_eta2))`. This is a frozen anchor summary, not a dynamically refreshed Jacobian.

Default residual centering models r=x-base_mean(y), falling back to context_x(y) as the center. The translation determinant is1 because the center is fixed given y. `--no-residual-center` gives the absolute x control. The full raw dx66 log density includes the TRAIN scaling and actual supplied frame determinant. There is no independent56-coordinate prior replacement. Gaussian/spline/triangular still generate all66 coordinates, and triangular models the entire strong conditional. Its default first16 coordinates are just a fixed order; a weak-geometric interpretation requires a legitimate supplied TRAIN frame.

`--frame` accepts a fixed TRAIN-only orthogonal matrix66 as NPY or NPZ key `frame`. Parent can compare identity, cached TRAIN Gram frame and fixed random orthogonal frame with identical training rules. No candidate-dependent frame or unaccounted determinant is allowed.

## Inference

```python
from p3r4.anchor_density import load_anchor_density
p = load_anchor_density('BEST.pt', device='cpu')
x = p.sample(y_raw, context_x, n_samples=128, seed=20261025031,
             batch_size=16, base_mean=None, context_J=None) # [B,128,66]
log_q = p.log_prob(x_raw, y_raw, context_x, base_mean=None,
                  context_J=None, batch_size=16) # [B], raw dx66
```

For CFM the base is a fullrank Gaussian centered at the externally supplied base_mean or context_x, with supplied positive-diagonal lower-triangular base_chol or anchor_scale*I. It is a computational reference, not a posterior covariance claim. Mean and covariance transform consistently through translation/TRAIN normalization/frame. No point-mass base or endpoint smoothing is used. CFM density still requires `allow_approximate=True,nfe=64`, with complete trace; it is not the finite RK4 map density and is unsuitable for exact MH.

The independent probability evaluator supports both checkpoint formats. Anchor-model online files must explicitly contain `context_x`; ambiguous `x` is never read because it can mean a supervised GT label. Native cached Jacobian may be `context_J` or `J`. These allowed context fields are sealed with the online output. Original raw observation files do not automatically supply anchor context. Missing context fails explicitly. Synthetic anchor inference/scoring passed three parents and four method/NFE variants with forbidden label arrays left unread.

## Matched-capacity initial comparison

All three use the same width128 condition encoder, default578 inputs:

| Family | Parameters |
|---|---:|
| spline |691308|
| diagonal Gaussian |691391|
| CFM |691022|

The differences are83 and286 parameters. CFM retains four width128 residual state blocks, with feedforward hidden547 for capacity matching. `--velocity-hidden` controls it explicitly in this new trainer. The raw-y CFM constructor's optional hidden argument defaults to its previous width128 and changes neither its default parameter count nor computations. With J summaries all three gain8448 parameters. Fullcov Gaussian and triangular controls remain available with actual counts reported.

Recommended first comparison: spline/diag_gaussian/cfm, seeds2026102501/2026102502/2026102503, same1024 TRAIN/256 validation cached96+64 contexts, same residual/frame/condition policy, 6000-update cap, validation including update0. Parent owns scheduling, cost receipts and model selection; no fit is launched implicitly. Further frame/metric/LR controls require their own fresh outputs and budget registration.

```text
python -m p3r4.train_anchor_density --train ANCHOR_DATA/density_train.npz --validation ANCHOR_DATA/density_validation.npz --out FRESH_RUN --kind spline --seed 2026102501 --updates 6000 --device mps --width 128 --lr 0.0002 --max-wall 900
```

Use PYTHONPATH=src for `-m p3r4...`, or invoke `-m src.p3r4.train_anchor_density` from repository root. All outputs and update0 selection mechanics mirror the raw trainer, with a separate anchor-specific contract. Existing outputs are refused. No SIGALRM handler is installed; external parent budget remains authoritative. A wall stop waits for the current update/validation to finish. Empirical data/capacity/LR checks are parent-run experiments; the trainer records finite scales, actual capacity and explicit LR without silently enlarging a failed run.

## Optional bounded velocity-error metric loader

`--metric none|maxwell|decoder|random` is implemented for CFM; this worker deliberately does not construct new physical metrics or tune geometry. Non-none requires per-parent NPZ `velocity_metric[N,66,66]` and `--metric-contract JSON` declaring:

```json
{
  "kind":"maxwell",
  "coordinate":"network_standardized_full66_velocity",
  "endpoint_label_dependent":false,
  "TRAIN_only_tuning":true,
  "bounds":[0.5,2.0]
}
```

Every matrix is checked finite, symmetric, and SPD with eigenvalues in[.5,2]. This only weights `(v-v_target)^T W (v-v_target)/66`; velocity is never hard-projected or separately penalized. Matrix values must be computed only from frozen observation/anchor context and fixed TRAIN constants, never hidden current GT endpoints. Metadata cannot prove data provenance; parent-owned generation receipts remain necessary.

The parent can build `.5I+1.5*lambda*(H+lambda I)^-1` with lambda frozen to a positive TRAIN spectral scale, and a fixed orthogonal conjugate preserving eigenvalues. Metrics must use the *same network-standardized coordinates* as the saved normalizer/frame, not raw units. For z=(r F-mu)/scale, column derivative Jz=Jx F^{-T} diag(scale); similarly for real material derivative Dz. H=Jz^T Jz or Dz^T Dz provides the declared normalized geometry. `fit_normalization` and `context_features` are reusable pure helpers for preparation. Zero positive TRAIN eigenvalues or missing provenance must be handled by parent rather than silently inventing lambda. All metric variants retain the exact same network capacity and seed rules.

`ALGEBRA_REPLAY.json` checks both centered/absolute models, with/without J summary, fixed-frame raw-density Jacobians, all five normalized families, fullrank CFM Gaussian reference, sample shapes, identity-W equality to plain FM error, bounded-SPD rejection and finite gradients. Fixtures are interface tests and make no new-parent calibration or posterior qualification claim. Resource receipts include wall/process CPU and peak RSS; physics calls and optimizer updates remain0.
