# Untrained normalized conditional-density pilots

Only `src/p3l/neural.py` and this audit directory were written for this task. No campaign data, observations, GT, teacher, physical arrays, ledger, network service or MPS device was read/used. These are independently supervised joint-simulation density architectures, with no distillation or posterior-certification claim. No optimizer, fitting loop or actual training run was executed. Random forward/backward computations validate differentiability only.

## Interface

```python
from p3l.neural import build_model, parameter_count
model = build_model('mlp_full')
logp = model.log_prob(values, scalar, image)
samples = model.sample(count, scalar, image)
size = parameter_count(model)  # model.parameter_count() is also available
```

Defaults: scalar context `(batch,1154)` = measurement512 + residual512 + anchor66 + spectrum64; image `(batch,4,24,24)`; values `(batch,10)` or `(batch,count,10)` =8 supplied fixed Gaussian weak coordinates +2 supplied probit-property coordinates. `log_prob` returns `(batch,)` or `(batch,count)`. Samples return `(batch,count,10)`. Input shape mistakes raise; no label clipping, standardization, probability-coordinate conversion or hidden GT lookup is performed. Caller owns dtype/device and input preprocessing.

The densities are normalized with respect to the exact label coordinates supplied to the module. Parent-owned affine standardization must add the exact inverse-scale log Jacobian when reporting original-coordinate density. Probit properties require the appropriate normal-CDF Jacobian when changing back to bounded-property density; a normalized probit density is not automatically a bounded-coordinate density. A10D weak/property head does not supply the missing58 strong Gaussian conditional coordinates or certify a full66D posterior.

## Models and parameter counts

| Model | Density head | Trainable parameters |
|---|---|---:|
| mlp_diag | diagonal Gaussian |461,844|
| mlp_full | full Gaussian |467,649|
| mlp_mdn |4-component full Gaussian mixture |493,320|
| mlp_flow |4 conditional affine RealNVP couplings |606,800|
| cnn_full | full Gaussian |467,441|

All heads use128 context features. MLP-only receives both scalar context and flattened image through `(1154+2304)→128→128` SiLU layers. CNN uses4→16→32→64 stride2 convolutions,576→128 image projection; scalar1154→256→128; fusion256→128. Full-Gaussian CNN/MLP counts differ by208, approximately0.0445%; both see the same inputs, so the distinction concerns architecture rather than removing image information. Other head families have their explicit additional parameter costs.

Gaussian diagonal scales and packed full-Cholesky diagonal entries use `softplus(raw)+.02`. Off-diagonal Cholesky entries are unconstrained. Density evaluation uses normalized multivariate normal formulas. Mixture weights use log-softmax/log-sum-exp; sampling draws a component per context/draw, then its Gaussian. Components carry no physical-mode semantics.

The RealNVP base is normalized standard Gaussian. Four alternating parity masks cover all coordinates. Each affine coupling conditions on the unchanged masked coordinates and context; transformed coordinates use `log_scale=2*tanh(raw_scale)` and unrestricted shift. Forward log determinant is the sum of active log scales; inverse uses the same unchanged coordinates and negative log determinant in reverse layer order. `log_prob=log p_base(inverse(values))+inverse_logdet`. The module exposes `forward_from_base` and `inverse_to_base` for audit; no learned normalizer or approximate ODE density is involved. Individual coupling scales lie in `(exp(-2),exp(2))`.

## CPU validation

`check_neural.py`, `results.json`, `check.stdout`, and empty `check.stderr` record14 controls. Seed2026102252, Torch2.14.0, CPU float64, one thread. Actual training runs=0, campaign data reads=0, physical calls=0.

- All five models: sample shape3×7×10, correct log-density shape, finite density and finite loss gradients.
- Independent SciPy normal/mixture log-density comparisons in2D: maximum errors `1.42e-14` (diagonal/full), `2.84e-14` (mixture). Direct201×201 grid integration on[-10,10]² gives `1.0000000000000002` in each case.
- Random conditional10D flow: inverse error `4.44e-16`, forward/inverse logdet cancellation `4.44e-16`.
- Full automatic-differentiation Jacobian determinant versus coupling formula: error `8.88e-16`.
- Change-of-variables density identity: error `3.55e-15`.
- Nonidentity2D affine flow versus independent SciPy Gaussian: maximum log-density error `8.53e-14`, finite-grid integral `0.9999999999999165`.
- Extreme negative Cholesky logits preserve the `.02` floor.

Normalization integrals are low-dimensional analytic controls; the10D normalized construction follows its Gaussian/mixture/bijective architecture. No check establishes useful conditional learning, calibration, material risk improvement, production tail adequacy, generalization or scientific posterior quality. In particular, a small positive Cholesky floor does not prove coverage.

Reproduce from worktree root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ../../.venv_nn/bin/python \
  research/p3l_local_posterior_imaging/AUDIT/neural/check_neural.py
```

Parent owns independent-parent splits, context and label standardization, all Jacobians, joint-simulation supervision, train/validation decisions, budgets, training and final scientific interpretation.
