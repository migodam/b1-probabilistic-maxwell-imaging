# Genuine raw512 field mean / joint Re-Im moment controls

Implemented only new files `src/p3r4/field_moments.py` and `src/p3r4/train_field_moments.py`. Existing active risk/density code and historical checkpoints were not changed. No trainer was invoked, no optimizer update occurred, and no physics was executed. `NUMERICAL_AUDIT.json` records model/label/PSD and historical adapter checks only; these are not new-parent experimental results.

Every NN receives all512 raw source-major measurements. No anchor, Jacobian, residual, inferred attributes, decoder derivative, local feature compression or GT-derived feature enters inference. The dense head jointly outputs all256 ROI locations. Outside the ROI the original complex host is fixed. This is a genuine observation-to-full-field model, not the existing LocalRisk/direct_moment alias.

| Encoder | Mean only | Mean + pointwise joint moments |
|---|---:|---:|
| MLP64, two hidden layers |70272|120192|
| MLP128, two hidden layers |148224|247296|
| Measurement CNN small + raw512 skip |67800|117720|
| Measurement CNN medium + raw512 skip |150160|249232|

CNN input is Re/Im ×32 receivers ×8 sources; three measurement convolutions, a pooled embedding and raw512-to32 skip feed a64/128-wide merge/head. CNN and MLP retain the same raw information and comparable actual capacity. No exact parameter equality is claimed. All output heads initialize to zero residual about the **TRAIN mean**, making update0 a legitimate mean/moment baseline.

CLI (parent runner only):

```text
PYTHONPATH=src python -m p3r4.train_field_moments --train research/p3l_r4_conditional_transport/inputs/joint/train32768.npz --validation research/p3l_r4_conditional_transport/inputs/joint/validation.npz --out FRESH_RUN --kind mlp_small --seed 2026102501 --updates 6000 --device mps --max-wall 900
```

Kinds: `train_mean`, `mlp_small`, `mlp_medium`, `cnn_small`, `cnn_medium`, their `_moments` variants, and `r2_reference`. `--second-moments` also enables the auxiliary heads. `r2_reference` requires `--reference-checkpoint` and performs no new optimization; `--seed` does not change historical weights. The helper loads an old R2 raw-measurement MLP directly, or a self-contained reference wrapper. The existing n4096/mlp_small_s11 reference was replayed **bitwise identically** through the packed-y adapter; old weights/statistics were not refitted or overwritten.

Input accepts current real `y[N,512]`+complex `material[N,576]`, or original R2 complex `observation[N,32,8]`+`chi[N,576]`, with unique disjoint parent IDs. Only if stored materials are absent may Gaussianized x labels be decoded through the unchanged pure P1 decoder; that label source is explicitly recorded. This does not launch a simulation. Default target/geometry come from the existing original P1 configuration and mesh.

TRAIN-only input column normalization and spatial target means are saved. Re/Im mean targets share one scale, so their MSE remains proportional to actual unweighted full576 complex material MSE. Auxiliary Re²/Im²/ReIm channels have individually declared global TRAIN scales. Their weight is fixed by `--moment-loss-weight`; inference/selection does not tune it from test truth. Validation chooses **raw full576 complex mean MSE**, including update0. Early stopping counts validation events (`--patience`); `--max-wall` is a soft between-update guard, and the parent hard budget remains authoritative. LAST optimizer/RNG state, BEST, complete curves, validation predictions and resource totals are saved. Existing outputs are refused.

```python
from p3r4.field_moments import load_field_moments
p = load_field_moments('BEST.pt', device='cpu')
mean = p.predict(y_raw, batch_size=16)             # [N,576], complex; alias predict_full
moments = p.predict_moments(y_raw, batch_size=16)  # optional moment/covariance actions
```

The optional labels are the same parent's Reχ, Imχ, (Reχ)², (Imχ)², Reχ·Imχ. The cross moment is not a product of separately predicted means. These are valid conditional-expectation regression targets even with one joint draw per observation; that single draw does not provide a per-observation posterior variance label.

Moment output saves raw E[Re²], E[Im²], E[ReIm], raw per-pixel real2×2 covariance, its separately named PSD projection, complex marginal variance, total material variance (mean pixel trace), and E[full material squared norm]. PSD projection clips negative eigenvalues of each pointwise block and leaves the mean unchanged. The raw covariance may be indefinite because independently regressed moments and plug-in squared means are inconsistent at finite sample/optimization accuracy. Projection changes the estimate and supplies no unbiasedness/calibration certificate. `COVARIANCE_WARNINGS.json` preserves negative-eigenvalue counts, correction size, label/plug-in bias warnings and mean excursions outside the legal convex rectangle. Raw means remain primary; legal-rectangle clipping is saved as a separate action only.

There is **no cross-pixel covariance estimate, spatial joint law, material density, ambient Gaussian posterior, joint material sampler or credible-interval qualification**. Mean-only and historical checkpoints explicitly report variance NOT_ESTIMATED. Per-pixel PSD is not a full P1 posterior certificate.

Numerical audit: original pure decoder labels; TRAIN statistics and exact mean initialization for all8 NN variants; full576 MSE normalization error2.3e-16; prior sample pointwise moment/covariance identity error3.9e-16; indefinite covariance PSD repair; finite CPU gradients; historical R2 MLP bitwise adapter equivalence. No training CLI was launched. Actual learning/generalization/risk and covariance calibration remain for the parent's bounded experiments.

## Independent evaluator integration

`tools/evaluate_r4_risk.py` now accepts `loader="field"`. An all-field specification must declare `fallback="none"`; mixed specifications retain the existing anchor fallback for their context methods, but field failures never receive an anchor substitute. Field-only inference reads only allowlisted observation/y and packed parent identity through `read_online_measurement`; it does not open q0/chi0/J/residual/GT fields, and does not need an anchor file. It checks frozen checkpoint identities and seals field prediction file identities before offline GT scoring in another PID. Saved `field__*` arrays include raw second moments and raw/PSD variances where available. Offline raw second-moment MSE is a moment-regression score, not a conditional-variance label or uncertainty qualification. Failed parents remain explicit; primary means and paired CIs are NOT_EVALUATED when required predictions are missing.

The existing r2/r4/anchorbin/anchor predictor branches are preserved, including the risk agent's format-aware anchorbin compatibility loader. `EVALUATOR_AUDIT.json` records successful first/second-moment online saving, forbidden object anchor/J/GT fields left unread in field-only mode, separate-PID scoring, explicit no-anchor failure handling, and a mixed anchor branch replayed bitwise identically. Field predictions were identical with and without the synthetic anchor context. No training CLI was run.
