# Full66 conditional density implementation

Status: numerical implementation audit PASS; posterior fidelity NOT_ESTABLISHED.
This worker ran zero Maxwell calls, zero optimizer updates, and no campaign fit.
The saved synthetic update0 checkpoints are interface audits, not fitted research models.

## Files and runtime

- `src/p3r4/density.py`: self-contained PyTorch RQ spline, full covariance Gaussian, matched-capacity diagonal Gaussian, triangular spline, conditional flow matching and saved-model inference.
- `src/p3r4/train_density.py`: one-fit CLI with TRAIN-only normalization, validation-selected checkpoint, update0 eligibility and bounded batches.
- Existing runtime: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`, PyTorch 2.14.0. No dependency added; no nflows/zuko/normflows installed.
- Read execution prompt, main probability algorithms and theoretical conditions, all Appendix B, R2 density/data conventions, model comparison and non-Gaussian negative results. Core scientific interpretation remains parent-owned.

## Input and inference contract

NPZ contains raw real Gaussianized `x[N,66]` and raw real `y[N,512]`. The observation packing is source-major: each source's Re32, then Im32. Optional `parent`, `parents`, or `parent_id` arrays enable uniqueness/disjointness checks; missing IDs are explicitly NOT_CHECKABLE. No test/locked file may be passed to this trainer. Caller owns provenance of supplied NPZ contents.

TRAIN fits all normalizers. Optional `--frame` is a frozen TRAIN-only orthogonal matrix66, supplied as NPY or NPZ key `frame`. Coordinates are `(x @ frame - train_mean) / train_scale`. Raw-coordinate density includes both the normalizer determinant and actual frame determinant. Sampling uses the stored inverse frame. No frame may depend on a candidate x. The CLI receives all raw y512 channels, through a shared architecture for the condition encoder. There is no PCA or omitted complementary latent block.

`load_density('BEST.pt', device='cpu')` returns a predictor:

```python
predictor.sample(y_raw, n_samples=32, seed=11, batch_size=128) # [B,32,66]
predictor.log_prob(x_raw, y_raw, batch_size=128)             # [B], raw dx66
```

Explicit density kinds are `gaussian`, `diag_gaussian`, `spline`, `triangular`. All provide full66 sampling and normalized log probability. Spline uses eight alternating-mask RQ coupling layers, ten bins, two width128 hidden layers plus a rank32 output bottleneck, and a conditional affine final map. The triangular model uses the same principle for `q(a|y) q(b|a,y)`, default a16/b50, retaining stochastic full-rank innovation in every strong coordinate. It is a fixed coordinate split; a geometry meaning requires a supplied legitimate frozen frame. The full covariance Gaussian uses a positive-diagonal dense Cholesky factor.

Actual default parameter counts:

| Kind | Parameters |
|---|---:|
| full covariance Gaussian | 375909 |
| diagonal Gaussian | 682943 |
| spline | 682860 |
| triangular spline | 549432 |
| CFM | 251842 |

The diagonal control retains the same width128 raw-y encoder as spline; its two width655 head layers match spline capacity to 83 parameters. Other model counts differ and must be reported. `--width`, `--diag-hidden`, `--layers`, `--output-rank` and a 2M default hard cap make capacity comparisons explicit; they do not perform an automatic capacity search.

## CFM limits

CFM uses full66 independent Gaussian noise given the same y and linear interpolation with no endpoint smoothing or projection. Four residual velocity blocks, Fourier time embedding and the raw-y encoder are implemented. A standard raw-coordinate Gaussian base is the default. `--base anchor` requires external y-only `base_mean[N,66]` in both NPZs, with optional positive-diagonal full-rank `base_chol[N,66,66]`. Without a supplied Cholesky, `--anchor-scale * I` is used. No anchor physics is launched here. Base provenance cannot be inferred from array values and remains the parent's responsibility.

The sampler uses RK4 with exactly 16/32/64 velocity evaluations when those NFE values are requested. CFM `log_prob` refuses by default. `allow_approximate=True, nfe=64` explicitly enables complete divergence, using all66 state derivatives and no stochastic trace. It is a numerical continuous-flow approximation. It is **not** the determinant of the finite RK4 sampling map, and must not be used as an exact-MH proposal density. The predictor accepts raw-coordinate `base_mean` and `base_chol` for anchor-base sampling/log density, transforms them consistently and requires them for an anchor-trained checkpoint.

CFM checkpoint selection uses heldout **raw full66 latent energy score**, four fixed independent draws, distinct-pair U-statistic and16 NFE by default. Validation FM MSE is a diagnostic only. This is not material energy score or a posterior qualification result. Frozen sampling noise, t and validation batching are reused across checkpoints.

## Training, LR and resource contract

Single fit:

```text
python -m src.p3r4.train_density --train TRAIN.npz --validation VAL.npz --out NEW_RUN --kind spline --seed 11 --updates 6000 --device mps --lr 0.0002 --max-wall 900
```

Direct script invocation is also supported. `--max-wall`, `--max-wall-seconds` and `--max-seconds` are aliases. A positive cap checks elapsed wall time between updates; zero delegates the hard wall budget to the parent runner. It installs no SIGALRM handler. SIGINT/SIGTERM request a stop and final validation/checkpoint. A currently running validation finishes before the stop takes effect. The external runner must enforce any truly hard wall limit. No internal campaign scheduler or implicit extra fit is installed.

Output: `RUN_CONTRACT.json`, `NORMALIZATION.npz`, `BEST.pt`, `LAST.pt`, `loss_curve.jsonl`, `loss_curve.npz`, `LOSS_CURVE.png`, `validation_metric.npz`, `RESULT.json`. `LAST.pt` includes optimizer/RNG state; automatic resume is not implemented. Existing BEST/RESULT destinations are refused. Initial validation is always saved and can remain the selected checkpoint. Nonfinite training/gradient/validation values fail rather than silently repairing the target or selecting omitted rows.

Data-scale checks record the TRAIN subset and fitted finite column scales (constant columns use scale1). Capacity checks enforce the declared parameter cap and report parameters per independent TRAIN parent. LR is explicit, positive finite, and recorded per curve point. Empirical LR/data/capacity comparisons have **not** been run by this worker. `RECOMMENDED_TRAINING_COMMANDS.json` contains 45 proposed argv lists: five kinds, three seeds (11,22,33), three existing TRAIN sizes (4096,16384,32768), 6000 updates each, fresh output paths, with900-second perfit ceiling. Parent must queue them within the shared budget; the manifest is not authorization to run all45. Separate LR controls at1e-4/5e-4 need separately budgeted fresh fits.

## Numerical evidence and observed repair

`NUMERICAL_CORRECTNESS.json` and executable `check_density.py` record inverse, normalization, analytic/autodiff Jacobian, full covariance reference, strong conditional dependence, raw coordinate Jacobian, Gaussian-base CFM density and RK4 convergence, CPU/MPS gradients and perturbed inverse checks. The random1D spline integrates to0.999998983 over[-10,10]. Full66 float64 inverse/logdet tests pass near machine precision; perturbed float32 CPU/MPS tests pass near1e-6. The analytic linear CFM audit improves with16/32/64 NFE, and a full-trace MPS check is finite.

The first audit found zero autograd slope at exact spline knots/tail seams caused by PyTorch2.14 clamp's endpoint derivative convention. Removing unnecessary theta clamps restored the spline's analytic derivative; the explicit seam tests now pass. This was an implementation repair, with no change to the statistical target.

`CHECKPOINT_IO.json` records CPU Gaussian and MPS CFM update0-only CLI runs, saved-model reload, sampling, finite raw density and selection of initialization. Their optimizer update count is0. The loss plot was visually inspected. Resource receipts are per process, using wall/process CPU and macOS peak RSS; the parent budget runner owns end-to-end accounting. No SHA256, provider calls, prior/history edits, or physical model execution occurred.
