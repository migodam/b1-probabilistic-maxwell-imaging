# Caller-array density trainer

Only `src/p3l/neural_training.py` and this audit directory were written. The trainer does not construct splits, generate observations, inspect campaign/held-out/truth paths, use physics, invoke an external service or select a device automatically. Parent supplies all arrays and owns independent-parent provenance. Actual campaign training runs=0. Validation below uses three optimizer updates on explicit synthetic arrays only.

## API

```python
from p3l.neural_training import train_density, load_density_checkpoint
result = train_density(
    'mlp_full', train_scalar, train_image, train_labels,
    validation_scalar, validation_image, validation_labels,
    output_dir=new_run_directory,
    label_coordinates='8 fixed Gaussian weak coordinates +2 probit properties',
    seed=2026102260, device='cpu',
)
restored = load_density_checkpoint(result['checkpoint_path'], device='cpu')
raw_logp = restored['model'].log_prob(raw_labels, raw_scalar, image)
raw_samples = restored['model'].sample(count, raw_scalar, image)
```

The wrapper uses all completed neural.py architectures without modifying them. Parameter matching, Gaussian/MDN normalization, coupling inverse/logJacobian and Cholesky floor remain their existing contracts. Shapes are inferred from explicit scalar/label arrays; images must be4×24×24. Default dimensions are inherited when caller arrays have1154scalar/10label dimensions.

Default single-run configuration: Adam1e-3, batch128 sampled with replacement by a local seeded NumPy generator, at most3000updates, validation every100updates, patience10 nonimproving complete evaluations, pre-step gradient norm clipping at5, cooperative wall cap1200seconds. Initial validation at update0 is eligible for best selection; final update is evaluated even if not a multiple of100. Smaller explicit caps/intervals support bounded controls. Requests above3000updates or1200seconds are rejected. No automatic hyperparameter search, three-seed expansion or training retry occurs.

Declared seeds:2026102260,2026102261,2026102262. Exactly one supplied seed runs per call. A different explicitly supplied seed is recorded as an override, rather than silently treated as one of the declared seeds.

## Scaling and density coordinates

Scalar and label mean/population SD are fitted on training arrays only. Exactly zero empirical SD uses scale1 and is recorded. Scales/means are stored as actual float32 values used by the model; nonpositive/nonfinite/unrepresentable scales fail. Validation arrays never fit scaling or enter optimizer updates. Images pass through unchanged, with no GT-dependent transform.

For raw supplied label y and standardized z=(y-mean)/scale:

```
log p_raw(y|context) = log p_standard(z|scaled_context) - sum(log(label_scale))
NLL_raw = NLL_standard + sum(log(label_scale)).
```

`ScaledConditionalDensity.log_prob` includes this exact affine change-of-variables formula, evaluated in the declared float32 runtime. Scalar conditioning normalization requires no label-density Jacobian. Sampling applies the inverse affine label map. No label clipping, probit conversion, bounded-property Jacobian, whitening transform or posterior likelihood factor is inferred. Parent must supply/record every non-affine label transform and add the corresponding Jacobians if changing reported coordinates.

`label_coordinates` is mandatory metadata. The reported raw-target proper NLL refers only to those exact caller-supplied label coordinates. It is not automatically an NLL in bounded properties, materials,66D full posterior or any unspecified target measure.

Default supervision is `JOINT_SIMULATION`. DISTILLATION is rejected unless `qualified_distillation=True` is explicitly supplied; even then checkpoint metadata says `CALLER_ASSERTED_NOT_VERIFIED`. The trainer cannot qualify a teacher or verify data independence. Both paths retain `NORMALIZED_ARCHITECTURE_NOT_POSTERIOR_CERTIFIED` status.

## Checkpoints, limits and billing

Lowest complete validation raw-target NLL selects `BEST_DENSITY.pt`, saved atomically with model weights, model type/options, all training scaling, seed/configuration, coordinate metadata, selected update and NLL. At completion the selected best state is restored, rather than returning the final optimizer iterate. `TRAINING_REPORT.json` records every completed or wall-interrupted validation, batch raw-target NLL, preclip gradient norm, elapsed wall, best selection and stop reason. Unhandled errors preserve `TRAINING_FAILURE.json` plus any prior best checkpoint.

Return value includes restored wrapper model, CPU best state dictionary, train scaling, checkpoint path, metadata and logs. Restore uses `weights_only=True` for this format and strict state loading, without data lookup or retraining. Existing completed run output is rejected; explicit checkpoint restore is separate from a new run. No implicit resume occurs.

Device must explicitly be cpu or mps. CPU runs set Torch to one thread within the call and restore the previous thread count afterward. MPS unavailability raises without CPU fallback. Data lives in explicit CPU arrays/tensors and transfers only individual minibatches to the requested local device. Initialization uses a preserved CPU RNG context; deterministic MPS behavior is not claimed.

On MPS, `gpu_wall_callback(seconds)` is called once in finalization, including exception paths. It receives synchronized allocated-job wall including host overhead, not isolated device-kernel occupancy; parent bills the callback value. CPU does not invoke a GPU billing callback. The callback should record billing without raising. No GPU budget is inferred from a single-run wall cap. Timing is cooperative between batches/validation batches: a running operation or checkpoint may overrun the cap. Partial validation is never accepted as a new best. Final checkpoint/report serialization adds overhead after the recorded pre-finalization wall value; the MPS callback includes finalization overhead up to billing.

## Dummy validation

`check_neural_training.py`, `results.json`, `check.stdout`, empty `check.stderr`, and the named dummy run directory preserve10 checks with16train/8validation synthetic observations, scalar dimension6, label dimension2, zero images, mlp_diag and3Adam updates on CPU.

- Training-only scaling remains unchanged despite deliberately shifted validation contexts.
- Four complete evaluations: update0 and each of three updates.
- Raw-coordinate log-density difference matches the label-scale Jacobian within `2.38e-7` in float32.
- Update2 has the best raw validation NLL `4.383848667144775`, and restored checkpoint densities are bitwise identical.
- CPU is not billed as GPU; coordinate/supervision metadata persists.
- Unqualified distillation and implicit device auto-selection reject; completed output cannot be overwritten.

No campaign arrays were used; no full training was performed. The small dummy fit validates loop/checkpoint/Jacobian plumbing, not useful inference, learned calibration, generalization, amortization or posterior correctness. MPS execution and billing remain NOT_RUN here.

Reproduce from worktree root (a fresh dummy directory is created each time):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ../../.venv_nn/bin/python \
  research/p3l_local_posterior_imaging/AUDIT/neural_training/check_neural_training.py
```
