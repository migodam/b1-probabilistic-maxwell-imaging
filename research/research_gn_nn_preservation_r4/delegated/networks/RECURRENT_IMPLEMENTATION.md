# Recurrent current/material/physics feedback — prepared, not trained

Only `src/recurrent_feedback.py`, `experiments/train_recurrent.py` and this document are added. Existing six-arm files, checkpoints, budgets and active processes are unchanged. Inclusion in a finite VAL ablation and the fixed rho remain parent decisions.

## Primitive and recurrence

`feedback_step(operators, current, material, public_y, anchor_current, rho, ...)` implements the supplied quadratic feedback using only fixed G/S products and their Hermitian adjoints. Material and rho are fixed across its inner steps. Public Y and original j0 squared norms remain fixed denominators. For exactly zero gradient or curvature, the per-observation step is zero; division uses a safe denominator before masking, so no numerical 0/0 is evaluated. Zero public normalization norms use the same 1e-20 floor as the existing network algebra and are exposed in diagnostics.

The returned field and receiver residual are updated by linearity and reused by a second feedback step. No LU, inverse, optimizer, truth/support input, or finite-difference physics occurs inside this layer. Autograd differentiates the entire unrolled computation, including the upstream rowLS material; “fixed material” means held constant during the inner feedback updates, not detached from outer training.

`RecurrentFeedback(mode, rho)` contains one `ResidualUNet('j_full')`, shared across exactly two passes. Each pass predicts an additive cumulative-current increment, computes the existing shared-real boxed rowLS material, performs zero/one/two fixed-material feedback steps, then refreshes rowLS after feedback. Material and current both propagate to the next pass.

All four modes have identical 121,528 network parameters and two network passes:

- `no_feedback`: original fixed-anchor Q channels; current/rowLS cascade only.
- `refreshed_q`: Q refreshed from the current but no physics feedback.
- `feedback1`: refreshed Q plus one feedback step per pass.
- `feedback2`: refreshed Q plus two feedback steps per pass, holding the first rowLS material fixed between these steps.

Q refresh uses the original TRAIN means/stds for channels 2–25. The material channel follows the same fixed .4 scale. gb/noise remain anchored. An available public anchor Q0 supplies the first refreshed input without another G multiplication. No new statistics are fit per layer.

## API, costs and artifacts

The forward call is `model(operators, features, m0, j0, public_y, stats, anchor_q=None)`. Results include existing candidate array names `m_proposal`, `j_aux`, `q_aux`, `delta_current`, plus `raw_delta`, `stage_materials`, `stage_currents`, `stage_fields`, per-step differentiable diagnostics and `preparation_counts`.

Costs expose separate batched G/G-adjoint/S/S-adjoint call counts, network passes, rowLS calls, feedback steps, complex matrix columns and forward FLOPs. They explicitly report zero LU, zero solve-RHS and zero dense-physics inverses. Matrix-multiplication column counts are separate from solve-RHS accounting. Setup, scalar work, normalization, BN/ReLU and backward arithmetic are excluded from the forward FLOP estimate; training records measured epoch wall time and does not invent backward counts. All arrays remain complex64 or MPS float32 real/imag pairs.

`train_recurrent(..., rho=..., authorized=False)` is a separate opt-in driver. It preserves AdamW .001, batch64, max100 epochs, patience15, normalized material L1+.1edge and the existing optional const/internal terms. Modes start from scratch with the same seed and shuffled batches. Full VAL material loss selects checkpoints; no TEST truth is opened. Existing R4 normalization may be supplied through `--stats-checkpoint`; only its statistics are used; its weights are never assigned to a recurrent model. Recurrent checkpoints use a separate format so existing inference interfaces are unchanged.

No training, prediction cohort, model allocation or extra experimental budget is scheduled by this implementation. The parent owns physical qualification, mathematical scope and scientific interpretation.

## Necessary wiring check

PASS on local MPS: one deterministic synthetic two-pass/two-feedback-step autograd connection check and a zero-gradient branch check. First update: head gradient 0.172445, stem gradient 0 (zero-initialized head). Second update: head gradient 0.144451, stem gradient 0.019760. Stage-material shape was (2,2,576); the zero-gradient sample returned tau=0 and finite current. No cohort or truth was read, and diagnostic weights were discarded. This is not an accuracy, convergence or random error-floor test.

For the checked batch of two observations, feedback2 used 6 forward G products, 4 G-adjoint products, 6 forward S products and 4 S-adjoint products; 2 shared network passes; 4 rowLS calls; 0 LU and 0 solve-RHS. Forward fixed-operator estimate: 690,094,080 real FLOPs; network convolution estimate: 130,572,288 real FLOPs. These are forward arithmetic counts, not measured speed or training-cost claims.


## Public recurrent inference supplement — prepared, not run

`src/predict_recurrent.py` provides `predict_recurrent_anchors(checkpoint_path, explicit_public_anchor_paths, operators_path, output_dir, device="mps", operators=shared, stats=shared)`. It opens the recurrent checkpoint, only the six public anchor inputs and public cost receipts; it never opens a manifest, evaluator, class labels or truth. The shared immutable G/E/S upload and TRAIN statistics are reused. Candidate files contain the existing evaluation arrays, both stage arrays, optional feedback diagnostics, per-case inference timing, `preparation_cost_json` and fixed-operator cost fields. Global setup and base-anchor solving remain separate, following the existing preparation-cost convention. Missing feature receipts remain explicitly incomplete.

`train_recurrent --predict-val` calls this module after each new best model is complete, writing `predictions/recurrent_{mode}_seed7_cascade2/<obs_id>.npz`. It uses exactly the full156 VAL records established by the existing central launch-time cohort check. Default new training scope is now only `refreshed_q` and `feedback1`, rho=1, seed7. Existing Jfull baseline is reused; no no-feedback or feedback2 model is scheduled. The original six-arm interfaces and active processes are unchanged. Parent authorization and the existing6h total budget still control whether these models actually run. No additional wiring check or inference/training execution was performed for this supplement.
