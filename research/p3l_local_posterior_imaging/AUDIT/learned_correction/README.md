# Learned proposal Full-MH comparison

The runner is authored only. No paid comparison, training, real checkpoint
inference, GT read, or offline risk evaluation was executed while authoring it.

`tools/run_p3l_learned_correction.py` freezes one Book-D job with exposed
DEV_000/028/030/026, all five declared neural families, checkpoint seed
2026102270, and 256 independent-MH transitions per pair. It never chooses a
model or parent by GT. Missing and incompatible checkpoints receive explicit
NOT_RUN statuses; the script never trains a replacement. CPU is explicit.

The saved A anchor supplies the exact scalar and image feature definitions
from the context runner: whitened y, whitened anchor residual, native anchor
q, log1p squared field-J singular values; image channels are real/imaginary
anchor material, A@z, and the norm of the field derivative. Training used
cheap 16-state contexts. The change to A-anchor contexts is explicitly saved
as a context-distribution shift. Existing train-only affine scaling and the
normalized neural wrapper are restored from the checkpoint without refitting.

The learned weak8/probit2 density is completed with the normalized correlated
strong56 Gaussian reference. A 10% standard Gaussian66 defensive mixture
gives full support. Every proposed complete x66 gets a paid unchanged Full
target evaluation. MH includes logQ(current)-logQ(proposed) using both
normalized mixture densities. Every rejection repeats q, x, Phi, prediction,
full density, and the physical endpoint. Saved finite A probit modes are
preferred deterministically; absent modes use an explicitly inward anchor
start. All initial Full states are freshly paid. Common per-parent acceptance
uniforms are fixed across models, with separate NumPy proposal and CPU Torch
streams. No model adapts during the comparison.

Each proposal intent is written before its Full call. Full outputs, endpoint
path, ledger sequence, and costs are saved before the MH decision. Every
transition then checkpoints raw proposals, the repeated accepted-state chain,
current cache, and RNG states. Unweighted model-proposal material moments
are separately labeled as proposal diagnostics, never calibrated posterior
moments. Failed jobs retain paid endpoints and the last accepted checkpoint.
The script refuses to overwrite an existing job protocol or comparison.

If all pairs run, the frozen bound is 5,140 fresh Full states and at most
four fresh receiver-J actions: 41,248 RHS, excluding no hidden gradient or
optimizer actions. Cache equation checks reuse paid LU arrays arithmetically.
The shared global Book lock enforces one paid worker. The wall cap is bounded
by 7,200 seconds; unused budgets are never transferred from another stage.

The synthetic control checks normalized independence-MH detailed balance and
shows that dropping the proposal ratio fails the same check. It independently
checks the scalar/image feature API, finite saved-mode preference, and explicit
inward fallback. Passing these controls does not establish mixing, scientific
benefit, posterior accuracy, or a qualified reference.
