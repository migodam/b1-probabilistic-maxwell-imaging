# Learned proposal Full-MH comparison

The runner is authored only. No paid comparison, training, real checkpoint
inference, GT read, or offline risk evaluation was executed while authoring it.

`tools/run_p3l_learned_correction.py` freezes one Book-D job with exposed
DEV_000/028/030/026, all five declared neural families, checkpoint seed
2026102270, and 256 independent-MH transitions per pair. It never chooses a
model or parent by GT. Missing and incompatible checkpoints receive explicit
NOT_RUN statuses; the script never trains a replacement. CPU is explicit.

The actual completed training audit now shows all 15 family/seed runs
selected best_update=0 and restored base weights equal their exact seeded
initializations. This is valid checkpoint selection, not a missing or invalid
checkpoint. The prospective preflight records
**NOT_RUN_VALIDATION_SELECTED_INITIALIZATION_NO_TRAINED_GAIN** and avoids
paying Full-MH for these untrained neural weights. If every family receives
this status, no Target, context, or Full query is constructed. This efficiency
decision follows the actual completed validation results and weight audit;
it adds no arbitrary numerical threshold and does not show that neural
architectures cannot learn. `STAGE_D_DECISION.json` records its evidence.

The separate prospective `tools/run_p3l_simple_correction.py` instead compares
the genuinely fitted selected PCA Gaussian and unconditional Gaussian from
`D/simple_density/RESULTS.json`. Its Book-D job is
`simple_correction_exposed_DEV`. Fresh matched contexts are shared by these
two controls; its conservative total bound is **19,008 RHS** (2,560 context
plus 16,448 MH). It reuses the exact normalized Full-MH/persistence code and
shared acceptance uniforms, with deterministic per-parent/family NumPy model
RNGs. Both model proposal moments and corrected chains remain explicitly
distinct; this is a strong simple-density baseline, not posterior calibration.

For each parent, a fresh independent coarse_TRF reconstruction from the
original prior start (zeros64, .5, .5), with a 16-fresh-state cap and no polish,
supplies the primary context. This exactly matches the training-context law.
The bounded recent/best Point cache in the current reconstruction is used.
All context physics, full Jacobians, endpoint records, and final paid cache
are billed to the same Book-D job and persisted. All five families share the
same per-parent context. Features exactly follow the context runner:
whitened y, whitened context residual, native context q, log1p squared field-J
singular values; image channels are real/imaginary context material, A@z,
and the norm of the field derivative. Existing train-only affine scaling
and the normalized neural wrapper are restored without refitting.

The earlier authored A96-anchor-context design has status
**NOT_RUN_SUPERSEDED_BEFORE_PHYSICS**. It is an unrun context-distribution-shift
stress design, not the primary comparison and not a physical failure. No
protocol or resource ledger is reset or rewritten by this source revision.

The learned weak8/probit2 density is completed with the normalized correlated
strong56 Gaussian reference. A 10% standard Gaussian66 defensive mixture
gives full support. Every proposed complete x66 gets a paid unchanged Full
target evaluation. MH includes logQ(current)-logQ(proposed) using both
normalized mixture densities. Every rejection repeats q, x, Phi, prediction,
full density, and the physical endpoint. Initialization first prefers
`B/joint_gaussian/{parent}/JOINT_MODE_GEOMETRY.npz:x_mode` from the completed
joint-mode Book-B job, then old saved finite A probit modes, then an explicitly
inward A-anchor fallback. All initialization sources and coordinates are
frozen before new context queries. Previous B/A setup is disclosed as shared
initialization cost; the context TRF never starts from these modes. All five
families share the same initial material, and all initial Full states are
freshly paid D. Common per-parent acceptance
uniforms are fixed across models, with separate NumPy proposal and CPU Torch
streams. No model adapts during the comparison.

Each proposal intent is written before its Full call. Full outputs, endpoint
path, ledger sequence, and costs are saved before the MH decision. Every
transition then checkpoints raw proposals, the repeated accepted-state chain,
current cache, and RNG states. Unweighted model-proposal material moments
are separately labeled as proposal diagnostics, never calibrated posterior
moments. Failed jobs retain paid endpoints and the last accepted checkpoint.
The script refuses to overwrite an existing job protocol or comparison.

If all pairs run, the conservative frozen budget is 4×16×40 RHS for fresh
matched contexts plus 41,248 RHS for MH: **43,808 RHS**. There are at most
64 new context states and 5,140 new MH initialization/proposal states. Context
setup cost is recorded once per parent and shared across the five families;
it is not silently omitted or billed five times. All actual actions remain
in the same Book-D ledger, with no hidden gradient or optimizer actions.
The shared global Book lock enforces one paid worker. The wall cap is bounded
by 7,200 seconds; unused budgets are never transferred from another stage.

The synthetic control checks normalized independence-MH detailed balance and
shows that dropping the proposal ratio fails the same check. It independently
checks the scalar/image feature API, completed-B then saved-A mode preference,
explicit inward fallback, and the frozen matched context call options.
Passing these controls does not establish mixing, scientific
benefit, posterior accuracy, or a qualified reference.
