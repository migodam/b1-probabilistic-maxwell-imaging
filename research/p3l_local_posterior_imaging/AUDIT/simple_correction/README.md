# Fitted simple-density Full-MH comparison

The prospective runner is code only. No paid context or Full-MH comparison,
real fitted-model inference, training, GT read, or risk scoring was executed
while authoring it.

`tools/run_p3l_simple_correction.py` uses the actual completed selection in
`D/simple_density/RESULTS.json`: the selected PCA conditional Gaussian and
`unconditional_selected`. Both are genuine fitted train512/validation128
normalized controls. Selection is reused without DEV/GT selection or refitting.
Each is loaded with `SimpleConditionalGaussian.load` and combined with the
normalized correlated strong56 reference plus a 10% Gaussian66 defensive
mixture through `LearnedFullProposal`.

The Book-D job is `simple_correction_exposed_DEV`, with four fixed exposed
parents and 256 transitions per control. Each parent gets a fresh independent
matched coarse_TRF16/no-polish context from the original prior start. Its
physics, full Jacobians, context arrays, best cache, and costs are paid D once
and shared by both controls. The initial x is the same for both controls,
preferring the completed B joint probit mode; previous B/A setup cost is
disclosed, and every initial Full state is freshly paid D.

The runner reuses the exact `run_pair` Full-MH implementation. Both normalized
mixture proposal directions enter the MH ratio. Every rejected transition
repeats x/q, prediction, Phi, density, and physical endpoint. Proposal intents,
raw Full outputs, endpoint paths, costs, transitions, accepted-state chains,
and RNG states are persisted immediately. The NumPy density model RNG is
deterministic for each parent/family (`[2026102270,parent_index,2,family_index]`);
the common MH uniforms remain shared (`[2026102270,parent_index,1]`).

The conservative bound is 2,560 context RHS plus 16,448 Full-MH RHS, totaling
**19,008 RHS**, under the global one-paid-worker Book lock and bounded wall
reservation. Existing experiment files and schedules are not overwritten.
Unweighted model-proposal material moments remain distinct from corrected
chains. This is a strong fitted simple-density baseline, not a posterior
calibration, mixing, or scientific-efficiency certificate.

Synthetic controls check the raw-coordinate normalized affine Gaussian,
deterministic model sampling, defensive full66 MH detailed balance, exact
reuse of the shared transition implementation, matched context options, and
the separate valid-but-untrained neural update0 status.
