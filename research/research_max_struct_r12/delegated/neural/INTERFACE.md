# Frozen NN interface study

`run_nn_interface.py prepare` verifies the four exact IDs from the actual R8
primary public manifest, freezes the comparison, loads the existing physics
weights, computes original ten-channel public features, and saves four NN
initializations. The supplied checkpoints and historical datasets are read-only.
No training is performed. All four new initializations reproduced their stored
historical initializations exactly.

`run_nn_interface.py run` executes one local CPU/one-BLAS-thread queue:

- `uniform_direct`: uniform .1 to original physical L-BFGS-B.
- `nn_direct`: frozen NN material to that same physical solver.
- `nn_dr`: complete physical lift of NN material to DR (scale6, omega1,
  grouped ports, 3000 steps, audit50), then the same physical solver.
- `nn_cie_single`: original R8 DPCS/Fourier CIE interface and 3000 auxiliary
  steps, then that same physical solver. As in historical R8 this initializes
  only material; `initialize_physical_current=False` is preserved.

All physical solves retain box [.02,20], the original MaxwellModel backend,
the original data normalization and smooth/TV prior. The shared corrector uses
the original numerical settings, at most3000 iterations/120CPU seconds each.
The cap is **corrector-only**. Proposal, projection setup, frozen features,
network load/inference and serialization are explicitly charged separately.
This is a shared-corrector comparison, not an equal-total-cost claim. CIE source
normalizations inside its original auxiliary objective are retained; online
selection always uses original physical Phi including the prior.

The exact R8 IDs are 0001,0004,0022,0027, selected before the new outcomes. These
are different datasets from identically numbered R10/R11 cases. 0001 has six
illuminations; the others have twelve. Public per-case acquisition contracts
and source aperture are taken from the actual R8 tasks. The historical0027
counterexample is deliberately included, so this is exposed diagnostic work.

`ALL_PREDICTIONS_SEALED.json` is written only after all16 outputs are terminal
and present. Only then may `score_nn_interface.py` open the four original R8
sealed truths. Scoring reuses R7/R8's pure material metric definitions, performs
no physics calls and cannot change predictions. No hashes are used.

Outputs:

- `FROZEN_PROTOCOL.json`, `PREPROCESS_RECEIPT.json`, `public/`: declared inputs,
  actual feature/inference costs and frozen NN initializations.
- `runs/{case}__{arm}/`: original-Phi paid history, every paid material,
  initial/precorrector/final predictions, physical operation counts, corrector
  termination and proposal costs. DR also retains full trajectory audit states;
  CIE retains its auxiliary history/candidates/counters.
- `OFFLINE_SCORES.json/.csv`, `OFFLINE_CPU_PREFIXES.json`, `RESULTS_ZH.md`,
  `materials.png`, `objective_cpu.png`: post-seal metrics and figures.

Paid CPU prefixes choose the minimum original-Phi eligible candidate, never the
best material error. A method that terminates early remains available at later
budgets without inventing extra expenditure. Nested operation counters and
phase timers overlap and must not be summed twice. Feature/network costs are
shared across arms, with per-object warm/cold attribution reported separately.
