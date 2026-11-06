# Six-case confirmation data contract — preparation only

Status: **PREPARED, NOT FROZEN, NOT GENERATED.** No materials, observations, features, NN outputs or inversions were generated during code preparation. The parent will freeze the generator and method contracts after current DEV evidence. `configs/confirmation_generator.json` contains the offline family/amplitude/noise schedule and root seed2026110609. It is not an inverse input.

The six fixed cases are ellipse(2.5,.001,full), disks(12,.01,half), annulus(18,.03,full), comb(12,.03,half), star(18,.01,full), perforated(12,.03,half). Each uses an independent child PRNG stream from the fixed seed and exactly one realization; no reroll or success filtering. This is a small diverse confirmation of the frozen noise-rule/generalization question, not a population-level statistical claim.

All cases reuse the exact copied R8 public operators: 24×24 generation and inversion, 12 exterior sources, 48 receivers, inherited scalar-TM disk-pulse VIE, exp(−iωt), radius.45h. Half/full denotes the inherited transmit aperture. There is no new operator build/backend, independent discretization, continuum validation or extra observation. Material geometry/antialiasing and background.1 are the existing R8 `parameters/material` functions; forward uses the existing `MaxwellModel.state`. The full/half public operator copies are in `data_public/operators/`; their preparation receipt records provenance without SHA.

Portable paths below are relative to the repository root, never the local worktree mount. `configs/confirmation_public_manifest.json` is a **planned observation list**, not an executable inverse queue or a claim that data exist. Actual generation will write the same record structure to `data_public/confirmation/manifest.json`. Parent-frozen inverse task manifests must list only each record's `public_files` (operator plus observation) and their separately frozen numerical solver config.

Every public observation NPZ has exactly:

- `Y`: 48×12 noisy complex data.
- `sigma`: one raw circular-complex noise standard deviation, with E|n|²=sigma² and sigma=relativeRMS·||Y_clean||F/√576.
- `noise_radius`: dimensionless √(chi2(1152,.99)/2), matching the R8 convention; multiply by sigma for an absolute full-noise radius. No noise draw is repeated to force it inside the ball.
- `grid_shape`: [24,24].
- `material_bounds`: [.02,20].

There are no public family, amplitude, parameter, seed, support, clean-data or feature fields. Opaque IDs are `r9_confirmation_0000` through `0005`. Operator files contain only inherited public physics/acquisition fields; they do not contain shape truth.

Generator-only outputs are physically separated:

- `sealed_truth/confirmation/<case>.npz`: m, offline evaluation mask and Y_clean.
- `sealed_truth/confirmation/score_metadata.json`: case_id and amplitude denominator only.
- `sealed_truth/confirmation/generator_metadata.json`: family/parameters/PRNG provenance for offline use only.
- `sealed_truth/confirmation/FROZEN_GENERATOR_CONFIG.json`: exact generator contract used, for offline reproducibility.

An inverse runner must never receive these three paths or `confirmation_generator.json`. GT scoring occurs only after parent-frozen predictions are sealed. No NN/features are generated or supplied. The generator imports the existing R8 material helpers without invoking its main/feature-generation path.

Future execution, after the parent freezes `status` to `FROZEN_FOR_GENERATION`, is `python research/research_max_struct_r9/experiments/generate_confirmation.py --execute`. The script refuses non-frozen contracts and existing public/truth/generation-receipt outputs. Generation failures preserve partial outputs/cost; there is no automatic retry. `audit/generation_confirmation.json` will retain actual per-object and shared setup costs, LU, complex RHS, CPU/wall and failed work. Expected successful work is one576-dimensional forward per case (6LU/72source RHS), but the receipt uses actual counters rather than assuming those counts; original copied-operator creation is reused external cost.
