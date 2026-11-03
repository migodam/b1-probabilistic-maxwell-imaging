# Confirmation preparation and remote source allowlist

Local preparation artifacts (not uploaded as solver code):

- `configs/confirmation_draft_frozen.json`
- `configs/confirmation_manifest_draft.json` (18 anonymous slots, no task IDs)
- `experiments/generate_confirmation.py`
- `experiments/package_confirmation.py`
- `CONFIRMATION_PROTOCOL.md`
- `CONFIRMATION_FILELIST.md`

Remote code copied explicitly from the current R4 worktree:

- `research/research_max_struct_r4/src/confirmation_solver.py`
- `research/research_max_struct_r4/src/cie_block_absorption.py`
- `research/research_max_struct_r4/src/cie_onehop.py`
- `research/research_max_struct_r4/src/cooperative_geometry.py`
- `research/research_max_struct_r4/experiments/cooperative_paths.py`
- `research/research_max_struct_r2/src/cie_global.py`
- `research/research_max_struct_r2/src/profile_current.py`
- `research/research_max_struct_r2/experiments/run_blind.py`
- `research/research_fbg_r3_experiments/src/physics.py`
- `research/research_fbg_r3_experiments/src/cie_tsom_baseline.py`
- `research/research_fbg_r3_experiments/src/full_maxwell_corrector.py`
- `research/research_z_alg_r4_diversity/src/physics.py`
- `research/research_z_alg_r4_basin/src/rational_rays.py`

The existing `cooperative_paths.refine` import pulls its module's public-only
R2 helper/current/corrector and R4 geometry/ray dependencies above. They are
included unchanged; no generator or scoring module is among these imports.

Copied unchanged from publication repository's
`research/research_z_alg_r4_diversity/`:

- `remote_runner.py`
- `vendor/threadpoolctl.py`
- `vendor/threadpoolctl-3.6.0.dist-info/` regular metadata/license files only

Generated later, after parent freeze and write-once generation:

- `data_public/confirmation/operators.npz`
- `data_public/confirmation/manifest.json`
- `data_public/confirmation/obs_<opaque16hex>.npz` ×6
- `configs/confirmation_frozen_PUBLIC.json` (private family-list key removed)
- `configs/confirmation_manifest.json` (18 actual opaque tasks, execution false)
- `PACKAGE_FILELIST.json` (concrete copy paths/sizes/source provenance)

Runtime outputs are new `results/confirmation/` receipts and
`predictions/<opaqueID>/<method>.npz` files under the isolated remote root. No
runtime output exists yet. Exclude `sealed_truth`, private mapping, seeds,
scoring/generation code, old results, old observations, credentials, Git metadata,
`__pycache__`, and any dependency installation or environment directory.
