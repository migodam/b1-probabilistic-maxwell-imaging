# R9 inherited baseline runner

This is a small policy/path/receipt adapter over the actual R8 `run_case.execute` and its inherited benchmark solver. It adds no physical operator, optimizer, learned model, training, GT reader or scorer. Nothing executes on import. No solver has been launched by this delivery.

Files:

- `research/research_max_struct_r9/src/baseline_adapter.py`
- `research/research_max_struct_r9/experiments/run_baseline.py`

## Parent-frozen inputs

Use an R8-shaped manifest with a unique task ID, `case_id`, `arm`, `public_files`, solver `config` and fresh R9 `output_dir`. The parent must choose the public development cases and review the rule before running. The frozen full R8 data remain in sibling `r8_worktree`.

Path rules are explicit:

- CLI manifest/config/output paths are absolute or relative to the R9 repository root.
- `public_files` and `initialization_file` are absolute or relative to task/manifest `input_base`; its default is the actual sibling R8 repository.
- `--r8-repo` can name another actual copy of the inherited R8 repository. Dependencies must come from that copy; an already imported other research checkout is rejected.
- Output and its sibling `__adapter_inputs` directory must be inside `research/research_max_struct_r9`; neither may already exist. Historical R8 output cannot be overwritten.
- Set the config to an explicit R9 file or an absolute inherited R8 solver config. There is no silent config search across checkouts.

Example schema only; these development identifiers are placeholders, not selected cases:

```json
{
  "config": "/absolute/r9_worktree/research/research_max_struct_r9/configs/dev_solver.json",
  "input_base": "/absolute/r8_worktree",
  "tasks": [{
    "task_id": "parent_selected_dev__a_plus__uniform__noise_ball_99",
    "case_id": "parent_selected_public_dev",
    "arm": "a_plus",
    "initialization": "uniform",
    "dpcs_rule": "noise_ball_99",
    "public_files": ["research/research_max_struct_r8/data_public/operators/PARENT_OPERATOR.npz", "research/research_max_struct_r8/data_public/PARENT_OBSERVATION.npz"],
    "output_dir": "research/research_max_struct_r9/RESULTS/baseline_dev/PARENT_TASK"
  }]
}
```

The config retains the inherited solver structure (`run_id`, `common`, continuation stages, A+ options, corrector limits and physical-direct limit). Its input `common.dpcs_rtol` must remain `.001` with no prefilled rank. This adapter changes only that rule in a generated effective config when `noise_ball_99` is selected. It does not shrink the schedule or modify scientific budgets.

Canonical inherited arms:

|Arm|Meaning|
|---|---|
|`cie_single`|Single-beta sourcewise CIE|
|`cie_dual`|Dual-beta sourcewise CIE|
|`global_dual`|Global dual-beta current update|
|`a_material`|Inherited original A readout|
|`a_plus`|Inherited global-current A+|
|`a_plus_sourcewise`|Same A+ readout with sourcewise current update|
|`physical_lbfgs`|Original physical L-BFGS-B, no auxiliary current solve|

Aliases `single`, `dual`, `global`, `A`, `A+global`, `A+sourcewise`, `direct` are recorded and mapped to those canonical arms.

Initialization is `uniform` (default inherited host, optionally explicit scalar `initial_material`) or `nn`/`frozen_nn` with a supplied public `initialization_file` containing `m_init`. NN material must be finite, real, 576-dimensional and in the common box; it is never regenerated or rejected by initial Phi. `initialize_physical_current` remains an explicit inherited option: when requested, R8 performs its legal DPCS/Fourier affine projection, rather than treating physical current as automatically DPCS-feasible.

## DPCS controls

`fixed_001` is the unchanged separate historical control. The R8 solver computes its own receiver SVD and selects `singular_i > .001*singular_0`.

`noise_ball_99` requires a public, finite, real scalar `sigma` common to every raw complex observation, with **E|n|² = sigma²** and independent circular complex Gaussian noise. It uses unnormalized public S and Y. Source-normalized auxiliary objective weights remain the inherited solver's weights; they are not used to transform this noise model.

For Y of shape `(nr, ns)`, set

```
R² = sigma² * chi2.ppf(.99, 2*nr*ns) / 2.
```

For S=UΣV*, enumerate contiguous prefixes r=0,…,r_supported and select the smallest r satisfying

```
||(I-U_r U_r*) Y||F² <= R².
```

A singular mode is numerically supported only when `singular_i > eps(float64)*max(S.shape)*singular_0`, matching the existing R8 numerical floor. Orthogonal-coordinate tail energies are accumulated from nonnegative terms, including the component outside the thin U range when needed, rather than subtracting large nearly equal energies. When U spans the complete observation space, its outside component is algebraically zero. Rank zero is legal. A zero-noise or under-supported case that cannot meet the criterion fails explicitly; it is not rescued by an arbitrary rank.

The event `||N||F² <= R²` has probability .99 and implies `||P N||F² <= R²` simultaneously for every orthogonal projector P. Thus one full noise event supports the rank choice, including a data-dependent choice, without per-r significance tests or rank-dependent degrees of freedom. This does not guarantee recovery, identifiability or absence of model mismatch.

Non-scalar sigma, nonidentical noise, covariance/precision, supplied whitening and normalized-observation flags are rejected for this optional rule. Unsupported noise distributions are not silently converted. `noise_radius` is not substituted for the declared sigma convention. For a direct-physics arm, DPCS is explicitly not applicable and no rank-selection SVD is performed.

## Reuse and accounting

The adapter generates immutable absolute-path `MANIFEST.json`, `EFFECTIVE_CONFIG.json`, `DPCS_RULE.json` and `ADAPTER_PREPARATION.json` under `<output>__adapter_inputs`, then calls unchanged R8 `execute` in the same process with one BLAS thread. Provenance identifies the actual inherited R8 wrapper, benchmark solver and CIE source paths and their bytes/mtime, without SHA checks. The old hard-coded R8 git metadata is labeled as inherited metadata, not an R9 commit or new physical implementation.

The inherited outputs are retained: stage joint m/J/Q/Fourier coefficients; initial and stage physical states; original-Phi candidates; original paid incumbent budgets; complete stop/corrector/FFT/material/LU/RHS costs and final receipt. No finalization forward is introduced.

The optional rank rule performs **one additional public receiver SVD** before R8 execution. It is separately charged in `R9_DPCS_RULE.json` and the preparation CPU/wall. R8 still performs its own receiver SVD, which its own counters retain. These are two actual SVDs, not falsely reported as a shared cache. Selection performs zero physical LU/RHS.

`BUDGET_STATES.npz` is preserved unchanged. `R9_BUDGET_STATES.npz` contains the same paid original-Phi incumbents with the measured R9 preparation CPU prefix added to the CPU column; LU/RHS are unchanged. `RECEIPT.json` retains inherited nested costs and adds explicit R9 total/preparation fields. Do not add an inclusive total to its nested components. R9 finalization/serialization after the last paid budget point is a separate tail, not additional physics. A supplied frozen NN file is only consumed: its earlier feature/inference/training creation costs are not rerun or included in these solver totals, and must be reported separately for a full deployment comparison.

Preparation failures, including no supported discrepancy rank, receive a write-once preparation failure receipt and do not invoke the inverse solver. Inherited solver failures preserve R8's failed receipt and paid costs.

## Run after parent review

From the R9 repository root, using the existing project Python:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" \
  research/research_max_struct_r9/experiments/run_baseline.py \
  --manifest research/research_max_struct_r9/configs/PARENT_MANIFEST.json \
  --task-id PARENT_TASK --dpcs-rule noise_ball_99
```

Use `--dpcs-rule fixed_001` as the separately frozen control. Config/output overrides are `--config FILE --output-dir DIR`; actual inherited source override is `--r8-repo DIR`.

`--prepare-only` reads public inputs and may compute the public noise-rule SVD, but never calls the inverse solver. It reserves immutable adapter-input files; use a fresh output directory for a later actual solve rather than overwriting that preparation. No NN training, GT read, scoring, launch, publication or dependency installation was performed for this implementation delivery.
