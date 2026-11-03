# R6 public single-task runner

The parent freezes the public manifest and scientific config. The runner executes
one task, on one CPU thread, and never opens GT or computes truth metrics. A finite
parent-owned process queue can invoke it independently for each task.

```sh
/Volumes/migodam\'s-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python \
  research/research_max_struct_r6/experiments/run_case.py \
  --manifest research/research_max_struct_r6/configs/primary_final_manifest.json \
  --task-id obs_b2449bfb61f845bf__cie_single \
  --output-dir /absolute/path/to/new_reproduction
```

`--public-list` is an alias input mode for the same JSON task-list contract.
`--config` overrides the manifest's config reference. `--output-dir` overrides the
task output directory. Relative paths resolve against the worktree repository.
The output directory must not exist; historical results are never overwritten.

Manifest shape:

```json
{
  "config": "research/research_max_struct_r6/configs/primary.json",
  "tasks": [
    {
      "task_id": "obs_b2449bfb61f845bf__cie_single",
      "case_id": "obs_b2449bfb61f845bf",
      "cohort": "DEV_ALL_DISCLOSED",
      "arm": "cie_single",
      "public_files": [
        "research/research_fbg_r2_adaptive_feedback/data_contract/stage4_public/operators/uniform12.npz",
        "research/research_fbg_r2_adaptive_feedback/data_contract/stage4_public/obs_b2449bfb61f845bf.npz"
      ],
      "output_dir": "research/research_max_struct_r6/RESULTS/primary/obs_b2449bfb61f845bf__cie_single"
    }
  ]
}
```

`tasks` may alternatively be a dictionary keyed by task ID. A top-level task list
is supported when `--config` is supplied. A task may set `config` and `output_dir`.
All public NPZ files are loaded with pickle disabled. Private payload keys are
rejected before material arrays are used. Repeated physics/geometry keys must agree. Public observation material_bounds overrides legacy operator metadata; numerical bounds still come from the frozen common config.
The fixed geometry is 24×24 real material cells, 12 sources, 48 receivers.

Shared config shape for the parent-authorized expanded main comparison:

```json
{
  "run_id": "PARENT_FROZEN_RUN_ID",
  "common": {
    "host": 0.1,
    "material_box": [0.02, 20.0],
    "dpcs_rtol": 0.001
  },
  "continuation_stages": [{"mf": 12, "maxiter": 1500}],
  "corrector_maxiter": 3000,
  "corrector_accepted_checkpoints": [100, 500, 1000, 2000, 3000],
  "corrector_evaluation_checkpoints": [100, 500, 1000, 2000, 3000]
}
```

The common settings may explicitly include coefficient_tol=0, miniter=50,
material_inner_maxiter=60, material_inner_ftol=1e-12 and
material_inner_gtol=1e-8; mismatches are rejected. Use either dpcs_rank or
dpcs_rtol. No wall cap is added by the runner. Any additional explicit budget
setting remains the parent's responsibility.

Four fixed arms:

| Arm | Numerical implementation | PR | beta_current | beta_material |
| --- | --- | --- | --- | --- |
| cie_single | FBG-R3 mature baseline, observed copy | sourcewise | 1,.2,.2,.2 | 1,.2,.2,.2 |
| cie_dual | FBG-R3 mature baseline, observed copy | sourcewise | 1,1,1,1 | 1,.2,.2,.2 |
| global_dual | Original R2 cie_global.run_cie | global | 1,1,1,1 | 1,.2,.2,.2 |
| a_material | Original R5 factorial material-only | global | 1,1,1,1 | 1,.2,.2,.2 |

The first three stages are fixed mf3/mf6/mf10, each 500 updates. Explicit mf12
continuation is appended inside that same function invocation: no restart of
material or J, and the original Fourier coefficients are transported. The source
data/state normalizations are identical across arms in the original public
source coordinates. Initial material is the common host .1; initial joint J is
the DPCS pseudoinverse current with zero APCS coefficients. The physical host
prediction is a separate original-Maxwell validation.

The physical corrector is an observation-only copy of R4 refine, retaining
scale=1e6, maxfun=4*maxiter, ftol=1e-12, gtol=.01, maxls=30 and the original
real material box. Its limit must be an explicit positive integer. An optimizer
may stop earlier under its original stopping rules; the receipt records why.
The corrector starts from the minimum original-Phi candidate among the host and
all true auxiliary stage endpoints. Selection includes its original trajectory
incumbent, including evaluated trial points, as in R4.

Outputs:

- `FROZEN_CONFIG.json`, `PUBLIC_TASK.json`, `NUMERICAL_SCHEDULE.json` and
  `SOURCE_METRICS.npz` identify the actual configuration and fixed source metrics.
- `HISTORY.json` retains every auxiliary stage/PR scalar record, including
  per-source PR/alpha arrays for sourcewise CIE and scalar values for global PR.
  Per-iteration observations include actual CPU/wall and physical model counts;
  unavailable internal operation counts are null rather than estimated.
- `initial.npz` and `stage_NN_physical.npz` contain material and physical
  j/Q/predicted data already computed by original validation.
- `stage_NN_joint.npz` contains joint material/J/Q/A at each true endpoint,
  with basis metadata and observed cost in its sidecar JSON. These retain the
  original native Fourier coefficient coordinates.
- `AUXILIARY_CANDIDATES.npz`, `AUXILIARY_COST.json`, `CORRECTOR.json` and
  `RECEIPT.json` retain scores, stopping information, complete model operations,
  LU/RHS counts, phase/process CPU/wall and physics timing. Nested costs overlap.
- Configured corrector evaluation and accepted-iteration snapshots copy already
  paid material, gradient, current and data. No extra forward/adjoint is issued.
  Scalar evaluation/accepted histories cover every optimizer point. Evaluation
  history CPU/wall starts at observer creation (thus includes auxiliary time);
  accepted history CPU/wall starts at corrector entry.
- `prediction.npz` is the retained original-Phi incumbent, exported from its
  already paid state. Finalization adds no physical forward call.

Provenance and implementation check: `audit/mature_cie_source.py` and
`audit/mature_cie_observation.diff` preserve the source and localized observation
changes. `audit/r4_refine_source.py` and `audit/r4_refine_observation.diff` do the
same for the corrector. `audit/OBSERVATION_EQUIVALENCE.json` records a bounded
implementation check using an arbitrary public observation and no material GT:
original and observed mature candidates were bitwise identical; original and
observed corrector materials/value histories were bitwise identical; all
physical operation counts matched. This is an instrumentation check, not
scientific qualification or a material-recovery result.
