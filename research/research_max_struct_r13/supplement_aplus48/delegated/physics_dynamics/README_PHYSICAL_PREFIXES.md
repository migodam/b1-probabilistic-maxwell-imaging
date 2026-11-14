# Source-only physical-prefix diagnostic

`SUP/src/analyze_physical_prefixes.py` is implemented and statically parsed. It has **not been executed**. No new physical solve, inversion, material candidate, public-array/GT read or plot was performed during preparation. The parent must first complete the original/new seals, authorize comparison, and finish `compare48.py` so its completed receipt, metric tables and `RUN_INDEX.json` exist.

## Fixed scope and guard

The coordinated `Aplus48_comparison_run_index_v1` schema has240 records:48 existing objects ×CIE/Dual-CIE/G19/DR/A+. All source paths are absolute. The script consumes only raw/corrected/primary plus15/30/60/90/120s candidate IDs: **1920 retained references**. Exact material duplicates are evaluated once per case/operator/prior/box, including across methods/views/prefixes. Missing or failing evidence remains a row, never excluded by quality.

The explicit parent `authorized_comparison_scoring=true`, scope `FOLLOWUP_EXISTING48`, root/supplement references and all three authorized terminal seals are checked through the existing stdlib-only `compare48.authorize`: core48, mechanism16 and extension32. The original every-third-ID16/remainder32 split, observation references and shared numerical protocol are retained. Index case/method sets, source folders, candidate files and public input paths must match the original/extension manifests. The completed comparison receipt must identify those same authorization/seal sources. These guards precede NumPy, model imports and every numerical array read. Original immutable view IDs and paid-min-Phi prefix rules are checked again against sealed receipts/public HISTORY; this is source-integrity checking, not a new scientific selection.

Additional RUN_INDEX `common_actual_CPU`/`common_physical_LU` selections are intentionally ignored. They do not enlarge the physical-evaluation scope without a separate parent request.

After authorization:

```text
EXISTING_GAUSSIAN_PYTHON SUP/src/analyze_physical_prefixes.py \
  --root ROOT --supplement-root SUP --legacy-repo EXISTING_R13_WORKTREE \
  --authorization PARENT_COMPARISON_AUTH.json \
  --run-index COMPLETED_COMPARISON_OUTPUT/RUN_INDEX.json \
  --out SUP/delegated/physics_dynamics/offline_prefixes
```

Use the existing Gaussian/.venv_nn interpreter with BLAS1; the source also fixes the BLAS thread environment. No installation. `--no-figures` is optional; plots are enabled by default. Output must be a new directory outside frozen data, source/config, prediction/queue and truth directories. No numerical request is dispatched until the guard passes.

## Exact physical quantities

The implementation follows `ROOT/src/analyze_stationarity.py::evaluate_material` without calling its fixed400-request driver. It uses the original FBG `MaxwellModel.data_value_gradient(m,Y)` with beta omitted, the exact original public G/E0/S/Y and frozen scalar prior/box. Execution input arrays must exactly equal the canonical original observation arrays. The shared regularizer is reconstructed with the original smooth/tv/epsilon parameters.

For each unique frozen material it records recomputed Phi,D,R; unscaled physical/data/prior gradient L2/infinity norms; the data/prior gradient cosine; the exact unit-step mapping `m−clip(m−gPhi,lower,upper)` and its infinity/L2 norms; exact box activity/violations; original physical state residual; and numerical agreement with saved Phi/D/R. Zero-gradient cosines remain undefined/null. Saved-Phi differences are numbers, not a pass tolerance or success gate. No optimizer1e6 scale, Hessian, SVD, new correction or alternative objective is introduced.

One retained forward LU serves all source RHS and its original analytic source adjoint. Actual counter deltas distinguish forward and adjoint RHS and aggregate complex RHS; unique costs, setup, CPU/wall, serialization and plotting scopes are explicit. Duplicate references share one actual evaluation; per-reference CPU shares are accounting allocations, not additional solves. Exact gradients are retained in per-case compact `GRADIENTS/<case>.npz` banks keyed by evaluation ID, so auditing direction calculations does not require another forward.

## Changes are selected-incumbent differences, not native steps

For each of the four adjacent paid-prefix intervals and the raw-to-corrected pair, compare the actual frozen material difference with minus the physical gradient at the old estimate. Record its norm/cosine, actual deltaPhi/deltaD/deltaR, `gPhi_old dot delta_m`, data/prior first-order contributions, predicted endPhi and nonlinear remainder. There are1200 retained transition rows.

A paid incumbent can remain unchanged while the algorithm takes many native steps. Conversely two saved incumbents may be separated by many iterations or line-search trials. These records are therefore explicitly **NOT native material/optimizer steps**, accepted-iteration traces, Hessian estimates or a basin/local-minimum certificate. Native-state instrumentation remains a separate later task.

No GT NPZ is needed or opened. Error annotations are attached only from the already-authorized comparison `ALL48_SCORES` and `PAID_PREFIX_SCORES` tables after matching the same frozen candidate ID; they never choose candidates. Missing annotations remain missing. Same-discretization original-Phi evaluation is **not independent physical validation**.

## Outputs after authorized execution

- `PHYSICAL_REFERENCES.{json,csv}`:1920 expected rows, exact source indices and saved/recomputed metrics, sharing references and comparison annotations.
- `UNIQUE_PHYSICAL_EVALUATIONS.{json,csv}` and `MODEL_SETUP_COSTS`: measured integer LU/RHS, CPU/wall and all source references. Repeated views/prefixes must not be summed as extra evaluations.
- `INCUMBENT_TRANSITIONS.{json,csv}`:1200 selected-estimate differences and first-order measurements.
- `PREFIX_SUMMARY`, `VIEW_SUMMARY`, `TRANSITION_SUMMARY`: explicit48-object denominators/available counts and descriptive statistics; no new success threshold or inferential gate.
- Three compact PNG/PDF plots: error/Phi/physical-gradient/box-mapping versus paid CPU; raw-versus-corrected estimates; incumbent direction cosine/first-order prediction versus actual change. Complete backing tables retain failed/missing rows.
- `RECEIPT.json`: parent authorization/seals/index snapshot, total/reference/unique counts, actual forward/adjoint costs, CPU/wall, zero GT-array reads and zero new optimization/candidate/Hessian/SVD operations.

All results remain `OFFLINE_DIAGNOSTIC` and `FOLLOWUP_EXISTING48_NOT_FRESH_BLIND`. Parent owns scientific interpretation and any later native-state or longer-schedule work.
