# Results support handoff

Status: **CODE_READY_NOT_EXECUTED**. No plotting command, tests, numerical experiment, physical assembly, or solve was run by this worker. Only `tools/g0_figures.py` and this delegated support directory were edited. The root reports, gates, and root `COST_SUMMARY.json` were not edited.

- `tools/g0_figures.py` reads the existing `results/g0_qualification/RESULTS.json` and registered `G0_CASES.json`. It reads the frozen independent-error limit from `GATES.json` (currently 0.1). Absent/empty results or individual missing cases remain `NOT_RUN`; an old NPZ alone cannot promote a case.
- Default output: `figures/g0_qualification/G0_FIXED_VS_FITTED_WHITE.{png,pdf,csv}`, `G0_DOFS_AND_PORTS.{png,pdf,csv}`, and three `G0_RESPONSE_<CASE>.{png,pdf,csv}` groups. Every response CSV preserves original real/imaginary components plus magnitude and principal phase. The default cases are `weak_center`, `strong_center`, and `low_loss_high_contrast`; no selection by measured discrepancy occurs.
- The response plots read only `fixed_data` and `reference_data` from the existing NPZs. They use stored TX/RX indices and principal phase in radians, without interpolation. The existing public mesh PNG/PDF and its NPZ are referenced for reuse.
- Runtime summary artifacts go here: `FIGURE_SUMMARY.json`, `COST_SNAPSHOT.json`, `WRAPPER_COSTS.csv`, and `PHYSICS_RESERVATIONS.csv`. The root cost summary is read only to preserve historical-bank and accounting-gap context.
- Accounting sums FE/RHS from each unique reservation and CPU/wall from each unique `job_completion`. Kernel timers, observed RHS fields, and duplicate receipt copies are not added. Conflicting repeated completions fail visibly. An in-child snapshot lists unfinished wrappers and therefore excludes the current plot wrapper until its final receipt exists.
- Scientific wording is bounded to the geometry prerequisite. Weak contrast and low loss/high contrast do not establish resonance. Port counts do not establish a material dimension or speed saving. The script does not make a whole-G0 decision.

Current source evidence, inspected without execution:

- `RESULTS.json` contains all 12 registered geometry cases. Stored whitened fixed-versus-fitted discrepancy ranges from **0.006924982802734923** (`weak_center`) to **0.020328110689468968** (`low_loss_high_contrast`). All 12 producer flags `numerical_qualification` are true. The interpretation remains limited to the recorded discrete comparison.
- Fixed mesh `ndofs=28903`; interface-fitted references span `50257–50815`; public target ports `2257`; touched ports span `43–124`.
- `PORT_PLAN.json`: setup `NOT_RUN`; whole G0 remains `PARTIAL: finite perturbation/Jacobian/12-class stress pending`.
- `results/jobs/implementation-tests-32/RECEIPT.json` reports exit code 0; its log reports `Ran 32 tests` and `OK`. These are existing implementation/algebra/accounting tests; the new plotting script has not been executed or syntax-tested.
- `results/jobs/g0-qualification-001/RECEIPT.json`: exit code 0, wrapper process CPU `29.532977 s`, wall `31.699708624975756 s`, peak child RSS `1183580160 bytes`. Kernel times must not be added to this wrapper total.

Parent action: run `tools/g0_figures.py` once through the existing serial wrapper (a `PREFLIGHT` plotting job is sufficient because the script imports no physical backend), inspect the emitted PNGs, and use that wrapper's final receipt when finalizing campaign costs. Do not use `report_preflight.py` after physical records exist: its guard intentionally refuses that overwrite. Final report/publication wording and stage interpretation remain with Codex parent.

Verification by this worker: static source review and `git diff --check` only; no script execution. Runtime/render inspection is **NOT_RUN**.
