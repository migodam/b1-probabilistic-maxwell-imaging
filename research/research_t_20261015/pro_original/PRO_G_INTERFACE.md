# Research T <-> Pro G contract, version 1.0

This contract describes candidate exchange, not a probability posterior. No probability weight is required or inferred from residuals.

## Bundle layout
- manifest.json: schema version, run/candidate identities, physical conventions, costs and validation status.
- material.npz: chi (shared across illuminations), grid geometry/mass, optional background and declared dispersion parameters.
- observation.npz: original y or immutable observation reference, W/whitening definition, noise ID and channel packing.
- prediction.npz: full predicted channels for exactly that material and nuisance history; per-block residuals.
- optional_state.npz: total fields/current arrays when exposed by the backend.
- structure.npz: absolute-material, air, protected-core and guard masks with their physical thresholds.
- lineage.json: parent candidate, proposal, corrector, accepted step history and verified path segments.
- cost.json: physical state attempts, factorizations, forward/tangent/adjoint RHS, failures, setup and wall time.

Never write invented file references or zero arrays for unavailable fields. Use null plus an explicit reason.

## Material and measurement semantics
State time-harmonic sign, susceptibility/permittivity convention, units, mesh, quadrature/self-term, dimensions and vector/scalar model. State real/imaginary packing and whether sigma is a complex RMS or a per-real-component standard deviation.

Each candidate has ONE material field for all source experiments and a declared frequency-dispersion law. Each world has ONE nuisance vector for the entire acquisition history, subject to a frozen legal set. Different worlds need not share the same fitted nuisance.

An abs(chi-background) threshold is an increment-support diagnostic, never an automatic air/cavity definition. In 3D distinguish connected components beta_0, tunnels beta_1, and enclosed cavities beta_2.

## Validation levels
- NOT_VALIDATED: no matching full endpoint.
- REJECTED: illegal material, nuisance, structure, or data mismatch.
- FULL_VALIDATED_NUMERICAL: matching full endpoint, declared tolerance and physics consistency checks passed; no rigorous output enclosure is claimed.
- FULL_VALIDATED_WITH_BOUND: additionally, a justified output error upper bound epsilon is available and residual_norm + epsilon <= tau.

A backward linear-system residual alone does not give a forward/output error bound without a stability/resolvent estimate. Analytic Mie can legitimately omit local constitutive-state arrays; record this explicitly, not as a zero defect.

## Validation operation
validate(candidate, frozen_observation, frozen_physics, budget) ->
  status, residual_norm, tau, error_upper_bound_or_null,
  checks, reasons, updated_cost_receipt

A validation request that needs new physics must receive an explicit budget. Reading a bundle must not launch a solver.

## Common corrector operation
correct(candidate, protected_geometry, legal_domain, observation_history,
        fixed_corrector_config, budget) ->
  attempted_endpoints, best_valid_endpoint_or_null, path_receipts,
  KKT_and_line_search_diagnostics, costs

Trust-region bounds and active constraints are enforced in the shared real material coordinates. Retain protected-core and guard geometry, so the corrector cannot succeed merely by filling the intended cavity or opening it to the exterior. If the guard makes the problem infeasible, report that local restriction rather than asserting global impossibility.

A fixed-coordinate chart failure can trigger only a preregistered equal-budget enrichment. No hidden full-space fallback for one arm alone.

## Distances
material_distance: mass-weighted norm or L1 with fixed units/normalization.
structure_distance: a declared combination of absolute-phase mask symmetric difference, localized cavity matching, and persistence/boundary descriptors. Report each component, not only one opaque score.
prediction_distance: whitened mean-response difference under explicitly fixed or jointly profiled nuisances; distinguish from each candidate's residual to y.

Closeness is not an equivalence relation. Clustering by a threshold is an algorithmic archive operation, not proof of a connected feasible component. A path is only recorded as feasible when every claimed segment has a justified tube bound or sufficient validated steps with its limitation stated.

## Status vocabulary
NO_QUALIFIED_BASE
NO_NEW_STRUCTURAL_CANDIDATE_FOUND
LOCAL_CHART_EXCLUDED_WITH_BOUND
SEARCH_BUDGET_EXHAUSTED
UNKNOWN_NUMERICAL_ERROR
MODEL_DISCREPANCY_SUSPECTED
VALID_NUMERICAL
VALID_WITH_BOUND

No finite search failure becomes GLOBAL_INFEASIBLE. No candidate count becomes posterior mass. No claim of complete coverage is allowed without an explicit enclosing set-inversion certificate.
