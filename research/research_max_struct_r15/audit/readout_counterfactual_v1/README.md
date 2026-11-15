# Fixed-state material-readout counterfactual

{
  "status": "COMPLETE",
  "expected_cases": 4,
  "expected_readouts": 16,
  "completed_readouts": 16,
  "source_policy": "queue_interventions.source_checkpoint: requested512, latest complete numeric<=512 else earliest complete; metadata-only availability selection",
  "CPU_seconds": 0.8509829999999999,
  "wall_seconds": 1.1133513753302395,
  "whole_process_CPU_seconds_at_receipt": 0.879723,
  "truth_access": "NONE",
  "current_PR_updates": 0,
  "main_candidate_pool_writes": 0,
  "corrector_calls": 0,
  "online_selection": false,
  "physical_LU_total": 20,
  "actual_complex_RHS_total": 288,
  "cost_scope": "One case physical gradient (retained initial LU+adjoint); four new-material forwards+prior only (exact cache hits retained). Case/arm/internal timers nested, do not sum; failed work charged. Auxiliary setup does not perform a new SVD.",
  "projected_gradient_definitions": "descent residual=m-clip(m-g,lo,hi); KKT projected gradient removes blocked descent at box-active tolerance1e-8. These are fixed diagnostics, not changed optimizer stopping rules.",
  "issues": 0
}

Each arm starts from the same frozen full m,J,Q,beta,Os,Ls and original prior. Only one inherited matrix_material_readout is applied per arm, with original seed/bounds/LBFGS8/acceptance. No current/PR state is changed or continued. No new result is inserted into a main candidate pool. No GT is read.

Physical-gradient alignment and linear prediction use the single original full physical Phi gradient at the shared source. New-material Phi uses only a full physical forward and prior, no adjoint/corrector. Nonlinear finite changes are not guarantees or online selectors. Parent owns interpretation.

Checkpoint inherited auxiliary/projector counters are preserved as provenance; every counterfactual operator starts a fresh Counter. Initial physical-gradient costs are shared once per case, per-arm physics deltas and full case totals are distinct; do not add these nested views. Optional real material sector directions use only saved geometry, no SVD. Failed case prerequisites still produce four failed readout rows and receipts.
