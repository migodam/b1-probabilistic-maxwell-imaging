# Paid primary endpoint validation (no GT)

Fixed denominator: 4 cases ×11 main methods =44. Complete 44; missing 0; failed 0. Missing/failed endpoints remain in ENDPOINTS.csv and are not replaced.

Only the new parent-authorized main264 prediction seal is accepted. Each endpoint uses the existing primary.npz real m. No candidate, optimizer, finite step experiment, GT score, or selection feedback is produced.

Original MaxwellModel.state(m) supplies J/Q/Y_pred and a retained LU. One repeated identical-m state call checks object/LU identity and zero added physical counters; data_value_gradient then reuses this LU for one ns-RHS adjoint. All actual counters, including a failed operation or unexpected extra factorization, remain charged. No full576 material Jacobian is requested.

D=.5||Y_pred−Y||_F²/||Y||_F²; R is the inherited smooth/TV prior; Phi and total gradient come from original data_value_gradient. g_prior comes from the unchanged model.reg.evaluate; g_data=g_total−g_prior. Unit-step box gradient mapping is m−clip(m−g_total,lower,upper); it is a diagnostic with declared scale1, not a new optimizer or candidate. Box activity uses inherited1e-8 reporting tolerance.

Raw sigma, delta and residual ratios use the original public fields. Missing/nonpositive denominators yield null, never an invented noise threshold. Unknown termination/cost fields stay null. Original endpoint selection CPU and new validation CPU are separate; OS child queue cost is retained as provenance and not refunded.

A small gradient or projected gradient does not establish a strict/local minimum, global optimality, information nullspace, identifiability or recovery. Scientific interpretation belongs to the parent.

Arrays: states/<case>__<method>.npz contains m/J/Q/Y_pred/Y_observed/residual, g_data/g_prior/g_total, box gradient mapping and active masks; no truth arrays. Reports: ENDPOINTS.csv/json, COST_EVENTS.csv, CASE_COSTS.json, FAILURES.json, RECEIPT.json.

Actual validation CPU at receipt: 1.142173000s; whole-process CPU at receipt 1.167117000s; wall 1.263929667s. Physical counters: {"DATA_ADJOINT": 44, "FULL_FORWARD": 44, "LU": 44, "actual_complex_RHS": 1056, "adjoint_RHS": 528, "cached_operator_imports": 4, "full_factorization_dimension": 25344, "source_RHS": 528}
