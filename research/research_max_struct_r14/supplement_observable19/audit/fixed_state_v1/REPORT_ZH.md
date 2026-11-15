# 固定同状态解释性诊断

固定7案例、internal_galerkin readout1/128/512实际RETURNED输入；缺失/中断不补跑。所有结果OFFLINE_DIAGNOSTIC，不进入主候选池，不读取GT。

原完整144终态父授权seal与逐项冻结manifest/inventory先核验，之后才创建此新输出。六个K直接从原档U/B载入，无新basis/DPCS SVD或G@U setup。每共同m一次原物理state/gradient，保留LU后一次nr-adjoint RHS构造O_m；各原LBFGS8提案的物理Phi另外付费，不能拿这些validation费用当在线算法成本。

||O_m(I-P)e||/||O_m e||可能大于1，包含相位抵消，不能标为材料信息覆盖百分比。正交current能量、Green输入遗漏、PGP额外输出裁剪与DPCS重合均是current空间描述；没有D_Q材料切向就不作材料Jacobian/可辨识性结论。

小19阶T最小奇异值/条件数只是离线诊断，不据此改变rank、停止或拒绝原readout。所有失败/警告/零分母原因与费用保留。GT若以后用于解释，需父线程另在本诊断输出封存后授权附加；此脚本没有GT读取路径。

{
  "status": "COMPLETE_WITH_RETAINED_MISSING_OR_FAILURES",
  "diagnostic_kind": "OFFLINE_DIAGNOSTIC",
  "parent_authorized_seal": "D:\\AI\\MAX_STRUCT_R14_DUALMETRIC_20261115\\supplement_observable19\\predictions\\PREDICTION_SEAL_APPROVED.json",
  "same_frozen_config": true,
  "truth_access": "NONE",
  "fixed_cases": [
    "r14_new_0000",
    "r14_new_0002",
    "r14_new_0009",
    "r14_new_0012",
    "r14_new_0015",
    "r14_new_0020",
    "r14_new_0021"
  ],
  "fixed_readout_calls": [
    1,
    128,
    512
  ],
  "expected_source_slots": 21,
  "actual_common_source_rows": 21,
  "expected_arm_rows": 126,
  "actual_arm_rows": 126,
  "arm_status_counts": {
    "PARTIAL_DIAGNOSTIC_FAILURE": 21,
    "COMPLETE": 105
  },
  "common_states_with_physics": 21,
  "new_original_physical_cost": {
    "LU": 126,
    "DATA_ADJOINT": 21,
    "cached_operator_imports": 21,
    "FIXED_STATE_PORT_TRANSPORT": 21,
    "adjoint_RHS": 1116,
    "actual_complex_RHS": 2628,
    "source_RHS": 1512,
    "FULL_FORWARD": 126,
    "full_factorization_dimension": 72576,
    "diagnostic_port_transport_adjoint_RHS": 864
  },
  "new_auxiliary_cost": {
    "diagnostic_original_Phi_gradient_attempts": 21,
    "diagnostic_original_Phi_gradient_completed": 21,
    "diagnostic_original_physical_prior_evaluations": 21,
    "diagnostic_Om_dense_actions": 147,
    "diagnostic_Om_dense_complex_RHS": 1764,
    "diagnostic_dense_G_actions": 273,
    "diagnostic_dense_G_complex_RHS": 3276,
    "diagnostic_orthogonal_projection_actions": 252,
    "diagnostic_orthogonal_projection_complex_RHS": 3024,
    "lowrank_Green_actions": 1663,
    "lowrank_Green_complex_RHS": 19956,
    "diagnostic_DPCS_projection_actions": 378,
    "diagnostic_DPCS_projection_complex_RHS": 4536,
    "diagnostic_basis_DPCS_overlap_products": 126,
    "diagnostic_T19_assemblies": 126,
    "diagnostic_T19_SVD_attempts": 126,
    "diagnostic_T19_SVD": 126,
    "diagnostic_gradient_parts_attempts": 126,
    "denominator_small_LU": 1537,
    "denominator_small_LU_dimension_sum": 29203,
    "denominator_small_triangular_complex_RHS": 36888,
    "diagnostic_original_regularizer_evaluations": 1642,
    "diagnostic_gradient_parts_evaluations": 126,
    "diagnostic_gradient_parts_source_metric_products": 126,
    "diagnostic_readout_source_whitening_dense_products": 252,
    "diagnostic_LBFGS8_readout_attempts": 126,
    "material_readouts": 126,
    "material_inner_evaluations": 1411,
    "material_prior_evaluations": 1411,
    "source_metric_products": 1411,
    "matrix_material_readouts": 126,
    "matrix_material_LBFGS_iterations": 1008,
    "matrix_material_readout_cpu_microseconds": 6453125,
    "diagnostic_LBFGS8_readout_completed": 126
  },
  "missing_or_failed_components": 21,
  "new_basis_SVD": 0,
  "new_G_times_U_setup": 0,
  "T19_SVD_is_paid_diagnostic_only": true,
  "original_LBFGS8_is_offline_only": true,
  "main_predictions_or_candidate_pool_modified": false,
  "CPU_seconds": 77.921875,
  "whole_process_CPU_seconds_at_receipt": 77.96875,
  "wall_seconds": 86.5227201001253,
  "CPU_scope": "Imports, source reads, original physical state/gradient+nr transport, all lowrank diagnostics/readouts and proposal forwards, copies and output; eventCPU subscopes are not added again. Final receipt/seal serialization follows timestamp.",
  "interpretation": "Current defect/projector diagnostics, not material information percentage, material-Jacobian identifiability, recovery predictor, new algorithm or novelty judgment",
  "finished_utc": "2026-11-15T14:33:12.164015+00:00"
}
