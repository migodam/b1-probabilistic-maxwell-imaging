# Research A 原生交接与独立 common exchange v1

`src/research_a/contracts.py` 只做数组写入、接口转换和一致性验证，**不会调用 Maxwell solver**。调用方必须先完成真实物理验证；本目录的测试 fixture 不是实验发现，也不会被导出为正式候选。

## 原生格式与用户 no-hash 指令的具体冲突

收到的 `MaterialHypothesisPackage.schema.json` 强制要求每个 asset 的 `sha256` 为64位十六进制字符串。该要求不能靠“省略一个可选字段”解决。本实现：

- 在 `original/` 原样保留收到的 schema、validator、原测试和原测试结果；没有改写附件。
- 使用独立 `MaterialHypothesisPackage.execution.schema.json`，schema 变化包括：asset.sha256 额外允许 null；analytic-Mie 未暴露内部场时 constitutive_relative_defect 允许 null，并要求显式状态。输出保留原生1.0.0字段布局，但明确标为 `research-a-native-1.0.0-execution-no-digest-v1` profile。
- 不计算、不核验，也不编造摘要字符串。`manifest.json` 的 asset.sha256 全为 null；`CONTRACT_AUDIT.json` 清楚记录与原 strict schema 的差异。
- `physics_target.hash` 沿用原生字段名，但值是调用方冻结的语义 target identity，不是假造的文件摘要。

原 `validate_package.py` 依赖当前运行环境中未安装的 jsonschema；本实现没有为此安装依赖。stdlib+NumPy validator 明确只覆盖保存 schema 实际使用的关键词，不冒充完整通用 JSON-Schema 引擎。16项接口/语义检查的结果在 `ZERO_PHYSICS_CONTRACT_CHECKS.json`，物理调用为0。

## Factory

```python
from research_a.contracts import write_packages
result = write_packages(output_dir, context, accepted_candidate_rows,
                        base=baseline_row, pair_distances=pair_rows)
```

`output_dir` 建议在本轮 `MaterialHypothesisPackage/<observation>/<method>/`。返回值含 manifest/exchange 的文件路径及对象。**只允许真实已接受候选**；空列表、超阈值残差、缺少已验证标志均报错。不会用 fixture 填充空 archive。

### context 必需字段

- `package_id, target_id, solver_id, backend_commit, material_law_id, time_convention`：真实版本与物理身份。时间 convention 必须与 solver 一致。
- `acquisition={id, frequencies:[...], frequency_units, n_sources}`。`n_real_measurements` 可省略，由观测长度产生。其余描述放单独 metadata，不往严格 acquisition 字段里任意增加键。
- `observation`：已声明 Re 后接 Im 的一维实数组；`noise_covariance` 与 `whitener` 是同维方阵；验证 `W Σ Wᵀ≈I`。
- `noise_id, realification_order`（默认 `all-real-then-all-imaginary`）, `noise_units`（可选）, `model_discrepancy_policy`。
- `grid, material_support, nuisance_support, prior_manifest`：JSON字典。`prior_manifest` 可为 None，因为 candidate-only 不强求概率先验。
- `self_term_convention, polarizability_map, boundary_conditions, dispersion_model`：具体物理约定，不用不适用的框架术语覆盖真实后端。
- `fit_threshold, requested_k`（后者可省略，默认输入候选数）, `distance_definition_id, candidate_ordering_rule`。
- `unexplored_regions:[str], search_restrictions:[str]`：本次实际未探索的支持与结构。完全覆盖或后验剩余质量均不由 adapter 推断。
- `total_cost`：本轮已计费累计成本字典；`upstream_manifest`：真实资产和版本 provenance；`failed_attempts`：失败记录 JSON列表或字典。
- `new_action_interface={evaluator_entrypoint, environment, allowed_actions:[{id,...}], cache_manifest}`。必须指向实际可调用的后端封装或明确的 caller-supplied callable 接口。新增action不可自动使用本轮余额；需要单独budget token和真正记账backend。

### 每个 accepted candidate row 必需字段

- `id`：文件安全的唯一字符串；`validated_full_maxwell=True`；`physics_target_id` 可省略（默认context target）。
- `material`：完整材料表示对应的 complex NumPy数组；`material_parameterization_id`。对同心球可以是**足以唯一重建全材料的径向分层数组**，并在grid/coordinates中存完整边界、每层ε；中心切面只是展示，不能单独冒称任意三维完整材料场。
- `prediction`：**已应用本候选共享gain**的旧 acquisition full-Maxwell 实向量；不可传 surrogate prediction 冒充 full。adapter 计算 `||W(prediction-observation)||`，必须重放 `whitened_residual_l2`。
- `whitened_residual_l2`，可选覆盖 `fit_threshold`；`nuisance_policy` 默认 bounded_profile；`nuisance_parameters` 如 `{gain: scalar}`。一个材料世界的一项gain必须跨所声明频率/源共享，不能给每个通道独立gain。
- `validation={constitutive_relative_defect, constitutive_relative_defect_status, solver_relative_residual, material_constraints_passed:True, all_illuminations_same_material:True, numerical_error}`。solver_relative_residual 必须来自真实求解器；Mie 的边界匹配 backward_residual 可以填写此项并注明规范。原Mie接口没有内部场/current，故 constitutive_relative_defect 必须为 null，且 status=`NOT_EXPOSED_BY_ANALYTIC_MIE_INTERFACE`。不得把边界残差映射为本构残差，也不允许编造0。共享scalar ε属于材料构造约束，不是实测内部场缺陷。实际别的backend暴露并测量时才可填非负值并注明status=`measured`。`numerical_error` 见下；独立网格核验缺省 false，不能由格式验证产生true。
- `structure={support_volume, cavity_count, component_count, feature_scale, descriptor_definition_id, cluster_ids:[str], event_ids:[str], topology_thresholds:[number], topology_persistence_rule}`；可加 `descriptor={...}` 原始几何描述，其余键按原生 schema 限制。
- `compensation={edit_family, edit_parameters:{}, material_compensation_parameters:{}, nuisance_compensation_parameters:{}, coordinates: ndarray}`；可加 `parent_candidate_id`（缺省null）, `protected_structure_constraints:[str]`。只有真正保存且逐点full验证的连接才能给pair声明 verified_path。
- `material_distance_to_base`：按冻结材料度量真正计算；measurement distance 则由gain后完整预测自动计算。
- `cost`：本候选及产生它的预算归属，不能把失败和拒绝漏记；`history:[...]`：生成历史，进入common exchange。

可选 `proper_log_prior` 与 `prior_score_is_normalized_density` 只在正确已知时提供；`proposal_log_density`、`posterior_weight` 在 candidate-only 一律null。对数likelihood也不在adapter中从残差临时捏造。

### cost 与 error 字典

`zero_cost()` 只用于确实未发生的工作。真实调用必须填写：

```
new_full_state_frequency_attempts, vector_rhs, forward_rhs, adjoint_rhs,
tangent_rhs, cache_hits, proposal_count, rejected_count, failed_count,
wall_seconds, peak_memory_bytes (可null), setup_allocation_seconds,
historical_cost_ledger (可null)
```

原生schema将所有计数字段设为必需；factory会为未指定的项目补0，调用方应明确确认这些项目不适用/没有发生。共享费用计在total，候选分摊规则必须在ledger中记录，不能把逐候选费用与total重复累加。

```
numerical_error = {
  'value': measured_value_or_None,
  'norm': '明确规范',
  'status': 'empirical_estimate' / 'solver_tolerance_only' / 'unknown',
  'method': '明确测量/计算方法'
}
```

不要将采样/网格检查称为 certified_bound。

### pairs

可提供完整pair列表，每行必需 `candidate_i, candidate_j, material_distance, material_metric_id, measurement_distance, measurement_metric_id`。adapter默认 `connection_status=not_tested`、未声明 hard_to_distinguish。候选距离是数值尺度，不能直接冒充后验概率或连接路径证明。

## 与 Research B 的独立共同交换

`common_exchange_v1.json` 由 `native_to_exchange` 生成，遵循本目录 `CommonMaterialExchange.v1.schema.json`。它**不声称**符合未读取的 Research B 原生 schema；B 需要显式读取此版本或自己写版本适配器。本任务不访问B工作区、不发送消息。

该格式保留 target/acquisition/noise、完整材料/旧观测预测、验证、结构与pair距离、补偿参数、共享gain规则、生成历史、每候选成本。所有概率/质量均null，未探索及模型外世界unquantified，不能进行 posterior information gain。

### 新测量动作接口

`evaluate_new_action(exchange, root, material_world_id, action_id, budget_authorization, evaluator)` 是显式调用接口：

```
budget_authorization = {
  'authorized': True,
  'remaining_attempts': 明确剩余额度,
  'ledger_id': '单独授权的账本ID',
  'target_identity': exchange['target_identity']
}

evaluator(material_array, action_dict, budget_authorization) -> {
  'target_identity', 'action_id', 'realification_order',
  'prediction': real_packed_array,
  'cost_receipt': actual_receipt,
  'validation_status': actual_status
}
```

序列化、验证、adapter均不会执行该接口。接口必须由真正受budget约束的physics backend实现预约和计费；上层budget token并非可绕过后端预算的安全凭证。未知action、无授权/无账本/无余额、target不符会在调用前拒绝；预测非有限、错packing或回执超额在返回后拒绝。候选已缓存的旧观测不需要新增求解。

没有实际新增测量实验时，接口的物理求值状态必须标NOT_RUN，不能用synthetic callback测试当物理验证。

## Analytic-Mie execution overlay补充

原strict schema没有表达“内部本构缺陷未暴露”的状态。本轮profile诚实保留null+显式原因，因此不宣称与原strict schema完全兼容。boundary matching残差与 constitutive defect是不同物理量；16项0物理检查包含将未暴露缺陷硬写0的负控制。

## 最小 Mie helper

`make_mie_context(...)` 生成当前B1球对称Mie固定物理约定、对角协方差和白化矩阵，减轻重复填表；**不生成prior、不生成预算、不求解物理**。必传参数为：

```
package_id, target_id, backend_commit, observation, sigma,
fit_threshold, acquisition, grid, material_support, nuisance_support,
total_cost, upstream_manifest, failed_attempts, new_action_interface,
distance_definition_id, noise_id, requested_k
```

`sigma` 是冻结协议的每个实通道SD，可为标量或与观测等长向量；不能传候选响应RMS。`grid` 必须声明完整径向层界面与材料数组语义；供画图的中心切面另外保存。固定 helper 只适用非磁性被动球层Mie，不可给A30二维VIE混用。其他字段可通过关键字覆盖，但必须对应真实冻结协议。
