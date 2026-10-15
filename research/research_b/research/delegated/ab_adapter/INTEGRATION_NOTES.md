# Research A → B adapter：已实现的合同与限制

本项实现跨项目数据导入，没有运行 Maxwell、训练或主动测量。最初没有真实 A 包，随后只读定位了 A 的 32 个真实 confirm 包（193 候选），并按固定字典顺序导出前 4 包（27 候选）；E3 已具备真实输入，不能继续写 `WAITING_FOR_A_REAL_PACKAGE`。实际付费闭环由主会话单独执行。原 15 项合同测试中的材料、solver 和响应仍是明确标注的合成 fixture，与后续真实资产检查分开记录。

## 实际读取的来源

- B 原始 schema：`research/research_b/source_package/schemas/MaterialHypothesisPackage.schema.json`，版本 `research-b-mhp-1.0`。
- B 原始 recommendation schema：同目录 `AcquisitionRecommendation.schema.json`，版本 `research-b-acquisition-1.0`。
- A 原生 schema：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/Research_A/package_received_20261015/Research_A_Maxwell_Compensation/MaterialHypothesisPackage.schema.json`，版本 `1.0.0`。
- A `validate_package.py` 已完整只读审阅。没有执行其资产摘要选项，也没有计算或核对任何文件摘要。
- A/B 原 schema 和 A 工作目录均未修改。

## 不能直接兼容的内容

| 项目 | A 原生 | B 原生 | 本 adapter 的明确处理 |
|---|---|---|---|
| 版本字段 | `version=1.0.0` | `schema_version=research-b-mhp-1.0` | 固定版本入口，其他版本拒绝 |
| 资产 | `path/format/key/shape/dtype` | `uri/format/shape/dtype/units/provenance_id` | 只读 A 本地资产，导出独立 B 资产 |
| 单位 | noise 有 units，材料 asset 无单位 | 每个 asset 必须声明 units | 必须有独立 import context 的材料 quantity、units 和 unit evidence；不猜单位 |
| target | `physics_target.hash`、solver、backend 等 | `target_id` 和 model | context 显式给出预期 physics identity；逐项比较不计算摘要 |
| history | acquisition manifest，未强制逐 observation ID | acquisition history + data history | context 必须给出 observation/noise/provenance IDs 和不重叠 real channel indices |
| gain | 候选各自 fit 参数；schema 未给跨 history sharing graph | nuisance sharing graph 和 conditional model | 一世界一组共享参数，跨历史保持；不同世界允许不同值；未来 action 不得重置 gain |
| event | 字符串 event IDs | 非负整数 event IDs | 确定性映射保存在 contract；partition version 绑定 source package，不跨包误拼接 |
| weights | candidate_only / learned_proposal / corrected / reference | none / proposal / importance_corrected / independent_reference 等 | 不从 residual 或候选数构造概率；reference 必须额外证明 target/history 对齐和独立来源 |

## 可运行 API

代码：`research/research_b/code/ab_adapter.py`。

1. `adapt_a_v1(package, package_root, import_context)` 返回 `{package, contract, assets}`；此步只读输入并检查数组，不写文件。
2. `make_recommendation(bundle, action_ids, criterion=..., ...)` 生成原 B schema 的 recommendation；不替 acquisition 算法排名，也不制造未计算的风险。
3. `updated_data_history(bundle, new_record, updated_history_id=...)` 接收独立模拟/硬件进程已生成的新观测，返回用于 A/corrector 的完整新历史；不计算后验。
4. `export_bundle(bundle, new_output_dir, recommendation=None, new_history=None)` 导出独立核心数组和元数据目录；拒绝覆盖已有目录。历史源文件和账本不变。context 中未来 action 的 geometry/noise 外部 URI 保留原值，由 context 生产方提供可解析资产，不冒称已经审计。
5. `waiting_for_a_real_package()` 返回真实状态，不能用 fixture 代替 E3。

当前环境可运行：从仓库根使用已有 `../../.venv_nn/bin/python -B -m unittest discover -s research/research_b/tests -p test_ab_adapter.py -v`。没有安装依赖。

`import_context.adapter_version` 必须为 `research-a-1.0.0-to-b-1.0-context-1`。完整例子在 `IMPORT_CONTEXT_FIXTURE_EXAMPLE.json`；它只是接口例子，不能作为真实 A 输出。context 必须显式给出：

- `expected_physics`：源 target 的 hash（只作 opaque identity）、solver_id、backend_commit、material_law_id、time_convention、dispersion_model、self_term_convention、polarizability_map、boundary_conditions。
- `material_quantity/material_units/material_shape/material_domain/material_metric/unit_evidence`、observation_units、frequency_units、discretization_id。
- `noise_packing`：当前支持 `real_only`、`all_real_then_all_imag`、`interleaved_real_imag`。原生声明必须完全一致；不自动猜别名。复杂测量按指定规则 realify，联合 covariance 必须是相应的实协方差。
- `history`：逐观测 ID、noise realization ID、provenance ID、完整且不重复的 real channel indices、原 B action 结构。`native_acquisition_id` 必须匹配 A 原 acquisition ID，历史频率必须覆盖源文件全部已条件化频率。
- `nuisance`：参数名、每参数共享 group、合法 bounds、proper prior 定义和来源证据。`sharing_rule=one_parameter_per_world_across_history`，`different_worlds_may_use_different_values=true`。
- `available_actions` 使用 B 原 action 结构；旧 action 及已条件化 observation 的 constituent 不能冒充新数据。
- `weight_base_measure` 和后定义 event 的 `event_definitions`。可靠 reference 还需 `reference_evidence` 的 target/history/provenance 和 `independent_reference=true`。

candidate-only 可以报告显式 `minimax_finite_candidates` 或 `declared_design_prior_not_posterior` 风险；不允许把它叫 posterior EIG。只有 source 权重声明 full support、总质量及 represented mass 为 1、unrepresented mass 为 0、独立来源匹配，且 A 原接口明确允许 posterior information gain 时，才开放 Bayesian risk。candidate-pool 条件权重和 partial region mass 会被保留在 sidecar，但不会自动重归一化成全后验。

`new_record` 必须提供 target/solver、units/packing、唯一 observation/noise/provenance ID、合法 action ID、原共享 nuisance group、source、noise evidence、新 values 和完整 old+new joint covariance。旧 covariance block 必须原样保留；若行动声明独立新噪声，其交叉 covariance 必须为零。更新使用完整联合 history，每个原始观测一次；不能把包含旧数据的 posterior 再乘一次旧 likelihood。

## 科学范围和未完成事项

- 检查源 candidate 的 full-discrete-Maxwell 记录、solver/target/material-law、跨源同一材料声明及保存预测与 whitened residual 是否一致。此为记录和数值一致性检查，不是新的 full-Maxwell recertification。
- 检查材料数组 finiteness、shape、相对介电常数/电极化率/绝对介电常数单位，及复数被动性符号。没有把 B1 的 13-cell 模型硬编码成 A 的材料支持。
- A 原 schema 本身无法独立证明材料单位或 gain sharing；缺失这些外部事实时，必须补交 context，不能把字段完整当成物理证明。
- 当前只读 `npy/npz/json`；hdf5/zarr、未知 packing 和未知版本明确不支持。没有通用 schema/solver 框架，没有为 A 重写 forward。
- covariance 的完整矩阵只使用一次，B 全局 noise 的 covariance 和 cross-history 引用的是同一完整资产，不是两份独立证据。
- 保存基础重建材料与参数用于后续 base-only 对照。更新标 `PENDING_INDEPENDENT_RECONSTRUCTION`，并明确固定候选权重更新不改变生成器 Best-of-K。
- 当前完成 15 项有实质意义的合同测试；各轮结果与准确运行成本见 `TEST_RECEIPT.json`。这些不算真实候选成像成功。

## 真实 A execution profile 与 B 独立菜单

- `adapt_a_execution_v1(package_root, expected_source=..., expected_observation=...)` 显式识别 `research-a-native-1.0.0-execution-no-digest-v1` 与 `research-a-common-exchange/1.0.0`，保留 source target 原 ID。
- 实际 solver 为 `B1-vector-Mie-tangential-boundary-matching`，backend commit `7ed4faef9bc3dd214d4de689c2203cc383f7495f`。仅作明确命名 alias；不凭 target 名称相似就套用 B posterior。
- `all-real-then-all-imaginary` 经过已声明的字面 alias 变成 `all_real_then_all_imag`，没有改动数组顺序。source/native/CommonExchange 的原观测与 A `CONFIRM_OBSERVATIONS.json` 的 online 观测逐元素完全相等。
- 材料为逐包不同的 exact radial union grid、complex relative epsilon；体积权重为 `diff(radius**3)`。绝不强制 13-cell。
- `prediction_i.npy` 是 A source-gain-applied 均值。`source_fitted_gains.npy` 和 `unscaled_old_fields.npy` 显式提供 `F=prediction/g_source`；未来要用 F 联合新旧数据重新 profile 同一个 gain，不能把 prediction 再乘 gain。
- 额外保存原物理坐标 `(r,h,er,ei)`、P2 latent 和原 decoder 的 effective geometry。真实边界案例的极小 h 会在原 `from_geometry→geometry` 浮点往返中成为 inactive atom；这不等于发现一个可分辨的小孔。原 exact material array 为准，诊断保存在 `PROFILE_COORDINATE_ROUNDTRIP_DIAGNOSTIC.json`。
- 原 A menu 只有 old1.6/new2.4。`extend_authorized_b_menu(bundle, frequency_covariances, authorization_id)` 读取 B 独立 overlay 和冻结 6-frequency menu；同源 solver/material/dispersion 才允许创建 B-only action IDs。旧1.6只可新增独立噪声 repeat，旧观测不会再次作为证据。没有借用 A 或旧 P2 余额。
- B-only snapshot 为 `research/research_b/A_IMPORT/confirm_00_B_menu_v3` 至 `confirm_03_B_menu_v3`。原 A menu IDs 与 B extension 声明均在 sidecar；六 covariance 使用 B 已保存 fixed-reference noise anchors。
- `B_MENU_EXTENSION_RECEIPT.json` 的真实资产检查同时拒绝错误 authorization、未冻结频率及错误原观测。所有这些导入和扩展为 **0 次物理**，不是一次主动测量收益实验。

## 结果图

`research/research_b/code/plot_results.py` 只读已保存结果，输出 E0/E1/E2 的 action/risk matrix、分选择池 Selected-MAE 和测量 block 成本曲线。每幅提供 PNG/PDF，无新增物理/观测模拟/训练。

- E0 图始终标明七个 exposed development cases、finite proxy reference、不能做盲人口推断。
- E1/E2 的两份噪声先在 parent 内平均。Mode 1 图指未加 posterior 权重的选材池；带星号 EIG/risk acquisition 仍属于 reference-assisted action selection，不能说整个流程纯 candidate-only。
- E2 复用 E1 parents，明示 exposed followup。初始 ranking 每个 parent/noise 只取首个有 utility 的记录，不把每策略重复 ranking 当独立样本。
- 成本曲线按每个预算上限保留每条 episode 的最后已保存结果；停止后仍纳入平均，x 为真正平均已付 measurement blocks。none 保持 0，two-step bundle 保持实际 2，不能虚增为 3 次测量。完整 Maxwell/setup/CPU 成本仍单列在主资源账本。
