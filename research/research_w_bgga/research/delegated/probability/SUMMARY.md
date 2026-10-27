# 有限 prior 区间概率核：实现交接

状态：代码与标准库 unittest 测试已写；静态语法解析通过。数值测试 **NOT_RUN**，由主线程统一执行、计费并保存 receipt。没有运行物理求解、完整 posterior reference、数据重放或新实验；本交付不构成 G1 通过。

## 实现依据与范围

已阅读 `source_review/PRO_REPORT_EXTRACTED.txt` 第 11 节（885–973 行）及附录 G1（1633–1720 行），并阅读相邻 `research_w_slam/src/wslam/finite_bayes.py`、`factors.py`。旧实现只作为已有有限 prior 与唯一观测 factor 语义参考；本任务没有修改旧文件。固定旧版本 `52351aee6ddd41b09851fe5e52d1a8861891ef13` 的整体一致性由主线程维护，本 worker 未另做版本或文件哈希核验。

新增代码仅依赖 NumPy 和标准库；概率核没有文件 I/O、reference、GT、未来观测或随机 oracle 参数。它处理用户声明的有限 joint prior；不是连续 prior 求积或 Maxwell/PDE 证书。

## 接口

实现：`src/bgga/intervals.py`。

| 接口 | 作用 |
| --- | --- |
| `LogLikelihoodBounds(log_lower, log_upper, labels=EMPIRICAL_INTERVAL, queried=None, provenance=...)` | 原子 likelihood 的 log 上下界；数组复制并设为只读；显式误差界默认经验标签 |
| `LogLikelihoodBounds.unqueried(n)` | 默认 likelihood `[0,1]`，即 log `[-inf,0]`；`queried=False`，保留全部 prior 原子 |
| `LogLikelihoodBounds.exact(log_likelihood, label=MATHEMATICAL_INTERVAL)` | 把调用者提供的 likelihood 视作声明的有限模型中的精确值；可保存 `-inf` 零支持 |
| `bounds.with_updates(indices, replacement)` | 替换指定原子，不改其余原子；拒绝重复/越界索引 |
| `likelihood_bounds_from_residual(r, e, label=EMPIRICAL_INTERVAL)` | 使用白化残差三角不等式产生 Gaussian likelihood 界；`e=inf` 恢复 `[0,1]` |
| `likelihood_bounds_from_potential(phi, eta, label=EMPIRICAL_INTERVAL)` | 使用 `abs(Phi-phi)<=eta`，并由 `Phi>=0` 把 likelihood 上界截为 1 |
| `combine_likelihood_bounds(factors, size=None)` | 不同观测因子在 log 域相加；身份唯一性需调用者保证，通常应使用下面的 registry |
| `finite_posterior_bounds(prior_weight, bounds)` | 返回 evidence log 界、逐 joint 原子的 posterior 界、标签；所有受支持 likelihood 精确时另保存 joint/log-joint 权重 |
| `posterior.event(boolean_mask)` / `event_posterior_bounds(...)` | 用事件与补集的分式界处理归一化分母 |
| `query_log_scores(prior_weight, bounds, cost=1, eligible=None)` / `select_query(...)` | prior × likelihood 区间宽度 / 预声明成本的固定选择规则；同分选首个索引，无正宽度时返回 `None` |
| `FinitePriorIntervals(prior_weight)` | 保存原始 log prior 与每个观测的当前 factor 版本 |
| `state.set_factor(id, bounds=None, observation_id=None, version=1)` | 新观测注册未知界；同一 factor 的替换必须保留 observation ID 并提高版本；重复 observation ID 报错 |
| `state.refine_factor(id, indices, bounds, version=None)` | 更新某些原子并递增 factor 版本 |
| `state.combined_bounds()` / `state.posterior()` / `state.query_log_scores(...)` / `state.select_query(...)` | 从 base prior 与当前唯一 factors 重建；不会重复乘入旧版本，保留 normalized 浮点权重无法表示的 prior log 尾部 |
| `joint_pushforward(joint_weights, cluster_ids)` | 保存原 joint/log-joint 权重、cluster mass、逐原子簇内条件权重；返回对象 `.lift()` 重建原联合分布 |
| `canonical_lift(cluster_mass, prior_weight, cluster_ids)` | 以原 prior 条件核 lifting quotient mass；返回 `joint_error_label=NOT_CERTIFIED`，不会把完整 joint TV 宣称为 0 |
| `task_cluster_diameters(task_values, cluster_ids, task_contract_id=...)` | 在声明的任务坐标中计算每簇全部点对的 Euclidean diameter，保存任务合同 ID、成员数与有限点范围 |

cluster mass 的顺序是 `np.unique(cluster_ids)` 的排序；ID 可为整数或字符串。`task_values` 第一个轴必须与原 joint 原子对应，其余轴展平为任务坐标；单位、白化和任务合同由调用者声明。

## 数学与数值约定

对每个原子，用 `u_i^- = p0_i L_i^-`、`u_i^+ = p0_i L_i^+`。事件 A 的界为

```
P(A)^- = U_A^- / (U_A^- + U_Ac^+)
P(A)^+ = U_A^+ / (U_A^+ + U_Ac^-)
```

零 evidence 下界可以保留；若 evidence 上界也为零则报错，不能虚构 posterior。完全不可能/必然的事件分别返回 `[0,0]`/`[1,1]`，其余未查询事件通常得到 `[0,1]`。未查询集合的上界因此自动包含报告中的 `U_num / (L_computed + U_num)`。

不在概率域提前计算 `exp(-Phi)`；似然 log 可小于 `-10000`。posterior 计算先减去共同 likelihood log offset，再加入 log prior，避免共同 log likelihood 为 `-1e100` 时抹掉 prior 比例。为了保留这一信息，`FinitePosteriorBounds.log_mass_lower/upper` 和 `EventPosteriorBounds.log_numerator_* / log_complement_*` **省略共同的 `log_mass_offset`**；还原绝对 log mass 时加回该 offset。`log_evidence_lower/upper` 已加回 offset。极端绝对 evidence log 的小量仍受浮点表示精度限制，posterior 比例使用平移后的量计算。

补集 log 总质量用 prefix/suffix `logaddexp` 计算，避免从几乎相等的归一化总量相减。`query_log_scores` 同样省略一个共同 upper-likelihood log offset，仅排序有意义。选择策略是区间贡献启发式，不是最优风险查询定理。实际使用 registry 时宜调用 `state.select_query()`，以保留 `state.log_prior_weight` 中的极小尾部。

`joint_weights` 仍可能在浮点显示时下溢为 0；同时保存的 `log_joint_weights` 保留有限尾部。组合因子的 upper log 若因超出浮点表示范围而溢为 `-inf`，代码报错，避免错误宣称 likelihood 精确为零。

## 标签与未完成证书

- `MATHEMATICAL_INTERVAL`：调用者声明的数学输入界及其有限模型代数后果；默认未知 Gaussian `[0,1]` 属于此类。
- `EMPIRICAL_INTERVAL`：任何 empirical forward/potential error；乘积、归一化和事件界不会自动将其提升为数学标签。
- `NOT_CERTIFIED`：输入界本身未获声明认证，或 canonical lift 尚未证明完整 joint fidelity。标签按最弱输入传播，与因子顺序无关。

posterior/event/任务直径结果显式保存 `floating_point_roundoff_certified=False`。本实现未提供 directed rounding，也未验证 Maxwell 稳定常数、物理 residual 到 forward error 的推导、误差界有效域或未来任务一致性。因此即便标签为数学区间，也不可将本文件称为严格有限精度物理证书。

联合 pushforward 保留原 posterior 条件核时可重建原分布；canonical lift 只使用 prior 条件核，一般不会恢复 posterior 条件相关。task diameter 逐簇检查所有点对，不把相邻距离阈值的链式连通当作 diameter 合格。任务合同改变后须由外层重新检查/重开；本模块不设新的实验科学 gate。

## 测试交接

`tests/test_intervals.py` 提供 18 个 `unittest.TestCase` 方法（无 pytest 依赖），覆盖：

- 小复数白化矩阵的 forward 三角界与经验标签；potential 上界截断和无限误差 fallback。
- 未查询 prior 支持、遗漏质量分母、16 个 interval-box 角点的事件分式紧界。
- 高负 log likelihood、极小事件、共同巨大负 log、零 prior、零上 evidence 及结构事件。
- prior/区间宽度/成本选择、确定性同分、浮点下溢尾部和精确原子不重查。
- 多 observation batch 等价、同数据版本替换、重复观测拒绝与旧零支持恢复。
- 经验/未认证标签混合顺序、atom 标签替换。
- 保留墙目标条件相关的 joint pushforward、prior canonical lift 与零 prior 簇拒绝。
- 三点邻接链的全 diameter 与任务合同保存。

本 worker 仅通过标准库 `ast.parse` 检查两个 Python 文件，未导入概率模块，未执行上述数值断言。静态 receipt：`research/delegated/probability/STATIC_CHECK.json`。

建议主线程使用项目统一的 unittest discovery 命令；例如从本目录在已有 NumPy runtime 下运行 `PYTHONPATH=src python -m unittest discover -s tests -p test_intervals.py -v`。这条命令是交接说明，尚未由 worker 执行。

## 修改文件

1. `src/bgga/intervals.py`
2. `tests/test_intervals.py`
3. `research/delegated/probability/SUMMARY.md`
4. `research/delegated/probability/STATIC_CHECK.json`

主线程负责 `src/bgga/__init__.py`、统一测试运行、计费、真实数据/物理接口及科学裁决。

## 2026-10-27 追加：明确任务范数与主线程测试记录

前一交接之后，主线程已运行项目统一测试。已只读核对 `results/jobs/implementation-tests-final/OUTPUT.log` 的 `Ran 30 tests ... OK`，以及对应 `RECEIPT.json` 的 `exit_code=0`、`process_cpu_seconds=1.195705`、`new_pde_authorized=false`。这是**本次范数修正之前**的 30 项实现测试通过记录；不会替换或改写原 receipt。它支持此前实现的有限代数/接口检查，不是 G1 通过，也不是物理或统计证书。

本次仅修改 `task_cluster_diameters` 及其结果的范数 metadata、相关测试和本追加说明。其余概率逻辑与文件未改。此处接口替代上表旧的 task-diameter 签名：

```
task_cluster_diameters(task_values, cluster_ids, *, task_contract_id, norm)
```

`norm` 为必填关键字，避免把拼接 L2 默认用于 `TASKS.json` 声明的 `max_j whitened response norm`。

- `norm="max_probe_l2"`：输入必须为调用者已经白化的实数数组 `(n_atoms, n_probes, n_real_measurements)`。每对原子先计算各 probe 的 measurement L2，再取 probe 最大值；最后对簇内全部原子对取最大值。没有明确 probe 轴的 1D/2D 输入拒绝，复数输入也拒绝；函数不代替调用者执行白化或把复数响应展开为实数测量。
- `norm="euclidean"`：保留原先全部任务坐标拼接后的 Euclidean L2；当调用者预先把坐标乘入所需权重时，其含义才是这些坐标对应的 weightedL2。没有隐式权重或白化，也不把该结果标为 max-probe 范数。

返回结果新增/明确 `norm`、`metric`、`norm_definition`、`coordinate_assumption`、`probe_axis`、`measurement_axis` 和 `input_shape`。它仍保留任务合同 ID、原子簇全点对 diameter，以及 `floating_point_roundoff_certified=False` 与有限点作用域。

新增两项 unittest 方法：两 probe 各具有同一单位响应差时，`max_probe_l2` 的 diameter 应为 `1`，拼接 `euclidean` 应为 `sqrt(2)`；同时拒绝缺少 probe 轴和未经实数测量表示的 max-probe 输入。旧的邻接链 diameter 测试改为显式选择 `euclidean`。本文件现在有 20 个测试方法；本次新增/修改测试 **NOT_RUN**，worker 未执行测试或物理求解，由主线程再次统一运行并计费。此前 `STATIC_CHECK.json` 是原始交付的静态记录，本次未改写它，也不把此前 30 项通过记录当作本次修改后的测试结果。
