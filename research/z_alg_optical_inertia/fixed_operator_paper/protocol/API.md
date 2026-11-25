# 固定算子研究原型统一协议接口

本目录实现输入协议 26–40 的 15 个入口，并逐条登记 09–25 的数据字段。可从 `api.py` 调用，运行 `selfcheck.py` 会用小控制或已冻结输入实际执行每个入口。数值严格性针对每条记录声明的离散算子。

- [接口代码](api.py) · [字段定义与历史映射](schema.json) · [自检回执](selfcheck.json) · [自检正负 frame](selfcheck_frames.npz)
- `schema.json` 的 `$defs` 给出统一交换字段。`x-protocol-field-map` 给出协议行号、原始记录及覆盖边界。`mapped` 对应已有等效字段或证据，`scope_limited` 对应只实现已声明分支，`partial` 对应历史记录中仍有缺字段。历史实验记录按其原格式保存，不能直接当成统一交换记录验证。
- 向 `ModelContract` 传入固定 `G,S,Bdict`、网格与单位。`Acquisition` 持有系数 `V`、观测 `Y`、接收坐标标识、每束能量与噪声来源，可带已声明的入射功率 Gram。接收数据与 `S` 必须已处于共同白化坐标，通用接口接收仪器提供的白化与功率标定。材料节点持有逐格实数边界，已知行写成等端点。框架冻结普通数组副本以保持数据标识一致。
- 下表中「严」表示经 Arb 检查且只对当前证据边界成立，「提议」表示候选输入还需要严格复核。接口对失败或未覆盖的受支持分支返回 `UNRESOLVED`，无效输入则报错。

| 协议 | 入口和返回 | 复用及可用范围 |
| ---: | --- | --- |
| 26 | `assemble_border(model, box, weights, acquisition)` 返回 `K,L,C`、全边框及合同标识 | `src/fixed_trace.py:assemble_border`。`weights=(w, optical)`，实材料，精确存储矩阵，光学乘法次序及源块保留。 |
| 27 | `local_pair_rows(model, indices, acquisition)` 返回可解构的 `(T,E)` | `src/fixed_trace.py:local_pair_rows`，每格先电流后场。入射项依赖已测组合，因此显式给 `acquisition`。 |
| 28 | `build_fixed_kernel(base, rows, verification=None)` 返回核或未决 | `src/fixed_trace.py:build_fixed_kernel` 严格区间枢轴与正补计数。通用核有精确函数行证明；A 的缓存另有正补数值下界。近零不删除。 |
| 29 | `update_box(kernel, local_D, acquisition)` 返回小边框 | `src/fixed_trace.py:update_box` 用增广迹保留线性项和源角项。此底层接口接收调用者给定的局部二次型；用于节点认证时 `certify_box` 从实际端点构造它。 |
| 30 | `certify_inertia(border, (t,j,r), frame_proposal=None)` 返回回执与 frame | `B/certify.py:sign_floor`。浮点谱仅提出正负 frame，Arb 验证投影正性与自身 Gram。差额是未决符号，不能推断精确零。 |
| 31 | `certify_box(node, acquisition, kernel_bank, budget)` 返回节点状态和检查记录 | 通用版接收有限的 `trials=[{weights,t,j,kernel?},…]`。全部变化都在已有物理迹中才复用，否则回到原边框。可另给 `material` 做独立见证校验。批量自适应设计沿用 A/B 原入口。 |
| 32 | `extend_traces(kernel, new_rows)` 返回联合核及重建耗时 | 用 28 对完整联合行重新验证，共享交叉项，未实现一般增量复杂度结论。命名行重叠时保留旧核并报告待显式去重。 |
| 33 | `exact_row_oracle(node, cells, a_acquisition, a_trace_bank, budget)` 返回覆盖叶与 frame | `A/compressed.py:interval_search`。使用 `prepare_a_chart_kernel` 绑定 A 的 `TraceBank` 及采集几何，逐格对照迹外背景。该分支使用 A 登记的迹行和恒等枢轴图。通用仿射矩形等式可用下方的区间图接口构造；此处的完整材料盒排除仍使用 A 的已验证图。 |
| 34 | `build_task_cover(task_contract)` 返回根节点、`refine(node)` 与覆盖说明 | `B/run.py:Contract.cover_box,members,quadrants`。有限格点中心，根盒是所有叶盒的坐标逐格外包。通过 B 的 `arbitrary_high` 模式枚举指定窗口的高值位置。 |
| 35 | `search_witness(node, acquisition, fixed_budget)` 后接 `verify_witness` | B 的非负方形网格有限预算搜索仅作提议，随后以原状态方程严格检验逐格边界、可逆与原噪声。直接 `verify_witness` 支持原存储模型的全格实材料。 |
| 36 | `build_response_envelope(node, certificate, c_acquisition)` 返回有界椭球或未决 | `C/design.py:form_border,envelope`。仅 C 的解析 TM、固定两束旧照明和经证明的联合 Frobenius 半径。奇异描述形式与不能证明的 `R_y` 返回未决。 |
| 37 | `score_action(action, envelopes, noise_contract)` 返回方向、Arb 分离边与最坏存活上界 | `C/design.py:split_query,certify_pair`，共享多束源 Gram；独立查旧、新物理源联合秩。该物理入口是 C 两个材料类；一般多类组合层独立处理。 |
| 38 | `select_and_freeze_action(experiment_path=C, diagnostic=False)` 写入或读取核对 `proposal.json` | C 设计在独立进程中运行，只读模型和旧数据；未冻结先采集报错，现有读数须能对应冻结标识。采集仍由 C 隔离仪器进程另行调用。 |
| 39 | `append_observation(old, new, expected_old_hash)` 返回联合采集 | 最小通用适配。拒绝旧标识、模型、接收基底或入射 Gram 的不符，保持旧 `V/Y` 原样并合并不重叠列组的误差平方。 |
| 40 | `replay_certificate(record, original_inputs)` 返回独立全边框复核结果 | 用 A 的物理迹二次型重装原边框并验证原保存的 frame，复核模型、数据、盒、权重及半径。不具备一般压缩或分层证书的反序列化，见 A/B 独立重放入口。 |

38 使用 C 的 `registration.json`、`public_model.npz` 和 `old/`。36 接收 C 布局名、旧 `source_indices`、`per_beam_energy`、`Y` 与 `noise`。旧阈值只使用经证明的 \(r_{\rm old}^2 I\)。37 接收实际新角的索引和每束物理能量，按两个包络的传入次序检查区分方向。多类时可用 `clique_cover_bound(labels, pair_checks, areas=None)` 聚合独立验证过的边。物理区分保证仍由实际端点检查提供。通用 `append_observation` 要求调用者已经按声明的功率采到新数据，C 的实际冻结和采集另留逐束回执。 [多类标签控制](../experiments/C/receive_socp/README_ZH.md)

## 严格重放与运行

原实验的原始入口和回执继续有效。A 使用 `experiments/A/replay.py` 复核盒端点、原状态和局部图。B 使用 `experiments/B/replay.py` 验证覆盖树和原边框。C 使用 `experiments/C/menu_replay.py` 只读公开旧数据重放动作，`readout.py` 再执行冻结支持测试和原联合尾界。历史全链路回执和新统一适配自检分别计数。

从项目根用已安装 NumPy、SciPy、python-flint、threadpoolctl 的 Python 环境执行：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B research/z_alg_optical_inertia/fixed_operator_paper/protocol/selfcheck.py
```

自检覆盖 15 个入口，包括有与无压缩的同节点决定、源角更新、冻结前重新生成提案、旧数据追加身份拒绝以及独立全边框复核。小控制中的未决图层仍作为合法回执留存。运行只写 `protocol/` 下的两个自检文件，C 提案重建在临时副本内进行。

## 矩形等式图

`api.rectangular_equality_charts(AffineEqualities(A0, B0, A_coeff, B_coeff), lower, upper, ...)` 接收有限二进制系数与共享实参数盒。输入是 \(A(\xi)x=B(\xi)Z\)，两个矩阵都可随参数仿射变化。显式指定 `pivot_size` 可消去一部分矩形等式。优先尝试指定枢轴、QR 建议和有限的备选子块。浮点逆仅作提议，以 128 位测试中的 Arb 全盒误差 \(\eta=\sup\|I-CP(\xi)\|_\infty<1\) 证明枢轴均可逆。回执同时给出逆范数和条件数上界。

`evaluate_chart` 可以对同一输入的一个参数值重算回代。`pullback_border` 用 `[w;Z]` 到 `[x;Z]` 的完整映射拉回原带边形式和源阈值，同时返回每个非枢轴行对应的一致性矩阵。精确等式图要求施加这些剩余方程。结果里的区间矩阵外包逐点的有理映射，区间分量之间允许相关。枢轴失败按闭盒二分，深度、叶数或时间预算用尽就保存未决叶。时间预算在尝试间检查。`OUTWARD_CHART` 或 `OUTWARD_COVER` 只表示代数外包，材料排除仍需依合同验证带边正性和全覆盖。

使用 `python -B research/z_alg_optical_inertia/fixed_operator_paper/protocol/chart_check.py` 运行独立小控制。测试含矩形稠密坐标与恒等枢轴的直接解、两个源、剩余行拒绝错误候选、阈值保持、需要更换枢轴的闭覆盖、真实秩下降和预算耗尽。秩下降示例 \(\xi x=z\) 在 \(\xi=0\) 时允许零源有非零状态，同时禁止非零源。这两个结论不能用相邻盒的逆外推。[自检与区间数据](chart_check.json)

## 历史数据可补录范围

| 协议 | 可从保存记录重建 | 无法还原的历史事实 |
| ---: | --- | --- |
| 13 | 已保存的源、接收与记录的度量。C 可重建每束能量和已施加的接收归一化。 | 从未记录的实测接收协方差、仪器标定和实际能量遥测，以及缺失原度量时的物理 Gram。 |
| 20 | 保留的 frame、间隙及根据原输入独立复核的新证书。 | 未保存的原求解器 frame、尝试次序或实际舍入历史。新重放无法充作同一次原调用。 |
| 21 | 已记录的总计、分项耗时、精度和内存峰值。 | 未记录的每次调用耗时、RSS 峰值、临时内存或废弃尝试时序。 |
| 25 | A/B 在许可的离线评分阶段由原材料和冻结映射补齐身份。C 逐路径保存四项字段。 | 离线补齐的标识不能追溯证明历史在线行在采集时已经记录该标识。 |
