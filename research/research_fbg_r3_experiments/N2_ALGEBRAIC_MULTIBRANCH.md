# N2：真实空间小支持的代数多根提案

**状态：结构实现成功；本次完整观测／成像增量门为 FAIL，不升级大规模实验。**

本线已通读 supplied `report_sources/report_body.md`、`memos/algorithm_branch_final.md`、`checks/exact_two_cell.py`，并读取冻结 Z-ALG R2 的 `general_quotient_solver.py`、`quotient_solver.py` 及其双线性系统接口。没有修改冻结源文件，没有安装依赖，没有使用 DeepSeek、递归代理或 GT 支持。全部预测封存后才交由父线程评分，本线没有读取 private 标签。

## 1. 物理核与真正运行的算法

在公共材料锚点 m0，L0=I-diag(m0)G。对真实坐标支持 A，令 U 为对应坐标列，V=solve(L0,U)，C=(GV)[A,:]，E_A=Q0[A,:]，O=SV。

```
F_A(d) = F0 + O solve(I-diag(d) C, diag(d) E_A)
m(d) = m0 + U d
```

这里完整保留全部原始源和接收器。只有 Q6/Q20 的提案方程选择一个已有接收通道与 k/2 个已有源；实际源列由带主元 QR 选择，接收器按 O 的行范数选择。4-cell 使用原 6×6 商特征问题，6-cell 使用原 20×20 商特征问题，均经 import 调用历史实现。图册先决条件失败、复 lift、不实材料、节点、box、近极点、重复材料分别保留分类，不能把失败图册说成物理无解。

合法根在完整有限核中评价，再以原始物理分数和材料距离聚类选择最多三个 seed。其中一个槽保留最低分数；有额外槽时保留与首个 seed 材料相距较远者。最终调用父线程共享的 `full_maxwell_corrector.corrector(method='trf')`，释放所有 16 或 36 个实材料，使用同一 box、smooth/TV 和停止准则。不存在合法根时仅把原锚点交给同一 corrector；此回退明确记录为 `NO_PROPOSAL_ANCHOR_FALLBACK`。

API：`src/algebraic_multiroot.py::run_multibranch(model,Y,m0,config)` 返回 top1 材料／预测、全部提案、根过滤档案、校正结果、实际费用、阶段费用与作用域。`src/finite_response.py` 提供精确有限响应、公共连续支持字典和完整数据 Jacobian。支持只来自真实 2×2、2×3／3×2 单元；没有稠密形状系数或 GT cell 列表。

## 2. 事前冻结的最小机制集

选择各 `(side,layout)` 组合中匿名 case_id 字典序最小的观测，选择规则不读形状或对比标签：

| grid | source layout | anonymous observation |
|---|---|---|
| 4×4 | jittered | obs_1bc8f5e504d740be |
| 4×4 | uniform | obs_1cc5bfaf8b33c1e0 |
| 6×6 | jittered | obs_5ec587d4bc0898c9 |
| 6×6 | uniform | obs_5a04b0ce1306e592 |

使用五个共享公共初始化。k=4、k=6 分开运行，每次三个支持：付费全数据梯度排名前两个，加上第一个尚未尝试的循环支持。每次最多一个实际接收通道。每个方法最多三个全材料校正 seed、每 seed 最多 40 次目标评价、每次运行 1 s 总 wall 上限。所有运行都在配额内完成，1 s 上限没有触发。

方法为 Q6/Q20、同核单起点 block-TRF、同核六起点 block-TRF、同核每支持 180 次普通随机提案、普通完整材料三起点，以及附加的严格 seed-count 对照。最后一项使用同核六起点 block-TRF，但全材料 corrector seed 数不得超过同 case/start/k 的代数路线实际校正数。四个观测 × 五个初始化 × 两个 k × 六个方法 = 240 次运行。

这些是相同预算**上限**，实际支出不同。不能把它说成每个方法消耗了完全相同 wall、完整 LU 或候选数。严格 seed-count 对照消除了“代数路线最多三个 seed、对照实际更多 seed”这一混淆，但求根与局部优化的提案成本本来就不同，必须按真实费用比较。40 个配置也不是 40 个独立材料：k 与初始化重复使用同四个观测。

## 3. 多根结构确实运行，完整数据门没有通过

| 项目 | k=4 | k=6 | 合计 |
|---|---:|---:|---:|
| 条件商图册 | 60 | 60 | 120 |
| 数值已解析有限端点 | 360 | 1200 | 1560 |
| 非物理共轭 lift | 296 | 1110 | 1406 |
| 实材料但出 box | 48 | 73 | 121 |
| 实合法子集根 | 16 | 17 | 33 |
| 完整数据 residual≤0.001 | 0 | 0 | 0 |
| 没有合法根的 case/start/k | 10/20 | 12/20 | 22/40 |

120/120 图册通过冻结数值门，所有期望端点均被数值解析；没有未解析图册、节点或近极点拒绝。**这不是严格代数认证或完整材料全根声明。** 条件完备性只涉及当前固定支持、固定 outside、选定测量系统；最终候选配额又进一步截断了提案集合。

33 个合法根的最大选定子集相对残差为 4.2710e-13；完整数据相对残差范围为 0.82212–1.25135，中位数 1.04527。也就是说，枚举和实材料恢复不是“没有根”，而是确实找到了完整观测不接受的子集根。

同一支持及相同测量图册中，实际出现 12 对不同合法材料根，完整材料 RMS 距离范围 0.87473–3.14267，中位数 1.32017。这是实际重建材料之差，不是复代数坐标差或仅响应差。两单元 supplied 有理数例仍作为已知端点、精确两根与通道障碍的演示保留；没有重跑精度扩展，也没有将其称为真实 Hankel 成像结果。

## 4. 共同 full-material corrector 后的结果与费用

| 方法 | data gate /40 | 中位完整相对残差 | 全阶 LU | 实际复 RHS | measured wall/s |
|---|---:|---:|---:|---:|---:|
| Q6/Q20 + corrector | 2 | 0.48712 | 2015 | 45138 | 2.193 |
| 同核 block-TRF | 4 | 0.42983 | 4374 | 96300 | 2.792 |
| 同核六起点 block-TRF | 4 | 0.38885 | 4481 | 98222 | 5.996 |
| 同核随机提案 | 0 | 0.40109 | 4506 | 98916 | 3.266 |
| 完整材料三起点 | 4 | 0.45777 | 4664 | 99312 | 2.407 |
| 同核六起点，匹配实际校正数 | 4 | 0.41293 | 1900 | 42192 | 4.788 |

**代数路线的两次 data gate 通过均来自同一个 6×6 uniform 观测、同一个初始化，在 k4/k6 都无合法根时的锚点回退；没有一个合法代数提案经过共同 corrector 后达到 data gate。** 因此不能把 2/40 写为多根算法的两次成功。

严格匹配实际校正数的 40 个成对配置中，代数路线的原始物理总分较好 15 次、相同 3 次、较差 22 次，比较容差 1e-8。原分数包含共同材料正则；不能仅凭较低 residual 解释材料准确性。它没有在当前门上给出可升级的独立收益。

代数路线支付 120 个商特征问题、600 个 full-system 支持 RHS、33 个合法提案的全源小核正演。自身每次没有对每个根重做全阶 LU；同核对照也享有相同权利。没有合法根也仍收取系数装配、QR、特征求解与过滤成本。全材料 corrector 全部收取原模型 LU 与 source／receiver-transpose RHS。详细小核、代数、失败、完整 RHS 等记录在 ledger；各阶段的 model cost 与 wall 分开保存。

费用边界：公共 Green 算子由父线程事前组装，原组装 wall 来自公共 operator 文件，summary 中按四套 geometry 单列；不假装每次免费重建 Green。每次公共 NPZ 读盘和 Model 构造 wall 未单独记录，缓存导入计数为一次/运行。最初五次可执行运行被直接纳入后续封存，未重复计算；这五次有总 wall 与完整费用，但 phase wall 未记录，其余 235 次有 phase split。不得把该缺口补成零费用。

## 5. 为何全数据校正收益不由小核定理保证

若完整数据的某个物理解只在 A 内不同，且条件商的简单根、非节点、非极点、rank、lift 等前提都满足，则精确枚举应包含它。这项包含性前提没有覆盖当前初始化：未知目标通常在 A 外也与 m0 不同。

对固定 outside，只用一个接收器与 s 个源相当于施加 2s=k 条实测量约束，所以可能有有限个合法实根。其余原测量施加了额外独立约束；同一 k 维受限材料集合一般无法同时消除 outside 错误引起的所有残差。精确解了子集方程并不创造其余数据相容性。本次 33 个子集根／零个全数据根是这个限制的实测表现。

全材料 corrector 随后解除 outside，故不是数学上“永远不能校正”；但局部 TRF 没有保证这些种子进入全数据好盆地，本次确实没有任何代数种子达到门。这也不能归咎于共轭根或箱约束检查“太严”：33 个实合法根已存在而全部观测残差仍大。

**决定：本线没有使用唯一允许的修复额度，也不新增物理数据、oracle 支持、更大 k、更多接收通道试探或更长校正来营救结论。** 一个未来可研究的、需重新预注册的改动是先支付共同完整锚点校正再构核，但本次未运行，不能写成已验证修复。当前结果足以判定本次固定初始化、小支持提案配置不支持升级。

## 6. 封存与复算

- `results/e2_mechanism_plan.json`：搜索前冻结规则。
- `results/e2_mechanism_results.json`：240 次完整观测指标、root／corrector 数量、实际费用；状态 `PREDICTIONS_CLOSED_LABELS_NOT_READ`。
- `results/e2_mechanism_materials.npz`、`_predictions.npz`：封存 top1。
- `results/e2_mechanism_candidate_materials.npz`：封存完整校正候选；离线 best@K 可另算，不能用于算法选择。
- `results/e2_mechanism_proposal_materials.npz`：封存所有原始提案材料。
- `results/e2_mechanism_ledger.jsonl.gz`：全部根、拒绝原因、校正来源、历史与明细费用。
- `results/e2_mechanism_summary.json`、`figures/e2_mechanism_summary.png`：无 GT 汇总与图。

现有运行环境：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`；CPU、单 BLAS 线程。执行 `experiments/e2_multibranch.py --tag mechanism` 可生成完整预测；`--summarize-prefix results/e2_mechanism` 只读取封存结果、不重跑物理。

父线程后续的材料 NRMSE／IoU 与结构类分析属于封存后的离线评分，应与这里的子集根、完整观测门、不同材料根距离和费用保持分离。它们不会把原始未通过全部观测的子集根追认为完整等响应物理解。
