# P3L-R2：实际迭代记录

本记录区分开发集研究、算法冻结及新测试。每个运行的实际费用以独立 `RESOURCE_LEDGER.jsonl` 和各作业 `CLOSED.json` 为准。没有借用上一轮余额，没有修改 P1 的 prior、decoder、noise 或 acquisition。

## 第一轮：把三个科学问题真正分开

**假设。** 上轮失败的 Gaussian teacher 不足以否定直接材料风险学习；Full posterior 推断、材料点估计和几何探索应该独立执行。

**实施。** 生成独立的 16,384 个 train 和 2,048 个 validation parent；训练 measurement-only mean/ridge/MLP/CNN、归一化 10 维坐标密度。另为 512 train 和 128 validation 对象建立不读取 GT 的 native48+probit48 物理上下文。重复使用历史 DEV 的合法物理缓存，只把明确的 TRAIN GT 用于监督奖励。所有新增 truth forwards、物理反演、拒绝及训练均记账。

**结果与解释。** 新网络取得了 update>0 的 validation 改善，因此不能沿用上一轮“全部 update0”的结论。较大的普通残差网络收益很小；5,266/5,842 参数的共享局部 CNN 明显更好。这是监督材料风险证据，尚不是 posterior 证据。32,768-parent 扩展保持同一 prior，继续改善 measurement-only MLP，但未据此删除强物理基线。

## 第一轮的独立概率反证

**实施。** 原生属性方框上的二维 Full-Maxwell SLICE，用完整方框分区积分对照四条 MH 链；同时运行真正的 66 维多中心三角 Student proposal、逐源同化 SMC、真实 Full-J metric MALA 与 frozen-metric source-adjoint MALA。归一化 proposal 和 Jacobian 正确并不自动使这些方法混合。

**结果。** 三角 proposal 的 4,096 次新 Full 提案全部拒绝。保存数组的诊断表明，linear latent 预测与真实有限响应有巨大差距；decoder 饱和转换与 Maxwell 余项均大，且存在显著负交叉项，不能只按各自范数归因。逐源同化到达八源 FULL target，但仍仅一个先验祖先。真实动态度量 MALA 没有在同 RHS 下胜过 frozen metric；两者均未取得 FULL reference 资格。

**决定。** 不扩大同一个失败 MCMC；监督风险和几何路线继续。属性 SLICE 的数值积分改善不被称为 66 维 marginal。

## 第二轮：检验局部 CNN 收益来自什么

**单一主要改动。** 加入输出空间、参数量与几何位置的控制，而非改变统计目标。

**固定对照。** 原 local-anchor CNN 与 local-geometry CNN；新增同样 5,842 参数但几何通道置零的 capacity control、固定空间置乱的 geometry control、相同输出维数/界限的随机材料 frame、只改属性、全坐标半径匹配的合法 decoder。三 seed 与 100 epochs 保持一致。另用 validation 选择普通 Gaussian blur/先验均值收缩；它选择恒等动作。

**执行失败。** 一次误传 200 epochs 的对照在参数检查时退出，未训练、未新增物理；保留 `A/jobs/residual_local_capacity_controls_r2/CLOSED.json` 及日志。随后新作业按共同 100 epochs 补齐缺失对照，没有重置费用或覆盖失败。

**决定。** 冻结全部已完成 checkpoint，包括负面方法；比较用独立 parent，不把三个训练 seed 当三倍测试对象。

## 第二轮：几何路径与风险动作

**实施。** 在相同 `.02/.04` 实际材料 RMS 任务、合法 decoder 与相同 Full evaluator 下，比较 frozen、dynamic、curvature threshold、普通随机方向、contextual bandit 与小型 model-based policy。训练和失败反馈来自明确的 TRAIN 轨迹；评价时不访问 GT。更新后的真实材料、有限响应、prior cost、signed progress、求解费用全部保存。

**解释。** 访问轨迹不是 posterior sampling。Arc-length 梯形、等节点和监督凸组合都是材料风险动作，其系数不是 posterior probabilities。数据相容、likelihood 竞争、任务尺度移动和最终材料风险分别报告。少量开发例的改善不替代 fresh test。

## 非 Gaussian 材料场控制

属性切片不能检验 saturation Jensen effect。因此在已暴露 DEV030 的 paid anchor，用局部风险算子选择两个 Gaussian field 方向，转换成正交绝对 Gaussian 坐标，冻结其余 62 维和两属性，运行真正非线性 field2 SLICE。GH 逐阶细化、独立 pCN、材料均值与 decoder-of-latent-mean 的差均保存；数值不收敛时保留 UNQUALIFIED，不能用节点数代替精度。

后续独立 prior-CDF 分区积分若执行，只精化同一条件 SLICE，不改变原 FULL P1，也不作为监督成像/轨迹测试的前置门槛。

## 全新测试冻结

`METHODS_FROZEN_FOR_LOCKED.json` 在第一份新观测产生之前保存：77 个单独 checkpoint、主模型/融合、原统计目标、256 个 parent IDs、三个训练 seed、普通额外16-state优化控制、轨迹覆盖、失败 fallback、2000 次 parent paired bootstrap 及两项共同主要比较的 97.5% 区间。

主物理基线使用 native96+probit64，比训练上下文的48+48更充分；这个 covariate shift 明确披露。全新测试的在线估计与轨迹全部完成后，独立 evaluator 才能打开 GT。后续报告以实际 completion receipts 为准；冻结不是实验成功证明。
