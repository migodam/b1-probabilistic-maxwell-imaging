# R5 Stage I：先观察，再决定扩展

2026-11-03。本轮独立工作树来源为 R4 发布 c55a4b46，普通 git 来源记录，不作 hash 核验。所有 R4 对象现均已公开，属于 DEVELOPMENT；新轨迹分析不得改写其确认统计。

## 固定样本与数据边界

旧开发弯条 obs_b20d5a55677d4eca；两个 R4 确认失败弯条 obs_42cb233e4f224bb3、obs_63abafa750fc483c；按公开 ID 取正常环 obs_100bd65498644322。前两类是用户指定的机制正反例，不能据此推断总体成功概率。每例 global/onehop 各一次辅助轨迹重放，不重复 corrector 或其余确认对象。旧 phase_pi/self 的三个阶段状态先直接复用。字段缺失保留，不编造旧 PR 向量。

重放调用原源码既有 callback，每阶段固定检查点 1、2、10、50、150、300、500；保留 m_before/m_after、更新后的辅助 J、系数和可用的上一步状态。额外 Q=Einc+GJ 只在检查点计算，单独收费。保持原初始化、源 metric、prior、box、三个500步阶段和停止规则；不因结果调整检查点。

## 先冻结的诊断

源 metric 记 Ω_s，数据噪声 whitening 记 W_y，材料切向记 C_mat，接收算子记 S，避免同名。

1. 同一辅助 J/Q 和材料 m，比较 ordinary CIE、self、正相及反相传播。计算本构缺陷、u、v=Lu、逐cell有符号 cross=2Re sum(conj(u) v Ω_s)，并保留其空间图。
2. 以原始 r=J−diag(m)Q 和 C_mat dm=diag(dm)Q 构造 frozen metric 的局部 LS 读出。C_mat 使用实材料、复场实化、全部源共同参数；不把每源解当合法材料。比较 transport-induced 读出差与实际 dm，而非仅看范数。
3. 在非零 r 处，将精确材料梯度分成 frozen-metric 部分与 (dP)r 部分；这是对根处公式适用性的必要检查。另把精确正反相梯度分成 self、交叉、二次项。prior 梯度另列。
4. 每个检查点原物理 m_before/m_after 的 Phi、相对残差、投影梯度、实际材料步、负梯度夹角、线性预测和真实 Phi 差。没有把物理下降等同 GT 改善。
5. 必要局部谱只在阶段端点求取：完整当前 Maxwell Jacobian 的 Gram/特征值，沿实际材料步的精确曲率与 GN 曲率；既有 LU 与多源 RHS 复用。其成本单列，不重算完整 Hessian 全空间。
6. 状态固定的完整材料 readout 反事实仅先在阶段末点做 ordinary/self/plus/minus，使用相同 J、Q、m、prior 和原有内层额度。它是开发诊断，不是在线 GT 选择的新算法。

上述指标设计、对象和检查点不读真实材料数组。待诊断产物封存后才作 OFFLINE_DIAGNOSTIC 的逐检查点材料误差分析，不能据此改变本批预测或规则。

## 下一步如何决定

没有新的数值 gate。先问冻结 metric 是否解释实际材料步，材料导数项是否不可忽略，正反相在相同状态是否出现有方向性的差别，以及优势出现于哪个阶段。若混淆仍在，用公式驱动的 current-only/material-only 或阶段交换消融区分路径因果；每次新增前记录它能排除什么。九数字、三维小模型和实测数据只在有明确待检验机制/反例后启动，不以空表冒充结果。

所有新执行保存状态、时间、LU/RHS与额外 Green/小系统计数。额外成本与 R4 历史成本分开，未测金额为 null；不扫描 β、分块权重或测试 GT。
