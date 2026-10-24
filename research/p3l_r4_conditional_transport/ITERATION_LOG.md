# R4 实质迭代与最终状态

## 风险路线

假设：完整属性–相事件矩会胜过直接回归。实际运行K8/16点阈值→精确h质量bin→anchor-induced bin初始化。依次改变求积误差与初始化机制，保持P1、标签/容量等共同条件。精确bin仍弱于强CNN，anchor-bin出现局部平台梯度抑制且风险更差，停止扩大该家族。

发现：首轮CNN与R2函数类别相同，收益先计为训练基线修复。随后扩大到1024context并保持decoder/random控制。另将全局model均值/方差作为四通道context，加入同容量Gaussian、直接field与空间打乱控制。新512/256独立确认主要风险改善13.10%/11.36%；Gaussian几乎相同，拒绝Flow特有风险归因。

## 完整概率路线

三seed Full66 Gaussian/Spline/三角/CFM以及数据规模、容量、lr真实比较。显式Spline和静态物理frame NLL有独立增量。满秩GN坐标在当前训练下失败，不强制使用。bounded geometry CFM改变训练误差权重而非硬投影，随机同谱控制、相同base与数据；两批energy score改善约21%。

进一步用真实Full source-adjoint做paired physics/supervised mixed-KL拟合：response改善但NLL/覆盖恶化，三族均未过复制门；停止seed扩张。CFM数值密度以双精度全trace自适应ODE独立核验，预算895秒后只得部分稳定性失败，其余NOT_RUN。不是重复Gaussianreference后结束其余路线。

## 物理计算路线

active-row精确复用因变化集合大而净wall失败，保留原结果。首轮完整区间链同Full决定却更慢；只修改逆范数包络与动作流程，在新DEV和预冻结8新对象quiet重放。新locked同核wall71.397→24.153秒、逐步决定/路径零差异；普通GMRES同时作为强控制。较大步长失效边界保留，未宣称posterior mixing改善。

六次localvalue/JVP和三次动作排序器真实训练，误差或相对普通fieldupdate控制不胜，停止扩大，不冒充完整LazyDINO复现。

## 动态几何路线

主线程发现四段Euler二阶自身加速度污染；保留旧实验，改Heun并做真实Full配对。小步核验不足以保证任务尺度，Heun288项仅9距离匹配、2全部支持；直接法94支持，停止四段路线升级。满秩冻结pCN四核实际运行，ESS仍UNQUALIFIED。共同archive严格RMS.02形成111新有限噪声候选，无后验质量声明。

## 优化诊断与工程修复

原SLSQP callback选点局限并非物理无解。零新物理best-query与新Full-J TR在GT开放前封存；λ0降Phi却增risk，λ1与原baseline相近。所有预算停和原失败保留，不事后改旧资格。

Closed fusion K元数据、PPC release状态与离线list矩阵评分故障使用独立只读compat或新输出恢复；原weights、partial、ERROR、预测和成本不覆盖。最终报告构建也保留首次schema assertion失败收据，仅修复对failure字典的读取，未改变科学数据。

本轮不再新训模型或付费物理。完整费用以RESOURCE_LEDGER和CLOSED收据为准，旧中期文本保留于AUDIT，最终主张由GO_NO_GO.md约束。
