# 实际恢复轨迹：理想谱之外发生了什么

本文件对应真实192条曝光轨迹、519199条完成native记录、170917个已付费候选，以及7个完整checkpoint的56条干预枝。所有候选的材料误差在封存后附注；prefix始终由公开minΦ决定。完整候选表保存在Release侧，Git中是紧凑任务/阶段汇总。

## 1. 正例需要空间传播，也需要合适的current路径

r13_fresh_0017与0045：完整01好，self-only01、global-Dual及从头11差；先01再晚切11保住，反向次序不能救回。0017换成所有四臂相同的χ-LBFGS60材料程序后仍只有01好。这排除了“仅Arc更准确”或“只要两块都使用更好norm”的简单解释。

但曝光0024在global-Dual上同样恢复；新r14_new_0021又只有11恢复，01及晚切11都差。早期不一致并非普适解答。四个标签代表具体metric函数和实际程序，不能从一个标签推导所有对象的可达性。

## 2. 同完整状态干预提供因果证据，短期Φ不是材料预测器

0002@128/512、0017@128/512/2500、0037@128/512的m/J/Q/A、PR梯度/方向、DPCS、stage及权重复制相同。每state运行四pair，分别保留与清空PR，观察1/8/32步。已有224快照含原Φ、D/R、构成缺陷、辅助port误差、box活动与168个真实材料步—负物理梯度角。

11能在若干状态更快降低短期Φ，但成功状态和失败状态都出现过这个现象。重置PR也会改变同一步材料方向，且不总改善Φ。这使“比较不同整轨迹的梯度范数”不足以独立归因；应先看同state干预，再看后续材料恢复。

## 3. signed block work的正确读法

每个stage分别计算current势的单步变化和material更新对该势的work。00/11在固定stage共用势；01/10一般不是共同下降。跨stage的β变化不能合并成一个全局Lyapunov函数。

本轮记录中current非下降方向重启全0，但PR+仍可截断且stage仍重置。相关性以48材料为单位、按pair分别计算，不把五十多万步当独立样本。正work或短时放大没有被证明是有用探索，更不是chaos证据。所有关联是曝光后的描述，不参与新24的算法选择。

公开minΦ incumbent的165467次变更中，有60580次NRMSE增加。这是依赖轨迹的现象计数，不是总体概率；它直接提示候选选择、corrector和真实材料质量必须分开。

## 4. 理想根谱、有限current误差和readout guard

8²理想01ρ=1.002712449。材料单独扰动的实际32步norm gain约.692，清空PR约.634；这不能反驳完整state的弱扩张。完整兼容current/material方向测试又发现01的128次读出全部nit=0，原gradient guard低于1e−8。这两项实验分别排除了“材料短期收缩就是全状态稳定”和“理想λ的32次幂就是实际程序响应”的误读。

实际方向加法缺陷在声明尺度上存在；current正规方程误差传入材料的项远大于材料内解误差。native非正交/冗余Fourier坐标、有限PR和容差都是实际算法的一部分，不能用精确BR替代后就把差异称作Maxwell材料动力学。

## 5. 非根metric导数力量级

在同56个first step，用原m、current更新后的J/Q和实际Ω/β，分解完整材料梯度为构成力、metric导数力和prior。直接调用原full derivative，不跑新forward或有限差分。metric/const范数比.1003–2.0673，各组中位约1；const·Δm>0有38/56，metric·Δm>0有6/56，full·Δm>0为0/56。

根处被O(||e||²)消去的项，在真实非根轨迹上并不小。P0也有该项，且其逐cell符号可解析写出。比值大可能有构成力抵消的原因；它不等于“更多材料信息”，也尚未成为恢复预测器。数组和方向work为下一轮Pro保留。

## 6. 五失败与corrector的边界

五固定严重失败的20个primary都有更低Φ的合法GT参考。28端点物理探针中，这五对象19/20在预定微小步下有下降，不能叫已证严格局部极小。小下降远不足以说明能跨过错误材料区域；预算/病态/表示路径仍需区分。

主48总体corrector增加联合恢复数，局部却能损伤材料。三维frame CIE是大幅结构损伤的明确例子；共享inner控制14/16的NRMSE增加则还需看幅度，不把所有微小增加都当灾难。raw/corrected/primary/secondary全部保留，任何保护策略都未用GT决定。

## 图与数据

- [恢复轨迹与signed work、same-state角度全图](analysis/recovery_dynamics_v1/report.html)。
- [192个轨迹汇总](analysis/recovery_dynamics_v1/TRAJECTORY_PRIMARY_ASSOCIATIONS.csv)。
- [56枝材料力分解](audit/material_force_decomposition_v1/report.html)。
- [根处有限响应](audit/local_closed_loop_v1/ANALYSIS_NOTES.md)、[完整兼容扰动](audit/coupled_mode_trajectory_v1/RESULTS_BRIEF.md)。
- [成对成本下的独立结果](analysis/pairwise_cost_independent24_v1/report.html)。

本轮最有解释力的预测反例是：**晚期一致化能保住0017，但只有从早期使用11的已测试配置救回新0021；根处ρ与短期Φ都不足以事先区分它们。**这不证明所有其他调度都不可能。下一轮应从完整state和非根力的方向关系提出可证伪规则，而不是从最终图倒推一个必胜标签。
