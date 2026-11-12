# 基线是否公平、充分运行，以及名称是否准确

源码逐项核对见 [BASELINE_SOURCE_AUDIT](delegated/full_current/BASELINE_SOURCE_AUDIT.md)。本文件给出父线程的比较范围判断；材料表现另见封存后的结果。

本轮保留成熟 sourcewise CIE、Dual-CIE，而没有用普通 GN 或单 beta 的弱配置代替它们。全部标量新对象保留576实材料未知量、共同 box/prior/观测、原物理 Phi 和付费候选规则。方法允许使用各自原生内部表示，转换缺陷不被隐藏。

## 同实际预算与充分日程分别回答

主 i9 确认集是120秒整进程 CPU 比较。在48对象中，CIE 有46份完成原生3000步、2份被 proposal CPU 截断；Dual-CIE 48份均完成；G19 只有1份完成、47份被截断；full APCS/full current 各2份完成、46份被截断。DR、整数事件、同范数随机各48份完成固定步数。程序写出 COMPLETE 不等于算法收敛。

因此，这组可以回答实际固定 CPU 下的性能与代价，却不能声称“G19在充分运行后没有搜索能力”。为避免这个混淆，在新 GT 尚未评分时另行冻结同一个16对象的 M4 参考：CIE、Dual-CIE、G19 均16/16完成原生日程；DR16/16完成；direct和physical GlobalDual均16/16由优化器返回。这个参考另加入16份原 Full-J TRF。各方法同机比较，其结果不与 i9 的时间直接拼接。优化器返回也不自动证明全局最小或严格局部最小。

## 容量与初始化的边界

24²上的 mf12 已是完整 APCS。G19 full APCS 从第一步解除早期 Fourier 截断；full current 还解除 DPCS 硬锁。它们保持 G19 本身的材料更新、PR、beta 和先验，但每步费用与 continuation 不同。不能把“空间更大”当作与计算成本无关的算法优越性。

G19 使用继承的物理 current LS warm start，CIE/A+ 使用继承的 cold coefficient start。相同 m0 并非同辅助状态。完整 current 减少初始化表示损失是接线事实，是否恢复更好另由实际材料指标决定。

`global_dual` 是 **global-PR Dual-CIE**。`physical_globaldual` 才是 beta=0 的 total-field source-Gram 原物理目标，保留 Gram 对材料的导数。这两个名字不互换。

## 不冒充的对照

`nfft_equations` 是 Eq3–6 的 joint L-BFGS-B/实材料/box/prior 适配。尚未完成原作者 NFFT-SOM 的全部 CG、材料步和停止细节复现，不将它作为“已击败成熟 NFFT-SOM”的证据。

三维使用 R8 的真正三分量 dyadic DDA/CIE；实测使用每个 source 自己的 receiver map、DPCS 与伴随。未实现兼容的 source-specific G19/A+ 或一般 vector G19，表中不填一个不合格移植来制造失败。3D/实测的 prior、归一化与标量新确认集不同，分别报告。

Austria 是完整物理合同下的公开文献对象适配，非只模仿轮廓；但逆网格32²而非原 NIE 的64²，另有共同先验/corrector和有限实材料盒。不存在原作者代码逐位复现的声明。

## 审计与成本

统一的 CPU 审计时钟仍允许不同候选数；原生阶段审计和已付费物理目标试探均保留，不用删候选制造表面公平。固定传播投影的分解、三角解、Green/端口动作与所有越限操作均计费。Physical LU 横轴单独展示，不能忽略 DR 的大量固定投影 RHS 后称等 LU 就等计算成本。

所有 raw/corrected、球内资格、Phi和材料误差分开保留。NN 的历史训练先验和训练费用未知单列。坏终点合法重启不冒称恢复了遗失的 PR 历史；实际状态保存范围见 [AUXILIARY_STATE_COVERAGE](delegated/full_current/AUXILIARY_STATE_COVERAGE.md)。
