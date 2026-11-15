# 方法、表示与消融

本轮不重设β、不调rank阈值、不训练NN、不按GT选择候选。所有在线未知材料为576维实contrast，原物理G/S从未被低秩近似替换。

| 标记 | 材料反馈 | 周期动作 |
|---|---|---|
| g19 | 原G右19的PGP | 只审计 |
| g19_s19 | S右19的PGP | 只审计 |
| twofold | 条件共同强U11的PGP | 只审计 |
| twofold_hop | 上项加一次N−PNP余项 | 只审计 |
| twofold_dr | 同一材料目标 | U01路由、保留记忆的反射 |
| twofold_ap | 同一材料目标 | U01路由AP |
| twofold_full_dr | 同一材料目标 | 完整空间周期反射 |
| twofold_reset | 同一材料目标 | 只重置PR |
| cie/dual_cie | 原成熟实现 | 原规则 |
| dr | 原完整场—current DR | 原规则，原corrector |

前八臂共用原global current PR、同stage/mask、完整移动材料导数的LBFGS8及同初始材料的物理current warm start。与原生CIE比较是完整工作流对照，不能把sourcewise/global等原有差别也归因于本轮分母。前八臂之间才严格隔离本轮改动。

## 周期动作的实际顺序

1. 完成原current及材料步。
2. 每16步及stage末，以g=(J/6,Q)加持久记忆构成反射状态。
3. 计算完整共享实材料投影和完整传播—观测球投影，按H=P01或I路由增量；AP对照替换反射动作。
4. 计算新的材料影子，付费full forward审计；不以立即Phi或GT拒绝程序更新。
5. shadow current投回原DPCS/Fourier仿射族，保留材料、重算Q、重置PR。
6. 在线primary始终只保留已完成合法候选的原Phi最小者；程序探索状态与被选输出分离。

phase、raw/corrected/primary/secondary、全gap/路由gap、注入误差、PR及state checkpoint均保留。仅zero routed gap不构成完整可行证书。新的全空间周期DR包含G19移动anchor和表示注入，与原独立DR不是一个算法。

## 为什么保留这些控制

S19检查收益是否只是观测方向替代内部方向；twofold检查删去S弱方向的效果；hop检查一次邻域是否增加恢复；reset检查PR重启；AP检查投影接口是否足以解释反射；full周期DR检查子空间路由是否减少全空间干预的损伤。各stateful臂周期审计机会匹配，额外计算全部收费。

原G19数值接线在一个公开对象前三步与冻结原程序逐项一致；该事实不表示额外计时/缓存开销不存在。完整原生CIE/Dual/DR的source与绝对预算改动记录在runtime和dispatch provenance，未把新方法给予更多未计费物理求解。

## 因果与结果边界

56fork共用完整状态，能检验有限程序干预；它们从中间状态开始，不是从最终错误局部极小重新逃逸。16读出反事实只固定材料块，不把单次局部方向当完整算法。44终点验证只提供投影梯度与活动面，不证明严格极小/局部Hessian正定。

原共同corrector是full-Maxwell L-BFGS-B，允许把unused proposal预算重新分配。原独立DR更早完成projection日程，确实获得较多corrector时间。这是实际效率的一部分，也是归因必须保留的因素。
