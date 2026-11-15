# R15 基线与比较范围

G19是用户从CIE发展出的内部Green反馈方法，不是外部算法。本轮新G19保留原global current PR，只改变材料残差变换的固定子空间和周期干预。`g19`、`g19_s19`分别为固定内部G右19与接收S右19的PGP/Galerkin分母；`twofold`按公开算子选U11，rank随孔径变为11或6。所有臂仍用完整576材料未知量和原Maxwell算子。

## 真实实现来源

- runtime/dualmetric_core.py是R14完整PR/current状态程序的逐字副本，引用实际R2/R4/R6/R7源；新接线通过独立函数全局绑定替换材料readout，没有重新写current PR。
- 原matrix_material_readout、固定低秩Woodbury分母及LBFGS8来自R10；原互不重叠4cell邻域动作来自R4。新代码补足一次残余邻域的完整材料导数。
- 原生CIE、Dual-CIE、完整场—current DR使用runtime/frozen_C2_runner.py的真实成熟路径；equal_budget_C2_runner.py只采用之前声明的绝对120秒corrector截止，dispatch在同一进程运行，不重启CPU计时。
- Baseline不被新方法的period16、U11或投影重启所替换。保留各自sourcewise/global、current表示及合法初值。因而新full-space周期DR不能被叫成原生DR，原生DR单列。

## 必要接线检验，不重验历史恒等式

audit/integration_v1包含7次实际验收：新G19对冻结原G19同config公开对象前三步的m/J/Q/A、PR旧梯度与方向、PR+和alpha逐项一致；五个周期臂各17步均实际触发事件，状态有限、材料合法；checkpoint恢复后下一步一致；周期前待处理DR事件也能精确重放。只证明已检查的短程接线，不保证所有硬件上同CPU终点逐位相同。

audit/denominator_wiring保留本轮新增移动材料导数、无邻域退化与伴随接线及全部尝试费用。没有重跑Atlas恒等式，也没有用GT做接线筛选。

## 公平性与限制

新八臂拥有相同周期审计机会、同current程序、同材料prior/box和公共初值。两个谱分母、邻域余项、PR-only reset、路由AP/DR及全空间周期DR逐项分开。所有新初始化投影、SVD、背景Cholesky、三角解、小Woodbury分解、CPU和序列化成本收费。

相同120秒上限不是相同迭代数。保存同CPU已完成候选前缀以及共同componentwise物理LU/RHS前缀；这些是付费候选池重放，不是免费早停重算。背景投影不计为材料Maxwell LU，也不免费。原生基线的额外审计和重启规则不被削弱。

本轮使用24个R14已公开对象，包含4个36到24跨网格采集、其余matched，噪声1%/2%、6/12 sources、全/半接收孔径。所有方法在同对象的同一核心顺序执行；不同对象跨P/E核心，不能把跨轮绝对CPU差当纯算法差。四对象在checkpoint512开始的五臂64-step同状态fork，以及36个16-step正负扰动另付验证成本，不混进120秒主对照；原checkpoint32设想已在运行前按V2合同改正。

R13的48对象、R8的NN对象、R11坏终点与本轮`r14_new_XXXX`不是同一分母，不能按尾号合并。本轮没有重训NN，没有加入新的3D/Fresnel数据，不能把历史三维资格和实测成果写成本轮twofold反射的成像证据。
