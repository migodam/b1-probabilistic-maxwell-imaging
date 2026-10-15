# Research T：共享 bounded material GN 实现审计

本子任务仅实现、读取源码并进行已声明的小矩阵概率外代数核对；**新增 Maxwell state=0、新 full-solver RHS=0**。实际物理执行和科学解释由主 Codex 负责。

## 已读证据

- Research T 包：`CODEX_RESEARCH_T_EXECUTION.md`、`PRO_G_INTERFACE.md`全文。
- 原 A30：`src/a30/physics.py`、`src/a30/representations.py`、`src/a30/budget.py`及 compensation 源码；实际 frozen adapter 来自既有 A29 后端。
- Research A：`src/research_a/e2.py`及 `results/E2/E2_FAILURE_ANALYSIS.md`。原方法使用16/32实参数 DCT/Haar，每轨迹仅两次修正，8条件只有3个原始base相容。invalid→valid轨迹为0，不能用这些结果声称 full-material corrector 不可能成功。
- A17 `CLOSURE_R1/manuscript/CORE_THEORY.tex`、`CONDITIONAL_EXCHANGE_THEORY.tex`、`A17_CLOSURE_REPORT_ZH.md`及`publication/START_HERE.md`。本任务不调用其旧自动哈希包装器，不改冻结稿。

## 真实物理接口与费用

`FrozenMaxwellAdapter.full_state`在求解前预约1个state，每个新24×24 TM端点有8照明forward RHS。`receiver_transfer(state)`首次通过32个接收器adjoint RHS构建 `T=S L^-1`，同一state的后续调用复用缓存。

材料变化与所有照明共享，

`J(q,s,pixel)=T(q,pixel)*state.field(pixel,s)`。

`corrector.linearize`按 source-major、每源先Re32再Im32实化，材料坐标是全部576个Re(chi)后全部576个Im(chi)。这是raw material坐标，不是A17/A30质量归一化坐标；没有遗漏的 `sqrt(volume)`，因为data和state.field沿用已归一化后端且方向本身是raw chi。通过已支付的receiver transfer构造1152列并不免费获得物理：必须支付相应32adjoint。

每步新线性化的典型成本是32adjoint RHS，第一有限trial为8forward RHS，最多两次回退另加每次8forward RHS；state尝试数由共享book记录。已付费起点不重复求解，旧state的Jacobian不允许用于新state。构造线性Jacobian不调用 full material JVP 的1152×8个独立solve。

## 普通修正器

文件：`src/research_t/corrector.py`。

`CorrectorConfig`冻结full或multiscale，trust radius仅允许`.3/.6`，默认LM damping=`1e-3`，line search为`1,.5,.25`，可以另设每arm full endpoint尝试上限。父线程拟定4更新/6端点、`.6`trust、`.95*tau`停止留余量；物理相容阈值`tau`本身没有放宽。

Full模式使用盒约束regularized LS与LSMR，逐像素限制Re `[-.5,4]`、Im `[0,.5]`，并与protected core/guard的显式上下界取交。固定core行从优化变量中消去，不用一个小的非法Im分量把所有方向缩成0。Multiscale使用相同真实材料行约束；DCT或A30原Haar表示与真值无关。Haar64实秩即32空间列的Re/Im块。受保护的等式先取nullspace消去，避免重复像素等式使低秩SLSQP假性失败。

不添加TV、不使用图像真值、不拟合自由gain、不修改sigma，不提供对某arm隐藏的full-space fallback。LM只是增量正则，不把它误称长期材料prior。

返回完整trial数组引用、接受路径、line-search residual、active bounds、KKT/数值收敛状态与失败原因。预算/物理异常保留已有结果，不无限重试。Full模式报告projected-gradient KKT诊断；chart模式报告SLSQP状态但不将无约束chart gradient冒称约束KKT证明。

## Compensation screen

`bounded_compensation_score`在实际编辑后的paid finite state上，解同一trust/材料/几何盒中的局部凸LS。返回attained primal upper U和合法dual lower L：令 `v=r+Jd`, `lambda=v/||v||`，则

`L=max(0, lambda^T r + sum_i min((J^T lambda)_i*l_i,(J^T lambda)_i*h_i))`。

这是**该线性盒图表**的界，即使优化未收敛也有dual形式；没有 uniform nonlinear remainder时，`nonlinear_feasibility=UNKNOWN`。不能把U<tau当作full-Maxwell接受，也不能把L>tau当全局排除。Final archive仍由实际full endpoint及几何检查决定。

## A17边界

A17 paired closure对指定Galerkin增广给出材料步差所在的双侧空间 `H0^-1 range[J0^T Y, N^T]`。它需要相同残差、先验/damping和正确real lift，不是对任意有限结构编辑的自动可补偿性定理。其材料injection必须是共享材料引起的 `dchi*field`，不能用任意receiver-dark current替代。

在当前24×24全材料GN中，完整paid receiver transfer已经包含Maxwell resolvent反馈；给普通方法移除反馈再与A17比较，会构造错误对照。合法额外消融应是相同paid成本和相同real-rank的paired reduced材料方向、material Krylov与multiscale子空间，且须建立相同端点的物理比较。A17旧结论已明确paired早期路线被等秩Krylov胜过、完整支付后更慢；当前任务不能预设paired占优。是否值得付费该消融由主Codex根据qualified base与机制证据裁决。

## 零物理核对

`ZERO_PHYSICS_CHECKS.json`保存确定性5材料、3接收、2照明的小矩阵链式法则：packing相对误差`1.5935e-16`；core精确固定、盒违反0；primal/dual gap `2.2482e-15`；stale cache被拒绝。这不代表真实Maxwell backend数值审计或实验成功。
