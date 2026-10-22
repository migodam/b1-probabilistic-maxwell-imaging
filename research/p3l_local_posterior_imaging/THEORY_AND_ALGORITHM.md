# 概率合同与已实现算法

## 物理与参考测度

保留原 P1 的 `q=(z64,u2)`，`z~N(0,I)`、`u~Uniform([0,1]^2)`。Gaussian field 经原 C² saturation decoder 形成 24×24 被动材料，256-cell ROI、8 sources 和 32 complex receivers；全部照明共用一个材料场。这里没有将 host phase 当成空气，也没有读取 GT 参与在线初始化。

Re/Im 噪声方差均为 `sigma²/2`。源主序 packing 为每个 source 的 Re32 后接 Im32；白化因子 `sqrt(2)/sigma`。目标为 `Phi=||F-y||²/sigma²`，`pi(q|y)∝exp(-Phi-||z||²/2)`。noise-ball 只用于材料资格诊断，不是 posterior 的硬支持。

在 `u=Phi_normal(eta)` 下，`p0(q)dq=N66(x)dx`，`x=(z,eta)`。这是精确测度变换。native uniform 属性在 q 中没有 Gaussian Hessian；在 x 中出现的 Gaussian 项来自 Jacobian，二者不能混用。

## 强弱分解及条件概率体积

采用 Gaussian64 的正交完整 frame `[Vw,Vs]`，`z=z0+Vw a+Vs b`，绝对体积 Jacobian 为 1。属性显式保留。generalized vectors 只用于选择 span，之后转成正确 prior 坐标。

完整弱 marginal 包含 `Zs(a,u)=∫p0(q) exp(-Phi(q)) db`。固定 b 是 SLICE，优化 b 是 PROFILE，有限区域条件化是 LOCAL-C；它们都不自动等于 FULL marginal。

在参考点 b0 附近，GN 二次参考为 `U0+g^T delta+.5 delta^T H delta`。局部积分有线性项修正 `.5 g^T H^-1 g`；只有 g 足够小，profile-logdet 公式才是相应中心近似。即便这个修正很小，非二次误差仍可能使 importance weights 崩溃，实际实验已看到这种情形。

## 材料矩与风险

Gaussian 强条件参考下，field 的每个像素边际是 Gaussian。饱和多项式 decoder 的一、二阶矩可由稳定 truncated-normal moments 精确计算，保留 latent 相关性对边际方差的贡献。两个属性固定或需在外层按正确测度积分；不能把相关属性与 field 的期望简单相乘。

这个解析模块只消除 Gaussian 参考内部的材料 Monte Carlo 费用。参考分布错误、forward 非二次误差、box 截断、外层积分和有限样本误差仍分别存在。真正 posterior mean 最小化 conditional MSE；它不一定属于单个 decoder 合法材料集合。complex MAE 的 Bayes action 是逐像素二维几何中位数，不能随意拆成 Re/Im median。

## 已实现的合法采样核

包括 native Gaussian field pCN、独立 property reflection block、未强制 bracket cap 的 Elliptical Slice、完整 probit MALA、变化 decoder mobility MALA、固定质量矩阵 HMC 及 adaptive tempered SMC。拒绝时保留重复状态。生产参数冻结，变化度量 Gaussian 提案保留两端 determinant 和 reverse drift。

变化 decoder mobility 使用已付费 `C_chi(anchor)` 与新 `Dg(q)`，但 proposal 的近似几何不替换真实 Maxwell likelihood 或 source-adjoint gradient。相同 full likelihood 的 MH 不要求 proposal geometry 本身准确；它仍需要真实混合诊断。

HMC 采用 inverse mass C：`K(p)=p^T C p/2`、`p~N(0,C^-1)`，完整 leapfrog、momentum flip 和 Hamiltonian MH。大但有限的能量误差仅作诊断，不任意提前停止轨迹来改变 kernel。

SMC 使用 `pi_beta∝p0 exp(-beta Phi)`、条件 ESS 温度、系统重采样和合法 MH mutation；保存全部祖先与权重。有限 N、自适应日程和模式遗漏风险明确保留，高 weight ESS 不是全局正确性证明。

## 归一化 learned proposal

Joint-supervised 模型学习固定 Gaussian weak8 与两个 probit 属性的归一化条件密度。训练尺度仅由 train 拟合，log probability 包含标签尺度 Jacobian。它不是 LOCAL-C teacher，也不是完整 66 维密度。

完整 proposal 可由该弱密度、显式相关 strong Gaussian 参考和 10% defensive Gaussian prior mixture 组成。正交 chart determinant 为 1，两个提案方向都可计算正规化密度；对原 FULL likelihood 使用真实 MH/IS correction。确定性 corrector 输出不会被当成 posterior 样本，短链或小 ESS 也不会当成可信参考。

## 几何与计算复用

动态弱子空间不必 Frobenius 可积；本轮测量的是具体投影场的法向 Lie bracket，并同时检查谱隙和有限差分收敛。单条轨迹没有横向概率质量。

primal/dual residual 和 stability lower bound 可以给出 output/likelihood 区间，但本实现的 floating-point inflation 仍有明确 IEEE/BLAS 假设，不是 directed interval certification。不可用的充分 bound 不意味着 full operator 奇异。固定 proposal、固定 uniform 的决策对照必须在看完整 reference prediction 之前作出。

A17 finite-state 表示与普通 Galerkin 在相同 current space 下没有自动的精度优势。本轮保留 optimized Full、exact T-matrix、普通 goal-ROM、Krylov 和相同空间 A17 对照。是否有成本增量由实际 setup、refresh、fallback 与所有 RHS 决定。

以上通用机制都有成熟先例。文献核验、429 及全文缺口在 `AUDIT/assets/literature/` 保留；没有检索到不等于原创。尚无完整 posterior/reference 与独立风险证据时，不提出 IEEE TAP 新颖性结论。

## 最终实验对理论的约束

128独立对象没有支持当前Gaussian条件积分成像收益。保留decoder的代理力显著减少响应误差，但未减少同RHS下material MCSE，且wall更高。数学合法的概率变换、解析矩与MH核不能替代混合与材料风险。详见最终报告、THEORY_RISK_INTERPRETATION与两臂实际对照。
