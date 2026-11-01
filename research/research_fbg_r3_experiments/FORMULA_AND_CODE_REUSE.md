# 公式、代码与实际复用边界

本表区分直接调用旧代码、在同一模型上实现等价公式、理论启发以及没有运行的模块。公式正确、源码可导入或 GPU 完成执行，均不自动证明成像优势或创新优先权。路径相对于 B1 仓库根；本轮代码位于 `research/research_fbg_r3_experiments/`。

## 物理、目标与成熟求解器

| 项目 | 原始资产／实际入口 | 本轮使用 | 范围限制 |
|---|---|---|---|
| Scalar-TM disk-pulse VIE | `research/research_z_alg_r4_diversity/src/physics.py::{Model,build_operators,Regularizer,material_basis}` | `src/physics.py` 直接 import，薄封装 state/forward/Jacobian/tangent；4×4／6×6重新组装 | 没有实现新 Born 正演器，也没有截取576子矩阵。粗材料参数仍送入相应完整物理。不是3D vector-DDA适配。 |
| 真实状态与导数 | L=I-XG，j=L^-1XE，Q=E+Gj，F=Sj；T_h=L^-1diag(h)Q | 保留一份 LU，多照明 RHS；Full-J 用 receiver-transpose solve；目标梯度用 source adjoint | 对实材料、固定 G/E/S；source→receiver→real/imag pack。复数符号不等于材料也允许复数。 |
| 原数据与共同先验 | Φ=||F-Y||²/(2||Y||²)+R(m)；原R4 smooth／smoothed-TV | `src/full_maxwell_corrector.py` 的解析 TRF／伴随 L-BFGS-B；`src/physical_budget.py` 只加计费门 | Λ 为共同 IRLS 局部曲率 surrogate，不称 exact smooth-TV Hessian。box交接投影不擦除历史失败。 |
| R2 fixed feedback／GlobalDual | `research/research_fbg_r2_adaptive_feedback/src/adaptive_solver.py::{FeedbackLane,solve}`；`src/feedback_fit.py::feedback_value_gradient` | CUDA E3 直接动态加载 R2；CPU小网格在 wrapper 中实现同物理源Gram及完整导数 | 固定β=.2；GlobalDual β=0；非零残差度量不是增加观测。一般自适应 B(m) 仍需其导数，不能移植成无条件任意反馈。 |
| Resident CUDA | `research/research_z_alg_cuda_pipeline/src/device.py::{SharedOperators,ObjectiveLane}` 与 `src/hybrid_svd.py` | E3 wrapper／`src/cie_tsom_cuda.py` 复用同一上下文和现成PyTorch库，缓存公共算子与接收/Fourier基 | 两lane有限队列、实际等待/CPU读出/传输计费。未写新LU/FFT核、未装新服务、未声称streams固定分配SM。 |

## CIE／TSOM：按原算法原则适配，不冒充逐字三维复刻

主要可核来源为2019原论文（DOI `10.3390/jimaging5020027`），本轮正文获取与公式范围见 `audit/CIE_SOURCE_AUDIT.md`、`CIE_TSOM_BASELINE_AUDIT.md`。核心入口为 `src/cie_tsom_baseline.py::{fourier_apcs,_readout,run_cie}`，GPU current port 为 `src/cie_tsom_cuda.py::{prepare_cie_shared,run_cie_cuda}`。

复用的原则包括：接收 SVD 的 DPCS／APCS 分解，j=j_D+U_A a；投影 Fourier 子空间；保留原生 Fourier 系数的 PR-CG 与二次 current 线搜索；modified contrast rβ=βm/(1+βm)；共享材料交替读出；分阶段继承；两种β经同一物理 m 正确运输。带共同 TV 的材料块使用有界内层优化，而不是把逐点 LS 当整个正则材料块的最优解。

这是**二维 scalar-TM、共享实材料**的算法适配。CPU 显式 IDFT 列和缓存矩阵作用不等于原文三维 CG-FFT 的计算复杂度；没有复现原三维采集或图像。原论文允许复材料，当前实读出是本问题的受限版本。mf 指频域 mask 的边长，不能称为总 current 未知量。成熟 CIE 的完整多阶段日程与原论文式停止校准保留，不能以早停弱原型承担成熟强基线上限。

当前 `native_tsom` 开发分支额外直接调用该 `fourier_apcs`：固定 j_D、mf2 的16个 APCS 空间方向，保留15个为 core，再从同一公共池选择一个替换方向；p=6时总192个实 current 自由度。其首候选恢复原mf2空间；DPCS 数据已经固定，不计为可自由优化的系数。该新分支与旧 physical-anchor N1 是两个 affine 空间，不能混并效果，最终结论等待对应新批次封存／评分。

## N1：共享实材料几何与 profile 释放

| 公式／操作 | 实现 | 复用判断与边界 |
|---|---|---|
| 本构缺陷 e=j-X(E+Gj)，D(Pe)[h,V]=PLV-Pdiag(h)Q+(DP[h])e | `src/joint_current_model.py::JointCurrentModel` | 从原 CSI/CIE 与实变量导数直接实现；保留非物理电流上的 DP·e。不是新的物理后端。 |
| 实法向 Π_Q^(⊥,R)LV；CIE有限偏置 β||e_CSI||²/(||q||²+β Re<q,j>)，一般2阶／真实切向4阶 | `src/current_normal.py`；supplied report与已存E0 | 理论机制与诊断，不是按normal energy最大化的成像保证；复材料投影不能替代实材料投影。 |
| Profile H=Bᵀ(I-P_A)B+Λ，g=Bᵀ(I-P_A)b+ell；释放减去DᵀD和Dᵀv | `src/profile_schur.py::{profile_schur,release_candidate}` | 标准实最小二乘／Schur代数在当前joint模型的实现。H是局部剖面曲率，不是Fisher信息。 |
| QR坐标 u=Ra_C，γ||u||²/2；α=(1+γ)^-1 | 同上；`current_normal.thin_basis` | 运输实际候选系数，保持换基／加旧空间的不变性；投影后rank参照原候选尺度，避免roundoff伪增广。不是靠缩放非零列渐释自由度。 |
| Signed material task gain，同时保留残差拉回与曲率 | `src/current_mode_selection.py::evaluate_modes`、`paid_physical_jvp` | **A17公式思想启发**。本轮没有实际调用原生A17 Galerkin current exchange、Schur误差证书或原A17 selector。即使某选择器有收益，也不能自动归因于A17原算法的独立增量。 |
| 有限±current seed、固定物理基hold、软释放、完整材料readback | `src/controlled_release.py` | 本轮实验机制；尺度／defect／current预算和实际可行尺度显式记录。两个已知根之间存在弦线通道，不等于未知另一端点已被生成。 |

旧 A17 的条件增量如 ΔJ=Y_C Γ_C^-1 N_C 与双通道分析用于理解“怎样影响材料步”；N1自由辅助电流扩展改变的是joint可行方向，符号与物理Jacobian近似问题不同。`C'=C+UK` 保持相同 span；dressing负控制不能算新增方向或新信息。正常数normal trace也不能独立提供有效模式排名。

## N2：精确有限核与冻结商代数

固定真实坐标支持 U，V=L0^-1U，C=U*GV，E_A=U*Q0，O=SV。`src/finite_response.py::ExactFiniteResponse` 实现

```
F(d)=F0+O (I-diag(d) C)^-1 diag(d) E_A.
```

这是固定 outside 的精确非线性消元，和 R10／A17-X 已有有限端口／Woodbury 类公式等价；**本模块没有导入 R10 现成 finite-kernel 源码**，不能把公式等价说成调用了其旧压缩 backend。k只可能是4／6个真实空间单元；稠密材料 decoder 的6个系数不能冒充rank6端口。

`src/algebraic_multiroot.py::_frozen_solvers` 直接 import：

- `research/research_z_alg_r2/src/quotient_solver.py::enumerate_quotient`：k4，1接收×2原源，条件Q6。
- `research/research_z_alg_r2/src/general_quotient_solver.py::enumerate_general_quotient`：k6，1接收×3原源，条件Q20。
- Q6 仍通过原 `branch_enumerator.py::BilinearDualSystem` 构造双线性兼容方程，冻结文件未改。

本轮新增的是适配真实小支持、全源／全接收评分、完整材料共同corrector与公平同核控制的实验接线。Q6/Q20代数本身是既有资产，非本轮新发现；只在原rank、leading-block、简单lift／分离式等条件下讨论。浮点6／20端点解析不是严格全根认证，更不是576材料全根算法。非实、出box、节点、极点、空合法根、图册未解析、子集合法而全数据失败分开记录。

两单元精确有理例保留在 supplied `checks/exact_two_cell.py`：它证明该已给定代数模型的两根与无正材料极点，不是Q6/Q20真实采集成像benchmark。真实小网格上的33个合法子集根与全部观测验收也必须分开；详见 `N2_ALGEBRAIC_MULTIBRANCH.md`。

## 什么没有被复用成新的科学结论

本轮未移植A17的3D压缩Cholesky/Woodbury backend或A18的向量FFT来替换scalar-TM；未运行原生A17 Galerkin交换来为N1提供免费Full-J误差证书；R10跨形状旧基的适用范围没有放宽。历史R2 8/12只是旧正面结果，新576成熟基线重跑仍用已披露对象。新native_tsom结果需要自己的封存与评分，本文不预写其成败，也不宣称已证明全局收敛、普适成像优势或理论首创。
