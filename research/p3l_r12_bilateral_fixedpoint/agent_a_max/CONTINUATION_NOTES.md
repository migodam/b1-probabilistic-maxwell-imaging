# R12 Sol-A 接续：公式与已有证据的解释

2026-10-30 UTC。接续范围仅为 R12 推导说明、living formula atlas 与已有文献身份记录的解释。没有调用 Maxwell、LU、正演、decoder 实验或新文献 API；没有打开 N 真值文件，也没有修改 physics、engine 或 agent_c 文件。最后阶段遵守用户“不做防御性编程、不重复比对检查校验”的指示，没有新增测试、数值恒等式重放、评分或独立审计。

## 交付

公式集追加 [第12章](../../maxwell_formula_atlas/12_BILATERAL_MATERIAL_MODES_AND_FIXED_POINTS.md)，九个稳定新编号 `F-BF-001`—`F-BF-009`，保留已有158条，共167条。章节给出尺寸、逐步推导、假设、资格与计算后果，并追加 registry、源码索引、总集和 changelog。文档渲染使用 R9 的文档 Reader 与中性外壳；不运行历史 builder 的 main，也不把 R11 的预算失败提示注入当前 atlas。

## 可供主线程直接采用的解释

- 当前完整双边 pair/Gram 与 Full-J 是同一局部数据因子的不同表示。top32 pair、单边或宿主因子会改变提案几何并丢掉信息；是否更有用须看原 decoder 下真实有限材料结果，不能从分解本身认定新颖性。[原始推导](../BILATERAL_MODE_DERIVATION.md)
- 原材料 RMS 分母保持576，Gaussianized66 的两个属性列保持 CDF 导数。先用材料实化 SVD 删除局部不可分辨方向，再分解显式增广因子 `[J_w T; T]`；下块是原 Gaussian 先验曲率，没有新增观测。[实现收据](../agent_a/IMPLEMENTATION_RECEIPT.json)
- 非根处 `DT_mu` 含残差 Hessian 项和逆矩阵微分项。精确根才可用 `mu/(sigma_i^2+mu)`；正残差驻点仍有残差曲率。近单位模态、decoder 零方向和人为固定坐标均不能直接作为连续材料根族维数。[既有谱推导](../FIXED_POINT_SPECTRUM_AUDIT.md)
- `shared/run_development.py` 为每臂设5400 RHS与240 LU的共同上限，并计入 COMMON；这是一致的上限，不表示实际各臂消耗相同。两锚点、三层、每层两次提案、每次五个 paid-J 状态、原先验感知几何和固定 latent radius 均为现有运行配置，本次没有改动。[运行源码](../shared/run_development.py)
- P1 是4个已暴露 DEV parent ×8臂，冻结 epsilon0.075 的经验覆盖全部0/4。H0 RML 平均离线最近材料 MSE0.0119253，低于 H1 full16的0.0130177、H2的0.0133596和HA的0.0131985；这是有限云的离线容量结果。32组在线输出仍全部是原锚点，平均所选 MSE0.0148154。因此当前证据没有显示分层/自适应刷新改善在线行动。[既有汇总](../agent_c/p1_score/SUMMARY.json)、[原锚点选择收据](../agent_c/P1_ADVERSARIAL_QA.json)
- 128条 finite secant 与96个投影对来自有限节点，采样包含被拒 trial 和分支切换，不能描述成连续接受路径或全线段界。[已有诊断收据](../agent_c/p1_diagnostics/ONLINE_DIAGNOSTIC_RECEIPT.json)、[secant记录](../agent_c/p1_diagnostics/FINITE_SECANTS.json)
- N1 的45组运行均未到达非锚点非线性机制比较，保留 `INITIALIZATION_BOUNDARY_FAILURE_NOT_PHYSICAL_TEST`。同锚点刷新不构成新非线性物理信息。N2 的修复、封存与真值访问顺序由主线程/agent_b继续负责，本接续不提前赋予结论。[N1既有收据](../agent_c/N1_BOUNDARY_AUDIT.json)

## 先例短名单及证据等级

本名单复用 agent_a 已完成的 ScholarQA 身份记录，2026-10-30。5/5 身份 resolved；检索中的429、distorted-Born 未覆盖与没有可用全文段落的缺口全部保留。以下均为 **Tier C 身份证据**；排序只是按题目与本轮主题的接近程度挑选待深入比较对象，属于相关性判断，不能支撑实质等价、优越性或新颖性结论。

| 待比较先例 | 年份与标识 | 本次允许的使用 |
|---|---|---|
| Zhong、Chen, Twofold subspace-based optimization method for solving inverse scattering problems | 2009，[DOI](https://doi.org/10.1088/0266-5611/25/8/085003) | twofold/subspace 主题的出版身份；具体投影对象与R12差异仍未实质比较 |
| Zhong、Chen, An FFT Twofold Subspace-Based Optimization Method for Solving Electromagnetic Inverse Scattering Problems | 2011，[DOI](https://doi.org/10.1109/TAP.2010.2103027) | FFT/twofold 主题身份；不能从题名推断与当前 pair 因子相同 |
| Xu、Zhong、Song、Chen、Ran, Multiplicative-Regularized FFT Twofold Subspace-Based Optimization Method for Inverse Scattering Problems | 2015，[DOI](https://doi.org/10.1109/TGRS.2014.2329032) | 正则化/FFT/TSOM 主题身份；出版年2015，DOI内2014不是出版年 |
| Zhang、Ma、Xu、Zhong, Wavelet-Based Subspace Regularization for Solving Highly Nonlinear Inverse Scattering Problems with Contraction Integral Equation | 2020，[DOI](https://doi.org/10.3390/electronics9111760) | 子空间正则化与非线性积分方程的主题身份 |
| Huang、Hu, Accelerating Uncertainty Quantification for Nonlinear Inverse Scattering Problems With High Contrast Media by Direct Envelope Methods and Krylov Subspace Iterative Integral Equation Solvers | 2021，[DOI](https://doi.org/10.1109/TGRS.2020.3035641) | 高对比度逆散射/UQ/求解器主题身份；不能据此认证R12有限云的posterior语义 |

完整既有依据：[CLOSEST_PRIOR_WORK_METADATA.md](../agent_a/CLOSEST_PRIOR_WORK_METADATA.md)、[SCHOLARQA_VERIFIED.json](../agent_a/SCHOLARQA_VERIFIED.json)。本接续没有刷新或扩大检索，也没有把身份核实升级成方法证据。

## 保留的科学边界

新增的是当前模型中可追溯的公式衔接、适用条件和失效解释。标准 SVD、Khatri–Rao、Hadamard Gram、链式法则、固定点微分及约束核空间均不作为新原理。原 P1 连续根流形、完整材料后验、总体覆盖保证与方法新颖性仍未建立；最终科学判断由主线程完成。
