# Research F — 来源审计与证据等级
日期：2026-10-16。正文使用 [R1]–[R9]、[L1]–[L12]；可复现定位见正文参考文献。
本目录没有打包第三方论文全文、用户原始二进制实验数据或字体文件。

## 已直接读取的项目内容
- B1 finite_feedback_fiber_diagnostic，固定 c76733676a7746627dfbf710c5f31b6375ae4335：
  START_HERE、中文完整报告、METRICS.json 核心前210行、REPRODUCE、
  LOCAL_PATH_AND_CURVATURE、THEORY_QUESTIONS_FOR_PRO。
  这是报告/部分机器可读数值核对，不是45个NPZ端点的独立重放。
- B1 LXZ，固定 dde0237eaa0d8c5507fedb07bb704c66a817028d：
  RESEARCH_REPORT_ZH.md，包括后验比较、TRAIN chart、失败和账本。
  其“完整读L/Z/X 41文本”是原项目的工作，不是本次已完成的读取范围。
- A26 decoder.py、training.py 的数据边界/归一化/损失有关代码、
  pilot.json、A26_RESEARCH_REPORT_ZH.md、REPRODUCING_A26.md。
  training.py 的后续定位使用固定2603f3bd16c48e48587fedc4363975779def17f5；
  blob e04ff8c44af75197848d13907de3bfdd0318cf70。
  原 loss 明确为复增量MSE + 实/虚分通道SSIM，不是幅值SSIM。
- A17 README、closure_final/START_HERE、PAPER_TAP_FINAL.tex、CORE_THEORY.tex。
  CORE_THEORY blob c31c5c6790820b153d2b6a738df94fda0c12165d。
  读取时默认分支；本次没有审计全部42个外部模块或重跑3D。
- A30 A30_E0_E1_DECISION.md，固定717d3d26a290248bd2b9aeb55f3d03694b407b84。
- C1 START_HERE.md 和B1 START_HERE.md 用于Joint P2与T/G/Prior-Constrained历史定位；
  不把入口摘要当成这些项目所有原始样本的独立复核。

## 跨会话
调用了可用个人上下文检索，返回相关研究的部分历史摘要/片段。
未取得用户列出的全部Session全文。未完整展开L/Z/X原始ZIP。
报告根据已取得仓库材料和明确的数学推导独立形成，没有把缺失内容写成已读。

## 文献读取深度
全文可解析并阅读了相关方法/数学段落：
1. e-GLASS，arXiv2207.00460，5页。
2. Deep Null Space Learning，arXiv1806.06137，相关正则化定理与证明。
3. Manifold Lifting，JRSSB2023，26页，尤其coarea/lifting/约束MCMC。
4. DINO，JMLR2025，68页，尤其目标、导数学习、几何提案和delayed acceptance。
“全文可解析”不意味着本次逐行独立复核了所有附录证明与所有实验表格。

部分正文/原论文公开方法说明或摘要级：
- Boink等非线性数据一致网络：出版页、原摘要与方法/投影说明；
  PDF正文获取未成功，因此未完成全部证明核对。
- Gabrié等NF增强MCMC：原摘要、方法与应用条件；全文入口遇到访问验证。
- DILI、CSI：原论文/作者来源的方法信息；未声称完整逐式新颖性排查。
- Maxwell Born/inverse Born 2017、Padé electromagnetic Born 2024：原出版摘要与记录。
- Magnus综述：原论文摘要/出版信息；本报告端点传播式是直接推导。
- 2026-08-18 inverse-Born neural operators：arXiv2608.18262摘要；
  HTML/PDF获取失败，列为尚未完成的最新近邻检查。

PDF截图接口多次返回内部错误；本报告不以未核对的图像、图轴或表格作为数值来源。
项目数字来自仓库文字/JSON，外部论文不提供本项目新增实验数字。

## 本次新增计算
仅运行 checks/algebra_checks.py：
六维随机矩阵恒等式、路径ODE、两单元Maxwell-form代数、坐标反例与软弱滤波界。
new_full_maxwell_states=0
new_complex_maxwell_rhs=0
new_mie_calls=0
gpu_training=0
historical_replays=0

“六维矩阵上的线性方程”不是既有物理实验后端的RHS，不计成新Maxwell evidence。
其作用是检查式子与代码一致性，不是性能benchmark。

## 结论等级
- 标准数学：常秩、coarea、MH、谱扰动、Grönwall、Tikhonov滤波。
- 本报告直接推导：有限响应递推、静态场端点传播、
  prior/coordinate/Maxwell曲率分解、材料fiber第二基本形式、
  多源current realizability消元、有限primal/dual误差接口。
- 条件性结果：定量局部chart、学习场稳定、余核曲率步长障碍、
  posterior小噪声极限、minorization覆盖界。
- 数值检查：小矩阵结果，无新的复杂材料/3D/学习收益。
- 待验证假设：现有P0 raw decoder曲率有多大比例是坐标效应；
  material-native retraction是否改善同成本有限位移；
  A17有限传输是否有增量；P1是否存在重要复杂多模式后验。

本报告不提供“全球首次”“已解决一般Maxwell多解成像”或“强TAP已成立”的保证。
