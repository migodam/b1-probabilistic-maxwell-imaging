# NIE公平基线资格：简短审计

结论：**现有 scalar CIE/FFT-TSOM 可作为明确披露的 NIE/CIE family adaptation；本轮没有核实到可独立拿来运行、许可明确的成熟“原2016 NIE”源码基线。** 保留其 `NOT_QUALIFIED_AS_INDEPENDENT_MATURE_NIE_IMPLEMENTATION` 状态，不把self-only、A或另一个短实现改名为该成熟基线。没有找到不证明不存在；本轮不是系统性软件检索。公开原作者代码也不是合法重实现的必要条件：可以依据论文另行实现并验证，但本轮没有完成一条与既有CIE区分开的成熟独立NIE复现。

原始 [Zhong–Lambert–Lesselier–Chen 2016](https://doi.org/10.1109/TAP.2016.2535492) 已从[作者Marc Lambert上传的完整论文](https://www.researchgate.net/publication/296475773_A_New_Integral_Equation_Method_to_Solve_Highly_Nonlinear_Inverse_Scattering_Problems)直接核查：上传声明与日期2016-05-15可见；Eq(3)是空间乘法beta构成的NIE，§III使用FFT-TSOM辅助电流子空间反演。本文已披露数值/实验测试，不能因没有本地成熟源码就称NIE仅是未落地的代数想法。作者全文支持这些身份和结构事实；本轮不复制其完整全文或图，不把作者上传正文等同于可重新分发软件许可证。

CIE2019作者出版XML已存本仓库，§2.1–2.4直接说明CIE-I、Fourier-TSOM、PR current/LS contrast alternation，以及current/material不同beta。当前baseline是这一结构的2D实材料与共同prior适配；新R6 FFT路径把显式IDFT列改成真正FFT action，并没有新创其科学表示。[Zhong–Xu 2019出版原文](https://doi.org/10.3390/jimaging5020027)。

本地差异审查（分析判断，核心科学结论仍归父线程）：

| 对象 | 实际含义 | 公平命名与资格 |
|---|---|---|
| scalar CIE | 对辅助本构残差作beta相关逐点尺度变换，搭配既有TSOM/交替求解 | 成熟结构的NIE/CIE family adaptation；须保留实材料/prior/数据归一化/停止与预算差异 |
| self-only | 用真实离散Gii进入分母，改变离开本构流形时的行权重 | self-term control；不能仅凭分母相似就等同原NIE全部目标和迭代 |
| A材料读出 | 在固定J,Q上使用self absorption后的一次块内offdiagonal residual exchange，原prior和原Phi候选验收 | CIE/FFT-TSOM上的residual/readout插件；不等同scalar NIE/self-only |
| 原2016独立baseline | 原作者FFT-TSOM/NIE的完整合格复现 | 本轮未核实其公开成熟源码、许可、原方案完整停止/验证；不冒名运行 |

[Bevacqua–Isernia作者综述](https://arxiv.org/pdf/2012.15252) §III.A/B区分CS-EB的域积分局部函数与NIE空间beta；域积分函数不能直接认作本离散Green的对角元。§V已有不同重写族的CSI混合先例。这些已知关系限制“新方程/新家族”的原创性表述；具体A效用要由主实验决定，不由名称推出。

软件核查范围：已读root/worktree AGENTS、CUDA catalog/REPRODUCE、FBG-R3 baseline audit/source、R5 factorial，以及R4 prior-art文件和FormulaAtlas固定算子/TSOM投影段。本地research下无独立 `nie*.py/NIE*.py` 实现；实际可运行路径是已审计CIE及其适配。Web检索仅两组NIE/FFT-TSOM源码定向式，未获得可核验作者代码链接；另一宽标题+code检索返回大量同名无关积分方程论文，明确排除，不作为NIE缺失的证据。

检索/验证回执：`audit/NIE_SCHOLARQA_COLLECT.json`是实际一次单并发HTTP429，未谎称Semantic Scholar成功；`NIE_SCHOLARQA_VERIFY.json`另记最终身份验证状态。HAL API返回原论文标题/DOI，但document请求404，见`NIE_PRIMARY_ACCESS.json`。作者ResearchGate全文可读及原论文Eq(3)/§III定位见`NIE_EVIDENCE_LEDGER.json`。因此旧R4“原2016全文未取得”的限制只对旧轮次成立，本轮已经补上作者全文读取；独立成熟代码资格仍未补上。

无NIE数值实验、无GT、无新拟合参数、无远端运行。所有资格判断只涉及所检索证据与本地接口，最终novelty与科学解释由Codex父线程确认。
