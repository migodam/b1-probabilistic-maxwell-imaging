Research P3L-R7-T — 交付说明
日期：2026-10-28（Asia/Singapore）
证据基线：d1bdc8aceae4cd6a843f873b61b019d5eee3cd1d

建议阅读顺序
1. Research_P3LR7T_Report_ZH.pdf：完整中文研究报告。
2. Codex_R7_Minimal_Experiment_Prompt.txt：可直接复制的下一轮任务书；默认 MODE=A。
3. CLAIM_LEDGER.json：逐主张的资格、证据与下一门槛。
4. REVIEW_AND_SCOPE.txt：本轮完成范围及必要限定。
5. source/：可编辑 LaTeX 正文、26 条原始文献索引和解析图源。

主结论
材料商和概率体积是理论深化；材料任务是有独立风险证据的分叉；减少 Native 当前属于工程优化。
优先检验非均匀材料上的真实补偿及其被新采集破坏的机制。暂不合写 A+B+C 的统一大论文。
Paper A 条件 GO、理论优先；Paper B 条件 GO、算法基础较成熟；Paper C 暂为计算模块。

下一轮任务
默认 E1+E2，上限 50,000 complex RHS，不启动大型网络。
MODE=B 的四维条件任务积分单独封账；总上限 1,600,000 RHS，先做 reference pilot。
两种预算均是建议的下一轮合同，本轮未执行任何新 Maxwell 求解或训练。

证据资格
已读上传 R6 全部指定报告、固定版本八份历史报告及 A17 两份最终理论文件。
R3 与 R5T 的完整独立 Pro 原报告未取得；相关历史采用用户本次说明及 R5/R6 引用，正文已标记。
未重新发布原始 R6 ZIP；请与本报告一同保留用户原始证据包。

重建 PDF
需 XeLaTeX（TeX Live）、Noto Serif CJK SC、常用 amsmath / fontspec / tcolorbox / tikz / needspace 包。
在 source/ 目录执行 xelatex report.tex，重复至交叉引用稳定。
若字体没有安装到系统，可以把 NotoSerifCJKsc-Regular.otf 放在自选目录并执行：
xelatex '\def\ReportFontPath{/absolute/font/directory/}\input{report.tex}'
本包不包含第三方字体。PDF 已嵌入显示所需字形。

解析图复现
在 source/ 执行 python3 make_figures.py，需要 numpy、scipy、matplotlib。
这三张图分别是共享补偿的代数反例、原子分支的精确积分、三点先验的任务信息。
它们均不调用 Maxwell，不构成新物理实验。

原文及来源
所有公式、条件、证明和引用保留在 source/01...09.tex。
正文内部证据与外部文献分别编号。2026 年新预印本明确标记，不把预印本性能当独立验证事实。

