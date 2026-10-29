# Research Z-ALG-R2

**限定问题上的代数求解成果，完整证据与失败一并归档。**

- [中文研究报告](RESEARCH_REPORT_ZH.md)
- [关键公式集：F01–F15](FORMULA_COMPENDIUM.md)
- [数学理论与证明](MATHEMATICAL_THEORY.md)
- [算法](ALGORITHM.md) / [可复现命令](REPRODUCE.md)
- [机器结果](RESULTS_SUMMARY.json) / [费用账](COST_LEDGER.md)
- [论文候选与不能声称的内容](PAPER_CANDIDATE.md)
- [独立科学审查](audit/FINAL_SCIENTIFIC_REVIEW.md)

本目录自包含四材料直接求根与48单元成像输入；`run_solver.py`只接受观测、采集几何和事前材料约束。`replay_holdout.py`将冻结的四方法比较写入新的输出目录。发布的私有GT仅在原预测封存后用于评分，现在作为复现证据公开；发布后重放不是新盲测试。

全部新实验支付53,951个完整复RHS。一般理论覆盖2s材料、单复接收、s独立源及明确非退化条件；四材料/六维商和六材料/20维商均已实际盲测；不声称连续Maxwell全域枚举、Full66恢复或严格数值全根认证。
