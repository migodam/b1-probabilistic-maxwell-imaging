# Research P3-L — 完整研究包

本包研究同一局部 Maxwell 材料解族中的概率、不确定性与材料成像，日期为 2026-10-19。

## 阅读顺序

1. `Research_P3L_Report_ZH.pdf`：完整中文报告，包含科学判断、六个核心条件命题、22 类算法、四个深入方案、五阶段实验 campaign 和 45 项研究文献。
2. `Research_P3L_Report_ZH.md`：报告的完整可编辑来源。
3. `technical_appendices/`：几何与风险、概率算法、Maxwell 与输运三个独立技术附录。局部文献和公式编号只在各自附录内有效。实验配额与最终判断以主报告为准。
4. `EXPERIMENT_PROTOCOL.yaml`：未来实验的机器可读建议，含维数、parent 划分、质量门槛与费用上限；不授权执行新物理任务。
5. `checks/`：已执行的数学控制实验、结果 JSON、科学图和脚本。所有检查都不调用 Maxwell。
6. `evidence/SOURCE_MANIFEST.json`：24 个已获取项目文件的固定 URL、字节数和 SHA256。原始文件不重复装入研究包。
7. `evidence/CLAIM_REGISTER.json`：关键主张的证据/推导/待验证状态。
8. `typesetting/`：PDF 构建程序和样式文件；字体获取与重建方法见其中说明。

## 已完成与尚未完成

已完成：固定项目证据审查、原始文献检索、独立条件推导、算法比较、零 Maxwell 数学核验、可执行实验设计。

未运行：新的 Maxwell 求解、posterior reference、网络训练、RL rollout、付费物理实验或大规模硬件测试。报告中的未来费用为建议 cap，不是运行记录或实际费用预测。

主要科研路线：强条件厚度与材料矩积分的局部 Bayesian imaging；并行研究可认证 Maxwell 后验决策，以及轻量条件密度/材料矩蒸馏。

## 重现数学控制

安装或使用已有 Python、NumPy、SciPy、Matplotlib，运行：

```bash
python checks/run_math_checks.py
python checks/check_frobenius.py
```

输出覆盖 profile 与真实边缘的体积反例、带线性项的截断 Gaussian、同随机数 MH 区间决策，以及 P1 饱和 decoder 的 Gaussian 矩。脚本写到自己的 `checks` 目录，不联网。

`checks/frobenius_check.json` 记录独立弱平面不可积反例的数值差分结果。它是普遍假设的反例核验，不是 P1 的实际 Maxwell 谱已经不可积的证据。

## 概率与报告约定

原 prior/likelihood、weak marginal、profile law、trajectory law 与 supervised conditional density 始终分开。局部 posterior 的逃逸质量未知时，不作完整 posterior 覆盖声明。P2 候选档案不能冒充 posterior teacher，oracle Best-of-K 不作部署风险成绩。

所有实验比较包括失败对象和 fallback。成像主终点为 complex material MSE/MAE、Bayes action、calibration 与材料函数 ESS/总费用。
