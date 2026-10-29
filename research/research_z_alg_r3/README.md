# Research Z-ALG-R3

从局部行列式商根搜索，转向全部观测的共同材料兼容性候选，并检验真实盲成像。完整非线性二维 scalar-TM VIE；12个实际材料单元全部未知。

**有限结果：** 最终四个新盲目标，保护＋兼容性候选达到3/4个预声明数据门，保护＋倒数Full-J达到2/4。新增案例在线材料RMS为2.12e-10对10.94。总wall没有优势；独立论文新颖性及一般优势仍为HOLD。原局部Q6路线产生远候选但未证明成像增益，失败原样保留。

- [中文研究报告](RESEARCH_REPORT_ZH.md)
- [核心数学与归因](MATHEMATICAL_THEORY.md) · [公式复用](FORMULA_REUSE_MAP.md)
- [latent与characteristic subspace如何接入](theory/LATENT_AND_CHARACTERISTIC_SUBSPACES.md)（仅推导，未运行Full66）
- [算法](ALGORITHM.md) · [直接复现](REPRODUCE.md) · [盲测合同](BLIND_EXPERIMENT.md)
- [机器结果](RESULTS_SUMMARY.json) · [完整费用](COST_LEDGER.md) · [独立科学审查](audit/FINAL_SCIENTIFIC_REVIEW.md)
- [论文候选与HOLD原因](PAPER_CANDIDATE.md)

```bash
python run_portfolio.py --public public_data/confirm/NEWCONFIRM12_00.npz \
  --method guarded_globaldual --output /tmp/zalg_r3_prediction.json
```

公开GT位于`truth_after_freeze/`，在预测和评分完成以后才导出，只含材料数组。生成seed保留本地，未发布。公开数据重放不是新的盲测。

研究沿用既有仓库，不新建repository：

- [网页报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r3/)
- [GitHub目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_z_alg_r3)
- [Release与完整/精简压缩包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r3-20261030)

所有物理费用、失败、终止超时、源码和预测冻结顺序见[实验技术记录](EXPERIMENT_TECHNICAL_NOTES.md)。本工作不重新运行、修改或覆盖R2/R9的冻结研究。
