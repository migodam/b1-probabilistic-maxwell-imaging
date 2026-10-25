# B：风险谱没有自动转成数据外增量

固定rank2/4/8/16与ridge网格，比较raw-y、decoder-only、physical-residual、same-spectrum随机四种表示。使用交叉矩的RRR，不用GT残差PCA冒充可修复风险。

实际18个baseline训练构造诚实嵌套二折残差。各最终修正均没有改善OOF/VAL风险，无法支持Maxwell约2%独立增量。几何仍条件于历史无监督frame，不能称完全fold-inductive新几何。

[RRR原始表](stage_b_material_risk/NESTED_RRR_SUMMARY.json) · [独立折分审查](audit/B_NESTED_OOF_REVIEW.json) · [完整B解释](stage_b_material_risk/REPORT_ZH.md)
