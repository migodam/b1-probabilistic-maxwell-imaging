# A参考资格

参考是明确归一化的二维截断子模型，不是完整66维边缘。原坐标采用全局/模式适配求积，独立响应坐标求逆积分。正常case有很好的证据/均值阶数收敛；宿主和复杂饱和case仍未普遍收敛，且有限网格不证明全局遗漏质量小于.01。

TV/KL只用于可归一化连续Gaussian近似对数值参考的比较，不将离散求积测度称作低TV样本分布。RTO为真实局部求逆尝试，全域归一化未建立。Full-66、普通Spline与manifold lifting NOT_RUN。

[逐case资格](stage_a_response_coordinates/OBSERVATION_SUMMARY.csv) · [独立A审查](audit/A_INDEPENDENT_REVIEW.md) · [混合秩审查](audit/A_V3_MIXED_RANK_REVIEW.json) · [最终阶段QA](stage_a_response_coordinates/A_FINAL_QA.json)
