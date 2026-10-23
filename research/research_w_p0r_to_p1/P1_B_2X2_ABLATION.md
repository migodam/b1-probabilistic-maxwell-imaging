# 四格消融：多假设收益与递进调度收益分开

以下均来自同一个新受控prior和16个parent各两份独立噪声实例，共32份新观测。统计单位是16个独立参数parent；同parent的两份噪声先平均。它们不代表16种独立墙体几何：墙体只有两个固定场。所有对象、包括概率CDF误差超标者均进入评分。

| 方法 | 目标x平方误差 | 目标x绝对误差 |
|---|---:|---:|
| direct_point | 0.0044293822 | 0.04873760 |
| progressive_point | 0.0044291620 | 0.04874046 |
| direct_SMC_mean | 0.0025451075 | 0.04739516 |
| progressive_SMC_mean | 0.0025466091 | 0.04736837 |
| continuous_quadrature_Bayes_mean | 0.0025566111 | 0.04749423 |
| zero_observation_prior_position_mean | 0.0023661589 | 0.04830233 |
| full_validated_continuous_point_control | 0.0040418124 | 0.04567858 |

主四格中，多假设direct相对direct格点point的平方损失下降42.54%，parent配对差为0.00188427，2000次percentile bootstrap 95%区间[0.00043418,0.00349844]。相对更强的actual-Full连续点控制，下降37.03%，差为0.00149670，区间[0.00029834,0.00279414]。这些是有限受控panel的描述性区间，未作多重比较调整；MAE没有相应确定改善。

递进SMC与direct SMC的平方损失差（direct减progressive）为−0.0000015016，区间[−0.0000178415,0.0000159221]；点方法差也接近零。没有独立调度风险收益。

零观测prior均值的平方误差为0.00236616，比本panel的direct SMC和连续求积Bayes均值还低；prior相对direct SMC的配对差为−0.00017895，区间[−0.00073936,0.00038532]。不能将这16个对象上的结果解释为正确posterior总体风险劣于prior：总体Bayes风险性质与有限样本波动不同。但也不能声称本panel证明了观测提供额外成像收益。

连续求积Bayes位置均值可位于两个位置区间之间的空隙。它是平方任务损失的动作，不是第三个合法目标位置，更不是第三个真实材料世界。方法未用GT选择分支或Best-of-K。机器表：results/P1B/offline_scores/。
