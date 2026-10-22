# 图像应当怎样读

[实部图](D/locked_evaluation/figures/fixed_first4_epsilon_real.png) 与 [损耗图](D/locked_evaluation/figures/fixed_first4_epsilon_imag.png) 使用冻结的前四个新对象，不按最后效果选例。Re 与 Im 各自共享色标。它们是二维 TM24 材料，不是三维矢量成像。

从左到右分别是：离线真值、16-state 便宜 context 图、独立 native96 点估计、probit64 坐标点估计、learned weak density × strong Gaussian 参考的材料平均、strong conditional-mode plugin 平均、原 prior 材料均值、128-draw Full-target SNIS 平均。

最后三类平均的概率资格都不能从图像外观推出。probit64 也不是 native MAP 或 MMSE；它是在另一坐标密度下、有限优化预算得到的点估计。

第一个与第三个对象中，参考材料平均出现很明显的弥散激活；强 conditional-mode plugin 也有类似现象。它们不是新发现的高概率材料世界，而是当前近似密度经过非线性 decoder 平均后的结果。完整128对象的 MSE 确实恶化，所以不能仅凭“看上去更平滑”称为 Bayesian 去噪成功。

原 prior 平均呈现均匀的圆形 ROI，因为此处没有使用观测。它是零物理信息的参考，不是另一个有数据支持的解。

Full-target SNIS 使用真实 Maxwell likelihood 计算权重，但128个对象的 ESS 都小于2。某一列看起来更锐利，可能只是有限权重集中在一两个 proposal，而非成功恢复 posterior。使用正确权重、做过 Full 求值、通过数据相容性与建立参考资格是不同事情。

所有这些图都是单个估计器的材料图，不是附带已认证质量的多个 posterior modes。真正的 posterior variance、credible region 与完整模式质量仍没有得到资格。图中真值只由隔离的离线 evaluator 在全部在线输出冻结和保存后读取。

[完整风险与配对区间](D/LOCKED_IMAGE_RISK.json) · [独立统计审计](AUDIT/locked_results_review/REVIEW.md) · [条件与近似的数学含义](THEORY_AND_ALGORITHM.md)
