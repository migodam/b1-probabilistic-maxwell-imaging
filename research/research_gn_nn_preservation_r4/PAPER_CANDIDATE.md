# Paper candidate — HOLD

**When Data-Dark Current Updates Become Observable: Material-Realizability Feedback in Learned Maxwell Reconstruction**

核心问题：固定Green／接收暗空间中的神经电流更新，经过共享实材料读出与完整物理恢复后，是否仍保护已准确的方向？与直接voxel补全相比，显式本构反馈能否带来可归因且成本合理的真实材料恢复？

|候选主张|当前证据|归因及限制|
|---|---|---|
|材料实现反馈解释辅助暗空间与真实强方向的偏离|精确实现误差与斜导数推导；新TEST真实电流阶段／Twofold／带符号交叉项|标准状态代数／链式法则在学习架构中的应用，不宣称新Maxwell定理；论文价值依赖机制的新实证与直接先例检索|
|受控表示比较决定哪些NN先验有效|504/78/108独立划分、六空间端点退化、三种子、未见4/7/字体、真实成像|voxel clean收益真实；noisy被高效普通停止对照显著缩小；监督先验是额外信息|
|无层内完整逆的共享反馈改善完整电流网络|quadRNN相对J-full部分noisy改善，完整物理验收与费用|只有seed7、未胜voxel、A=.4劣化、16步外推失败；SOM-Net／CSI／Landweber直接先例|

不允许声称：新的全局逆散射定理、概率多解集、严格根遍历、GT唯一性、high-contrast新成功、全闭环收敛、全系统加速或普遍current-domain优势。

关键未完成工作：用相同高效公开初始化，在新鲜holdout上区分精确物理场刷新与无逆近似刷新，验证一个具有独立增量的反馈机制，并让它跨seed胜过voxel及噪声感知普通对照。当前更适合作为机制／负面结果与后续方法的基础，尚不足独立宣称新的TAP成像算法。

直接先例：[SOM-Net, IEEE TGRS, DOI10.1109/TGRS.2022.3217162](https://arxiv.org/abs/2209.03567)；本轮采用的CIE/Fourier适配源码及来源见vendor/PROVENANCE.json。最终归因由主研究者裁决，不由“模块更多”决定。
