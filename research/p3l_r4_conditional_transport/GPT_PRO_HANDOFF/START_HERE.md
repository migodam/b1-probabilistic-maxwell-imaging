# 下一轮 GPT Pro：基于新实验证据，而非算法清单

先读[完整报告](../RESEARCH_REPORT_ZH.md)、[裁决](../GO_NO_GO.md)、[评分数学解释](../AUDIT/JOINT_SIMULATION_PROPER_SCORE_INTERPRETATION.md)与[机器指标](../RESULTS_AND_METRICS.json)。本轮真实完成多模型训练、768新对象独立评价及物理重放。

[交接阅读顺序、证据版本与立即复核入口](READING_GUIDE.md)区分阅读ZIP、完整原始分卷和封账后补充，避免把小包当成全部模型/物理证据。

已公开：[本轮Pro阅读ZIP](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-p3lr4-20261025/Research_P3LR4_Complete_Results_and_Handoff.zip)（约5.2MB）· [完整原始证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-p3lr4-20261025)（175GB分卷）。最终出版读回晚于阅读ZIP生成，回执另附，不改原包或科研账本。

## 十项明确答复

1. Spline有稳定显式full66 proper-NLL收益；CFM有样本energy-score收益，但其数值log-density未过稳定门，二者不能按一个虚构统一NLL排名。
2. 完整posterior失败具体表现为原噪声尺度的response失配、模型center和dispersion共同贡献、重要性权重退化及物理训练后的coverage恶化；不是只报ESS小。
3. 相阈值恒等式与精确bin正确，当前head不胜CNN；anchor mass初始化还出现平台梯度抑制。
4. 直接生成材料均值不胜CNN；作为全局context却有独立风险增益。Gaussian一样有效，不能归因于Flow概率正确。
5. 静态Maxwell frame/权重对一些proper scores与小幅风险有重复增益；动态weak/bracket没有匹配成本优势。
6. 当前J/动作学习弱于便宜field+exactdecoder更新，且有昂贵teacher/候选Full费用，不值得直接扩大。
7. active-row finite feedback没有净优势；残差/势区间驱动数值重放有66%wall收益，范围限定且普通Krylov也有效。
8. 非完整路径没有材料统计增量；直接法更容易达到同RMS。111有限噪声候选不等于modes或posterior mass。
9. 本轮没有新Maxwell定理已被证明。标准conditional-risk、proper-score、resolvent和retrospective判定必须正确引用，不能包装原创性。
10. 最值得继续的问题在下方，必须围绕现象提出新机制。

## 三个具体数学问题

**一：response calibration与forward-KL覆盖为何冲突？** 同容量Spline已显著改进joint-test NLL，固定draw Phi仍约8千；加入真实Maxwell反向KL降Phi却使NLL增加逾10nat。能否根据decoder饱和image、response共维数、条件厚度与有限容量，构造不坍缩的可逆概率输运或可量化的response-tube factorization？需要给出比再次加大likelihood权重更具体的条件。

**二：局部patch缺失的任务信息是什么？** 全局Spline均值很差，却帮助6,418参数局部head；Gaussian几乎相同，空间打乱失去收益，约99.9%增益落在基础图过渡相。能否把此现象写成可计算的transition-cell conditional-risk decomposition，区分真正Maxwell sensitivity信息与普通空间/属性context？标准条件期望正交性只是起点，不是新定理。

**三：何时只为概率决定付物理精度有利？** 相同Full核上已实现数值决定零差异和66%wall节省，大步长则fallback/拒绝上升。能否把preconditioner收缩、current injection、acceptance margin与material-function mixing一起优化，而非只优化同一坏混合核的求值？需包括缓存setup、有效统计移动、浮点包络和更大二维系统成本。

## 不应重复的路线

不再把Gaussian强补全当FULL posterior；不拿纯trajectory频率当概率；不重复Euler四段错误；不把有限Neumann恒等式当多解定理；不将teacher/anchor成本视为免费；不使用两批已暴露GT做新blind claim。

本轮GO-A成立，GO-B仅相对，GO-C仅限定数值核，GO-D动态部分未成立，GO-E尚缺关键机制。新理论可以推翻原假设，但必须解释实际数组与负结果。

归因缺口需保留：CFM没有专门的同谱decoder-only权重臂；none/随机控制的正面结果尚未分离DF·Dg中传播与decoder各自的贡献。已有风险CNN的decoder消融不能跨算法代替它。下一轮若要主张Maxwell特有机制，应先补上该辨别性对照，而不是只扩大模型。


## 封账后的数学解释补充（零新增物理）

固定样本的PIT尾事件还可通过标准数据处理不等式，给出平均条件材料TV/KL的保守下界，不需要MCMC teacher。完整26模型×2批校正后，rawSpline seed1的两批TV下界为0.953/0.939；这是原模型与理想独立随机流前提下的post-hoc解释，不是每观测保证、预注册主结果或新Maxwell定理。详见[数学证明与全部限定](../THEORY_CALIBRATION_ERROR_BOUND.md)。
