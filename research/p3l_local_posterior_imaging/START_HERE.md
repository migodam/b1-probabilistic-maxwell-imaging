# P3-L：已完成并封账 · REDIRECT

**当前近似概率成像没有改善重建。** 128个冻结新对象全部保留：材料均值MSE **0.0615021**，probit点估计 **0.0110309**，误差增加 **457.5%**；配对差95%区间 **[0.0427873,0.0582865]**。

这不是对真正posterior mean的否定。**11组四链Full参考全部未合格**，不能提供校准posterior、teacher或材料credible coverage。

实际完成800个独立parent的观测、32个development基础反演、128个冻结成像评估、44条链、SMC、40组强条件实验、解析材料矩、受限SLICE积分、15次神经训练、简单密度与Full校正、两臂物理代理力、合法投影及A17/普通求值对照。未执行的teacher distillation、RL与完整trajectory-tube概率研究均保留NOT_RUN，不冒充完成。

点估计的probit坐标修订有11.71%独立MSE改善，但它不是native MAP或Bayesian平均。15次神经训练均选中初始化；Full-SNIS的128个对象均ESS<2。保留decoder降低代理误差，却没有改善材料MCSE且更慢。

实际费用 **8,374,472 complex RHS、707,328 state attempts、3.805h受控wall，GPU0**。历史账本不变，中断/失败计费，剩余3,625,528RHS封账。

[中文完整报告](RESEARCH_REPORT_ZH.md) · [网页版报告](report.html) · [逐列图像说明](IMAGE_READING_GUIDE.md) · [实部图](D/locked_evaluation/figures/fixed_first4_epsilon_real.png) · [损耗图](D/locked_evaluation/figures/fixed_first4_epsilon_imag.png)

[原始128对象风险与CI](D/LOCKED_IMAGE_RISK.json) · [独立审计](AUDIT/locked_results_review/REVIEW.md) · [参考资格](POSTERIOR_REFERENCE_AUDIT.md) · [最终裁决](GO_NO_GO.md) · [成本](RESOURCE_SUMMARY.md) · [复现](REPRODUCTION.md)

[同一B1的证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3l-local-posterior-20261022)：大数组、原始chains、模型权重与完整追加ledger。上传/匿名读取状态单独记录，不以网页可达替代科学通过。

本轮是二维scalar-TM/VIE、原P1 prior；没有新Mie、三维vector-Maxwell或外部云计算。停止当前Gaussian成像与无收益NN/A17整合，转向可信Full参考及非Gaussian条件机制；当前不具备强IEEE TAP算法证据。
