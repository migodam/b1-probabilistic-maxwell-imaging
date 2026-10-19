# Research P2：真实材料探索已完成，整体REDIRECT

2026-10-19。原P1、24个fresh parents、A–E与补充强基线均实际执行。不是新的概率posterior论文。

[中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p2_material_fiber_20261019/report.html) · [材料图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p2_material_fiber_20261019/gallery.html) · [证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p2-material-fiber-20261019) · [GitHub目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/p2_material_fiber_20261019/)

- 15/24合格物理解，held-out 4/8；旧B由309.787恢复至23.69左右。
- 预冻结archive在held-out有29个材料去重替代、1个结构变化、0个拓扑改变；概率为null。
- 全8 held-out联合Best-of-4 MAE .06513030→.06424926，仅1.35%；不是单算法同成本胜出。
- 便宜hard-prior是强对照；去掉角向限制后普通Full-J同样能发现替代材料。
- 当前Gaussian/Flow没有增量；A17摊销约为full的11.1倍，停止集成优势主张。
- 最终2236/4000 states、46728/132000 RHS、受控wall319.243s，GPU0。旧账本未改。

阅读顺序：RESEARCH_REPORT_ZH.md → GO_NO_GO_DECISION.md → METHOD_COMPARISON.md → THEORY_REFINEMENTS.md → RESOURCE_SUMMARY_FINAL.json。

METRICS.json是补充对照前的预冻结主分析；SUPPLEMENT_RESULTS.json是评分后EXPOSED强基线审查，二者不混为新的盲测。MaterialHypothesisPackage/FINAL为最终交接包；不相容base留在failed_baseline，绝不作为合格假设。原始payload和所有失败保留。
