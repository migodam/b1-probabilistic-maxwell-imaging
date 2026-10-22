# Research W continuation：高精度前向通过，物理歧义尚未得到见证

**P0-Q PASS / QUALIFIED_DISCRETE_FULLWAVE；P0-B NOT_FOUND_WITHIN_ALLOCATED_SEARCH_BUDGET；P1 HOLD。**

- 原两墙完整数据细化差0.0065–0.0332，低于原0.1；高阶贴合FEM的线源圆柱相对误差0.0266%。
- 三个context四方法完成，最好cal白化距离19.91，未达到1；第四个context缺失。没有joint、posterior、NN或三维。
- 搜索还存在明确实施缺口：没有严格执行原0.1最小特征条件；候选不能用于正向物理见证，不据失败推断唯一性。
- 不把搜索失败当唯一性，也不把数值平台通过当成成像优势。

[完整中文报告](RESEARCH_W_CONTINUATION_REPORT_ZH.md) · [科学裁决](SCIENTIFIC_DECISION.md) · [机器判决](GATE_DECISION.json) · [复现](REPRODUCE.md)

[圆柱图](figures/01_cylinder_accuracy.png) · [八项收敛](figures/02_wall_convergence.png) · [真实参考世界/场](figures/04_known_world_and_field.png)

代码在src/wq，原始指标与小数组在results，大型复场数组在同仓库Release research-w-continuation-20261023，路径见ARRAY_MANIFEST.json。旧P0保持原目录/commit，费用分开。
