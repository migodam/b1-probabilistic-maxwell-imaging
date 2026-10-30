# N2 离线简报

9 个设计数字 × 6 条方法已在主线程实际盲读后，用原评分器离线评分一次，新增物理求解为 0。**六条方法 ε=.075 材料覆盖均 0/9。** 数字 6 每条方法的合格云均为空，保留在分母；所有方法在线输出真实 noise 门均 8/9。

|方法|平均在线材料 MSE，全部 9 个对象|选点 Q 通过|选点 noise 通过|材料/occupancy 拓扑 exact|
|---|---:|---:|---:|---|
|H0_GN|.02880951|8/9|8/9|5/9；3/9|
|H0_DATA_RML|.03036269|6/9|8/9|5/9；4/9|
|Fixed8、H1、H2、HA，各臂|.03171430|8/9|8/9|5/9；4/9|

Q 是 noise、合法和材料 RMS .002 贪心去重旗的合取，Q0 可能只是重复，不能等同 noise 失败。逐方法 CSV：[`N2_METHOD_SUMMARY.csv`](N2_METHOD_SUMMARY.csv)。

数字 8 的每臂 7 个 raw 候选，共 42 个，在 Reχ≥.75 和 occupancy≥.5、FG4/BG8 下都没有 GT 的双孔。主线程打开最终 GT 比较及六法材料 mask 后也记录了该失败。此结论只针对已评估点；最近 GT 候选未替代在线输出。

主线程原盲读中的真实 3 读作 `5or3`（2/5 分），真实 5 读作 `3`（3/5 分），真实 8 为 `0or8`（1/5 分）；这些分数及 A05/A06 材料视图近空白的原备注保持不变。原始盲读与真值/方法映射在 [`N2_BLIND_READ_MAPPING.csv`](N2_BLIND_READ_MAPPING.csv)。该 benchmark 此前已暴露，不能称为新的独立人类留出评测。

原 A05/A06 近空白属于原盲读显示限制：存储材料最大实部为 1.144–1.464，χ=host+δa 残差≤1.12e−16，不是材料幅度塌缩。最终图直接展示原封存 24×24 数组、统一色阶，无插值；原盲图与原记录保留。

主图：[`DIGITS_123456789_RECONSTRUCTION.png`](../figures/n2_online_and_clouds/DIGITS_123456789_RECONSTRUCTION.png)。数字 8 全候选图已按主线程要求最终调整行距；Q/noise/合法旗及失败候选保留。全部图入口：[`FIGURE_CONTENTS.json`](../figures/n2_online_and_clouds/FIGURE_CONTENTS.json)。详细数值说明：[`N2_NUMERICAL_REPORT_ZH.md`](N2_NUMERICAL_REPORT_ZH.md)。
