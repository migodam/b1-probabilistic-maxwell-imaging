# Research W P0 — Gate decision

|Gate|状态|理由|
|---|---|---|
|0 解析复现|PASS|完整一维参考面、PEC级联、幺正性及带宽负结果独立复现|
|1 二维资格|UNQUALIFIED_SOLVER|128→192修订版多频完整数据白化差0.951–1.229，未达0.1|
|2 Calibration机制|NOT_RUN_UNQUALIFIED_SOLVER|未启动候选搜索，没有合格有限墙体见证|
|3 Target-present joint|NOT_RUN_UNQUALIFIED_SOLVER|没有联合目标/仪器补偿实验|
|4 带宽/微扰鲁棒性|UNRESOLVED_2D|只有一维指定补偿负结果，不能外推二维；未做轻微扰动压力测试|
|5 算法/P1价值|HOLD|未有qualified joint见证；P1没有执行授权|

总体：**P0_B_UNRESOLVED_NUMERICAL_QUALIFICATION**。没有耗尽总预算；是资格门槛阻止搜索。不能把20项测试通过、258个离散状态或四块恒等式误写成物理见证成功。物理存在性与搜索算法价值均没有被否定或确认。

精确数值/统计位置：[results/GATE_DECISION.json](results/GATE_DECISION.json)、[RESULTS.csv](RESULTS.csv)。
