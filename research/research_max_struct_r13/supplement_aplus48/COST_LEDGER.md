# 成本账本

## 新增实际执行

下表仅计本补充新增消费，来自Windows外部 `GetProcessTimes` 全子进程寿命。所有导入、初始化、候选审计、原物理corrector、最终写回和预算越限原子操作均收费。内部/辅助/history累计计时有重叠，不叠加。

|新增队列|运行数|整子进程CPU秒总和|物理LU|complex solve RHS|队列wall秒|
|---|---:|---:|---:|---:|---:|
|原A+缺失32|32|3803.421875|26493|537072|543.221847|
|被动native状态重放|40|4480.296875|33197|729462|620.447893|
|固定Arc版本与NN接口|112|12373.671875|96749|1773918|1728.826853|
|合计新增|184|20657.390625|156439|3040452|各队列分别报告|

worker CPU总和不是并行历时，也不是从首次准备到发布的总人/机耗时。三个前台连接的实际历时分别564.826334、686.840767、1815.081006秒，另含队列/归档等过程；不能再把它们当子进程CPU加一次。队列controller CPU另为.546875/.71875/2.40625秒。

新数据对象数为0；补充使用原48公开材料。184包括40记录副本与112后续控制，不能称184个独立样本。旧基线、旧A+16、旧NN/G19/DR/direct引用成本可追溯，但不再次计入新增消耗。stack分析中的 `reused_original_run_reference_only` 也包含此前本补充已付的32；新增总账仅在此32队列计一次。

## 固定投影、FFT与求根没有免费

`analysis/costs/DECLARED_OPERATOR_COUNTERS.csv` 分开列物理、DR固定背景投影与辅助操作的所有实际声明字段、总量及报告该字段的回执数。DR的H Cholesky、setup/iteration RHS、port eigensolve、A/C投影、Green/receiver正伴随和multiplier均保留。CIE/A+的Fourier/FFT与条件求根计数留在原 `AUXILIARY_EXECUTION.json`；没记录的计数写未知/不适用，不能用CPU以外的零掩盖。

总RHS采用 `actual_complex_RHS`，已含源/伴随等角色；不再加角色字段。DR背景H三角解RHS与材料依赖Maxwell RHS分开。单次物理LU维度576，原 `full_factorization_dimension` 是跨调用累加计数，不是一次矩阵维数。不同维度的分解不能按“1次”当等价成本。

内存回执提供endpoint RSS，不是完整采样峰值；未测的峰值保持未知。新运行使用CPU、无新增CUDA移植；GPU使用率不被当恢复证据。模型token/金钱费用没有可核验账单，保持未知，不虚构美元数。

## 离线验证与评分

|过程|CPU秒|新增物理LU|新增complex RHS|范围|
|---|---:|---:|---:|---|
|五方法全48指标/图|33.558605|0|0|封存后的既有选择评分|
|已付费候选原Phi/梯度诊断|14.010751|1262|25320|1920引用去重；原离散，不是独立验证|
|native状态原物理与梯度诊断|22.952092|844|18924|561状态，844唯一材料；非新反演|

stack完整指标/NN图的实际CPU/wall见 `analysis/stack_comparison/COMPARISON_RECEIPT.json`。动态图初次/排版修订、历史提取、表标签修正、费用汇总及后续归档/上传均各自保留实际receipt，不另起一个混合“算法秒数”。分析内部重叠timers不加到whole-process总数。

新增在线＋上述离线物理求值共158545 LU、3084696 complex RHS；作总资源盘点时可相加，作在线性能表时必须分开。

## NN与历史成本

本次没有训练、没有重新推理。复用的R13原16推理准备为1.285002 Mac CPU秒，包含模型加载、特征/3接收SVD/144列动作和写盘；最终receipt/manifest写回及exit不在其计时内。R8主训练＋训练/验证生成为402.558889 Mac CPU秒，其他历史训练/开发另计，来源 `research_max_struct_r8/NEURAL_WARM_STARTS.md:55`。

同16 NN对照表主列是solver CPU；独立部署应另付适用的特征/推理与训练摊销。Mac与i9秒数不能直接相加制造端到端速度排名。共享setup只在相应原准备回执记一次，单次部署与批处理摊销须注明。

## 失败与资源诊断

首次脱离SSH启动失败、源码默认编码提取错误、一次源归档控制脚本转义语法错误均留在操作回执；后者发生于远端代码执行前，新增科学求值0。184终态科学任务无运行失败；材料恢复失败仍完整保留。

CPU低频审计是只读元数据，远端完整CPU未独立测量、不得填零。它没有更改系统设置；具体往返耗时和本机传输CPU见 `audit/REMOTE_METADATA_COSTS.jsonl` 与 `analysis/cpu_frequency/`。传输/下载/压缩与GitHub发布回执另存；没有将全部费用藏进setup或忽略失败。

## 同成本的具体含义

1. 相同120上限是部署合同，不等于每方法消费相同。
2. common actual CPU只选择当时已完成物理审计的候选，不从GT选前缀；它不是新的实测停止wall。
3. 同物理LU仅作物理预算诊断，不等额化固定投影/FFT/读出/setup。
4. NN uniform-six与NN-six的限制不同；同初值效应用所有12流水线共同限制，不能错配两个six预算。

明细：`analysis/costs/NEW_RUN_COSTS.csv`、`NEW_CAMPAIGN_TOTALS.csv`，以及两个比较目录的 `COST_LEDGER.csv`、`COST_TOTALS.csv`。
