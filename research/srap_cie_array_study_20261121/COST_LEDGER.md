# 实际计算费用

单臂实际墙钟合计包含原始轨迹与两条额外的 15/30 秒轻校正分支。两个 raw 前缀共享轨迹，仅计一次。状态构建、FFT、背景谱、物理审计、未入选状态和截止后才返回的计算均计入。各并行臂墙钟之和计作计算使用量。

| 数据组 | 主机 | 轨迹 | 累计臂墙钟 s | 累计单进程 CPU s | 无可用材料的前缀 | 超过预算返回最大 s |
|---|---|---:|---:|---:|---:|---:|
|arc_cap|XINAN|32|1266.5|811.0|0|0.93|
|arc_cap|migodamdeMac-mini.local|16|169.4|166.1|0|0.00|
|confirmation|XINAN|224|8774.0|8692.7|0|0.57|
|confirmation|migodamdeMac-mini.local|112|2458.6|1911.7|0|0.14|
|development|migodamdeMac-mini.local|32|739.2|550.9|0|0.10|
|followup|migodamdeMac-mini.local|96|2161.0|1847.8|0|0.31|
|lazy_fresnel|migodamdeMac-mini.local|24|853.0|708.5|0|1.75|
|legacy|XINAN|70|2758.9|1838.2|0|29.26|
|legacy_mac|migodamdeMac-mini.local|36|864.5|779.5|0|0.09|
|line_cap|XINAN|64|2544.7|1492.4|0|1.15|
|line_cap|migodamdeMac-mini.local|32|371.0|365.3|0|0.00|
|measured_extra|XINAN|4|155.6|153.1|0|1.27|

共 742 条二维轨迹，累计 5.366 CPU 小时。原始求解停止状态为 {'SCHEDULE_COMPLETE': 322, 'WALL_LIMIT': 420}。共同 corrector 的独立状态、截止前完成候选的可用时间见逐例表。所有原始求解臂只有各自 15、30 秒前已完成审计的候选能参加相应预算的比较。单次不可打断的操作可能越过截止，迟到结果无法入选。

每个任务读取观测包、导入运行库及最后序列化状态发生在臂计时以外；冷启动物理算子及分解发生在臂内。并发运行会影响可达的阶段，因此同一对象的两臂在同机顺序测量。不同硬件的胜率只在各机配对后计入。

[逐轨迹完整计数与超时](results/cost_per_trajectory.csv) · [同工作量跨主机测量](engine_dev/HARDWARE.md) · [二维测量生成计数](data/generation_cost.json) · [补集生成计数](data_followup/generation_cost.json) · [三维 CUDA 计数](vector/REPORT_ZH.md)

新 24 件二维材料共享 744 个细网前向右端。生成费用和线性方程残差逐对象存档；实测无需额外目标生成。矩阵自由 H 求解的 FFT 与合法 Fourier 回注的 FFT 分开记数。二维 32² 的原 Maxwell 审计仍用稠密材料相关 LU，60² 使用既有大网格物理正演。不同网格的求解计数应结合实际时间阅读。

Mac 工作进程的5秒间隔采样中，最大单进程驻留内存为 438.6 MiB；短于采样间隔的峰值未知。

三维 48 个正式冷启动臂合计 1679.8 臂墙钟秒，包含两条 corrector。中途纠正共面圆环为矩形面阵的费用及三网场误差另见三维报告。二维线程与三维作业同时运行；GPU 内核运行时间没有通过设备 profiler 单列。

评分耗时与开发开销分别列在各 cohort 的执行记录、同状态诊断与三维报告。部分早期开发诊断仅记录关键操作耗时，未记录启动到退出的总 CPU，整项研究的 CPU 总额高于二维反演小计。[公开开发诊断](engine_dev/REPORT_ZH.md)
