# 同窗口计时复测：独立只读审查

本次复测属于**评分后计时控制**，不是新模型、新风险测试、重新选择方法或重开原global8。每个noise condition固定metadata顺序first16，共32个已暴露parent；三个物理算法各独立初始化，按parent/condition循环顺序运行，共96路径。额外实际费用268,544 RHS，低于预留482,560，全部保留在D账本。

本审查不运行物理、不读GT。逐个重新比较96份saved context与原封存q0/chi0/prediction，最大绝对差均为0。逐个closed receipt实际forward+derivative RHS相加=268,544。两份20,000×16 parent-index bootstrap与saved means、三种节时率/95%区间均重放一致。详见 `TIMING_DRIFT_QA.json`。

|noise条件|算法比较|平均anchor节时|描述性双侧95%区间|
|---|---|---:|---:|
|DLOW ×.5|Low24/P64 vs Full96/P64|27.65%|[16.82%,36.13%]|
|DHIGH ×2|Low24/P64 vs Full96/P64|32.77%|[23.39%,39.58%]|
|DLOW ×.5|Direct prior-P128 vs Full96/P64|29.31%|[16.64%,42.83%]|
|DHIGH ×2|Direct prior-P128 vs Full96/P64|36.34%|[24.79%,46.82%]|
|DLOW ×.5|Direct prior-P128 vs Low24/P64|2.29%|[−15.41%,20.68%]|
|DHIGH ×2|Direct prior-P128 vs Low24/P64|5.30%|[−16.51%,26.08%]|

这些是**anchor wall**的同窗口比值：包含对应context持久化/cleanup，不是另一次完整NN系统benchmark。抽样单位是16个固定parent，不是96个独立对象；没有重新生成truth或重新评分。这里的95%描述区间不与原global8的99.375%单侧界混称，不宣称合并family-wise覆盖。

原主实验两个条件cost_GO的PASS保留为**原执行窗口**的冻结结论；目前Pareto图已明确标记original execution window。新窗口DLOW的节时下界仅16.82%，所以“跨窗口至少节省20%已获确认”不成立。DHIGH描述性区间虽整体高于20%，也不能把此小样本控制升级成原global8规格的新确认。Direct与Low的区间两档均跨0，不能据点值2.29%/5.30%宣布Direct是计时赢家。

端点完全一致使该差异可归于执行窗口/硬件负载/IO时序等计时环境，不能归于新的重建结果。这个解释不唯一定位到某种硬件原因，且不能仅由本复测估计所有未来窗口的性能分布。最终GO范围由root综合风险、成本与这些不确定性裁决。

源：`track_c/common_window_timing/PROTOCOL_PRE_TIMING.json`、`RESULTS.json`、两份 `*_TIMING_DRAWS.npz`、96份TIMING job receipts和原/重放contexts。此审查新增RHS=0。
