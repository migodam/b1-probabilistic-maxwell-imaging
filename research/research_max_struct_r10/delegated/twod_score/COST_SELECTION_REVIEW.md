# Fresh评分的成本与预算选择审查

**结论：三个预算clock已分开，选择依据为已付成本与原Phi-incumbent，不使用GT选时点。外部GetProcessTimes只有完整child寿命终点，不能生成真实逐步外部CPU曲线。** 本审查仅读源码并作机械修正/语法检查，未运行scorer、未打开fresh truth/score metadata、未增加forward或改变inverse预算。

## Clock与计数来源

1. `r8/src/benchmark_solver.py::PaidStateObserver`保存`self.best['m']`，其best只由原物理Phi严格改善更新。CPU为observer的process_time增量加原runner prefix；物理LU与complex RHS直接来自原Maxwell model计数。保存点为首次、每10次已付记录及auxiliary阶段，最后另存一次；因此不是完整逐步轨迹。材料快照为float32，原Phi仍是对应已付float64候选的值，不能说重新优化了量化图像。
2. `twod_adapter.py`只给预算CPU列加一个R10 preparation常量prefix，并保留原file。`R10_BUDGET_STATES`的第一列是**内部已付process CPU**，不含完整child启动/导入及退出尾部的所有开销。追加prefix不会改变快照之间的CPU增量。
3. `twod_queue.py::external_cpu`在child结束后使用Windows GetProcessTimes取得user+kernel累计CPU。这是**完整child寿命终点**，没有逐步外部采样。scorer只接受明确GetProcessTimes scope；没有child成功创建的launch failure可记已知child CPU=0，queue controller成本另列。

总成本将disjoint child lifetimes与两个queue controller各加一次。R10/R8/aux/corrector内部计时互相嵌套，只作分解视图，不相加。未测得的外部child CPU保留unknown；physical counters缺失也另列，已知和不冒充全量。

## 固定预算点的选择

内部CPU50/100/200/400与physical LU50/100/300/600均取**成本不超过阈值的最后一份已保存Phi-incumbent**。不按最小NRMSE、最高IoU、最佳GT阶段或曲线插值取点。同一cost出现多个保存点时保留最后记录，仍是公开Phi选出的incumbent。若已完成终点的对应clock/counter已不超过预算，可直接使用该已付终点；selection_kind会明确区别终点与保存快照。

物理LU预算只控制原576维Maxwell LU数量。它不等于总CPU，也不保证19/28维小分解、source/adjoint RHS及其它工作都相同；对应RHS和外部CPU同时保留。

保存快照稀疏，故“nearest”严格是nearest **saved** already-paid state，不能声称精确重建预算截止瞬间的在线incumbent。未覆盖预算、缺预算文件或错误schema保留null/原因，不能用之后的GT优点补齐。机械修正后每任务都保留12个预算单元，缺失不会从预算表中悄悄消失；完整168任务预期共2016预算单元。

## 外部CPU只有保守包络，不能伪造逐步观测

设内部快照clock为t_i、最后保存值为t_last，外部完整child CPU为T。若确认是同一child的process CPU计时、计数嵌套关系成立且差非负，则

`t_bound(i)=t_i+(T−t_last)`

将所有不能分配的导入/尾部CPU放在快照之前，给出保守上界。因为内部两点之差保留相同process CPU增量，这种常量平移不会把未知开销“测量”出来。它不是观察到的external prefix，更不是每步外部CPU的估计真值。通过该包络选择的状态可保证已付成本上界不超阈值，但可能过于保守，并不等于该full-child CPU预算下实际最后在线状态。

只有已完成整个child且T≤预算时，其final endpoint拥有准确的完整child成本。其余快照若无法证明包络条件，明确`EXTERNAL_PREFIX_CLOCK_NOT_PROVEN`，没有拿内部时间冒充外部时间。

每一行已加`selected_cost_semantics`区分：准确完整child终点、嵌套内部终点、观察到的内部快照、保守外部上界、原物理LU计数。`total_CPU_equal_time_claim=False`防止下游将不同clock混画为equal total CPU。

当前scorer绘制的cost/performance图仅使用**完整child GetProcessTimes终点**散点；没有把内部CPU预算曲线画成equal-total-CPU曲线。若交付时另画预算曲线，内部轴必须标“recorded internal paid CPU”；外部包络轴必须标“conservative full-child CPU upper-bound envelope”，同时保留coverage/缺失数，不能改称逐步实测external CPU。

## 原物理与候选小LU

`physical_model_cost`的LU为576维、source_RHS为已付forward列数、adjoint_RHS为伴随列数、actual_complex_RHS为原模型累计复RHS。预算文件LU列也只来自该model，不含候选小分母分解。

矩阵候选的`auxiliary_execution_receipt.operation_counts`单列小LU次数、LU维数和、small triangular complex RHS。score行保留这些原始计数及维数，不把19/28维分解混到物理LU中。基于原算法结构，每次候选value/gradient会有一次小LU及正/伴随小系统调用；计费使用实际counter，不用这一结构公式重造计数。先验、FFT、低秩乘法、setup等实际耗时已经进入child CPU，不将局部counter解释成完整FLOP计价。

## 本次机械修正与验证

补丁保存为`COST_SELECTION_MECHANICAL_DIFF.patch`，修改前源码快照为`score_before_cost_review.snapshot.txt`。只改评分程序：

- 为每预算行补充clock/上界语义与非equal-total-CPU标记。
- 验证完整四列schema，确保physical LU轴没有误用其它counter。
- 缺失/错误预算file仍保留12个null预算单元；若已知完成终点在预算内，可以只报告这个有准确完成成本的终点，并保留snapshot缺口标签。
- 明确完整child外部采样只有终点；记录包络条件。
- 单列576维physical LU与候选small-LU维数；launch failure已知零child CPU与未测CPU分开。

语法编译通过，scorer main没有执行。没有读取fresh GT、没有生成fresh评分结果、没有改变producer/runner或inverse预算。下一步仍须root在168任务全终态后明确postseal授权。
