# P0-R FWI轨迹：数据拟合、驻点与局部非凸分开

状态：**WORK_IN_PROGRESS / 保存结果说明**。本文件只读现有结果、审计与源码文字；新增 Maxwell、LU、RHS 均为0。科学判决与后续实验由主线程维护，本页不修改门或账本。

本页使用`RESULTS.json`与完成评价历史，不通过“进程结束”推断算法收敛。数据距离d小于1、投影梯度小、完整几何合法、独立端点数值资格是不同条件。

## 1. Oracle与在线完成结果

|seed|band|d_mean|projected_gradient|iterations|success|stop|
|---|---|---|---|---|---|---|
|2026102201|wide|0.967572217|0.467013702|4|False|objective_floor_reached (not gradient convergence)|
|2026102202|wide|0.987976116|0.305986096|7|False|objective_floor_reached (not gradient convergence)|

|seed|band|method|完整数据d|投影梯度|iterations|停止/状态|原始结果|
|---|---|---|---|---|---|---|---|
|2026102201|narrow|TR|12.3099031|1|14|max_iter reached|[JSON](results/R3/online_2026102201_narrow_TR/RESULTS.json)|
|2026102201|wide|FREQUENCY_CONTINUATION_TR|14.0792428|见stage04|8+6+6+10|ALL_FOUR_NUMERICAL_STAGES_COMPLETED|[JSON](results/R3/online_2026102201_wide_FREQUENCY_CONTINUATION/RESULTS.json)|
|2026102201|wide|LBFGSB|45.1492658|1|0|Requires passive 0<=f*Im(epsilon)<=0.15 (not gradient convergence)|[JSON](results/R3/online_2026102201_wide_LBFGSB/RESULTS.json)|
|2026102201|wide|LBFGSB_REPAIRED|17.6030426|1|17|max_calls reached (not gradient convergence)|[JSON](results/R3/online_2026102201_wide_LBFGSB_REPAIRED/RESULTS.json)|
|2026102201|wide|TR|19.6273841|0.995385742|14|max_calls reached|[JSON](results/R3/online_2026102201_wide_TR/RESULTS.json)|
|2026102202|narrow|TR|6.89265014|0.988820158|14|max_iter reached|[JSON](results/R3/online_2026102202_narrow_TR/RESULTS.json)|
|2026102202|wide|TR|24.9986688|0.999542236|13|max_calls reached|[JSON](results/R3/online_2026102202_wide_TR/RESULTS.json)|

在线表中的首个LBFGSB保留边界编码失败；REPAIRED是独立修复后轨迹，不覆盖旧行。TR与Oracle停止信息也原样保留。原始结果均无足以宣称在线全局收敛的终止证据；Oracle的objective floor更明确写着not gradient convergence。

## 2. frequency continuation不是把不同子集当作一条loss

[完整续接结果](results/R3/online_2026102201_wide_FREQUENCY_CONTINUATION/RESULTS.json)保存四阶段，当前状态`ALL_FOUR_NUMERICAL_STAGES_COMPLETED`，最后全五频d=14.07924278。

|stage|频率个数|λ|迭代数|子集d|投影梯度|stop|
|---|---|---|---|---|---|---|
|1|1|0.002|8|3.07927183|0.973271224|max_iter reached|
|2|2|0.002|6|5.51640243|0.977974631|max_iter reached|
|3|3|0.002|6|7.94543943|1|max_iter reached|
|4|5|0.0|10|14.0792428|1|max_iter reached|

频率依次为[.8]、[.8,.9]、[.8,.9,1]、[.8,.9,1,1.1,1.2]；前三阶段λ=.002，最后λ=0。各阶段残差属于不同数据维度，不能用3.079→14.079宣称同一目标变差。每阶段都因max_iter结束，支持始终为同一测量选择tile集合，没有外层支持交换。工具将四条`FWI_HISTORY.json`分别画在子图中。

## 3. 控制：已知形状也可以停在坏的局部区域

[控制原始结果](results/R3/controls/RESULTS.json)与[src/wr/controlled_fwi.py](src/wr/controlled_fwi.py)明确这不是未知墙体实验：几何已知，只求两个全局材料幅度，频率f=1，初始归一化参数[.1,.1]。

|kind|method|d_mean|投影梯度|归一化恢复参数|success|stop|
|---|---|---|---|---|---|---|
|vacuum|TR|6.03893257e-09|9.00855486e-11|9.00855486e-11, 1.73472348e-18|True|projected_gradient_converged|
|vacuum|LBFGSB|0|0|0, 0|True|projected_gradient_converged|
|block|TR|10.3082656|0.000161134741|0.0944350003, 1.97215226e-31|False|max_calls reached|
|block|LBFGSB|10.3082656|1.95510727e-06|0.0944349979, 0|False|max_calls reached (not gradient convergence)|

vacuum的两方法达到其记录的投影梯度停止条件。block的数据由同一离散模型在归一化参数[.375,.2]生成；两算法却保留约[.094435,0]，同一d≈10.308265566，投影梯度已较小，但都因call cap结束、success=false。这个受控现象与局部非凸坏盆地相符；它提醒我们小梯度不能取代数据拟合，也不能拿“同两算法停到同处”作为全局最优证明。这里没有做全局景观扫描或证明该点是严格局部极小值。

[R_SEALED.json](R_SEALED.json)列disk、smooth_block未完成。disk仅有部分TR评价历史，没有完成RESULTS，必须标为`UNCOMPLETED_R_SEALED`；smooth_block未产生结果，不能编造对照。

## 4. 可复查轨迹而不是补推成功

[生成残差表](tables/saved_results/completed_evaluations.csv)只读取`evaluation_complete`，call ID包括accepted refresh，不是物理RHS或独立材料数。对零中心先验，d由`sqrt(2*fun−λ||x||²)`恢复；固定在线λ=.002、Oracle/controls λ=0，续接阶段按保存λ。开始但未完成或失败的调用不生成残差点。

本页没有重跑任何优化、Maxwell或端点验证，没有扩充已封账R预算。后续B2为独立分支，不能用它替代在线拟合成功。
