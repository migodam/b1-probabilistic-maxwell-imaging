# P0-R 在线与Oracle：相同距离，不同信息条件

状态：**WORK_IN_PROGRESS / 保存结果说明**。本文件只读现有结果、审计与源码文字；新增 Maxwell、LU、RHS 均为0。科学判决与后续实验由主线程维护，本页不修改门或账本。

## Oracle实际使用了什么

[oracle_fit.py](src/wr/oracle_fit.py)使用原真值子像素投影决定非零支持，并以其95%作为初值。两个wide对象达到d=0.967572217/0.987976116；投影梯度0.467013702/0.305986096，均非梯度收敛。结果身份`ORACLE_DIAGNOSTIC`不能重命名为measurement-only online。

这提供的是带真值信息的局部数据拟合证据，仍未证明原最小特征语义下完整几何合法；更不能由两个端点构造后验模态权重、概率或分支覆盖率。

## 在线实际使用了什么

[online_inference.py](src/wr/online_inference.py)只读测量均值、公开天线/频率/网格，使用背景Born初始化和测量分数选取2×2 tiles。其支持是受限集合，不读取真值mask。接受迭代保持该集合不变，未实现外层支持交换。下面全部数值来自各自完成结果：

|seed|band|method|完整数据d|投影梯度|iterations|停止/状态|原始结果|
|---|---|---|---|---|---|---|---|
|2026102201|narrow|TR|12.3099031|1|14|max_iter reached|[JSON](results/R3/online_2026102201_narrow_TR/RESULTS.json)|
|2026102201|wide|FREQUENCY_CONTINUATION_TR|14.0792428|见stage04|8+6+6+10|ALL_FOUR_NUMERICAL_STAGES_COMPLETED|[JSON](results/R3/online_2026102201_wide_FREQUENCY_CONTINUATION/RESULTS.json)|
|2026102201|wide|LBFGSB|45.1492658|1|0|Requires passive 0<=f*Im(epsilon)<=0.15 (not gradient convergence)|[JSON](results/R3/online_2026102201_wide_LBFGSB/RESULTS.json)|
|2026102201|wide|LBFGSB_REPAIRED|17.6030426|1|17|max_calls reached (not gradient convergence)|[JSON](results/R3/online_2026102201_wide_LBFGSB_REPAIRED/RESULTS.json)|
|2026102201|wide|TR|19.6273841|0.995385742|14|max_calls reached|[JSON](results/R3/online_2026102201_wide_TR/RESULTS.json)|
|2026102202|narrow|TR|6.89265014|0.988820158|14|max_iter reached|[JSON](results/R3/online_2026102202_narrow_TR/RESULTS.json)|
|2026102202|wide|TR|24.9986688|0.999542236|13|max_calls reached|[JSON](results/R3/online_2026102202_wide_TR/RESULTS.json)|

目前保存的在线完整数据距离最低为6.892650137（seed2202 narrow TR），仍大于1。不同band的数据频率不同，两个narrow结果不能替代对应wide结果；同seed的多方法也不是更多独立真值样本。修复后LBFGSB与续接分别留下17.60304259/14.07924278的wide端点，均没有达到原均值距离1。

## 不能从差异直接得出的结论

Oracle与在线不仅初值不同，支持集合也不同；旧TR还受边界roundtrip异常影响。所以这些结果不能分离地估计“只换优化器”“只换初值”或“只换支持”的因果贡献，更不是等成本算法排名。可分离的历史数值与异常记录见[代码审计](P0_R_CODE_AUDIT.md)。

全域/全类可辨识性也不由当前在线失败决定。R封账为12000 RHS；它只限制本轮搜索，不证明不存在其他合格材料或其他测量选择支持。

## B2与采样身份继续分开

`results/B2/`保存的是`ORACLE_CONSTRUCTED_PHYSICAL_WITNESS`身份的连续世界分支候选：参考支持来自已知物理世界。它可以作为独立物理存在性探索；即便某个候选将来满足校准阈值，也不会追认P0-R在线成功。页面只呈现其实际d_cal、d_T、成本和待独立验证status，不在这里作见证资格判决。

C/P1A/P1B没有真实结果时保持NOT_RUN；优化轨迹不是采样链，不产生ESS、后验概率、模态权重或覆盖率。campaign页面持续WORK_IN_PROGRESS。

[数据页面](index.html) · [机器可读报告](RESULTS_REPORT.json) · [R封账凭据](R_SEALED.json)


## R阶段封账时点说明

上述保存说明中的WORK_IN_PROGRESS与“后续未运行”描述R封账时点，不是整轮当前状态。R原始数字与阶段Gate保持原样；最终已完成的B2/C/受控P1-A/P1-B及其限制见[整轮裁决](SCIENTIFIC_DECISION.md)和[阶段汇总](GATE_SUMMARY.json)。本补记没有重开R预算或将Online失败改成PASS。
