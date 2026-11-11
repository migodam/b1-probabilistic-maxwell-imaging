# Pro toy 成熟基线对照结果

40个反演均完成；scaled/clipped各20预测先独立封存，再离线读truth评分。两套都保留，不能按GT挑选主套。所有 auxiliary 满3000步，corrector在原规则下自然停止或达3000步。

每个材料100个独立实未知量、同G/E/S/Y、实box[0,12]、初值.1、无先验。实现与适配边界见 `INVENTORY_AND_REPRODUCE.md`。

|套|case|方法|corrector前材料误差|最终材料误差|前/后data误差|CPU秒|corrector步/评估|
|---|---|---|---:|---:|---:|---:|---|
|scaled|bridge_high|a_plus|0.907720|0.920613|0.736585/0.520419|3.589|128/138|
|scaled|bridge_high|cie_dual|0.912092|0.920613|0.815846/0.520419|2.410|127/135|
|scaled|bridge_high|cie_single|0.887138|0.920613|0.738730/0.520419|2.568|113/118|
|scaled|bridge_high|global_dual|0.910084|0.920613|0.795504/0.520419|2.468|133/148|
|scaled|bridge_high|green_strong_full|0.880527|0.920613|0.630383/0.520419|5.208|122/134|
|scaled|disk_low|a_plus|0.154887|0.307623|0.010009/0.009390|4.048|1325/1361|
|scaled|disk_low|cie_dual|0.167798|0.357254|0.009934/0.009388|3.184|2669/2704|
|scaled|disk_low|cie_single|0.171716|0.343487|0.010648/0.009388|3.241|2708/2738|
|scaled|disk_low|global_dual|0.164649|0.353834|0.009922/0.009388|3.193|3000/3034|
|scaled|disk_low|green_strong_full|0.293489|0.353868|0.010913/0.009388|5.018|1993/2022|
|scaled|ring_high|a_plus|1.000000|1.000000|1.000000/1.000000|2.549|0/1|
|scaled|ring_high|cie_dual|1.000002|1.000000|1.005556/1.000000|2.525|1/2|
|scaled|ring_high|cie_single|1.000002|1.000000|1.005925/1.000000|2.420|1/2|
|scaled|ring_high|global_dual|1.000000|1.000000|1.000000/1.000000|2.321|0/1|
|scaled|ring_high|green_strong_full|1.000000|1.000000|1.000000/1.000000|3.549|0/1|
|scaled|split_shell_high|a_plus|0.417425|0.050836|0.051316/0.009124|4.892|3000/3227|
|scaled|split_shell_high|cie_dual|0.404796|0.051777|0.031966/0.009127|3.368|3000/3219|
|scaled|split_shell_high|cie_single|0.357626|0.049245|0.023281/0.009126|3.382|3000/3216|
|scaled|split_shell_high|global_dual|0.405096|0.051193|0.032157/0.009128|3.308|3000/3207|
|scaled|split_shell_high|green_strong_full|0.537442|0.052576|0.074163/0.009126|7.317|3000/3232|
|clipped|bridge_high|a_plus|1.099197|1.060722|0.868592/0.577892|2.949|156/166|
|clipped|bridge_high|cie_dual|0.987668|0.920613|0.978758/0.520419|2.459|144/166|
|clipped|bridge_high|cie_single|0.987668|0.920613|0.978758/0.520419|2.419|144/166|
|clipped|bridge_high|global_dual|0.987668|0.920613|0.978758/0.520419|2.446|144/166|
|clipped|bridge_high|green_strong_full|1.069361|1.060722|0.924582/0.577892|3.740|139/156|
|clipped|disk_low|a_plus|0.120578|0.340156|0.009847/0.009388|6.386|2438/2462|
|clipped|disk_low|cie_dual|0.125255|0.370364|0.009832/0.009388|3.196|2383/2413|
|clipped|disk_low|cie_single|0.163217|0.338976|0.010638/0.009388|3.237|2801/2834|
|clipped|disk_low|global_dual|0.122174|0.346184|0.009818/0.009388|3.232|2971/2998|
|clipped|disk_low|green_strong_full|0.238886|0.353857|0.010764/0.009388|8.751|2732/2766|
|clipped|ring_high|a_plus|1.000000|1.000000|1.000000/1.000000|2.097|0/1|
|clipped|ring_high|cie_dual|1.000050|1.000000|1.015890/1.000000|2.421|1/2|
|clipped|ring_high|cie_single|1.000014|1.000000|1.011021/1.000000|2.383|1/2|
|clipped|ring_high|global_dual|1.000050|1.000000|1.015890/1.000000|2.361|1/2|
|clipped|ring_high|green_strong_full|1.000000|1.000000|1.000000/1.000000|3.032|0/1|
|clipped|split_shell_high|a_plus|0.383015|0.049743|0.044095/0.009128|4.647|3000/3217|
|clipped|split_shell_high|cie_dual|0.368466|0.051934|0.024548/0.009126|3.221|3000/3227|
|clipped|split_shell_high|cie_single|0.320648|0.051641|0.019238/0.009132|3.050|2423/2606|
|clipped|split_shell_high|global_dual|0.377723|0.052206|0.025546/0.009127|3.133|3000/3237|
|clipped|split_shell_high|green_strong_full|0.545317|0.053114|0.067530/0.009126|6.836|3000/3232|

环的所有方法最终为零材料，原物理data残差1；不是噪声球内错误解。桥在scaled套五臂由不同辅助端点进入相同corrector吸引域，clipped套A+/G19落入更差端点。因此“不同算法最终一致”不能误解为辅助阶段轨迹相同。壳最终误差约5%；低对比disk的data都降到噪声水平，但材料误差仍约31%–37%，再次说明data拟合不等于材料恢复。

共同CPU前缀结果分别保存在 `OFFLINE_SCORES.json` 和 `OFFLINE_SCORES_clipped.json`，规则为该case五臂最小实际CPU，只取此前已付费最低Φ快照。没有使用truth挑选前缀或快照，也不把不同实际开销叫作等成本。

成本下界：scaled 70.557396 CPU秒，28723 physical LU，344178 complex RHS；clipped 71.995712 CPU秒，29904 physical LU，358344 complex RHS。合计 142.553108 CPU秒。

这些计时包含每次run内导入/model/setup/科学计算与部分输出；不包含Python启动、全局NumPy/SciPy导入、launcher/prepare/seal/score/report的完整费用。因此全任务CPU是更大的未知量，不能称上述为完整端到端总成本。计数详细值及停止原因在 `COMPARISON_AND_COSTS.json`。Green分母small LU/解与dense Green RHS在每个receipt另列，不同维度操作不能合并成同一种LU。

科学边界：这证明这两套已声明的成熟算法适配在给定四个toy上的具体结果；不证明连续PDE精度、所有成熟超参数无解、DR/AP的普适优越性或新机制的新颖性。no-prior toy与历史有prior的24²合同不同；已保留两种current continuation敏感性，未扩展搜索救结论。
