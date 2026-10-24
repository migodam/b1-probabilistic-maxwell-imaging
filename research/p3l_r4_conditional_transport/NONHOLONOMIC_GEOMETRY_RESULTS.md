# 动态弱几何：有限材料可达性与非完整路径的区别

真实prior material metric与Jacobian按rank2/4/8等原配置计算，包含decoder核、谱旋转、缝跨越和prior代价。训练、GT离线评分与candidate archive分开。

四段Euler不是纯commutator：二阶含 [X,Y]−DX·X−DY·Y。发现后保留原错误实现证据，改用Heun并做固定小步检查，不把错误Euler结果当理论反证。

原任务尺度尝试及离线2%距离适用门：Euler4/288匹配、全部支持0；Heun9/288匹配、全部支持2；直接方向288/288匹配、全部支持94（原首轮95）；frozen144/144、支持59；dynamic144/144、支持83。未完成和失败保留在原计划分母。原teacher1%门不改成2%；新表只回答哪些比较真正处于同一材料距离。

因此未匹配Heun/direct平均GT误差不能用来排名。直接方向更便宜且更容易保持距离；当前四段非完整路径没有成本或材料统计增量。

八份开发观测的共同candidate pool按在线Phi/UID排序、完整χRMS≥.02严格去重，最终每份观测有15/12/7/16/28/16/15/10个（含anchor），共111个新的prior-legal、有限噪声竞争材料。48个仅略低于.02的pair仍被拒绝，不用epsilon放宽。其支持含noise-ball、ΔPhi≤2、prior_delta≤2、signed progress及backward residual条件。

这些是不同观测下的candidate集合，不能合并称119个同一世界解，也不能称posterior modes/质量。它们未要求两世界白化预测距离≤1，未构造exact-fiber证书。weights/probabilities均null。

另外实际运行四种冻结满秩pCN核、8对象×1024steps、262,144RHS。Q正定且Q<I，原Gaussian可逆；完整MH保留拒绝重复。有限材料函数ESS约2.5–4.1，均UNQUALIFIED；bracket-enriched未建立独立混合收益。静态几何在监督评分上的正面结果不能用来追认动态路径成功。
