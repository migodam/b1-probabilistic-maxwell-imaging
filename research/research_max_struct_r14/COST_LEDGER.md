# 本轮逆解阶段成本账本

已按真实终态回执核对 **576个成功solver child**，实际OS CPU合计 **64294.503713秒**。其中Windows为63043.828125秒，Apple为1250.675588秒。另有shared旧失败attempt的16个child，实际OS CPU为0.514267秒；成功与失败solver合计64295.017980秒。该计数表示进程/费用闭合，不表示收敛、恢复或科学验收。

|Campaign|设备|成功child|实际child OS CPU（秒）|该队列wall（秒）|物理LU|物理复RHS|
|---|---|---:|---:|---:|---:|---:|
|dualmetric_pilot12_v1|Windows|48|5729.765625|997.009203|43280|813768|
|dualmetric_remaining36_v1|Windows|144|17267.906250|2508.409878|127472|2382942|
|self_only16_v1|Windows|48|5804.234375|816.870982|44195|877932|
|stage_control16_v1|Windows|32|3872.390625|548.704795|27063|539394|
|independent24_v1|Windows|264|30369.531250|2885.279222|207339|3892614|
|shared_readout4_v1|Apple_local|16|952.745000|486.547589|43330|886170|
|vector4_v1|Apple_local|24|297.930588|156.623597|14331|1714800|

每个wall只描述对应队列，**不相加为整体并行wall**，也不作跨设备绝对速度结论。Windows五队列已逐child核对GetProcessTimes与队列合计。Apple shared/vector只有批次OS总量；按method的OS CPU未独立采集，因此在细表保留null，不按body时间比例分摊。各method的实测body计时另列，且不与外层OS CPU重复相加。

## 物理求解与辅助动作分别计数

|设备|物理矩阵维数|成功child|原始物理LU|原始物理复RHS|
|---|---:|---:|---:|---:|
|Apple_local|576×576|16|43330|886170|
|Apple_local|648×648|24|14331|1714800|
|Windows|576×576|536|449349|8506650|

scalar矩阵为576×576；vector矩阵为648×648，对应216个共享实材料cell、每cell三个场分量。两种物理与尺寸不混为同等LU工作量。

|DR任务数|固定背景H Cholesky|H setup RHS|H iteration RHS|
|---:|---:|---:|---:|
|24|24.0|960.0|673092.0|

DR的固定背景H分解与迭代RHS单独保存，既不加到材料相关物理LU，也不加到物理Maxwell RHS；其CPU已包含在实际child计时。material-force-decomposition的56分支含152064个block复RHS，属于缓存4cell辅助动作，**不是物理solver RHS**；其LAUNCH外部CPU仅计一次，内层receipt与plot body为嵌套计时。

## 费用与来源边界

当前RECEIPT/ACTUAL_CONFIG和终态队列足以核对576个solver费用。新independent264的完整NPZ原档已经按全部任务分卷上传到本轮GitHub Release；本地主要保留compact review及12个新机制checkpoint，完整下载范围以DOWNLOADS和公开readback为准。本地未重复下载所有NPZ不是科学缺失。shared旧失败root与成功attempt_02分开保留；vector使用attempt_01。

exposed48_v1的分析输出重复pair错误已保留其20.179077秒分析CPU，不计为物理失败；v1/v2数值比较由主线程验收。初始import、输出格式失败及未测传输费用的CPU仍保留null，不补0。controller、generation、validation、scoring与archive/transport在详细账本单列；不把solver、wrapper、body与嵌套事件重复相加。最终上传/出版费用由主线程后续追加，本表不等待该阶段、也不把目前已知费用称为最终全项目总成本。

证据文件：[cost_ledger_solver_complete_v1/REPORT.md](audit/cost_ledger_solver_complete_v1/REPORT.md) · [按campaign与method分组](audit/cost_ledger_solver_complete_v1/SOLVER_METHOD_TOTALS.csv) · [逐child费用](audit/cost_ledger_solver_complete_v1/SOLVER_TASK_COSTS.csv) · [DR背景费用](audit/cost_ledger_solver_complete_v1/DR_BACKGROUND_COSTS.csv) · [576进程核对](audit/cost_ledger_solver_complete_v1/SOLVER_COMPLETION_VERIFICATION.json)。本次只读取成本/状态metadata，未读GT数组、未运行新physics或改变算法。

## 出版、传输与失败操作（不并入solver费用）

原始档案的无损打包、离线评分、图表、上传和访问核对均为额外工作，不能把这些CPU合成新的反演吞吐量。可实测的远端末三组打包102.031250 CPU秒，新264打包70.640625 CPU秒，四个远端上传批次合计46.343750 CPU秒；Apple fullstates上传4.953722 child CPU秒（44.717459 wall秒）。完整原始回执在audit/publication。新12个Pro checkpoint原字节导出另0.031250 CPU秒，没有新物理计算。

归档计划曾误对已知大数值JSON重复作字符串扫描，停止前已知CPU至少473.875秒，未测尾部保留null；未进行任何新Maxwell步骤。另有partial-clone历史大小查询被停止，CPU/下载量未独立采集而保留null。两次操作与失败记录均保留，不从费用中隐去。最终本地source/diagnostics/Pro/离线真值打包、网页构建与发布核对费用继续记录在audit/publication的实际receipt；其中嵌套计时不重复求和为虚假的全项目总价。API/模型的真实货币账单在本环境不可得，标null。
