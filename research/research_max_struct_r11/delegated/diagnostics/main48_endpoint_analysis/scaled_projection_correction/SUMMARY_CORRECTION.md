# 缩放投影映射解释纠正（零FW）

必须纠正v1的解释：一般box内，`||m-clip(m-1e6*g)||inf<=.01`与`||m-clip(m-g)||inf<=1e-8`不等价。投影不与gradient缩放交换。v1的未缩放数值和费用全部保留，1e-8计数仅作未缩放端点描述，不能称作原solver停止判据本身。

这里对288个公开封存minPhi材料，重新逐项精确匹配已有m/g，然后只做数组运算；未调用物理、优化器、评分器或新版本probe。g仍是含同prior的未缩放原Phi梯度。

|方法|未缩放mapping<=1e-8|scaled mapping<=.01，在封存点|两判据不一致|scaled满足时最大未缩放mapping|原PGTOL stop字符串|
|---|---:|---:|---:|---:|---:|
|cie_single|27|45|18|1.6124642e-08|45|
|cie_dual|27|44|17|1.7999375e-08|44|
|global_dual|29|45|16|1.4151967e-08|45|
|fullG19|25|41|16|1.7196887e-08|41|
|freezeA|32|43|11|1.3149317e-08|43|
|contractedT|27|43|16|1.7198753e-08|43|

全部94个不一致均是scaled满足而unscaled不满足，反向0例；对应最大未缩放量为1.7999375e-8。所有v1未缩放数值逐项完全复现，证明本次改变的是解释及新增scaled量，不是换端点或改主评分。全部不一致case/method与两种数值见`CRITERION_DISAGREEMENT_CASES.csv`。

## 封存候选不自动等于优化器终止点

上表scaled量在**sealed minimum-originalPhi incumbent**计算。原corrector stop字符串来自optimizer result.message；二者来自不同记录语义，不能仅凭字符串把这份gradient当成terminal gradient。

所有288个封存m也与LAST_PAID_STATE逐项相同，但LAST仅是最后一次已付物理evaluation，不是单独导出的optimizer result.x。`refine_observed`只返回best材料及stop/历史，没有保存独立terminal x。因此`optimizer_terminal_point_identity_established=False`逐行保留；没有把LAST相同升级为终止点身份的证明。已有PGTOL字符串的261行在该封存点scaled量均<=.01；这只是两种记录的一致性观察，不重建历史停止事件。

## 源码与版本依据

原生产代码`research/research_max_struct_r6/src/benchmark_solver.py:171–179,190–198`在observer保存原gradient后给L-BFGS-B返回1e6*g，gtol=.01。Root已核对SciPy官方projgr参考实现：v1.14.1的lbfgsb.f:2834–2841以及v1.15.1的__lbfgsb.c。精确reference URLs及来源声明在`SOURCE_BASIS.json`；它们是已核对的源码参考版本，不是生产环境版本声明。

twod owner仅从已有R11/R10环境回执搜索，没有找到SciPy runtime version/package path。二者在本补充中均为null；只知道D:/python/python.exe命令和threadpoolctl vendor不够推断SciPy版本。没有新remote/local SciPy probe。

## 保留、费用及流程偏差

`v1_tables_preserved/`保留原全部表、脚本与费用。`JOINED_FINAL288_CORRECTED_PROJECTION_LABELS.csv`保留原材料/数据/成本字段，并新增scaled量与更正后的解释；不改主scoredCSV、runner、选择、预算或成功门。旧`FAIL_CN99_PG_CROSSCOUNTS.csv`的PG仍只是unscaled描述，不能作为solver终止集合来解释。新目录只补端点计算，不给信息不可辨识或局部极小标签。

新增完整childCPU0.279910秒，0FW/adj/反演/评分/bootstrap。v1的0.242041秒是以前已计费用，不再次收费；两部分累计已记0.521951秒仅作非加法参考。

在计算前查来源时发生一次过宽文本搜索，意外打印paired8私有生成metadata的一行；没有载入GT数组，但不能称本任务过程全无私有metadata接触。已单独披露，详见`../../private_metadata_search_incident/INCIDENT_AUDIT.md`。本次288数值计算的所有输入均为明确public48端点列表；流程偏差未进入online输入、算法、选根或paired8已冻规则，没有重抽/修正数据。
