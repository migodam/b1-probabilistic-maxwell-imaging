# 费用、时钟与计算范围

主反演OS整子进程CPU总17199.390625秒；12路并行队列实际历时1614.144759秒。队列控制器CPU4.25秒另列，不能把线程wall之和称并行历时。所有144次实际花费保留，无失败对象从分母删除。

|臂|OS整子进程CPU秒|576阶物理LU|实际物理复RHS|
|---|---:|---:|---:|
|G-PGP|2838.046875|16003|313860|
|S-PGP|2882.562500|16703|323748|
|O0-PGP|2889.609375|16578|322512|
|G-GP|2866.187500|16490|317958|
|S-GP|2853.328125|15845|310854|
|O0-GP|2869.656250|16768|324600|

物理合计98387 LU、1913532 complex RHS。RHS以aggregate actual_complex_RHS为准，不再把source/adjoint role计数加一次。FFT、Green动作、分母19阶小LU/RHS、SVD、prior及材料内迭代完整留在COST_LEDGER.json和原回执，不能折算成同一次576阶物理LU。O0额外接收伴随RHS已经计入物理aggregate，不双算。

每臂120秒是内部候选资格上限，proposal最多72秒；未用proposal时间给原物理corrector，绝对截止仍120。候选只在完整原物理审计后才有资格。最后操作和退出序列化不可总在时钟点截断；全部OS花费最小84.859375、最大120.84375秒，130次实际退出CPU略超过120。超时完成的候选仍不能进入预算内选择。不能说144次恰好消耗144×120秒。

## 额外研究与失败成本

|事项|已测CPU秒|物理LU / 复RHS|范围与回执|
|---|---:|---:|---|
|本地必要接口接线及import preflight，3次记录|4.550838|8 / 192|audit/VALIDATION_COST_TOTALS，含两个接线版本|
|远端import-only preflight|0.90625|0 / 0|audit/REMOTE_DEPLOY_AND_PREFLIGHT|
|远端部署的本机SSH控制CPU|0.009823|0 / 0|wall15.211510，不与远端计时相加|
|主队列控制器|4.25|0 / 0|queue receipt；不再加各snapshot|
|未授权预测封存|1.046875|0 / 0|PREDICTION_SEAL.receipt；之后父审核另存批准清单|
|正式离线评分与主图|87.390625|0 / 0|analysis/observable19_v1/RECEIPT，wall137.503906|
|fixed_state_v1（保留21条G-PGP输出缺陷）|77.96875|126 / 2628|实际整进程CPU；不能删掉或归零|
|fixed_state_v2（修复后全126条完成）|75.65625|147 / 2880|无GT、无新basis SVD；不进主候选池|
|旧源码/指标机械抽取，两次|0.104770|0 / 0|audit/history_evidence，未重算科学结果|
|实际结果机械汇总|0.112428|0 / 0|audit/result_summary|
|同状态图表汇总|2.260284|0 / 0|analysis/fixed_state_summary，wall11.608300|
|完整raw归档|23.15625|0 / 0|原字节885133907，ZIP787518418，wall76.364945|
|raw上传本机控制|见audit/RAW_RELEASE_LAUNCHER.json|0 / 0|wall331.301436；归档数据未经本机中转|

全物理新研究成本：主反演98387 LU/1913532 RHS，加接线8/192及两版诊断273/5508，总98668 LU/1919232 RHS；这些是相同576阶模型的实际聚合计数，不将19阶小系统混入。编辑、理论思考、交互式检查及子agent生成代码的硬件CPU没有完整仪表记录，写unknown，不伪造0或美元费用。

## setup、观测与预算归一化

各臂setup都已在整进程CPU内，不另加一次。G-PGP沿原prepare_full还包含其原projected-basis SVD；其余臂不作该未使用的SVD。O0复用已经付费的初始化LU，新增nr伴随RHS和O0 SVD。原R13名义setup字段与本补充actual_factory_setup分别保留，不能用原名义常数盖掉实际。

原同case六臂同一逻辑核心顺序，统一单BLAS与优先级。核列表沿既有8P+4E映射，但每条队列core_type为空时不凭逻辑ID重新推断每条记录的P/E性能。没有改变系统功率、CPU频率、驱动、CUDA后端或其他任务。

完整paid_prefix与common_CPU_prefix、common_physics_pool均由已经支付的候选池回放原Phi选择。共同LU/RHS预算按同对象六臂终态的逐分量最小值计算；实际所有开销仍收费，不退还额外计算，也不宣称运行了新的提前停止版算法。旧R14 CIE/Dual/A+/11/G19/DR是独立标注的冻结参考，不混作本次六臂重放成本。

归档/公开工程后续的计划、实际压缩、上传、匿名回读和页面检查都有新回执；它们属于publication，不计成online反演提速。所有失败attempt和原文件保留。最终公开入口为evidence.html的已核验附件及各PACKAGE_RECEIPT；没有新增科学运行。
