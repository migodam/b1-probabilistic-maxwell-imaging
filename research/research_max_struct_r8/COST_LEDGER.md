# R8实际计算与发布成本账本

成本状态：`FINAL_RECORDED_COSTS_WITH_EXPLICIT_GAPS`。回执截点：2026-11-06T00:22:43.721569+00:00。原主750/750、独立anchor126/126、Fresnel25/25实际outer已登记；32个anchor开发任务、48个原验证任务与21个历史开发任务分开记账。

**当前可审计的已记录CPU下界为 118,512.075778 秒，即 32.920021 个累计进程CPU小时。** 这不是精确完整总成本，也不是项目历时。已有11个加和行的CPU未记录；另外保留未仪表化的调度、失败启动、复核和共享准备缺口。NULL不得解释成0。

其中明确标记为工程中断/OS采样下界的已知CPU合计 1,070.930000 秒，其余有数值的记录为 117,441.145778 秒；这些数字并不补齐未知部分。不同机器/物理模型的进程CPU可作为已支付工作量累计，不能据此给出硬件速度比、能耗或材料恢复优势。

已记录的逐任务wall合计为 132,961.525193 秒，包含并发任务和中断下界，**不能当项目经过时间**；更不能把CPU合计称wall历时。真实项目历时没有由本账本推断。金额、计费单价、货币费用和能耗均未记录，不报告估算金额或零费用。

## 不重叠费用分类

|费用阶段|加和行数|已记录CPU秒|CPU未知行|
|---|---:|---:|---:|
|合成观测/训练数据及公开特征|14|70.739036|0|
|实际网络训练（含epoch验证选择）|3|1,406.920852|0|
|实际冻结初始化/装载/推理|6|0.453587|0|
|已明确分出的setup与调度控制|983|18.917280|2|
|在线反演主体（已扣单列setup）|1114|115,415.891857|0|
|独立诊断/验证/制图QA及部署计时实验|19|68.315690|1|
|独立离线评分/分析/可视化|8|71.522019|0|
|失败及工程中断尝试|17|1,071.523750|3|
|远程结果压缩/传输与已执行发布准备|16|387.791707|5|

每个warm outer先扣其明确setup，再把setup另列一次；auxiliary/current/corrector/physical validation、训练history累计timer、队列child总计不再加。论文分支使用原会计器已扣复用的executed-cost行；向量3D只计成功臂明细，不能再加包含这些臂的进程总计。LU/RHS/FFT/SVD及矩阵维度保留原物理模型口径，不给出跨标量/矢量或不同维度的“等价LU”。

## 确认、开发与硬件分开

|研究队列|硬件|实际任务|outer CPU秒（setup已合回）|任务wall和秒|
|---|---|---:|---:|---:|
|原主确认750|Mac|405|24,750.170775|31,135.550276|
|原主确认750|XINAN|345|53,928.734375|55,405.293371|
|独立anchor确认126|Mac|66|3,928.733645|4,508.382259|
|独立anchor确认126|XINAN|60|11,138.812500|11,407.323743|
|原验证开发48|Mac|48|2,985.544517|3,693.846923|
|anchor试点8|Mac|8|529.186418|721.789962|
|anchor补充验证24|Mac|24|1,499.351317|1,883.592460|
|暴露历史开发21|Mac|21|1,274.461107|1,583.460057|

最终实际主任务分配为Mac405/XINAN345，anchor为Mac66/XINAN60。远程主队列465与anchor队列81的dispatch记录包含141个`RELOCATED_BEFORE_EXECUTION`确认；该确认不是物理求解。迁移的7对象全部方法在Mac执行，逐份实际receipt平台优先；迁移控制CPU未知，仅记录已测控制wall。

## 实际研究花费与部署情境

上面的总账只计实际发生的独立过程。公开特征生成、原网络训练、实际共享初始化批次已分别收费。每个方法的warm/cold deployment值是同一求解硬件的部署情境：给对应solver加独立测得的特征/推理或冷装载前缀；**不能把这些重复的部署前缀逐行再加进实际研究总账**。XINAN情境与实际曾在Mac生成冻结初值的工作流分开，训练摊销不是新发生的CPU。下表按队列、cohort、硬件和方法列每任务中位数，仅用于费用口径展示，不是新的科学胜负判据。

|队列|cohort|硬件|方法|任务数|solver中位CPU秒|warm部署中位CPU秒|cold部署中位CPU秒|
|---|---|---|---|---:|---:|---:|---:|
|主750 NN|CROSS_GRID|Mac|a_plus|6|62.426440|62.427876|63.590434|
|主750 NN|CROSS_GRID|Mac|a_plus_sourcewise|6|62.698642|62.700002|63.862583|
|主750 NN|CROSS_GRID|Mac|cie_dual|6|53.519807|53.521244|54.683801|
|主750 NN|CROSS_GRID|Mac|physical_lbfgs|6|28.755606|28.757256|29.919536|
|主750 NN|HELDOUT_FAMILY|Mac|a_plus|6|75.995024|75.996529|77.158856|
|主750 NN|HELDOUT_FAMILY|Mac|a_plus_sourcewise|6|71.533704|71.535209|72.697536|
|主750 NN|HELDOUT_FAMILY|Mac|cie_dual|6|63.174207|63.175835|64.338115|
|主750 NN|HELDOUT_FAMILY|Mac|physical_lbfgs|6|37.084534|37.086161|38.248441|
|主750 NN|HELDOUT_FAMILY|XINAN|a_plus|6|158.234375|158.244263|159.994263|
|主750 NN|HELDOUT_FAMILY|XINAN|a_plus_sourcewise|6|162.898438|162.908325|164.658325|
|主750 NN|HELDOUT_FAMILY|XINAN|cie_dual|6|138.929688|138.939575|140.689575|
|主750 NN|HELDOUT_FAMILY|XINAN|physical_lbfgs|6|100.656250|100.666260|102.408447|
|主750 NN|IN_DISTRIBUTION|Mac|a_plus|10|63.627690|63.629299|64.792388|
|主750 NN|IN_DISTRIBUTION|Mac|a_plus_sourcewise|10|70.069579|70.071224|71.234313|
|主750 NN|IN_DISTRIBUTION|Mac|cie_dual|10|64.514187|64.515611|65.678699|
|主750 NN|IN_DISTRIBUTION|Mac|physical_lbfgs|10|30.638999|30.640766|31.803601|
|主750 NN|IN_DISTRIBUTION|XINAN|a_plus|14|134.851562|134.860474|136.602661|
|主750 NN|IN_DISTRIBUTION|XINAN|a_plus_sourcewise|14|147.023438|147.033813|148.783813|
|主750 NN|IN_DISTRIBUTION|XINAN|cie_dual|14|142.359375|142.370728|144.120728|
|主750 NN|IN_DISTRIBUTION|XINAN|physical_lbfgs|14|62.921875|62.931885|64.674072|
|独立anchor126|CROSS_GRID|Mac|a_plus|6|53.610891|53.611944|54.774271|
|独立anchor126|CROSS_GRID|Mac|a_plus_sourcewise|6|51.425811|51.426670|52.588997|
|独立anchor126|CROSS_GRID|Mac|cie_dual|6|50.524252|50.525305|51.687632|
|独立anchor126|HELDOUT_FAMILY|Mac|a_plus|6|64.545106|64.546734|65.709014|
|独立anchor126|HELDOUT_FAMILY|Mac|a_plus_sourcewise|6|59.683644|59.685508|60.847812|
|独立anchor126|HELDOUT_FAMILY|Mac|cie_dual|6|65.474569|65.476073|66.638400|
|独立anchor126|HELDOUT_FAMILY|XINAN|a_plus|6|216.414062|216.424194|218.166382|
|独立anchor126|HELDOUT_FAMILY|XINAN|a_plus_sourcewise|6|214.601562|214.611938|216.354126|
|独立anchor126|HELDOUT_FAMILY|XINAN|cie_dual|6|193.921875|193.931885|195.674072|
|独立anchor126|IN_DISTRIBUTION|Mac|a_plus|10|63.726514|63.728291|64.891126|
|独立anchor126|IN_DISTRIBUTION|Mac|a_plus_sourcewise|10|53.378977|53.380586|54.543674|
|独立anchor126|IN_DISTRIBUTION|Mac|cie_dual|10|57.271096|57.272575|58.435663|
|独立anchor126|IN_DISTRIBUTION|XINAN|a_plus|14|157.070312|157.080444|158.822632|
|独立anchor126|IN_DISTRIBUTION|XINAN|a_plus_sourcewise|14|170.031250|170.042603|171.792603|
|独立anchor126|IN_DISTRIBUTION|XINAN|cie_dual|14|165.054688|165.065674|166.807861|

## 其他物理分支与失败

Fresnel25个成功outer累计 10,603.740701 CPU秒，包含原20个比较/NIE、两个6+8GHz joint/control、一个air0开发诊断、两个corrected-air full-budget开发实验；原结果全部保留。添加6GHz观测属于额外数据，不据此声称A+算法收益。独立制图/QA费用另见阶段表，嵌套物理校验没有再次相加。

|失败/中断类别|条数|已记录CPU秒|CPU未知条|
|---|---:|---:|---:|
|fresnel_engineering|4|405.520000|1|
|paper_benchmarks|3|320.950000|1|
|vector3d_engineering|1|344.460000|0|
|vector3d_serialization|1|0.000000|1|
|warm_remote_dependency_failure|8|0.593750|0|

8³未缓存工程尝试的344.46秒OS CPU下界已经包含完成LSIE的227.39秒，不重复相加。早期PR48 numpy-bool序列化失败的首次A+数组被同合同重跑覆盖，初次快照NOT_RETAINED、该失败CPU未测；不能伪造完整历史或bitwise一致性。Fresnel首个layout失败CPU与paper首次64grid失败CPU未记，三个Fresnel中断CPU合计下界405.52秒；paper两个中断下界320.95秒。8个threadpoolctl缺失的真实失败进程单计0.593750 CPU秒，失败队列wrapper不再加。

新direct-budget分析首次启动因语法错误在解析阶段失败，没有读取输入或调用物理；该启动CPU未测。成功分析费用1.555670秒已单列，之后只读verification CPU未独立计时，不能回填为0。旧PILOT绘图及部分原调度/资格检查也有未测费用。完整缺口逐条列于COST_LEDGER_FINAL.json的gaps与failures，不用已完成的科学状态抹除这些缺口。

## 远程提取、复制来源与发布准备

5次实际incremental fetch（主3次、anchor2次）远程压缩共 263.515625 CPU秒、593.720487 任务wall秒；控制总wall剔除嵌套压缩wall后的SSH/传输/解压余项为 155.183650 秒，其本地CPU未计量。最新FETCH.json和audit/remote_imports下的复制审计只作来源镜像，不是新工作，不重复计费。

最终compact/ZIP执行费用已由PACKAGING_EXECUTION_RECEIPT独立登记：106.938462 CPU秒、131.197553 wall秒；没有新增物理调用。

当前含scored-truth的权威计划：Git候选 **87,628,690 bytes**，Release原成员 **16,334,522,276 bytes / 52 包**。这是stat得到的候选/非压缩源字节，尚非实际ZIP尺寸；已实际创建1126份独立compact；最终Git字节仍会随root文档/费用元数据更新。实际ZIP共5,308,028,340bytes、52包，全部路径见delegated/cost_package/ARCHIVE_PATHS.txt。包源上限480,000,000 bytes，实际ZIP强制≤500,000,000 bytes。

自有FIGURES PDF共14份，最大1,473,547 bytes；没有超过90MB阈值的自有PDF，当前无超大PDF警告。

完整原prediction/current/j/Q、history/budget/checkpoints、训练数据和大physics算子均保留Release原相对路径，compact不替代原证据。两套750/126封存评分门控已通过，计划中自有GT单独scored_truth分包；本writer和packager的plan mode都不打开GT数组。GT实际发布后42对象是公开开发集，不能再称未来盲测。

永远排除inputs原Chat、data_local第三方原始Y/测量、*.local.*全文下载、第三方PDF、SSH身份/credentials/.env、本地vendor/cache和临时传输ZIP。自有FIGURES PDF允许发布。恢复与root实际执行命令见delegated/cost_package/FINAL_BUILD_AND_RELEASE.md；真实上传和直接网页链接由root另行完成，不把计划当已发布。

数据来源：delegated/cost_package/COST_LEDGER_FINAL.json（实际outer回执账）、ARCHIVE_INDEX.json（权威stat成员计划）、RESULTS/warm_scoring/SCORES.json与delegated/warm_anchor_analysis/CONFIRMATION_SCORES.json（部署口径标量）。计费和部署表均未读取原始状态数组或测试GT。
