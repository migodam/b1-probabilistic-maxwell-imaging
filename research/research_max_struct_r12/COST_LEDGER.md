# R12 费用台账与计时边界

这是按实际目录与计时边界整理的台账，不是已结束整轮的精确总账。独立24的192次运行已封存，整组外部子进程费用已取得；所有准备/审阅/发布成本仍不齐全；**全局精确CPU合计=`null`，全局wall=`null`**。缺失不是0。下列数字来自现有回执，只读汇总，没有启动新实验。

## 计费优先级与去重

- 同一执行有external child lifetime时，采用它计费，保留inner/phase CPU为解释字段，不能再相加。whole-child包含imports、setup、serialization等；某些旧body计时不含启动，明确标为下界。
- `audit/remote_snapshot/RESULTS/basin_v1`是`RESULTS/basin_remote`的运输副本，不另收费。原v1完整端点与其240CPU前缀属于同一次已支付运行，前缀不能再加一次。
- 物理LU、source/adjoint RHS、固定背景H Cholesky、port eigensystem、projection RHS、FFT与接收/Green products分别记录。576-cell scalar-TM、648维矢量DDA及其他维数不能混成一个“等价LU总数”。
- 模型生成、算法、后处理、重试、transport为不同过程；同一clean forward被复用时只计一次。CPU并行累加不等于wall；序列wall之和也不等于并行队列wall。

## 已完成目录：可用计费行

|目录/范围|可用CPU秒|边界/缺失|来源|
|---|---:|---|---|
|Pro原始源码重放|`null`|只有每方法wall、LU/RHS，whole-process CPU/import未测；不能把seconds冒充CPU|[REPLAY_RECEIPT](experiments/pro_reproduction/REPLAY_RECEIPT.json)|
|成熟toy baseline，scaled|70.557396|inner计算下界；包含run内部setup/解算/部分序列化，外部启动与报告不全|[COMPARISON_AND_COSTS](delegated/baseline/COMPARISON_AND_COSTS.json)|
|成熟toy baseline，clipped|71.995712|与scaled不同实际运行，分别计费；同上边界|同上|
|本地basin_v1，3个实际运行|241.994978|三份body CPU回执之和；外部lifetime=`null`；不是远端12条的镜像|[本地目录](RESULTS/basin_v1/)|
|远端basin，12个完整运行|2263.406250|body CPU，不是12×240；导入不在body内；完整archive/corrector部分会超名义窗口|[结果](RESULTS/basin_remote/)、[预算更正](audit/BASIN_BUDGET_ACCOUNTING.md)|
|moderate vector：18原端点|402.267099|per-arm near-exit下界；含模型setup、审计归档/corrector，imports和最终receipt不全；不要再叠加PROJECTION_QUEUE外部子CPU|[vector ledger](delegated/vector/COST_LEDGER.json)|
|moderate vector：18统一long continuation|87.623642|新增独立执行；continuation queue external child另存且重叠，不重复加|同上、[CONTINUATION_QUEUE](delegated/vector/CONTINUATION_QUEUE.json)|
|moderate vector生成|.541862|生成wholeCPU；上述两个solver行不含它；vector目录已知下界490.432603|[GENERATION_RECEIPT](delegated/vector/GENERATION_RECEIPT.json)|
|vector_high：12算法运行|276.859628|setup/iteration/archive/corrector阶段包含在solver内，imports在时钟前；不再加phaseCPU|[high ledger](delegated/vector_high/COST_LEDGER.json)|
|vector_high生成|.503190|与solver独立；6³/8³物理生成计数另存|[generation](delegated/vector_high/GENERATION_RECEIPT.json)|
|vector_high后封存评分|1.311740|0新LU/RHS；receiver projections720columns|[scoring](delegated/vector_high/SCORING_RECEIPT.json)|
|vector_high保存状态分解|.130509|0新Maxwell/0新truth；receiver projections360columns；不能混入算法预算|[decomposition](delegated/vector_high/DECOMPOSITION_RECEIPT.json)|
|冻结NN接口16算法运行|265.601552|moduleimports在execute前，不全含；proposal/corrector嵌套phase不可再加|[neural summary](delegated/neural/COST_SUMMARY.json)|
|冻结NN预处理|.763496|含实际features/推断；没有新训练；初始NumPy import在时钟前|同上|
|冻结NN两次评分/图版|1.535992 + 1.860799|两个实际后处理执行，0新物理LU/RHS，不是算法收益费用|同上|
|external原12算法children|**402.666742**|采用外部lifetime；内层399.330165是重叠下界，不再加|[external ledger](delegated/external/COST_LEDGER.json)、[queue](delegated/external/QUEUE.json)|
|external后封存评分/图版|1.247248|0新物理；最终label更正/本报告准备CPU=`null`|同上与score输出记录|
|external_continued完整Dual，3children|**222.964183**|新增实际完整allowance控制，不覆盖原12；外部lifetime含imports/output|[full Dual ledger](delegated/external_continued/COST_LEDGER.json)|
|完整Dual保存状态分析/绘图|`null`|未独立完整计时；0新物理正演，不能写0CPU|[诊断](delegated/external_continued/AUSTRIA_OBJECTIVE_MATERIAL_DIAGNOSTIC.json)|
|external_selection保存池诊断|.186379|0新正演/0新反演；515/1027/630/1112/1011实际保存对；不同池不是独立统计样本|[selection结果](delegated/external_selection/RESULTS.json)|
|external_physics：GT+shadow诊断|.536772|恰1次1024维物理LU/20sourceRHS；另8G/8S block products，各160columns；计时含该次operator setup|[physics结果](delegated/external_physics/RESULTS.json)|
|external_matched32：3children|**146.655343**|新增匹配诊断；外部lifetime；原150→32结果仍单列|[matched ledger](delegated/external_matched32/COST_LEDGER.json)|
|matched新噪声生成|.001657|0新LU/RHS，复用上述已付clean32GT forward，不再次收取1LU|[generation](delegated/external_matched32/GENERATION_RECEIPT.json)|
|matched评分/图版|.335322|0新物理；后封存GT评分|[matched ledger](delegated/external_matched32/COST_LEDGER.json)|

现有R8 Austria观测生成、NN训练与历史R11终点的获取是继承资产，不假装本轮免费从零生成，也不把历史费用再计入本轮实际CPU。具体历史生成回执被external ledger保留；原NN训练回执位于[继承R8训练记录](../research_max_struct_r8/RESULTS/training/physics/RECEIPT.json)。本轮无新训练。

## 独立确认与工程失败：执行已结束，后处理另计

|范围|已核实CPU秒|物理计数/边界|状态与去重规则|
|---|---:|---|---|
|DEV6观测生成|.503048|6次576-cell LU、72sourceRHS；新operator build0|完成；[生成回执](delegated/confirmation/DEV6_GENERATION_RECEIPT.json)|
|CONFIRM24观测生成|.622058|24次576-cell LU、288sourceRHS；GT分离|该行只计生成；192反演成本另列，不与生成重复；[生成回执](delegated/confirmation/CONFIRM24_GENERATION_RECEIPT.json)|
|本地dev资格CIE子进程|27.347478|whole_child字段；其inner26.532494不另加|与远端重启/确认队列区分；[receipt](delegated/confirmation/RESULTS/dev6/r12_dev6_0000__cie_single__240CPU/RECEIPT.json)|
|本地dev资格DR子进程|18.768386|whole_child字段|独立实际本地执行；[receipt](delegated/confirmation/RESULTS/dev6/r12_dev6_0000__dr__240CPU/RECEIPT.json)|
|startup失败13children|**.906250**|全部在model setup前失败，非科学失败；原source与terminal receipts保留|绝不可删除/重写成一次成功；[failure summary](delegated/confirmation/engineering/startup_failure01/receipts/STARTUP_FAILURE_SUMMARY.json)|
|startup修复后一个有效remote CIE|128.484375|当前交接记录的外部GetProcessTimes；已经保留的有效dev预测|是远端dev cohort的一个子集；未来汇总dev总CPU时不得再加第二次；[handoff](delegated/confirmation/CURRENT_EXECUTION_HANDOFF.md)|
|远端DEV6有效48子进程|7189.687500|externalGetProcessTimes；包含上述128.484375修复运行，不另加它|48/48封存|
|远端CONFIRM24有效192子进程|28147.796875|externalGetProcessTimes；317576次576维物理LU、7598304复RHS，projection与小LU另计|192/192封存，0确认执行失败|
|确认后公开封存|51.562500|独立后处理近写出CPU，末尾未全计；0GT/0FW|不加到单算法收益预算|
|紧凑公开候选导出|7.453125|独立导出近结束CPU，0GT/0FW|完整状态另行归档|
|已记录transport32行快照|.489021 child + .062149 controller|只计本地传输/SSH过程CPU，不是远端科学CPU；后续JSONL可能增加|[transport JSONL](delegated/confirmation/CONFIRMATION_TRANSPORT_COSTS.jsonl)；最初ad-hoc probe未全计时=`null`|

startup失败的原因是浅Windows包路径上 eager `ROOT.parents[3]`默认表达式；修复为lazy显式路径fallback，未改变科学模型、预算或原算法。失败前后原bundle与源文件均保留。独立队列采用240 whole-child CPU、144proposal/96corrector以及atomic-overrun付费记录；这不追溯改变旧basin_v1的body-CPU语义。[当前执行交接](delegated/confirmation/CURRENT_EXECUTION_HANDOFF.md)。

## 仍缺失或不可加总的项

源码/文献阅读、调度者推理、部分prepare/seal/scoring/report、core projector工程基准与旧ad-hoc探测的完整进程CPU均为`null`或仅有局部计时。发布/打包/远端可访问性检查费用当前=`null`。全局GPU费用=`null`（本文已列本地组均CPU执行；不能由此推断整个项目历史无GPU费用）。

因此本初稿不生成一个混合body/external/phase的“精确全局合计”，也不把8臂240CPU预算乘法当成实际花费。确认attempt已按下表汇总；后续归档/发布费用独立登记，不重跑或补造科学结果。

## 确认组逐方法实际外部CPU与物理计数

主/次选择共用一次付费轨迹，费用不能加两遍。外部CPU包括当前实现的投影/FFT/小系统/写出；下表物理LU列不代表所有算子费用。

|方法|24运行CPU秒|平均CPU秒|576维物理LU|复物理RHS|
|---|---:|---:|---:|---:|
|cie_single|3278.625000|136.609375|32471|778140|
|cie_dual|3268.437500|136.184896|32175|771024|
|global_dual|3257.031250|135.709635|32293|773844|
|a_plus|3648.578125|152.024089|30647|734328|
|fullG19|4290.187500|178.757812|32265|773124|
|dr|3082.750000|128.447917|33764|792768|
|multistart|5262.718750|219.279948|89838|2156124|
|nn_direct|2059.468750|85.811198|34123|818952|

整campaign已结束attemptCPU35338.390625 = 有效开发7189.6875 + 有效确认28147.796875 + 原13启动失败0.90625。不能再将前表startuprepair子集另加。recovery02控制器wall3936.777206秒，包含最后4个dev及确认，不等于纯确认wall；controllerCPU12.46875单列。来源：[CONFIRMATION_EXECUTION_COSTS](delegated/confirmation/CONFIRMATION_EXECUTION_COSTS.json)。后续评分、归档、运输与发布见各独立回执；总体全覆盖CPU仍null。

## 封存后分析与恢复验收

确认评分及制图CPU为27.371442秒（独立后处理，不加入算法预算），见[评分费用](delegated/confirmation/RESULTS/offline_scoring_v1/SCORING_AND_FIGURES_COST.json)。附加raw/噪声诊断、3D图版及候选pool oracle各自保留receipt，均不新增反演。

2026-11-11外盘I/O错误导致归档中断，12日恢复；已完成科学结果没有重跑。中断传输、已下载partial、本地压缩和恢复费用见audit/release，发布运输回执继续追加。计时未覆盖的人工阅读/推理与旧准备部分仍为null，不提供虚假精确研究总额。
