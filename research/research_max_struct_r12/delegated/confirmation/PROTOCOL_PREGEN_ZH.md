# R12 独立6开发 / 24确认：生成前协议

后续执行更新：root已写入ROOT_CORE_FREEZE.json与CONFIRM24_GENERATOR_FROZEN.json，再授权生成24确认；24现已生成，数据未评分。八方法队列已按root授权在remote启动并完成工程启动修复。当前状态与成本历史见CURRENT_EXECUTION_HANDOFF.md。下文保留生成前约束及当时的待定事项，不应当作当前“尚未生成”的状态声明。

本文件在任何新材料生成前建立。后续root授权只生成6dev，已完成，见DEV6_GENERATION_RECEIPT.json；没有24确认生成，没有远程启动。root要求先以6开发对象评估实际运行与资源费用，不以开发GT调参。24确认配置明确为 `PENDING_ROOT_CORE_FREEZE_DO_NOT_GENERATE`，脚本拒绝执行；即使改成 FROZEN_FOR_GENERATION，还要求 root_core_freeze_receipt。root负责DR核心、统一目标/先验与最终科学冻结。

## 固定物理与成熟基线

读取已有 R10 `research/research_max_struct_r10/data_public/twod/operators.npz` 的真实 G/E0/S，24²材料、12sources、48receivers、k=2π、radius=.018、原source/receiver geometry。完整576实材料自由度，box[.02,20]，均匀.1初值。没有toy grid替换、Fourier映射、dense Green注入或新的正演模型。

基线直接重用既有R8 `experiments/run_case.py` 与 `src/benchmark_solver.py`、`src/cie_fft.py`，同R9 noise-ball DPCS策略；sourcewise CIE/Dual、global PR Global/A+/G19、原共享实材料读出，原 native Fourier mf3/6/10/12及500/500/500/1500步。CIE current/material β同为1/.2/.2/.2；其他四臂current β恒1，material β同前。PR+、固定材料精确二次线搜索、仿射DPCS offset与native系数transport均未改。没有GT或形状规则进入逆解。

A+直接调用原R7 Möbius/MM（arc candidates、两sweeps、inherited gtol）；G19直接调用R11 full_material_readout/prepare_full，最终转发原R10 .1公开Green谱阈值与Galerkin PGP、完整材料导数。实际rank从每次公开setup receipt取，不强迫19。它们仍用相同global current，不把G19称为A+别名。

先验是原 `Regularizer((24,24))` 默认：smooth=tv=1e−4，epsilon=.01，原 D 使用19.5差分尺度，按edge数平均；原物理Φ=.5||SJ(m)−Y||²/||Y||²+R(m)。本轮没有去掉或重新调整先验。DR也必须按root冻结满足相同先验合同；这部分目前未实现，不能把无先验可行性DR直接当成最终六臂结果。

原corrector完整继承R6数值body：scale1e6，ftol1e−12，scaled gtol.01，maxls30，3000 accepted iterations、12000 function上限、原box。统一每臂240**整child CPU秒**，含导入、DPCS、setup、辅助、true-physics audits、corrector。root后续冻结proposal最多144CPU（含setup），corrector再最多96CPU，并受整child240限制；不通过缩短成熟stage迭代数模拟预算。proposal在原数值轨迹的检查边界截断，handoff此前已经做原物理验证的阶段端点/初值；未完成stage的current不伪装成新免费physical候选。若proposal截断，完整loss history不再返回，薄prefix日志明确标注，operation counts保留。预算检查按R11 WholeChildBudget，atomic完成但越界的工作仍付费，越界状态不能进入cap预测；只保留cap内已经付费的最低原Φ候选，无额外final正演。30/60/120/240 CPU快照保留。外部Windows GetProcessTimes全生命周期为权威成本，内部timer不相加。

## 新分布与不可更改的抽样

`generate_observations.py::material` 是显式生成分布，不是GT支持输入给solver。

|家族|固定生成机制|mixed含义|
|---|---|---|
|levelsets|六个固定频率、随机phase/amplitude的三角场，固定阈值0，椭圆窗口，随机中心/旋转|support内由旋转线性分区分配high或2.5|
|branches|一个随机扰动junction及三侧端点/根部，四条带正弦弯曲的分段支路，随机固定范围width|按固定v>.015分配high或2.5|
|mixed_islands|五个分离椭圆分量，角度/半径/轴长在固定范围独立扰动|分量0/2/4为high，1/3为2.5|
|continuous_random|高斯平滑随机场乘紧支撑平滑椭圆envelope，整体归一到max1；材料.1+(high−.1)field|再乘固定低/高平滑partition后归一，得到连续材料，不伪装成二相GT|

dev6固定seed202611111206，六个独立geometry/noise子流；设计见DEV6_GENERATOR.json。独立confirm24 seed202611111224，每家族6个独立draw，合计high8/high18/mixed8/mixed18各6个；同名重复是独立draw，不是同geometry改变contrast。确认前不生成、不窥视、不依据开发恢复结果更换family、contrast或seed。

structured材料必须support非空，mixed必须high/low两部分均非空；如果固定draw违法，**显式失败并保留，不重抽**。不按forward residual、conditioning、幅值、重建质量或任何阈值筛选/reshuffle。连续族不要求二相support标签。主比较应报告连续材料误差与原数据目标；若root需要binary IoU或成功标签，必须在确认前定义连续族的合理口径，当前脚本没有偷偷增加成功阈值。

噪声为raw iid circular complex，E|n|²=σ²，σ=.01||Yclean||/sqrt(Yclean.size)；geometry与noise使用两个独立子流。公开文件只写Y、σ、噪声模型与全局grid/box元数据。truth、Yclean、family/regime、seed、支持、真实amplitude与生成参数全部单独放sealed_truth，solver manifest不含这些字段。生成整childCPU、每个forward的LU/RHS与setup/copy均记录；已有operator的原构建不重做。

## 实际产物与执行边界

- `DEV6_GENERATOR.json`：六开发生成设计，已按root授权生成，未评分。
- `CONFIRM24_GENERATOR_DRAFT.json`：24确认设计，等待root核心冻结，禁止运行。
- `DEV6_BASELINE_SOLVER.json`：五成熟基线及240CPU合同，直接继承历史24²设置。
- `DEV6_SIX_ARM_PLAN.json`：36个预期task，30个成熟基线准备就绪，6个DR明确PENDING；不是36个已执行结果。
- `generate_observations.py`：只生成public+separate truth；自动生成纯public baseline_tasks.json及DR pending清单。
- `run_baseline.py`：R11 wrapper的R12输出目录副本；目录、legacy repo路径、固定readout guard、receipt study与144/96观察性phase预算改变，原PR/readout/corrector数值body未动。变更清单见INHERITED_WRAPPER_CHANGES.json。继承的R11_输出文件前缀保留，receipt study为R12。

Windows既有只读代码根拟用 `D:/AI/MAX_STRUCT_R11_DEV_20261107`；新R12目录拟用 `D:/AI/MAX_STRUCT_R12_20261111/.../delegated/confirmation`。既有Python `D:/python/python.exe`，vendor路径仍用旧根vendor。**root负责资源协调/新目录同步与启动，本worker未连接或启动remote。** 不覆盖老结果，不修改旧代码，不安装服务/依赖。

直接多起点与NN→direct仅是root提出的后续可选控制，尚无冻结算法/初始化成本合同，不出现在已准备五基线清单；不得把占位写成已完成。
