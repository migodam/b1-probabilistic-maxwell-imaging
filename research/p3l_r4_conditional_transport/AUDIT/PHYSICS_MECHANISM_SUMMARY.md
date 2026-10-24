# 物理机制控制：封账事实表

本页汇总已完成的 6 次局部前向/JVP 拟合、3 次动作误差排序器拟合，以及前 8 个 DEV 的 4 种冻结几何核。此次提取没有物理查询、拟合或真值读取，也不作最终 GO 判断。逐种子/父对象机器事实与来源路径见 PHYSICS_MECHANISM_FACTS.json 和 physics_mechanism_tables/ 的四份 CSV。

共享几何教师拥有成本是 459,712 RHS（保存 Full-J 状态部分 414,040 RHS；10,063 个独特配对），种子之间不重复计费，复用也不意味着冷标签免费。公共在线锚点取得另外列出。若部署未缓存，单锚点 Full+接收器 J 需 8+32 RHS。

## 六次局部前向/JVP 拟合

|监督|种子尾号|更新/所选|最大白化 L2|95分位 L2（描述）|最大绝对 Phi 误差|原验证/算子测试门限|
|---|---:|---:|---:|---:|---:|---|
|value|501|1000/100|329.3|58.66|9.564e+04|False/False|
|value|502|1000/0|333.5|59.09|1.044e+05|False/False|
|value|503|1000/100|327.9|58|9.69e+04|False/False|
|jvp|501|1000/200|330.8|60.64|1.014e+05|False/False|
|jvp|502|1000/100|334.8|52.51|1.041e+05|False/False|
|jvp|503|1000/100|348.4|55.74|1.085e+05|False/False|

实际保存的通过条件是所有点最大白化 L2≤.05 且最大绝对 Phi 误差≤.05；六次原验证和算子测试门限全部未通过。这里没有增加“95%通过”规则，95 分位仅描述误差分布。每次算子测试为 8 个已暴露开发父对象、790 配对；热查询约 79–83 微秒，但失败门限仍要求每候选 8 RHS Full 回退。每模型反事实全回退 6,320 RHS 未在此执行或计费；不能用短推理时间宣称前向加速。精确概率用途仍需 Full 校正。

## 三次动作误差排序器拟合

|种子尾号|所选材料动作误差|frozenJ|decoder exact|field updated|Broyden|
|---|---:|---:|---:|---:|---:|
|501|0.008953|0.074738|0.010242|0.007242|0.022863|
|502|0.008446|0.074738|0.010242|0.007242|0.022863|
|503|0.008798|0.074738|0.010242|0.007242|0.022863|

每次完成 6,000 更新，误差先按测试父对象再平均。教师是 H=I+JᵀJ 对应的实际 GN 材料动作。排序器输入需要候选点 Full 响应/场特征，冷特征费为每候选 8 RHS；它不是不查询候选前向的响应替代器，经验误差或刷新选择不是部署证书。该表仅陈述与四个实际基线的数值关系。

## 冻结几何核：逐父对象有限 ESS

4 核×8 对象×1,024 步实际增加 262,144 RHS、32,768 Full 状态/LU；每个核/对象 8,192 RHS。实际整作业受控时间 879.573 秒。ESS 是 8 个预定材料函数的 initial-positive-pair 有限运行估计；常量函数为未评估，不是无限 ESS 或 MCSE=0。无共同后验参考，不跨父对象合并。下表为各核有限函数 ESS 的最小值，完整向量见 FROZEN_GEOMETRY_PER_PARENT_KERNEL.csv。

|父对象尾号|isotropic|weak|weak bracket|weak random|公共锚点 RHS|原几何+Heun 拥有 RHS|
|---|---:|---:|---:|---:|---:|---:|
|0|3.573|2.770|2.746|2.569|2776|30168|
|1|3.931|2.864|3.007|3.315|5664|30168|
|2|2.945|2.645|2.479|2.679|4000|30152|
|3|3.480|2.759|2.799|3.051|3088|30112|
|4|2.646|2.553|2.515|2.644|2944|30120|
|5|2.466|2.557|2.606|2.479|1680|29936|
|6|3.010|3.261|4.128|3.822|4440|30168|
|7|2.903|2.514|2.612|2.778|6136|30168|

公共上下文、原几何与 Heun 成本每父对象只拥有一次，四核共享；几何拥有 RHS 是同对象研究的全部配置，不是单冻结核的边际构造费。所用 rank4/66 中央 FD 的 480 RHS 已在拥有费用内，不再计第二遍。线性期望材料 RMS=.002；这些冻结核没有提供任务尺度 .02 的新证据。所有 ESS/MCSE 和核均保持有限运行 UNQUALIFIED。

## 前瞻 locked 四臂：已准备，未释放

入口 tools/run_r4_interval_chain_locked.py，未来名称 interval_chain_locked_v1、Book IV，R4L1_00000–00007，seed20261025171，8×1,024 步 pCN beta=.002。四臂 Full/Richardson/globalGMRES/固定精度 GMRES 共享原协议的缓存/存储权利与小 LS 收费；组合上限 300,000 RHS（65,536/78,154/78,155/78,155），各臂 1,200 秒、3GiB。目标与门限不变。

AUDIT/locked_samecache/FROZEN_PROTOCOL.json 已冻结；configs/INTERVAL_CHAIN_LOCKED_RELEASE_TEMPLATE.json 为 released:false。只有 FROZEN_LOCKED、在线 locked 上下文/缓存/原观测、已关闭 StageIV 上下文收据和精确匹配协议身份的父释放全部具备才可运行。没有等待后自动启动逻辑，缺失输入在 Book 前停止；当前 preflight 明确缺 FROZEN_LOCKED。

fast 引擎只添加默认仍为 StageII/DEV 的 stage/context/cache hooks；旧 global 入口与已封账源码快照不改。两项零物理包装器检查通过不意味着 locked 物理实验完成。paired physics VI 的执行源码没有改动。

事实提取首个分析作业在写正文时使用了错误的核展示键 bracket_enriched，报 KeyError；实际 CSV/JSON 提取与零物理测试已完成且保留，原 ERROR 收据不改。最终整理只读取这些事实，将展示键改为实际 weak_bracket/weak_random，没有重做任何物理或训练。相关成本见 II/jobs/physics_mechanism_summary_and_locked_prepare_v1/CLOSED.json 与 physics_mechanism_summary_finalize_v1/CLOSED.json。此次新增 RHS、状态、LU、神经优化均为 0。

费用字段补充：selected_rank4_66_centralFD_shared_RHS=480 指原三档 h 的整套中央 FD（每档4个Full-J状态×40RHS）；单独 h=.001 的四个状态是其中160RHS。两者均已在原教师拥有费用内，未新增或重复计费。见 physics_mechanism_tables/FD_COST_SCOPE.json。
