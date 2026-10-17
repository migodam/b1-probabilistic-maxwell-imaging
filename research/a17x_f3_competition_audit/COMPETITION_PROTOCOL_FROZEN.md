# 下一轮竞争协议 v1：设计已冻结，P1绑定与执行未授权

2026-10-17。**design_frozen=true；P1_bound=false；physics_authorized=false。** “冻结”仅指下面公平比较和停止规则，不表示尚未存在的P1已被偷偷定义。当前所有真实物理组NOT_RUN。若后续P1不适用，必须建立新版本，不能填入有利结果后继续沿用此版本。

## 进入条件与材料约束

先取得共同P1 prior sampler/decoder及其法测度、被动材料和空间尺度、固定acquisition/noise/gain、本构与mesh、合格physics-only anchor、原始state/current/field/prediction和相应费用。不得由GT选anchor/方向/端口。任务材料RMS s_P1、结构差异和likelihood参考必须在观测质量结果前由TRAIN或物理任务确定；不能拿F3的.00527555硬套到新问题。

优先两个独立material parents，每parent一个已合格anchor、一个GT-blind共同方向、一个任务半径，形成两个机制块。只改变背景不能算新parent。方向依据共同先验及可付费粗物理规则冻结，全部方法共享并付相同setup成本。没有合格anchor则该块NOT_RUN，不用GT替代。

本轮不建议先训练新P1 generator。对任意已提供P1，先检查其实际DeltaX支持和port谱，材料相位/混合/相关性/geometry约束全部保留。有理proposal不合法应拒绝或在原坐标中修正，不能clip为box材料后算成功。

## 轴一：proposal竞争

共同完整Maxwell evaluator、相同prior、材料RMS目标、起点、seed、进度约束和终止精度：

| 提案臂 | 具体允许机制 | 不允许的免费信息 |
|---|---|---|
| F3-native | 原生合法decoder、材料RMS retraction，必要时有限非线性修正 | GT、未付费Full-J或其他臂current |
| X-rational | exact-port或明确标注的近似H/O，shared-material有理proposal/约束修正 | 把自由复chart直接当P1 prior；将A17误差当0 |
| Ordinary Full-J | receiver-adjoint共享的完整J、标准有界GN/LM，保留非零进度 | 故意逐source/逐变量重复分解的弱baseline |

每臂每块最多2个intermediate完整state加1个独立endpoint audit，合计18states。回溯失败、拒绝和无进度均计费；距离不匹配单列FAIL，不缩小尺度救成功。所有臂最后评价共同原始完整材料，不用approximate端口材料替代它。只发现candidate，不叫posterior sampling。

## 轴二：固定提案的求值竞争

先冻结轴一产生的全部完整材料，保留所有失败与拒绝，再使用同一批材料比较：

1. Optimized full solver，复用同材料factorization和source packing。
2. Exact T-matrix/Woodbury，使用真实U及真实DeltaX。全支持时照实支付s≈N的setup和dense有限反馈成本；不能省略这个强对照。
3. Ordinary goal-oriented primal/dual ROM或residual Krylov，允许相同basis pool、增量Schur复用与同精度fallback。
4. A17 finite exchange，明确记录当前resolved exchange、full−new剩余tail、H/O误差及失效fallback。

基于本轮端口秩审计，A17评价默认采用MECHANISM_REVISION中的**直接有限状态RHS**版本：固定完整材料，在L_D及DeltaX E0上增量更新current空间，不要求先构建完整R0 U/H。低端口版本可作条件消融，不能把两个实现事后择优后只报告赢家。Exact-port仍是独立强对照，不能因默认版本不用端口就省略其实际setup费用。

普通ROM与A17若用相同[V,C]，在精确算术下有限预测必须相同；该配对是实现控制。额外比较只能来自不同、同成本action选择/复用/有效证书，不能把恒等式当方法胜出。

Full reference与GT位于隔离evaluator，不向selector免费泄漏。所有candidate modes的构造、被试但未选中的actions、G/G*作用、正交化、full/JVP/adjoint、端口列、Schur core、bounds、refresh及失败均记录。没有可靠尾界时不标certified early rejection；经验区间与完整fallback分开。theta≥1记BOUND_UNAVAILABLE，不等于物理奇异。

## 指标与费用规则

proposal主要终点为相同材料距离、prior合法、原noise及固定likelihood要求下的非重复进度数和cost/有效位移；结构差异单列。评价轴主要终点为固定提案full输出误差、likelihood/接受决定一致性、未解析tail及完整付费成本。没有posterior reference不报mass或校准。

共享physical调用仅在实验总账记一次；每个算法从零运行的cold-equivalent费用另列，不能只给A17免费cache。warm需明确复用内容、刷新条件和部署频次。两套臂不能因共享physical reference而都被宣称零成本。首pilot不安排重复timing；速度排名NOT_MEASURED，后续计时须单独计费。

## 最小后续预算建议：尚未授权

前提是**P1的两份合格付费anchor和prior已经交付**；本建议不包含从零构造P1、重建anchor、posterior reference或NN。两个anchor端口数为s1,s2时：

| 用途 | 新state cap | complex RHS cap |
|---|---:|---:|
| 3proposal臂×2parents×(2中间+1完整audit) | 18 | 144（8source前提） |
| 两anchor共享Full-J（32receiver） | 0 | 64 |
| Ordinary Full-J臂最多4次刷新J | 0 | 128 |
| 两anchor精确R0 U | 0 | s1+s2 |
| 共同primal/dual pool，每anchor最多16+16列 | 0 | 64 |
| 6失败/验证states及额外失败/probe预留 | 6 | 112 |
| **合计** | **24** | **512+s1+s2** |

对24×24标量全支持，s1,s2≤576，RHS最多1664；受控CPU physics wall建议600秒，GPU/Mie/新3D/NN仍0。这是可审计的预算请求公式，不是“已批准1664”。实际P1若source/receiver/本构不同，先重新计算协议及权限，不套用此表。小端口支路可更省，但不能事前假定s=7。

若实际目标不可比较、端口setup没有可接受摊销、prior law缺失或新预算未获授权，只交设计并停止；不为寻找正结果挪用旧额度。
