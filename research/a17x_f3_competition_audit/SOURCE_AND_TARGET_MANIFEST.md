# 来源、目标与当前权限

2026-10-17。本轮为独立零新增物理审计，输出目录与所有旧实验分离。当前源码分支research/a17x-f3-audit，起点B1 bf27a17；科学baseline严格固定a323ef223731870b390f875c55f041aa8cc855b4。

## 实际读取与执行

完整解压并读取最新版包README、F3_FIRST_GATE、F3 integration、1340行正式报告、执行提示、claim ledgers、原代数源码和既有结果。包内旧128/4096计划仅作历史建议；入口更新明确audit_only、0 states/0 RHS。未执行SHA256SUMS，也未进行文献再检索。原报告中基于旧F1的“F2未运行”由F3 integration及固定commit实际证据覆盖；原文保留、不回写。

直接核对固定F3的商先验、驻点、四臂、MH、路径、A17区分与账本；没有重跑P0四臂、40万次prior MC、旧45-state诊断或任何真实forward。新缓存代数只读取已付费chi/current/field/prediction、G/S及C/D/T，不读取或应用LU/pivots、不读取GT。

A17原始CORE_THEORY与CONDITIONAL_EXCHANGE_THEORY从公开仓库固定df61a1b533682a883cdaf65259b4528c354bd832取得，正文与本地CLOSURE_R1对应文件一致。A30原始报告固定717d3d26a290248bd2b9aeb55f3d03694b407b84；C1最新入口固定83056881e89329212b16f5be47a4d0dbc4465943。详细路径及远端固定链接见evidence/REMOTE_SOURCE_PINS.json。这里没有重新认证所有后端源文件。

## 不可混用的目标

| 对象 | 材料law及物理 | 当前作用 |
|---|---|---|
| F3 conditional smooth P0 | 24×24 scalar TM/VIE；同parent250073的ring/weak；8源、32复接收；原固定gain/noise | 已封账控制，n_parent=1 |
| 原A30 P0 | 原A25 family生成测度、固定复相位；全材料capacity不等于online恢复 | 历史负结果，不是新P1 |
| C1-P1 | 两连续材料坐标、受限vector-DDA intensity；旧密度修复 | 名称相同但不是新丰富材料P1 |
| Joint P2 proxy | 四维径向、缓存响应插值 | 与full-Mie分开，不用于本轮新发现 |
| Joint P2 full-Mie | 受限径向3D vector-Mie、多频计费 | 旧阶段，不借余额 |
| 待建物理P1 | 未取得冻结prior/decoder/acquisition/likelihood、合格anchor与任务尺度见证 | UNKNOWN；真实实验NOT_RUN |

P0的q=(六原shape均匀坐标,t)，p_q=2L p_rho((2t−1)L)。不改为7维均匀，不加材料Gram因子。合法闭边界与零点密度不同。source-chart的自由复Lebesgue Jacobian不能替代这个law。

## P1核查范围及结论

已核对本机B1所有worktree列表、主仓research/configs/src的P1相关文件、A30原始决策、C1-P1决策/共同接口和最新C1入口；仅查看当前research相关进程，不访问旧locked truth或无关私人目录。B1最新main仍是本轮前F3发布，无新增P1阶段提交；未见匹配的物理研究运行进程。该检查不证明其他设备或未连接会话不存在P1。

找到的是F3的P1建议与NOT_RUN，及不同目标的C1-P1。**新物理P1状态为UNKNOWN/未就绪，而不是COMPLETE。** 未找到可复用P1 Maxwell states、已注册prior law、独立parent竞争见证或本轮新增物理授权。不得从“高维”或“允许自己思考”推断已满足启动条件。

执行限制是用户要求的“P1目标与新预算同时准备好才做真实实验”。新physics预算为0；旧F3 55/712、旧F 24/960等全部原样封账，不借用任何余额。
