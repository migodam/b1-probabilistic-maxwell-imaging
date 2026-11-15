# R15 配置接线与继承字段说明

新manifest保留R14 protocol的全部旧字段以便逐字追溯，同时新增r15命名空间；不能把每个遗留描述字段都解释成R15正在运行的设置。主算法源码冻结于ba6d68a87，原始文件不在运行中改写。

- 本轮身份以顶层EXPOSED24_FOLLOWUP、tasks和EXPERIMENT_CONTRACT为准。遗留protocol.scope仍写原R14首次独立生成时的名称，不将本轮重新变成盲测。
- 新八个stateful臂由TwofoldG19显式采用物理材料初值的current warm state；这与原G19接线检查相同。继承的protocol.initialize_physical_current=false是旧metric-factorial字段，不控制该构造器。
- 新八臂材料读出使用原G19 matrix_material_readout的LBFGS8，完整移动分母导数、旧候选比较和box/prior保持。遗留material_inner_maxiter=60属于原scalar/factorial设置，不覆盖matrix_inner_maxiter=8。
- 当前步始终使用原global scalar current PR；beta_current=1，材料beta为1/.2/.2/.2。cycle的DRscale6与omega1见protocol.r15；没有AP准备。遗留protocol.AP_preparation=128仍供未改动的原生baseline自己的逻辑使用，不用于新路由算法。
- 所有stateful臂mf=3/6/10/12、maxiter=500/500/500/1500，coefficient_tol=0；原生CIE、Dual-CIE、DR保留各自历史合法设置。sourcewise/global差异不被悄悄改成同一个算法。
- 构造器的TWOFOLD_SETUP.json、每次ACTUAL_CONFIG.json、完整checkpoint frozen_config及原生WRAPPER_PROVENANCE共同提供实际接线。新旧G19短程数值等价仅在已检查公开对象上成立，不保证跨硬件同CPU终点相同。
- 120/72是所有臂实际绝对整进程CPU期限；外部child CPU还包含退出序列化尾部。任何超预算操作仍计费，其晚完成候选不进入预算内选择。

这份说明只消除元数据歧义，不改变冻结预测、参数或在线选根规则。
