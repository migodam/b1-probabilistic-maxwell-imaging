# E3：四份真实 Research A 材料候选交接

这里只导出已完成实验，没有重跑 Maxwell、优化或更新。每个 parent 提供 robust 新测量、protected epsilon corrector、三起点连续材料 corrector 三份独立目录。

每份 `HANDOFF_COMPLETION_RECEIPT.json` 索引两个已经完成的 corrector。原 robust 测量 history 中的通用 `PENDING_INDEPENDENT_RECONSTRUCTION` 字段按原样保留；它不表示后续 corrector 尚未执行。实际完成状态、完整候选目录与负端点日志见外层 receipt 和 `completed_corrector/`。

A 原材料与 gain-applied 预测、去 gain 的原频响应、旧观测和完整 old+new history 都保留。候选没有 posterior 权重；residual 只作共享有界 gain 的 profile 选择。四个 parent 仅支持描述性结论，不能与 B 的32个 parent 混作一个测试集。

`raw_e3/` 保留已支付 full fields、真实新增模拟观测和全部端点记录。复制与导出新增物理调用为0。
