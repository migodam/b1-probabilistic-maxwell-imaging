# 失败及不能宣称的事项

1. Geometry FD paired-row索引错误影响两轮机制解释；旧131072 proxy保留，单独repair131072。目标不变性并不意味着实现了正确geometry。
2. E2首任务35states/312RHS后continuation的GaussianCDF坐标达到数值端点而失败。这些35状态由于当时end-of-anchor保存设计丢失endpoint原数组；不能以ledger代替原始证据，也不重构成发现。后续改每endpoint先落盘，缩为每arm5步，56states完成；最后5states用于导数诊断，总96，未提高上限。
3. 原本material-normalized weak方向映射成1998倍prior参数方向；盲目有限步与一阶预测不相容。trust控制保护prior后多数计划幅度变成同一端点，field去重后不能把图数叫新结构数。
4. developmentN3无原P0合格anchor，不能把低MAE和noise球条件互换。已知TRAIN校准不能填补未知物体的inverse验证缺口。
5. Global MH胜过hybrid，gain-aware geometry只恢复普通local水准；当前无必须加入local/global复杂混合或A17的证据。
6. 多事件/ridge不是分离mode证明，E3不运行；fullMie budget未花。没有新训练、noise放宽、prior扩展或locked truth读取。

辅助记录写入一次stdin UTF8失败，physics/proxy费用均0。原异常日志与FAILURE_RECORDS.json保留。所有未测指标明确NOT_RUN，不以漂亮图库替代缺失证据。
