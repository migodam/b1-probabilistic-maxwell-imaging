# 完整成本与独立授权

用户Research T独立overlay：192full-state、6144complex RHS、7200秒受控CPU；GPU0、新3D0。不叠加Pro建议，不借A30/A/B/Joint P2余额。

|阶段|full-state|complex RHS|受控job wall|
|---|---:|---:|---:|
|资格Q，两臂8case|16|448|5.217s|
|机制v1，16编辑＋base重分解|68|2112|80.830s|
|机制v2，同16编辑|68|2112|83.836s|
|机制v3，8编辑|36|1088|42.550s|
|实际Jacobian/独立直接解核验|4|72|0.779s|
|**合计**|**192**|**5832**|**213.212603s**|

RHS独立核算：192×8=1536 forward；8个Q初始Jacobian＋40编辑×3刷新＋1核验Jacobian=4128adjoint；21次缓存LU核对×8=168，合计5832。没有免费adjoint或失败重试。

第一瓶颈state预算已耗尽，剩312RHS不可单独支付新fullstate。RUN_CLOSED.json及唯一Book在求解前拒绝超额；历史ledger原样只读。新账本COST_LEDGER.jsonl包含所有STATE/RHS/ACTION/GN试点/JOB_BEGIN/END，global snapshots不能逐bundle再求和。

独立零物理E0成功wall4.949秒；离线三个成功后处理累计19.250秒；导出/分析各receipt另列。一项失败后处理wall未记录（null），所以213.213s是完整计费solver jobs wall，不是全部人工实现与文件处理end-to-end wall。cached Green/S矩阵运算、保存背景setup/历史摊销未知额外列为null，不冒写0。

未运行同配置冷/暖重复物理计时；各新进程job是实际测得的流水线wall，已有效支付cache读取和refactor。旧bundle可零物理重放，但I/O快不能称warm correction加速。没有用不同轮算法差异冒充cold/warm对照。
