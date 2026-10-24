# 同权限 global GMRES 四臂：实际结果附录

闭合比较：PASS_SAME_KERNEL_FOUR_ARM_NUMERICS，全部判决、路径和中间势包络错误 0。额外逐项检查已存的候选 x/q/chi/logU（包含拒绝候选）与 logU，差异 0 / 0。

|方法|总复 RHS|全网格 RHS|小 LS RHS|Full 状态/LU|Full 回退|外部整进程秒|受控秒|物理秒|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|full|65536|65536|0|8192|0|95.336|94.893|67.612|
|richardson|4608|4608|0|0|0|23.718|23.353|0.126|
|global_gmres|14377|4608|9769|0|0|34.692|34.330|0.140|
|global_gmres_accuracy|20489|4608|15881|0|0|43.800|43.432|0.146|

组合费用 105010 / 300,000 RHS，所有方法按预分配子预算运行一次，没有重跑或增加迭代。每种近似方法各自冷付 8×576=4,608 逆列 RHS，同等享有 R0/G/S 八源矩阵运算；R0 矩阵乘法不是新 LU 求解。GMRES 是把 576×8 电流视作单个复向量的普通 left-preconditioned global GMRES，非八套独立或 block GMRES。小型复 LS 逐次先保守计 1 RHS，再作实际 ≤5×4 一右端解；小 LS 与全网格项均列出。

接受近似点没有变成新的精确锚点，原始 R0 保持不变；两个端点的当前精度/basis 可保留，跨材料 basis 不复用。所有尝试、拒绝、Full 回退、最终电流/场和版本/中间 proof 标量保留。共同 stored/canonical endpoint-once 存储权限适用于全部方法；外部整进程时间包含导入、冷 R0、写入、proof 和审计前封存。Book 受控时间与物理调用时间分列，小 LS/矩阵计算不属于 Full 物理时间。

普通 Krylov 同样具备已付缓存权限，因此结果不能自动归因于 A17/current-space 理论，也没有要求 Richardson 胜出。decision 方法与固定精度方法比较的是不同停止准则的实际代价，不能将固定精度更贵泛化为所有科学任务均无需精度。该固定 pCN 核的数值一致性不证明混合、独立样本、MCSE、后验参考或 float64 定向区间证书。

原协议 PROTOCOL_ZH.md 和全部封账结果未改写。机器证据 EXECUTED_AUDIT.json、EXECUTED_COSTS.csv；正式逐项比较 results/II/interval_chain_global_gmres_v1/COMPARISON.json。float32 材料框只读补充见 AUDIT/raw_risk_actions_DEV_float32/REPORT_ZH.md。最终科学获益和普通 Krylov 是否解释观察收益由父线程判断。

实际数值路径补充：所有臂接受 3,286/8,192 步，三个近似臂 Full 回退均为 0。两种 GMRES 实际复 LS 分别 9,769 与 15,881 次，shape 仅 2×1、3×2、4×3，没有 happy breakdown；其纯数值 LS 总时间 .2263/.3717 秒，含付费 ledger 的 callback 总时间 5.6040/8.9699 秒，二者不能混淆。Richardson、decision GMRES、accuracy GMRES 的第三次 refinement proof 分别 19、8、8 条，未用满预设第四次。

缓存峰值说明：初版 EXECUTED_AUDIT 对每父对象 matrix_actions 统一求和，故其中名为 basis_cache_peak_bytes 的聚合值是各父对象峰值之和，不是单进程峰值；真正的臂峰值应取 max，见 BASIS_PEAK_CLARIFICATION.json。两个 GMRES 臂实际 basis 缓存峰值分别为 1180288 / 1180288 字节。原子收据峰值 RSS 不受此聚合命名问题影响，未来汇总代码已改用 max；没有改写已关闭结果。
