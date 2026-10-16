# A17与高阶反馈的实际增量

本轮D实际用了原2D full material injection→内部current response→receiver transfer：一次8RHS的一阶current方向，再用已付费receiver transfer得到二阶Maxwell响应。没有重复计算旧resolvent恒等式，没有将低维材料变量当作低秩operator更新。

D在weak任务尺度的ΔΦ1.6327劣于普通native的1.1406；ring稍优于native的5.1964（D4.5023），但都未达到competitive阈值。计完整56RHS预算，D在weak只有1个合格端点，普通native前7个有2个；ring均0。现有小pilot不支持为此增加高阶网络。

**独立A17 Schur/paired-update/ROM加速消融NOT_RUN。** D是普通full-solver方向Hessian加几何修正，不能直接称A17算法。其没有新增收益既不证明A17所有新接口失败，也不影响原A17冻结GN理论。新Hessian和补偿数据保留，未用于训练。

没有严格误差区间、certified early rejection、full-vector 3D或跨物体结论。不宣称普通full/Krylov/A17的全成本优劣已测。
