# A17增量：NOT_RUN / NOT_ESTABLISHED

此次有效局部几何基线使用原full material derivative/receiver adjoint，radial几何是普通proxy Jacobian及gain Schur投影。它们不属于A17增量证据。没有同proposal、同精度且计完整setup的A17 likelihood实现，也没有output uniform bound，故没有certified early rejection。未重新跑frozenGN gap，也未修改A17冻结稿。

原P0 smooth校准规模full-LU单state很便宜；优先保留ordinary full。只有未来存在更大可行chart且具体reduced-likelihood有实测总成本优势，才值得追加A17。A17原局部GN结论保持独立。

## 最新独立N3应用结果交叉核对

本地A17/N3_APPLICATION_R1正式入口、PHYSICS_ADAPTATION、MATCHED_COMPARISON及FAILURES已查。原8暴露条件N3MAE .041756，anchored fullGN约.043917仍未改善；同rank反馈交换的frozen gap下降24.1%/40.5%没有转化为MAE优势，ordinary materialCG平均更好。这是已有独立实验，不重复求解、不并入本轮预算，也不当作本轮A17 likelihood消融。其材料box/anchored正则不是原P0_EXACT先验，不能填补本轮合法anchor缺口。
