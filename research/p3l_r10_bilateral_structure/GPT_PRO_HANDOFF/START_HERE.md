# R10交接：一个真实条件计算收益，一个清楚的材料覆盖负结果

先读 [完整中文报告](../RESEARCH_REPORT_ZH.md)、[GO/NO-GO](../GO_NO_GO.md)、[独立QA](../INDEPENDENT_QA.md)，再读 [新18条推导](../../maxwell_formula_atlas/09_BILATERAL_SHAPE_AMPLITUDE.md)。不要把旧R9或本轮候选负结果从论文主线中删去。

**保留的机制：** 原P1 $D=\delta A(z)$ 精确形状—幅度分离；固定形状传递 $\delta L(I-\delta B)^{-1}R$；标准block-Krylov降阶；由材料损耗与Hermitian虚部下界得到非正规resolvent界；完整residual控制观测与likelihood误差，必要时full fallback。

**实测：** 10形状×32新属性节点，固定r64全部过门；实际M32冷形状成本比直接LU节省51.82%，共享缓存读取也计入则48.54%；M1/M4亏。一个DEV条件积分836节点同规则重放节省59.22%。这些是条件任务收益，不是完整形状成像加速。r64不是最小实现证书，数值稳定性界不是区间认证。

**否定的具体假设：** 宿主双边Gram可通过更好的固定形状目标改善材料GT覆盖。本轮Fresh10：bilateral2/10，Full66随机目标3/10，shape目标3/10，anchor2/10。三新方法同预算；R9对照同K但不同预算。配对bootstrap不支持新提案优越性。

**概率资格：** 3个固定形状Uniform2数值参考通过；初版求积失败与修复完整保留。没有Full66 posterior，没有候选softmax质量，没有连续exact fiber。真正非线性属性profile未运行。

**原创性：** T-operator、eigenpermittivity/material-independent modes、SOM/TSOM、block-Krylov、Rao–Blackwell、RML均为已有基础。当前可写的是有误差门的结构组合与实际成本，强TAP资格HOLD，不能把标准公式重新命名成新理论。

**唯一下一优先方向：** 形状更新时可复用响应的范围及风险控制刷新。验证它能否进入完整非线性形状搜索，而不仅加速同一形状的属性扫描。先作有限、强基线、可证伪的路径实验。详见 [待解决问题](OPEN_PROBLEMS.md)。

完整付费账本323024跨阶数保守RHS、9273LU256；原数据与源码在同一B1 Release，复现见 [REPRODUCE](../REPRODUCE.md)。

## 原预算内追加更新

此前推荐的形状更新复用裁决已执行：旧加权状态V64在24形状192查询中0通过，实际全LU回退使费用达直接LU的2.249倍；每形状重建也因8查询未摊销而为1.531倍。不能继续把这一具体方案列为尚未验证的积极前景。固定形状多查询的正面证据保留，材料覆盖和强TAP仍未突破。下一项仅保留数学上不同的固定总场/电流TSOM共享材料约束对照，未执行；不得再次简单冻结同一加权Krylov基。见研究根目录followup_shape_reuse/REPORT_ZH.md。
