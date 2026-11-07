# Fresh24：冻结比较、真实成本与负例

**新矩阵metric没有得到总体优势；它在0001产生大幅正例，在0004/0022产生大幅反例。多数其它对象差别小，不能用11/24对8/24把所有退化都描述为大失败。** 本文仅解释168份完整封存端点及已付保存状态，不改变方法、预算、阈值或选点。

第一GT读取在`ALL168_PUBLIC_SEAL.json`写出后，使用root明确的postseal授权。168任务均有效，0任务失败/缺失、0评分错误。主门严格沿用冻结合同：RMSE除sealed amplitude、预测/GT都按chi>1计算IoU、NRMSE≤.15且IoU≥.7且材料合法。legacy半幅度Otsu与full-support Otsu独立；原值未替换主门。全24材料图按low/high/mixed三页展示，没有筛对象。

**Fresh Twofold实际rank30（S强11＋补投影G强19）；历史DEV才是28。** 已仅修正fresh图标签，未改数值、阈值.1、方法或LU计数。原误标签图保存在`delegated/twod_score/label_refresh_rank30/previous_rank28_label_archive/`，重绘独立收费。

## 主结果与完整child CPU

| 方法 | 主success/24 | low/high/mixed各8成功数 | 平均NRMSE | 中位NRMSE | 完整child CPU/对象均值(s) | 全24 CPU总和相对CIE |
|---|---:|---|---:|---:|---:|---:|
| CIE-single | 11 | 4/3/4 | .162871 | .138007 | 157.705 | 1.0000 |
| Dual-CIE | 10 | 3/3/4 | .163460 | .138285 | 153.982 | .9764 |
| A+global | 8 | 2/2/4 | .180444 | .145312 | 169.493 | 1.0747 |
| A+source | 10 | 3/2/5 | .180453 | .144790 | 174.979 | 1.1095 |
| Twofold r30 | 8 | 2/2/4 | .180994 | .145439 | 192.527 | 1.2208 |
| Green19 | 8 | 2/2/4 | .176687 | .146383 | 187.656 | 1.1899 |
| GlobalDual，内部消融 | 10 | 2/3/5 | .163693 | .138944 | 150.135 | .9520 |

完整child CPU来自远端GetProcessTimes，168个disjoint child合计28475.421875秒；queue controller、generation/transport/archive另在原账本，不与nested R10/R8/aux/corrector计时重复相加。主物理576维LU与候选19/30维小LU分列。这是实际端点支出，不是同迭代就等总算量。

对象paired bootstrap固定seed20261107、10000次。平均NRMSE差“新−CIE”的描述性95%区间：Twofold为+.018123，CI[−.026356,+.068553]；Green19为+.013816，CI[−.029851,+.061470]。两者均只在9/24对象NRMSE低于CIE。Green19相对Twofold平均差−.004308，CI[−.012848,+.000146]，15胜/9负；不足以称接收Twofold有特有增益，也不确立G19的总体优势。区间只描述既定对象，未作广泛总体保证或多重比较校正。

## 二值门不能掩盖差异尺度

Green19相对CIE有21/24对象的绝对NRMSE差≤.005；Twofold为20/24，中位绝对差分别约.00100/.00085。主要大分叉是0001、0004、0022，详见`THREE_CASE_INTERPRETATION.md`。

0000：CIE .149114、Twofold .155165、G19 .154003；0021：CIE .148489、Twofold .152298、G19 .152648。IoU均过门，success变化主要因NRMSE跨越.15，绝对增量约.0038–.0061，不应称灾难性材料崩溃。它们与0004由.111到.515、0022由.097到.337/.437的退化性质不同。全部103个material FAIL（含内部GlobalDual）保存在`ALL_MATERIAL_FAILURES.csv`，不是进程/评分失败。

## 固定已付预算切片

下表为“平均NRMSE / success数”，每格24对象均有可用保存状态。外部CPU中除已完成终点外，是**保守完整child CPU上界包络**，不是逐步外部CPU采样；内部CPU轴独立保留。physical LU只约束一种576维计数，不包含小LU和其它work。

| 预算 | CIE | Dual | A+global | A+source | Twofold r30 | G19 |
|---|---|---|---|---|---|---|
| 外部保守100s | .1672 /9 | .1665 /10 | .1836 /6 | .1868 /7 | .1868 /5 | .1850 /5 |
| 外部保守200s | .1628 /11 | .1635 /10 | .1803 /7 | .1804 /9 | .1785 /6 | .1773 /7 |
| 外部400s，均已完成 | .1629 /11 | .1635 /10 | .1804 /8 | .1805 /10 | .1810 /8 | .1767 /8 |
| 原物理LU300 | .1782 /7 | .1708 /7 | .1847 /7 | .1899 /6 | .1809 /9 | .1844 /6 |
| 原物理LU600 | .1699 /8 | .1670 /10 | .1820 /9 | .1838 /10 | .1781 /9 | .1809 /7 |

2016个预算单元、原选中index/成本/clock语义及coverage均保留在`BUDGET_SCORES.csv`；汇总和图为`BUDGET_SUMMARY.csv`、`BUDGET_CURVES_SCOPE_EXPLICIT.png`。只取最后已付原Phi-incumbent，未按GT选时间或插值。某些原物理LU预算下Twofold success较多，但它的材料/小系统费用不受此单一counter控制；完整CPU切片并无对应总体优势。更久/更低Phi也不保证材料门单调改善。

## FFT/current能量：不等于材料可辨识性

Fresh同case七臂共享实际DPCS rank10–14及mf3/6/10/12的affine族。实际系数range秩为36/144/400/562–566。每case oracle只算一次，20个rank/mf空间cache支撑96份case/rank/mf oracle，没有为七arm重做SVD。末mf12不可表示GT current能量占GT current的.000225%–.014336%（中位.003732%），远小于许多实际current误差；它不能独自解释大的trajectory失败，也不保证物理一致的材料可恢复。floor只约束辅助J，final物理j不受该affine限制。

固定tau.1汇总的材料误差平方能量：CIE80791.37、A+128003.28、Twofold147423.77、G19133077.81；对应物理current误差为6074.13、8664.79、10004.32、10625.11。新候选没有在总体绝对能量上隐藏改善。材料Gob补空间比例约66.65%–69.54%；实际current则约99.64%–99.92%在Gob补空间，但这不是材料Jacobian或唯一性证书。

反过来，Gob强带虽然只含很小的current输入误差能量，却承担CIE约99.980%、Twofold约99.924%、G19约99.929%的接收场误差平方能量。输入分量占比与传播后场误差占比不能混用，不能从current补空间巨大就推出数据残差也必然小。

辅助J与物理j必须分开。末stage按24对象合计，G19的aux误差11023.16、匹配物理误差12954.50、gap2919.25、signed cross−4850.59；Twofold为10463.61/11826.18/1764.82/−3127.39。不能丢cross或将辅助误差当物理恢复误差。三tau全部表已保留，paper图固定.1；Twofold串联、交集/remainder、overlapping span依然分列，传播孤立能量不可强归一化。

## 数据合格与材料合格是两件事

CIE只在0001/0016未通过单case CN99，其余22例数据合格；主材料成功仍只有11。利用已缓存GT current和公开S/Y/sigma，GT24/24都通过各自CN99。CIE“数据过门但材料失败”的11例中，GT与估计均满足硬box并同在实际y的单case CN99球：0003/0007/0008/0009/0010/0011/0012/0013/0014/0018/0020。

这是**本噪声允许集合里的两份合法材料实例**，不是精确多解、材料信息论下界或24个对象联合99%事件。smooth/TV是soft prior，不是硬可行域；它仍可对两份材料给出不同偏好/Phi。不能把这些noise-compatible材料偏差都称为wrong-data basin，或期待单纯放宽优化预算必然恢复GT。

## 明确标记的posthoc部署诊断

额外只对保存的CIE和G19端点按最低原Phi选择，GT完全不参与selector，所有24对象均支付两臂费用。结果选G19 11例/CIE13例，平均NRMSE .149599、中位.127391；但主success为10/24，低于CIE11/24：新增0001，同时丢掉0000/0021两个阈值附近成功。

费用是8288.65625完整child CPU秒（2.1899×单CIE）、79400物理LU/1903224复RHS，G19小系统另计，未选臂不免费。该portfolio是看到结果后请求的补偿性部署诊断，不是新确认算法或数学创新，不能把降低平均NRMSE包装为保住所有CIE成功。脚本及完整24选择表为`posthoc_portfolio.py`、`POSTHOC_PORTFOLIO_OBJECTS.csv`。

## 成本与交付入口

原scoring CPU12.864411秒，无新physical calls；energy CPU44.827705秒，恰24个OFFLINE_SCORING GT forward LU/288 RHS、0 adjoint、0端点/阶段补forward。后续标签重绘、保存结果汇总和posthoc/noise诊断各有独立receipt，不混入原inverse成本。所有fresh图已目检，包含全部24对象；源码/receipt保留算法冻结与标签更正的区别。

这些已完成分析窗口统一索引为`POSTSEAL_ANALYSIS_COSTS.json`；geometry及逐GT solve子计时为nested，未重复加到energy总计。新增right3和conditional4的inverse/评分费用不在这里。

关键数据：`SCORES_FINAL.csv`、`OBJECT_WIDE_COMPARISON.csv`、`REGIME_METHOD_SUMMARY.csv`、`PAIRED_BOOTSTRAP_SUMMARY.csv`、`BUDGET_SCORES.csv`、`CIE_GT_CN99_FEASIBLE_PAIRS.csv`；能量完整表位于`delegated/twod_energy/scored/`。此处没有纳入root之后新增的right修正或独立conditional补偿确认，它们不能改写原168比较。
