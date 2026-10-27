# R6独立收尾：D冻结协议与输出定义审查

2026-10-27。此为**评分前只读协议审查**，不是D结果收尾或最终科学裁决。未导入/运行D评分器、未运行物理/训练/模型推断、未打开任何D真值。D仍在完成DHIGH在线contexts；不存在由本审查给出的D PASS。只写本audit目录，冻结源码/模型/协议未修改。

## 结论

当前 `track_d/score.py` 与冻结global8/20000-bootstrap合同一致，未发现阻止继续当前物理实验的评分定义错误。已把恢复计时的方向性限定单独告知root；其他资格限制必须随最终结果保留。

## 主统计合同逐项核对

|项目|实际实现|判断|
|---|---|---|
|独立单位|每condition256 parents，三seed loss先算术平均；不重采样seed|一致；不是768独立对象|
|主loss|每parent、seed，576 cells平均complex ∥prediction−truth∥²|一致；不是ROI分母或图像ensemble|
|noise两条件|DLOW sigma×.5，DHIGH×2分别评分、分别bootstrap|一致；不混合后报平均噪声收益|
|全局family|2条件×4比较=8，global alpha .05，单项 .00625|一致；单侧99.375% upperCI|
|bootstrap|20,000次，parent有放回；每个condition全部4对比共用同一index矩阵|一致；seed202610270699/+1，保存indices/draws|
|主门1|selected q-Low MSE−1.01 Full-originalDet MSE，upperCI<0|一致；线性NI contrast，不以ratio点值代替|
|主门2|selected wall−.8 Full wall，upperCI<0|一致；包含各自真实anchor与learning chain|
|机制门|selected q-Low MSE−same-low actual-gradient MSE，upperCI<0|一致；是独立增量门，未用通过/未通过NI之差代替|
|次门|selected q-Low−FullDistill MSE|一致；标SECONDARY，不参加champion selection|
|cost_GO/overall|仅前两门都通过|一致|
|mechanism_GO|前三门都通过|一致；成本成功不能升级机制结论|
|区间资格|percentile bootstrap，linear quantiles|有限样本近似，不是确定性证书或严格覆盖率保证|

## 输出loss、seed、phase与成本

- `PARENT_LOSSES.npz`保存[N,3]原始full-grid MSE；compact CSV保存每parent/seed raw与fixed-passive-clip的full/ROI MSE+MAE。seed-mean为主，先平均预测图再算loss的ensemble另列且仅描述。
- 固定clip沿用原 `p3r2.fusion.passive_box`：ROI内实/虚分别投影到host与host+原属性上界的矩形；ROI外固定host。其范围包含所有合法truth，但不保证单一decoder世界或响应拟合。clip不是新的选择规则；raw主门未被替换。
- phase由**Full96/64 anchor**的场/半宽在±1分类，不用truth；ROI外代码−1。固定ROI=256；phase contribution分母576；per-cell MSE分母是该phase跨parent总cell数。先对seed误差和取平均，再除cell counts，未把空phase计为零风险。
- Full/Low物理点基线直接用各自 `chi0` 材料计算，同一Full-anchor phase划分；它们不是新的学习模型或新增主比较。
- optimizer success按原Native/Probit报告统计；finite-best返回进入总体。pipeline任何缺失/非有限输出都中止，未按GT fallback。`pipeline_fallback_count=0`只指外层pipeline，没有自动把非收敛有限端点当作成功。
- phi从保存context读取，使用该noise condition的真实sigma；平均/中位/极值分别报告。stationarity、global optimum、posterior资格不由phi小或SciPy success推出。
- timing每parent=对应Full/Low独立付费anchor wall+该action/seed实测batch64 inference总秒/256，再对seed平均。head/provider加载时间另列；warm推断包含predict_contexts的packing/stacking/feature/head，不含事先从磁盘加载contexts。GPU训练、数据生成、两种路径实验总预算与部署单路径wall分别解释。
- anchor运行顺序按冻结parent index奇偶交替Full→Low/Low→Full，各自独立初始化，无初值/J跨路径传递；GC清理计入anchor wall，Book close/admin receipt wall另在账本保留。
- 推断只有每seed一次batch64时间测量。parent bootstrap主要量化对象组成差异，不能估计重复硬件计时噪声；不得称为跨设备、batch1、冷加载或部署总体加速证明。

## 隔离与epoch

`validate_global`先验证两condition的全部方法×seed预测、phase、Full/Low context和REPORT身份，再允许第一次GT读取；每report的freeze/context身份、GT_read=false和path_complete均检查。root另外提供精确global seal绑定的GT_UNLOCK。online在GT unlock后停止，global seal不重写。当前无自动调用评分的online路径。

初始freeze与serializer恢复后freeze逐项独立比较：groups、seed、manifests、selected config、六arms、四reference IDs、统计字典与checkpoint identities完全相同。仅增加行政修订记录并更新有真实bug修复的source identities。当前全部冻结source BLAKE2b一致；原source/初始freeze、TypeError失败与CLOSED FAILED、零RHS恢复sidecar均保留。

## 异常与需要随报告说明的限制

1. **首Full计时恢复的“保守”方向。** CLOSED=4.1316412499s用于账单是保守上限；原ANCHOR_REPORT run wall=4.1085700828s。增加Full基线却有利于低成本相对节时，不能称它对speed gate保守。丢失区间≤.02307117s，DLOW均值影响≤9.0122e−5s，wall contrast影响≤7.2097e−5s。保留冻结主结果；若接近20%门，root应加预声明描述性lower-Full计时敏感性/区间资格，不虚构精确cleanup或重跑求解。
2. `freeze()`身份验证完整性仍取决于root所列dependencies。当前source列表覆盖实际关键loader/feature/closure和原物理优化入口；不能由“检查函数通过”推断所有未知外部环境也被锁住。
3. `score.py`是冻结、未执行的待运行评价器；本次源码审阅不替代结束后对实际NPZ/CSV/draws的机械重放。root所有512预测封存前不得以本审查为评分授权。
4. 全部范围仅离散TM24、原prior、两个已知noise倍率、冻结模型。不给Full66 posterior、拓扑校准、采集变化、连续Maxwell或全球原创性结论。

## 评分后追加：更强同闭合Full-q参考（0 RHS，仅事后资格）

root另建 `track_d/FULL_Q_REFERENCE_POSTHOC.json`，比较已封存的selected q-Low24/P64与Full q96/P64。两者使用相同q闭合特征算法；这个参考比原primary FullActual更强，不能由原global8中的FullActual非劣性自动替代。本审计只读取已有[N,3]loss与原20,000份同parent bootstrap indices，独立重算全部数字，无GT读取、无模型推断、无训练、无物理。

|condition|selected相对Full-q均值变化|1% NI线性contrast|原pointwise水平单侧99.375%上界|资格|
|---|---:|---:|---:|---|
|DLOW|+0.360564%|−3.3992033e−5|+1.2674483e−5|未确立1% NI；HOLD|
|DHIGH|+0.025439%|−9.3100532e−5|−7.9709250e−5|仅事后pointwise支持指示，不新增确认声明|

重算与root文件一致（`FULL_Q_POSTHOC_INDEPENDENT_QA.json`，无差异）。该额外两比较未在原global8 family内，借用其pointwise置信水平**不保留扩充family的原整体覆盖**；不得把高噪结果升级为新confirmatory PASS。没有重开原family、修改门、调参或重训。

原global8 PASS继续保留，但最终成本/风险GO必须明确限于“对原FullActual参考、原执行窗口”。不能概括为“对所有Full估计器”或“对最强Full-q已确认1%NI”。低噪对强Full-q仍HOLD；叠加共同窗口低噪20%节时下界16.82%的结果，跨窗口/强基线的更广成本优势仍未完整建立。

检查时 `RESULTS_SUMMARY.json` 的GO_C/GO_D文句仍较宽（仅写registered noise experiments / two noise scales）；已报告root，由root收窄到具体参考与窗口，独立审计不改最终科学裁决或summary。
