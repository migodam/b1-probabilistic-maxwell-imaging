# B2实际结果检查点

状态：**SAVED_COMPLETED_EVIDENCE / ORACLE_EXISTENCE_NOT_ONLINE**。本次构建只读已完成保存记录，新增Maxwell/LU/RHS为0；没有优化或预算操作。P1状态另按其实际保存结果显示。

## 1. 信息条件与几何限制

native分支从已知参考世界构造108参数固定支持chart，属于Oracle物理存在性探索，不是在线重建。参考与FULL_J_FINAL实际材料见[双世界图](figures/saved_results/B2_joint_world_material_gallery.png)。目标沿x平移0.1，圆盘半径/介电对比/质量不变，target_mass_relative_change=0。原minimum_feature标量语义保持**UNRESOLVED**；conditional fixed support不认证全部原几何材料类。科学解释以父线程摘要及原始gate为准。

## 2. 已保存native候选

|来源|方法|d_cal|d_T|direct-joint distance|status|
|---|---|---:|---:|---:|---|
|results/B2/native_cal_2026102201_narrow/RESULTS.json|TASK_AWARE|0.23753082048631677|1.118000639718423|—|NOT_YET_INDEPENDENTLY_VALIDATED|
|results/B2/native_cal_2026102201_narrow/RESULTS.json|RANDOM_MATCHED|0.3791567383950729|1.22673993006909|—|NOT_YET_INDEPENDENTLY_VALIDATED|
|results/B2/native_cal_2026102201_narrow/RESULTS.json|FULL_J_WEAK|—|—|—|NOT_RUN_NO_BOX_FEASIBLE_MATCHED_DIRECTION|
|results/B2/native_cal_2026102202_narrow/RESULTS.json|TASK_AWARE|0.241497554939435|1.022869606671758|—|NOT_YET_INDEPENDENTLY_VALIDATED|
|results/B2/native_cal_2026102202_narrow/RESULTS.json|RANDOM_MATCHED|0.6104296053755314|1.0079276257618743|—|NOT_YET_INDEPENDENTLY_VALIDATED|
|results/B2/native_cal_2026102202_narrow/RESULTS.json|FULL_J_WEAK|—|—|—|NOT_RUN_NO_BOX_FEASIBLE_MATCHED_DIRECTION|
|results/B2/native_cavity_direct_2026102202_narrow/RESULTS.json|FULL_J|—|—|0.5214006441236753|NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION|
|results/B2/native_cavity_direct_2026102202_narrow/RESULTS.json|RANDOM_MULTISTART|—|—|0.5214006458788107|NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION|
|results/B2/native_joint_direct_2026102201_narrow/RESULTS.json|FULL_J|—|—|1.112671461067532|NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION|
|results/B2/native_joint_direct_2026102201_narrow/RESULTS.json|RANDOM_MULTISTART|—|—|1.1126717849390928|NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION|
|results/B2/native_joint_direct_2026102202_narrow/RESULTS.json|FULL_J|—|—|0.7977170125472409|NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION|
|results/B2/native_joint_direct_2026102202_narrow/RESULTS.json|RANDOM_MULTISTART|—|—|0.7977166354692127|NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION|

## 3. 独立coarse/fine完成项

|验证ID / 条件|coarse两世界距离|fine两世界距离|world0网格差|world1网格差|
|---|---:|---:|---:|---:|
|cavity_seed2_narrow_bestcal / cal|0.521325486513|0.520776484651|0.00810962925561|0.00746378611354|
|cavity_seed2_narrow_bestcal / joint8|0.921362828122|0.921876198547|0.00567462310677|0.00693900749933|
|joint_seed2_narrow_fullj / cal|3.56306363214|3.56051191121|0.00810962925561|0.0178573100151|
|joint_seed2_narrow_fullj / joint8|0.797717012547|0.798716101969|0.00567462310677|0.0152014532531|

同一保留材料对在无目标校准条件与平移目标joint条件下具有不同数据距离：不能把joint小于1写成cal小于1。两网格差是经验一致性检查，不是连续PDE误差上界。实际复数预测差见[coarse/fine图](figures/saved_results/B2_joint_predictions_joint_seed2_narrow_fullj.png)，来自joint8_PAIR四数组；没有重新求解或重新拟合。

## 4. 任务类进度只按已完成记录公示

* cavity_seed2_narrow_bestcal：当前保存24/24个probe结果；已完成子集最大距离1.37906100106。24项保存齐全；资格解释仍归主线程。
* joint_seed2_narrow_fullj：当前保存24/24个probe结果；已完成子集最大距离4.81819338104。24项保存齐全；资格解释仍归主线程。

原始来源：[native cal](results/B2/native_cal_2026102202_narrow/RESULTS.json)、[native joint](results/B2/native_joint_direct_2026102202_narrow/RESULTS.json)、[GRID_COMPARISON](results/B2/validation/joint_seed2_narrow_fullj/GRID_COMPARISON.json)、[ALL_TASKS](results/B2/validation/joint_seed2_narrow_fullj/ALL_TASKS.json)。其他seed/分支不由本seed2检查点推定完成；不据Oracle候选数量产生概率。

生成UTC：2026-10-23T18:49:24.318228+00:00。本页对应已封账研究。数值精度、科学门槛和适用先验分别列示；报告发布本身不等于科学通过。 不修改科学门、物理源码或账本。
