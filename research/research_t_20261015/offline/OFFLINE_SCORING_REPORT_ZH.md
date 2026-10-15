# Research T：最终离线材料评分

状态：FINAL_FROZEN_ONLINE_PHASES_OFFLINE_SCORE。读取原已暴露A30真值，未读取locked truth，未运行新forward/adjoint/tangent；新物理/RHS均为0。所有MAE只在独立离线进程产生，不进入候选选择或在线优化。字段、数组与图像可复现于本目录。

## 基础解：数据相容不等于材料准确

|规则/端点|数据与数值合格|平均白化残差|complex MAE|Re MAE|Im MAE|
|---|---:|---:|---:|---:|---:|
|original_N3|2/8|29.139042|0.041756|0.039688|0.009997|
|full_final|8/8|20.148396|0.167560|0.141362|0.062463|
|multiscale_haar_final|8/8|21.162010|0.054811|0.048762|0.017934|
|qualified_base|8/8|20.136726|0.147656|0.123295|0.057465|
|qualified_minchange_base|8/8|21.927005|0.050757|0.045796|0.015730|

原最低观测残差选择改善相容性，但明显增加full-grid材料误差。随后最小N3位移规则是受这轮exposed离线汇总启发的**post-hoc算法修订**，不是预注册盲验证。该规则在每个case仅用已合格端点与N3间的材料距离，不用真值MAE选择。它减少了材料偏离，仍不能宣称比原N3更接近真值。

## 有限编辑与补偿：逐轮保留

|轮次|显著合格编辑|不合格→合格|raw平均残差|实际选中端点平均残差|全部端点complex MAE|core外平均相变化细胞|
|---|---:|---:|---:|---:|---:|---:|
|mechanism_v1|15/16|15|331.962875|22.293967|0.504148|132.06|
|mechanism_v2|10/16|10|331.962875|23.240037|0.506147|113.69|
|mechanism_v3|7/8|7|321.953730|21.919004|0.495841|133.75|

“实际端点”优先使用保存的best_valid；没有合法端点时保留final失败端点。全部端点均值没有悄悄去除失败；合法子集均值另在AGGREGATE_METRICS.json。raw→corrected通过数是物理相容与保护结构的机制见证，**不是**重建真值、posterior mass、独立材料族数、或搜索方法优越性。

v1与v2采用同一最低残差base和相同raw编辑；v3换成minchange base，只比较固定编辑indices0/2，必须与原base分开。MATCHED_ITERATION_COMPARISON保留24个配对，包含base身份及raw residual差。跨轮32个通过记录不能直接当32个已去重材料世界。

## 图、掩膜和边界

所有epsilon实部共用[.5,5]、损耗共用[0,.5]。原增量support与绝对近空气掩膜分开：primary近空气为|chi|≤.1，另列.075/.125敏感性，严格真空为|chi|≤1e-10，不能把这些定义混用。空气IoU在两空集时取1，并附empty_both标志，不能因此宣称发现空腔。边界距离是24×24像素中心栅格量，不是连续几何或网格收敛证书；任一边界空集时为null。

图库index.html含8case基础反演、40条机制编辑、实际best_valid、全8源散射拟合、保护core/guard与材料差别；原始full-material/prediction保存在OFFLINE_MATERIAL_AND_PREDICTION_ARRAYS.npz。真值图仅离线展示。

## 费用与失败

物理求解费用由父阶段唯一ledger计费，离线脚本不写该ledger。独立零物理后处理费用存EXECUTION_RECEIPTS.jsonl。首次最终刷新在图已生成后出现重复字段TypeError，FAILED_REFRESH_RECEIPT.json保留异常；其wall未被外层计时，明示null/UNRECORDED，不冒写0。修复后复用同一冻结输入已生成的图，重算全部指标，无新物理。所有概率字段为null/candidate_only。
