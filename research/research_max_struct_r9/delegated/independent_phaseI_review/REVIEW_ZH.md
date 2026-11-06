# PhaseI独立证据复核

结论：在本次限定范围内，**未发现推翻三个关键表述的致命错误**。发现两处容易误读的记号/费用措辞，建议root在最终审阅时澄清；没有修改报告、运行算法或原始回执。核心数学有效性与最终科学判断仍由root负责。

本复核只读已完成PhaseI报告§1–7、THEORY、确认评分/汇总、失败诊断和必要代码/费用元数据。没有打开任何私有GT数组，没有新forward、重评分、3D自审或在途classic/original设置/GT访问。可复核算术和48行配对表见`EVIDENCE.json`、`PAIRED_SUCCESS_MATRIX.csv`。

## 1. 未发现问题：四类方法2/6→4/6且无A+独有成功

从`RESULTS/confirmation_scoring/SCORES.json`的66条原始评分中筛选四方法×两规则×六例，共48条；每格均有六个不同case。直接重算成功集合、平均NRMSE和平均CPU，与`RESULTS/summary/CONFIRMATION_SUMMARY.json`一致至1e−12：

|方法|fixed_001成功case|noise_ball_99成功case|
|---|---|---|
|single CIE|0000,0001|0000,0001,0003,0005|
|Dual CIE|0000,0001|0000,0001,0003,0005|
|Global Dual|0000,0001|0000,0001,0003,0005|
|A+|0000,0001|0000,0001,0003,0005|

A+独有成功集合为空。这里材料成功使用已冻结的amplitude-normalized NRMSE≤.15、prediction-only Otsu IoU≥.7及box合法；实际四个修正成功对象也都通过独立CN99噪声检查，故此处没有因二种gate区别改变4/6计数。报告关于“六个对象而非66个独立样本”的分母限制正确。A+相对修正single-CIE平均CPU约高19.4%，报告“约19%”准确；不能由六例推总体非劣或显著提升，报告已保留该限制。

原文位置：`RESEARCH_REPORT_ZH.md:7`、`:23`、`:72`–`:96`。证据：`RESULTS/confirmation_scoring/SCORES.json`、`RESULTS/summary/CONFIRMATION_SUMMARY.json`与本目录配对CSV。

## 2. 未发现重复收费：40ac弱阶段交接

原始已完成评分给出NRMSE **.0583512451458**、IoU **.969230769231**。成本逐项为：

|项|CPU秒|来源/范围|
|---|---:|---|
|本次新物理continuation|16.352524|weak-handoff原RECEIPT的process_CPU_seconds|
|原profile启动前缀|.452537|profile_homotopy的pre_solver_CPU_seconds|
|已付低罚权首阶段前缀|2.828670|BUDGET唯一profile_0_0.01_endpoint|
|方法standalone合计|19.633731|16.352524 + .452537 + 2.828670|
|原完整profile实际运行|28.096411|包含原首阶段，已独立计入实际研究总账|

`score_weak_handoff.py:31`–`:48`取精确已存前缀，并校验public初值等于该阶段的保存材料；`:71`后把它加到standalone方法成本。`audit/ALL_COSTS.json`的实际inverse_receipts对本次只收费16.352524，19.633731放在明确的`standalone_weak_recipes_NOT_ADDITIONAL_ACTUAL_COST`视图，未把3.281207再加进实验总账。standalone同时补计5个physical LU/108 complex RHS与105次current Cholesky/1260 RHS，不只补秒数。

**措辞问题W2（非运行错误）：**`RESEARCH_REPORT_ZH.md:64`的“独立成本19.63CPU秒”可能被误读为一次从冷启动完整重跑并独立计时19.63秒。实际是测得的新continuation加已付首阶段前缀的**standalone方法归属成本**。建议写成“含已付初始化前缀的standalone成本19.6337秒；本次新增执行16.3525秒”。不应称这条开发控制为未暴露确认对象上的独立泛化证据。

## 3. 未发现问题：GT可行Phi见证的逻辑边界

诊断六例均标记GT合法、同一24×24离散、继承同G/E/S。见证由已存Y_clean对public Y的数据项加相同prior构成；无新全波调用。六例`GT_data_term + GT_prior = GT_feasible_Phi`的元数据算术一致。两共同失败的最小端点比值为：

|共同失败对象|GT可行Phi|最小已存端点Phi|最小端点/GT|
|---|---:|---:|---:|
|annulus0002|.0004452215931|.004522012838|10.1567689|
|star0004|.00005638817101|.08834832937|1566.788349|

这排除的只是“这些当前端点没有更低已知可行点”/“已达到这个见证值”的辩解；不证明GT全局最优、唯一性、非唯一性、全临界锥局部极小、整个PDE的盆地结构或新算法保证。`RESEARCH_REPORT_ZH.md:104`与诊断`qualification`都正确保留这一限制。

还需保留selection区别：失败诊断的optimizer stop通常属于corrector最后iterate，未必就是最终按最小原Phi选出的incumbent；诊断明确说明没有对每个被选incumbent重算梯度。不能把REPORTED_PGTOL自动升级成该保存预测的驻点或局部极小证书。报告目前没有这样升级，因此这不是发现的错误。

## 4. 推导记号问题W1：局部Goff应明确

`THEORY_DERIVATIONS.md:15`–`:19`写`R_A=β(I+D_z G_off)(J−D_z H)`，但首次出现未说明`G_off`是局部算子。实际`research_max_struct_r4/src/cie_onehop.py:26`–`:35`只抽取不相交2×2空间块的4×4子矩阵，并把每块对角置零；R8/R9 onehop读出调用这一实现。

这不是数值运行错误，但把它理解成完整物理`G−diag(G)`会改变空间耦合范围、费用和方法含义。建议root在首次出现时明确`G_off=⊕_B(G_BB−diag G_BB)`或采用`H_local`记号，并显式定义`Q=E+GJ`。完整物理forward中的G仍然是全算子；不能把二者混写。此建议仅为合同/记号澄清，不替root判定核心理论。

## 5. 另外核对的分母与prior

报告§2的NRMSE是逐像素RMSE除声明contrast幅度，不是`||mhat−m||/||m||`。例如0000 single-CIE fixed的RMSE=.1946453154、幅度2.5，报告NRMSE=.0778581262；同记录的relative-L2为.2575760273。两数不同且均正确；表格外引时需保留“amplitude-normalized”前缀。未发现报告把二者互换。

`score_primary.fixed_prior`使用24×24相邻差除19.5，smooth=TV=1e−4、平滑epsilon=.01，与实际MaxwellModel继承的Regularizer一致；19.5是冻结prior尺度，不是每例contrast幅度。66条确认预测从已存Y重建的Phi与回执最大差约2.78e−17。这里只核对已有结果和代码，没有重跑求解、重算GT或改变prior。

## 分类与限制

- 致命错误：本次范围内未发现。
- 措辞/记号：W1局部Goff未明确、W2 standalone费用不等于独立冷运行计时。
- 核对无问题：2/6→4/6与无独有成功、19.633731费用归属且不重复收费、可行见证不推唯一性。
- 未做：总体统计显著性、核心定理的最终裁决、全临界锥证明、在途classic/original或任何3D结果审计。
