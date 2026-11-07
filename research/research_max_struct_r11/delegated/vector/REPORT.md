# R11 真实三维有限确认：完成报告

本轮50个冻结反演全部完成：8个独立材料×5方法=40个matched12任务，另2个预选材料的18→12跨网格观测×5方法=10个任务。两组分开统计；交叉观测不增加独立材料数。GT在50个预测与400个公开cost-prefix逐项封存并经root独立审查后才离线评分。未新增物理调用、算法、对象、预算或3D成功门。

8例matched的材料主轮廓通常可辨认；这与冻结IoU数值、分相幅值是否准确是不同问题。full-r120相对普通Dual在5/8例的χ相对L2较小，但八例平均仍略高于Dual；有限预算下总体接近，不应说所有方法“未成像”，也不足以据此把Green分母设为新的普适默认。局部改善与退化均保留。

full-r120与contracted-r120仅改材料块，current始终ordinary CIE；二者使用同一公开背景源基U=orth([Einc,G Einc])、rank120。它们不是二维G19、不是BothBlocks，也不声称rank120覆盖全dyadic G。full指材料依赖分母，非全秩G。full/contracted与Global匹配global-PR current，三个阶段的β_current均为[6,4,3]；Dual采用sourcewise-PR且β_current同为[6,4,3]，CIE采用sourcewise-PR且β_current=[6,3,1]。所有方法β_material=[6,3,1]。因此full对Global是材料块匹配对照；对Dual的差异还包含PR聚合方式，不能只归因于材料分母。

## 冻结合同与指标

12³、1728个实各向同性χ、5184×60复电流、180×60观测；CM+RR dyadic物理不变，FP64/complex128，m0=.25、box0..8、R0、3%独立圆复高斯噪声。MF2/3/4各max500/min20，corrector320 accepted/800requests；50/50最终均达到320 accepted上限，PG∞仍高于1e−8，不标作收敛。

χ-relative L2与amplitude-normalized RMSE同时列出，避免空背景分母美化误差。fixed IoU是原先固定阈值的连续重叠指标，不补设3D成败阈值。Phi3是source-normalized拟合目标，CN99是另一个按σ定义的全观测噪声球一致性指标；CN99内不等于材料正确。

## 8个matched材料终点

|方法|mean χ-relative L2|median χ-relative L2|mean amp-NRMSE|mean fixed IoU|CN99内|
|---|---:|---:|---:|---:|---:|
|CIE|0.228436|0.180327|0.052930|0.808690|8/8|
|Dual|0.226079|0.180061|0.052185|0.812481|8/8|
|Global|0.231274|0.189406|0.053312|0.803822|8/8|
|full-r120|0.227803|0.184074|0.052675|0.813349|8/8|
|contracted-r120|0.227738|0.181933|0.052625|0.812590|8/8|

下面列全部8例，不作GT选择删例。

|对象|几何|bulk对比度|主机|
|---|---|---|---|
|c11 high|非球形包覆|[4, 8]|MAC_CPU|
|c12 mixed|非球形包覆|[1, 8]|WINDOWS_CUDA|
|c13 high|开口矩形腔|[8]|WINDOWS_CUDA|
|c14 mixed|开口矩形腔|[1, 8]|MAC_CPU|
|c15 high|相邻分离双体|[8]|MAC_CPU|
|c16 mixed|相邻分离双体|[1, 8]|WINDOWS_CUDA|
|c17 high|有颈桥|[8]|WINDOWS_CUDA|
|c18 mixed|有颈桥|[1, 8]|MAC_CPU|

### χ-relative L2（越低越好）

|对象|CIE|Dual|Global|full-r120|contracted-r120|
|---|---:|---:|---:|---:|---:|
|c11|0.1212|0.1164|0.1184|0.1211|0.1223|
|c12|0.5580|0.5601|0.5633|0.5598|0.5627|
|c13|0.1122|0.1116|0.1105|0.1091|0.1074|
|c14|0.1497|0.1522|0.1449|0.1466|0.1518|
|c15|0.1793|0.1785|0.2083|0.1928|0.1846|
|c16|0.1813|0.1816|0.1857|0.1753|0.1793|
|c17|0.3223|0.3161|0.3259|0.3152|0.3195|
|c18|0.2036|0.1920|0.1931|0.2024|0.1943|

### fixed IoU（连续值，越高越好）

|对象|CIE|Dual|Global|full-r120|contracted-r120|
|---|---:|---:|---:|---:|---:|
|c11|0.9252|0.9307|0.9202|0.9164|0.9191|
|c12|0.6489|0.6424|0.6359|0.6506|0.6382|
|c13|0.8836|0.8886|0.8896|0.8913|0.8970|
|c14|0.7910|0.7894|0.7917|0.7970|0.7868|
|c15|0.8945|0.8977|0.8676|0.9070|0.9028|
|c16|0.7500|0.7552|0.7368|0.7371|0.7539|
|c17|0.7244|0.7302|0.7323|0.7294|0.7412|
|c18|0.8519|0.8657|0.8565|0.8779|0.8618|

full相对Dual的局部χ误差改善出现在c12/c13/c14/c16/c17；c11/c15/c18方向相反。c16的χ误差full更小，但IoU低于Dual；c18的full IoU较高，χ误差反而更大。不同指标不能互相替代。c13开腔、c15双体、c17桥在固定切片/体图仍有可辨认主结构；c12混合包覆的内高值区明显展宽、幅度偏低，五方法χ误差都约.56。这里是具体误差描述，不加新的二元标签。

## 分相偏差

直接复制原评分表的分箱：background χ=0，low 0<χ≤1，middle 1<χ≤3，high 3<χ≤8。边界/有效体素也可能落入这些区间；它们不自动等于设计的bulk相。c11的4与8都归high，不能声称已经单独测得χ=4相偏差。空分箱保留缺失，不补0。全部50×4分箱见[PHASE_BIAS_AND_RMSE.csv](completed_analysis/PHASE_BIAS_AND_RMSE.csv)。

以下为matched低分箱的bias（估计−真值）；完整background/middle/high bias与RMSE在原表中。

|对象|CIE|Dual|Global|full-r120|contracted-r120|
|---|---:|---:|---:|---:|---:|
|c11|-0.0593|-0.0573|-0.0638|-0.0526|-0.0508|
|c12|-0.0793|-0.0821|-0.0797|-0.0804|-0.0792|
|c13|0.0167|0.0077|0.0221|0.0027|0.0029|
|c14|-0.0200|-0.0291|-0.0238|-0.0053|-0.0181|
|c15|0.0081|-0.0127|0.0042|0.0045|-0.0208|
|c16|0.0187|0.0142|-0.0013|0.0034|0.0078|
|c17|0.3244|0.3149|0.3340|0.2879|0.3078|
|c18|-0.1034|-0.0501|-0.0516|0.0120|-0.1071|

## 同成本前缀与工程敏感性

400个原先封存的选择均为已付状态、公开minPhi，未按GT重选。c17 Dual原attempt0在315 accepted后因nvidia-smi查询超时丢进程，额外320个5184-LU已计；独立attempt1从同初值按同320/800完整重跑。该臂重跑prefix明确是attempt-only，不能当端到端成本。主材料表仍保留8例；成本敏感性整例排除c17，使用7例共同对象集。这个工程排除不是由GT决定，也不替代原8例。

以下是7例共同集的mean χ-relative L2。CPU/wall是每对象在其固定主机上的预算视图，不把Mac秒与Windows秒相加称算法速度；LU为同5184阶、all60物理factor预算，额外120阶材料分母工作另有记账。

|预算|CIE|Dual|Global|full-r120|contracted-r120|
|---|---:|---:|---:|---:|---:|
|CPU_seconds 480|0.272986|0.275537|0.262353|0.292668|0.282006|
|CPU_seconds 960|0.224142|0.236028|0.219623|0.240806|0.219703|
|wall_seconds 480|0.274544|0.277612|0.263768|0.294243|0.284299|
|wall_seconds 960|0.225699|0.241099|0.220379|0.242764|0.220561|
|physical_LU 80|0.329389|0.323762|0.325901|0.329138|0.329128|
|physical_LU 160|0.244364|0.245251|0.247621|0.250212|0.245178|
|physical_LU 240|0.220065|0.215676|0.220529|0.222754|0.217194|
|physical_LU 320|0.214921|0.213179|0.217818|0.215365|0.214734|

在这个共同7例视图，full在480/960CPU和wall预算下均值较Global更差；contracted在960CPU/wall附近接近Global；320-LU视图的差异仍小，Dual均值最低。它们是冻结前缀描述，不是新停止规则或显著性结论。全部8例attempt视图、7例共同集和2cross分开的精确数值见[COST_PREFIX_SUMMARY.csv](completed_analysis/COST_PREFIX_SUMMARY.csv)。

## 2个跨网格观测单列

|对象|clean18与clean12相对差|CN99|解释范围|
|---|---:|---|---|
|c13_cross18to12|17.6827%|五方法均在球外|固定12³反演对18³生成的模型/离散风险，不统称分母失败|
|c18_cross18to12|3.5287%|五方法均在球外|固定12³反演对18³生成的模型/离散风险，不统称分母失败|

|方法|mean χ-relative L2|mean amp-NRMSE|mean fixed IoU|mean CN99 ratio|
|---|---:|---:|---:|---:|
|CIE|0.464833|0.158622|0.752628|1.157470|
|Dual|0.454391|0.155923|0.753932|1.155290|
|Global|0.433562|0.150384|0.762925|1.159018|
|full-r120|0.467138|0.158298|0.763631|1.157607|
|contracted-r120|0.463195|0.157487|0.762708|1.156837|

跨格Global的两例平均χ误差最低，Green两臂没有改善这个平均；不能将这两例并入matched成功比例。c13的17.68% clean forward差远高于3%噪声，c18的3.53%也非零；观测模型差、有限预算与算法误差在这里并存。c18跨格图出现亮颈/空心化等材料分配偏差，IoU相近仍可能有显著χ误差。

## 成本与失败保留

统一总账[COSTS.json](COSTS.json)覆盖本worker负责的全部vector生产者：公开basis/导数/固定状态、DEV12、一次6LU库选择、10观测生成、50确认、原失败与重跑、controller/helper/transport、评分与图、这份机械汇总；不是仅fresh。root环境清单、root理论/发布与2D生产者明确排除。

已知CPU下界按平台分开：Mac 37715.173199 s，Windows 26000.843750 s。全vector 5184阶物理LU20507次；17496阶生成为2个all60逻辑批/120RHS/6392 GMRES迭代，未从RHS别名重复相加。小32/120阶分母LU独立列示，不混入5184物理LU。wall不可跨并行lane相加为总历时，CUDA event区间也不是活跃kernel时间。

原GPU failed phase包含15个已封存臂+第16臂315步，替代这些嵌套arm计时；独立重跑phase额外计费。原315目录归archived attempt0，sealed重跑attempt1入canonical，下一臂mkdir sentinel受控0FW停止，后9新进程用30s有效空闲资源快照缓存。查询等待改善仅为工程资源协作，不算算法速度。

查询、导入/退出、部分旧remote pack/extract等未测量部分在null gaps中保留。所有历史日志/非零exit都保留，受控controller drain/c18 waiting取消不误列科学失败。Mac20个临时spool均在完整原数组归档后删除，当前临时目录为空。[ATTEMPT_PATH_REGISTRY.json](ATTEMPT_PATH_REGISTRY.json)给出路径语义。

## 本轮收束建议

保留三维ordinary基线与两种rank120材料分母的局部差异证据；本轮平均值与同成本前缀不支持把Green材料分母直接替换成统一默认。跨网格观测按明确模型/离散风险单列。这个判断只覆盖规定的8个材料与冻结预算，不推广到未知总体、全G近似、BothBlocks或收敛极限。本轮到此结束，不追加算法/对象/物理试验。理论与原创性及综合论文结论由root审查。

## 全30图入口

每个体图为真3D体素外表面χ>.5、统一χ=0..8与相机；内部材料分相须结合固定中心正交切片和区域误差，体图遮挡不证明内核准确。所有8matched与2cross、五方法及GT均保留。
- [c11_matched12_orthogonal](confirmation_postseal_scores/figures/c11_matched12_orthogonal.png)
- [c11_matched12_volume](confirmation_postseal_scores/figures/c11_matched12_volume.png)
- [c12_matched12_orthogonal](confirmation_postseal_scores/figures/c12_matched12_orthogonal.png)
- [c12_matched12_volume](confirmation_postseal_scores/figures/c12_matched12_volume.png)
- [c13_cross18to12_orthogonal](confirmation_postseal_scores/figures/c13_cross18to12_orthogonal.png)
- [c13_cross18to12_volume](confirmation_postseal_scores/figures/c13_cross18to12_volume.png)
- [c13_matched12_orthogonal](confirmation_postseal_scores/figures/c13_matched12_orthogonal.png)
- [c13_matched12_volume](confirmation_postseal_scores/figures/c13_matched12_volume.png)
- [c14_matched12_orthogonal](confirmation_postseal_scores/figures/c14_matched12_orthogonal.png)
- [c14_matched12_volume](confirmation_postseal_scores/figures/c14_matched12_volume.png)
- [c15_matched12_orthogonal](confirmation_postseal_scores/figures/c15_matched12_orthogonal.png)
- [c15_matched12_volume](confirmation_postseal_scores/figures/c15_matched12_volume.png)
- [c16_matched12_orthogonal](confirmation_postseal_scores/figures/c16_matched12_orthogonal.png)
- [c16_matched12_volume](confirmation_postseal_scores/figures/c16_matched12_volume.png)
- [c17_matched12_orthogonal](confirmation_postseal_scores/figures/c17_matched12_orthogonal.png)
- [c17_matched12_volume](confirmation_postseal_scores/figures/c17_matched12_volume.png)
- [c18_cross18to12_orthogonal](confirmation_postseal_scores/figures/c18_cross18to12_orthogonal.png)
- [c18_cross18to12_volume](confirmation_postseal_scores/figures/c18_cross18to12_volume.png)
- [c18_matched12_orthogonal](confirmation_postseal_scores/figures/c18_matched12_orthogonal.png)
- [c18_matched12_volume](confirmation_postseal_scores/figures/c18_matched12_volume.png)
- [cross18to12_Phi3_CPU_seconds](confirmation_postseal_scores/figures/cross18to12_Phi3_CPU_seconds.png)
- [cross18to12_Phi3_and_sigma_CN99_separate](confirmation_postseal_scores/figures/cross18to12_Phi3_and_sigma_CN99_separate.png)
- [cross18to12_Phi3_physical_LU](confirmation_postseal_scores/figures/cross18to12_Phi3_physical_LU.png)
- [cross18to12_Phi3_wall_seconds](confirmation_postseal_scores/figures/cross18to12_Phi3_wall_seconds.png)
- [cross18to12_all_material_metrics](confirmation_postseal_scores/figures/cross18to12_all_material_metrics.png)
- [matched12_Phi3_CPU_seconds](confirmation_postseal_scores/figures/matched12_Phi3_CPU_seconds.png)
- [matched12_Phi3_and_sigma_CN99_separate](confirmation_postseal_scores/figures/matched12_Phi3_and_sigma_CN99_separate.png)
- [matched12_Phi3_physical_LU](confirmation_postseal_scores/figures/matched12_Phi3_physical_LU.png)
- [matched12_Phi3_wall_seconds](confirmation_postseal_scores/figures/matched12_Phi3_wall_seconds.png)
- [matched12_all_material_metrics](confirmation_postseal_scores/figures/matched12_all_material_metrics.png)

原始评分[主表CSV](confirmation_postseal_scores/MATERIAL_AND_PHYSICS.csv)、[400前缀CSV](confirmation_postseal_scores/COST_PREFIX_MATERIAL_SCORES.csv)、[失败/预算类](confirmation_postseal_scores/FAILURE_CLASSES.csv)、[public SELECTIONS](confirmation_host_union/SELECTIONS.json)、[完成复现说明](COMPLETED_REPRODUCE.md)。
