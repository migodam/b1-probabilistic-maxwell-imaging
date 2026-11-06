# 冻结结果数字备忘录

只摘录既有评分/汇总/置信区间；未重跑评分、实验、训练或bootstrap，未读新GT，也未覆盖图。精确原表行号、所有硬件层及逐对象公开条件见同名JSON。以下均不作科学胜负判断。

原主比较630=42×15；额外120仅固定8对象（ID4/OOD2/cross2）；anchor126=42×3独立扩展。876个回执均为SEALED_PUBLIC_PREDICTION，不能据此称材料恢复成功。

## 最终材料：均值NRMSE；继承复合成功数/分母

|方法|ID24|OOD12|cross6|
|---|---:|---:|---:|
|physical_lbfgs / uniform|0.24966; 2/24|0.40474; 1/12|0.19593; 2/6|
|physical_lbfgs / nn|0.09647; 9/24|0.16245; 6/12|0.08431; 2/6|
|cie_dual / uniform|0.36843; 7/24|0.63825; 1/12|0.32821; 2/6|
|cie_dual / nn|0.32822; 7/24|0.56779; 2/12|0.32835; 2/6|
|global_dual / uniform|0.35689; 8/24|0.62538; 2/12|0.32874; 2/6|
|global_dual / nn|0.28919; 8/24|0.60481; 1/12|0.32834; 2/6|
|a_plus / uniform|0.34062; 8/24|0.63748; 1/12|0.32833; 2/6|
|a_plus / nn|0.29609; 7/24|0.59393; 2/12|0.32933; 2/6|
|a_plus_sourcewise / uniform|0.34368; 8/24|0.63096; 2/12|0.32831; 2/6|
|a_plus_sourcewise / nn|0.32182; 7/24|0.59164; 2/12|0.32811; 2/6|
|cie_dual / nn_affine_anchor|0.29007; 7/24|0.52299; 4/12|0.25988; 2/6|
|a_plus / nn_affine_anchor|0.31395; 5/24|0.51831; 4/12|0.22800; 2/6|
|a_plus_sourcewise / nn_affine_anchor|0.31293; 6/24|0.50160; 4/12|0.22833; 2/6|

NRMSE/IoU定义及成功门槛完全继承；表不合并三个总体。源：PRIMARY_630_cohort_summary.csv和CONFIRMATION_cohort_material.csv；JSON逐单元保留原CSV物理行号。

## NNdirect vs uniform：相同原物理算法、同硬件实际CPU

|cohort/硬件(n)|solver CPU均值 uniform→NN|warm CPU均值 uniform→NN|warm NN/uniform|
|---|---:|---:|---:|
|cross6/Mac(6)|16.719→30.821|16.719→30.823|1.844|
|OOD12/Mac(6)|19.650→39.122|19.650→39.123|1.991|
|OOD12/XINAN(6)|46.852→128.508|46.852→128.519|2.743|
|ID24/Mac(10)|30.213→32.032|30.213→32.034|1.060|
|ID24/XINAN(14)|51.999→69.757|51.999→69.765|1.342|

CPU比为配对同对象的**均值之比**。原feature/inference字段为Mac实际计时；deployment_prefix为独立同机计时；actual_workflow可包含Mac预处理+XINAN求解。所有25个核心NN/uniform硬件层均值/比值见JSON，未混硬件平均。

## NNdirect相同at-most预算端点：NRMSE差 NN−uniform [95%区间]

|cohort/硬件(n)|CPU allowance|LU allowance|RHS allowance|
|---|---|---|---|
|cross6/Mac(6)|-0.1132 [-0.2168,-0.0285]|-0.1128 [-0.2185,-0.0287]|-0.1128 [-0.2185,-0.0287]|
|OOD12/Mac(6)|-0.2284 [-0.3328,-0.1281]|-0.2271 [-0.3301,-0.1275]|-0.2271 [-0.3301,-0.1272]|
|OOD12/XINAN(6)|-0.2631 [-0.4720,-0.0844]|-0.2627 [-0.4711,-0.0845]|-0.2627 [-0.4710,-0.0843]|
|ID24/Mac(10)|-0.2331 [-0.3620,-0.1152]|-0.2359 [-0.3634,-0.1190]|-0.2359 [-0.3641,-0.1158]|
|ID24/XINAN(14)|-0.0818 [-0.1395,-0.0303]|-0.0819 [-0.1392,-0.0299]|-0.0819 [-0.1399,-0.0297]|

预算是每对象min(两方法最终实付成本)，按该预算内原Phi incumbent取已保存状态；停止方法不续算。不是精确等支出；unused和checkpoint间距来自现有端点表。所有原NN/均匀及anchor三方法对原NN/NNdirect/global nn_joint的CPU/LU/RHS区间已逐行摘录JSON。例：ID/XINAN global-anchor−globalNN在CPU预算NRMSE差+0.11044 [+.00077,+.23583]；同组anchor−nn_joint为−.00828 [−.07402,+.06107]。不把区间包装成显著性或普遍机制。

## 原physics网络训练摊销

固定分子402.558889 Mac CPU秒=生成1536训练11.396809 + 192验证3.472078 + 全120epoch训练387.690002；生成总体已含feature/setup，未重复加。Born/domain费用不计入此部署分子。NNdirect五个cohort/硬件层平均CPU节省全≤0，均无有限CPU摊销优势。其他支路的全部25条均值摊销记录在JSON；Mac为CPU记账，XINAN仅给(cMac/cXINAN)乘系数，未指定价比不能称经济盈亏平衡。

## 两完整支路portfolio：最终均值NRMSE；成功数

|NN physical+辅助支路|ID24|OOD12|cross6|
|---|---:|---:|---:|
|cie_single|0.11768; 11/24|0.37895; 6/12|0.22790; 2/6|
|cie_dual|0.11653; 11/24|0.30903; 6/12|0.22811; 2/6|
|global_dual|0.07726; 13/24|0.34577; 6/12|0.22809; 2/6|
|a_plus|0.07754; 13/24|0.34498; 6/12|0.22882; 2/6|
|a_plus_sourcewise|0.07744; 12/24|0.34531; 6/12|0.22789; 2/6|

例如global A+ portfolio相对单NNdirect均值差ID −.01893、OOD +.18253、cross +.14450；相对等候选数uniform portfolio均值差ID −.17071、OOD −.29260、cross −.09927。不是GT最优候选，而是已封存原Phi选择。全部5种组合的材料差和旧CI在JSON。

|global A+ portfolio硬件层|NNdirect warm CPU|NN两支路warm CPU|uniform两支路warm CPU|NN两支路/NNdirect|
|---|---:|---:|---:|---:|
|cross6/Mac|30.823|97.659|83.235|3.168|
|OOD12/Mac|39.123|113.755|100.423|2.908|
|OOD12/XINAN|128.519|306.438|230.526|2.384|
|ID24/Mac|32.034|95.933|97.554|2.995|
|ID24/XINAN|69.765|214.986|220.811|3.082|

两条完整solver成本均支付，共享NN prefix只计一次。没有冻结部分支路执行调度，因此没有等预算portfolio轨迹。

## 固定8对象初始化/Born-only：仅global A+示例，全部3支路×初始化在JSON

|初值/网络|ID4 均值;成功|OOD2 均值;成功|cross2 均值;成功|
|---|---:|---:|---:|
|uniform|0.24599; 3/4|0.88100; 0/2|0.07030; 2/2|
|uniform1|0.24359; 3/4|0.87086; 0/2|0.06886; 2/2|
|born|1.23642; 2/4|0.87527; 0/2|2.08652; 1/2|
|dpcs|0.24387; 3/4|0.86759; 0/2|0.06944; 2/2|
|random|0.24410; 3/4|0.87674; 0/2|0.07065; 2/2|
|nn|0.24427; 3/4|0.64886; 0/2|0.07222; 2/2|
|nn_joint|0.46267; 1/4|0.89184; 0/2|0.07190; 2/2|
|nn_born_features|0.06421; 4/4|0.87823; 0/2|0.07117; 2/2|

Born-only网络是另一个训练控制，不能据ID4的结果宣称42对象优势。解析Born和Born-only NN是不同初值。

## 初值进入corrector及“更低prePhi、但finalPhi与材料皆劣于NNdirect”

|方法|ID stage−1 /双劣|OOD stage−1 /双劣|cross stage−1 /双劣|
|---|---:|---:|---:|
|cie_dual/nn|4/24 ; 10/24|3/12 ; 5/12|0/6 ; 2/6|
|global_dual/nn|5/24 ; 11/24|2/12 ; 6/12|0/6 ; 1/6|
|a_plus/nn|4/24 ; 10/24|2/12 ; 5/12|0/6 ; 2/6|
|a_plus_sourcewise/nn|3/24 ; 10/24|2/12 ; 5/12|0/6 ; 2/6|
|a_plus/nn_joint|3/24 ; 11/24|1/12 ; 6/12|0/6 ; 2/6|
|cie_dual/nn_affine_anchor|5/24 ; 6/24|2/12 ; 6/12|0/6 ; 2/6|
|a_plus/nn_affine_anchor|4/24 ; 7/24|2/12 ; 5/12|0/6 ; 1/6|
|a_plus_sourcewise/nn_affine_anchor|4/24 ; 8/24|2/12 ; 4/12|0/6 ; 1/6|

stage−1只指初始候选进入校正，非最终图像保持不变。“双劣”为严格符号比较，会包含浮点微差，不能称每一计数都是实质basin失败。NNdirect上限6000 vs辅助corrector3000也使计数不是统一固定C映射的证明；不另设差值阈值。

### 更具体的继承成功标签分歧：direct成功/辅助失败；反向增益

|辅助方法|ID24 丢失/新增|OOD12 丢失/新增|cross6 丢失/新增|
|---|---:|---:|---:|
|cie_dual / nn|6/4|4/0|0/0|
|global_dual / nn|6/5|5/0|0/0|
|a_plus / nn|6/4|4/0|0/0|
|a_plus_sourcewise / nn|6/4|4/0|0/0|
|a_plus / nn_joint|8/2|5/0|0/0|
|cie_dual / nn_affine_anchor|5/3|2/0|0/0|
|a_plus / nn_affine_anchor|6/2|2/0|0/0|
|a_plus_sourcewise / nn_affine_anchor|6/3|2/0|0/0|

每个丢失/新增对象ID集合已列入JSON；完全继承原材料成功标签，无新门槛。

### 预定对象0027：停止原因与实际accepted次数

|方法|prePhi→finalPhi|NRMSE / IoU|accepted / evaluations|corrector stop|
|---|---:|---:|---:|---|
|cie_dual / nn|0.0251650727→0.000462581142|1.013703 / 0.235795|3000 / 3053|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|a_plus / nn|0.0156508798→0.000463503877|1.013697 / 0.228412|3000 / 3054|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|physical_lbfgs / nn|0.309652839→0.000428658662|0.120818 / 0.849673|1814 / 1834|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|a_plus / nn_joint|0.153560364→0.000494317379|1.064779 / 0.224090|3000 / 3052|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|a_plus_sourcewise / nn|0.01425595→0.000463775596|1.015332 / 0.226519|3000 / 3040|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|global_dual / nn|0.0367780764→0.000462607469|1.013547 / 0.237288|3000 / 3041|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|cie_dual / nn_affine_anchor|0.083849416→0.000492833808|1.064148 / 0.232493|2789 / 2825|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|a_plus / nn_affine_anchor|0.0915991718→0.000492844495|1.064394 / 0.226891|3000 / 3033|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|a_plus_sourcewise / nn_affine_anchor|0.160004199→0.000493941242|1.064775 / 0.219444|3000 / 3047|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|

0027的NNdirect在1814 accepted步按PGTOL正常停止（低于3000）；原四个NN辅助方法均达到3000 accepted上限且材料标签失败。它比微差符号计数更具体，但辅助若继续迭代会怎样仍未运行；不能据此断言不可逃离盆地。

## NNdirect未达继承材料复合门槛：公开条件及组件

|cohort|未达门槛|仅NRMSE|仅IoU|两者|
|---|---:|---:|---:|---:|
|ID24|15|0|12|3|
|OOD12|6|1|1|4|
|cross6|4|0|3|1|

全部25个NNdirect未达门槛对象的case_id、contrast、源数、孔径、噪声、生成网格、NRMSE/IoU、corrector stop和PG norm在JSON逐行列出；另含全42×核心方法结果与公开条件。不从图像推断未公开的具体家族。

## 缺口与解释边界

- Final endpoints have unequal actual costs; NNdirect/uniform are the same algorithm with6000-step ceiling, auxiliary flows have3000 auxiliary plus3000 corrector ceilings.
- Matched-budget endpoints use min(finalA,finalB) independently per object and cost axis, sparse paid Phi incumbents, unused allowance and checkpoint gaps; they are at-most allowance, not exact equal expenditure or a GT iterate oracle.
- Hardware partitions are different fixed cases: ID Mac10/XINAN14, OOD Mac6/XINAN6, cross-grid Mac6. Absolute CPU is never pooled across hardware or interpreted as monetary cost.
- Global A+ anchor vs nn_joint is the available truncated physical-current control. Dual/sourcewise have no analogous nn_joint control; anchor changes their initialization and affine family together. No42-object uniform-anchor arm exists.
- Fixed8 has ID4/OOD2/cross2 and120 additional rows; its extra initializers/Born-only network cannot be assigned a42-object denominator. Training costs of Born/domain repair are not the original physics-network numerator.
- All876 solver predictions were sealed; a solver can terminate successfully while failing material metrics. NRMSE/IoU and noise discrepancy are separate, and smooth support topology is diagnostic.
- Cross-grid36→24 remains within the same discrete scalar-TM physical family; it is not independent continuum-PDE validation.
- Strict signed pre/post flags include possible floating-point micro-differences; not every counted case is a substantive material/basin failure. Direct ceiling6000 vs auxiliary corrector3000 also prevents interpreting these counts as one common fixed correction map.
- Stage−1 means the initial candidate entered correction, not that the final image was unchanged. Lower prePhi followed by worse endpoints is descriptive; it does not identify a basin or establish causality.
- Portfolios charge both complete solvers and share measured NN prefixes once. They were selected by public Phi, not material error; no partial-branch budget schedule was frozen, so no invented equal-budget portfolio trajectory is claimed.
- Existing bootstrap intervals are descriptive object-paired percentile95% intervals; no significance claim, threshold search or new resampling was performed here.

证据行：每个JSON摘录原CSV `_evidence.artifact/line`；derived均值差/比值列出两个或三个原表行。原score逐case证据以rows索引和唯一task_id定位。所有CI均来自已有分析输出，未重新采样。
