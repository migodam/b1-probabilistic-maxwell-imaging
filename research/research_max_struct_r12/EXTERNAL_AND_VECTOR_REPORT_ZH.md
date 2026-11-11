# R12 外部实测与矢量挑战：结果及边界

外部与3D结果不支持把DR的路径能力直接升级为普遍材料恢复优势。公开实测中共同物理校正后的留出预测接近，但没有材料pixelGT；Austria原150²→32²存在显著离散模型差异，匹配32²后仍有有限预算失败和校正器损害材料；两个3D组中较低Phi同样不保证较好材料。独立24确认已另行完成封存评分，见delegated/confirmation/CONFIRMATION_RESULTS_ZH.md；不与本开发/外部组混合。

## 1. Fresnel：预测拟合不等于绝对材料标定

实际读取继承的2001 `dielTM_dec8f.exp`、`twodielTM_8f.exp`，2GHz、24²、0.16m方域。每源使用真实接收行；27训练源×39行共1053训练通道，另711通道留出。源多极展开只拟合空目标场；没有用目标散射留出重标定增益。原始测量未复制进交付包。完整规则与来源见[合同](delegated/external/CONTRACT.json)、[源清单](delegated/external/SOURCE_INVENTORY.json)。

|实测对象|direct留出残差|Dual|AP|DR|完整原生Dual补充|
|---|---:|---:|---:|---:|---:|
|single|.09150|.09150|.09150|.09152|.0914988|
|twin|.09857|.09856|.09868|.09855|.0985633|

初轮Dual300总PR步之外，另执行原完整500/500/500/1500调度：每对象恰3000PR/材料读出、15000inner评估、27SVD、5次物理stage验证。原Phi选出的auxiliary仍是初始uniform.1（stage−1），共同1000步corrector再产生最终拟合；所以接近的最终图像不能归功于所有proposal均有同样能力。两次完整Dual分别84.799/84.819CPU秒；细节见[完整预算控制](delegated/external_continued/REVIEW_READOUT.md)。

空场验证中位数为.0251475723/.0240806883；独立圆柱界面参考与等面积VIE在24²的差异.0243024323/.0189421309，36²降为.0129375430/.0120232929。这些是校准/离散误差诊断，不是内部场或材料真值证明。最初示意图坐标轴解释曾产生>100%测量误差，失败保留；后来的训练-only名义方位诊断误差.1889438838/.3060280655，已知参考增益转移的留出.1607348202/.2961524435未用于反演。[精确R8记录与限定](delegated/external_continued/FRESNEL_CALIBRATION_FORWARD_EVIDENCE.json)。

噪声协方差、绝对场单位和接收喇叭响应仍未知；.15||Ytrain||是预先固定的MODEL_DISCREPANCY_DIAGNOSTIC，不是CN99证书。拟合峰值chi约9–13；没有pixelGT，不能把低留出残差写成正确介电常数。[所有raw/校正材料的完整0–20显示](delegated/external/EXTERNAL_MATERIALS_FULL_RANGE.png)、[不拟合增益的留出预测](delegated/external/FRESNEL_HELDOUT_PREDICTIONS.png)。

## 2. Austria：原150→32、保存池选择与匹配32必须分开

实际复用数据的PUBLIC_DATA和DATA_RECEIPT均写明150²生成、32²反演。早期报告误写64²，现以独立[元数据修正](delegated/external/METADATA_AMENDMENT.json)更正，旧冻结合同不覆盖；代码始终读取同一个实际Y。epsilon4对应chi3，不是chi4。原生Dual另补原paper runner每stage300 allowance、mf5/8/10、beta8/2/.5/currentbeta8、literal DPCS15；未按已评分结果换成R9另一noise-adaptive策略。

|原150→32|aux材料NRMSE|校正材料NRMSE|校正数据残差|
|---|---:|---:|---:|
|direct|.986564|.986564|.545680|
|Dual，初轮300总PR|.728915|1.419644|.107260|
|AP5000|.955725|.989414|.545675|
|DR5000|2.034990|2.131982|.122674|
|full native Dual，673实际PR|.613774|1.032829|.093899|

完整Dual改善了初轮欠预算的raw材料，但同一物理corrector仍以材料恶化换取目标下降：Phi .1202229→.00441375，R从1.030512e-6升到5.245214e-6，数据项从.1202218502降到.0044085014，最终有3个chi20 cell。这里是**支付更高先验代价换数据残差下降**，并非只怪DR找路失败。此prior只占最终Phi约.1188%；这不证明一个通用替代prior，只证明当前Phi选优不能当成材料正确性。[分解](delegated/external_continued/AUSTRIA_OBJECTIVE_MATERIAL_DIAGNOSTIC.json)。

对实际保存池做一次明确事后选择诊断：固定原.12球，eligible中取最小继承R。Dual v1/full分别从注册最终NRMSE1.41964/1.03283改善为1.17860/.79797；这恰是各自eligible池oracle，但仍比球外的raw形状.72891/.61377差。direct/AP/DR池无eligible点，DR最好的已注册物理残差.12267未被四舍五入成可行。未保存的native-Dual其他stage与未审计inner迭代没有补造或重放。该结果是POSTHOC_EXPOSED_SELECTION_DIAGNOSTIC，不是新算法或确认。[保存池报告](delegated/external_selection/REVIEW_READOUT.md)。

随后仅做一次离线粗GT正演，1LU/20RHS：真实32²GT相对实际观测残差**.2252565**，已在.12球外。没有可用的fine cleanY，因此没有新增fine正演或反推噪声。原DR的A–C相对gap从50步.12782降到5000步.01335，但C传播缺陷仍约.179 incident norm；A数据残差.12不代表两个shadow相交或材料满足完整物理。见[GT与shadow诊断](delegated/external_physics/REPORT.md)。

为隔离上述模型差异，另获准生成匹配32²、seed2026111232的新10%固定范数噪声，复用已经支付的clean32GT预测，不重复正演。真实目标残差.0995131在同一.12球内；其余contrast/box/prior/scale/budgets不变：

|匹配32诊断|aux材料NRMSE|校正材料NRMSE|数据残差|
|---|---:|---:|---:|
|direct|.970767|.970767|.566303|
|full native Dual|.606380|1.081429|.093104|
|DR5000|2.190625|2.253906|.134903|

匹配控制消除了原网格差异，仍保留DR未进球与Dual校正损害材料；不能把原结果全部归因于离散误差，也不能把这一次事后诊断称为新确认。原150→32全部坏结果不变。[匹配32报告](delegated/external_matched32/REVIEW_READOUT.md)、[完整幅度图片](delegated/external_matched32/MATCHED32_MATERIALS_FULL_RANGE.png)。

## 3. moderate true3D：短corrector优势在较长共同预算后消失

复用R8 dyadic CM+RR DDA，216实chi、648矢量current、60偶极源、180接收分量；两个非球形几何加一个8³→6³配对，chi box[0,4]，不是高chi8组。120原PR cycles对3000AP/DR步并非相同费用。初轮各法共同120corrector均碰迭代上限，之后在读GT前统一增加1000iter/60CPU物理continuation，18+18预测一起封存。

|观测|初轮CIE/Dual/AP/DR NRMSE|较长共同corrector后六方法范围|
|---|---|---|
|matched包覆立方体|.10565/.10410/.09252/.09456|.09763–.09766|
|matched连接夹杂|.13614/.13234/.08948/.08146|.08231–.08243|
|包覆体8³→6³|.22367/.22309/.21922/.22462|.22554–.22558|

因此初轮改善不能写成稳健最终盆地收益；包覆体AP继续降低Phi却从.09252恶化到.09763。较长continuation都以relative-F停止，不提供局部极小/唯一性证书。[原始表与费用](delegated/vector/REPORT.md)、[长corrector图](delegated/vector/MATERIAL_SLICES_LONG.png)。

## 4. vector_high：matched raw优势、prior交换与cross-grid负例

此组由projection worker独立记录，本文只汇总，不修改其成果。seed2026111213冻结在随机几何前；两独立非球形材料分区chi8/mixed2–8，另首对象8³→6³配对。相界用解析体积交得到cell平均，逆问题仍是216/648、60/180，chi[0,8]，最大chi8对应epsilon9。原R8 projected-IFFT默认SVD range白化、统一prior和原source-normalized Phi均披露，不能称逐字2019实现。[原报告](delegated/vector_high/RESULTS_ZH.md)。

|观测|CIE raw→final|Dual raw→final|AP raw→final|DR raw→final|
|---|---|---|---|---|
|matched chi8|.457843→.129785|.454820→.129783|.103885→.129786|.061631→.129783|
|matched mixed2/8|.618636→.227237|.608627→.227232|.216416→.227260|.078315→.227239|
|8³→6³配对|.397404→.432218|.395368→.432220|.424022→.432217|.491405→.434481|

两matched对象DR raw进全局CN99球（ratio .99017/.98395），但共同corrector后ratio约1.14389/1.12570，反而出球。其chi8例D从.000102102升到.000135767、R从.000717317降到.000591582；mixed例D从.000100340升到.000130326、R从.000358556降到.000245045。这里不是“更好拟合换材料恶化”，而是**降低prior代价、容忍更差data和材料，仍降低总Phi**。与CIE/Dual提高R换D大幅下降的路径必须分开写。[保存状态分解](delegated/vector_high/RAW_CORRECTOR_DECOMPOSITION_ZH.md)。

matched两对象最终均值CIE/Dual/DR/AP约.178511/.178508/.178511/.178523，DR没有最终材料收益且总CPU更高。cross-grid中DR Phi稍低（.003415070对其他约.003415922），材料却稍差（.434481对约.43222）。其8³/6³ clean数据差异.261567，noise-only CN99球没有覆盖这项模型差异。baseline raw精确全局残差因未存对应current如实缺失，只用已保存source-normalized D的范数界判定球外，不补造数值。

此高chi离散挑战的最大内部cell相位约1.885（6³）/1.414（8³）；它不构成连续Maxwell精度资格。配对不是第三个独立随机对象。所有raw/full-proposal费用、corrector阶段和已冻结10/30/60/120CPU候选分别保留；raw选中时刻不等于跑完整proposal花费。[费用](delegated/vector_high/COST_LEDGER.json)、[材料图](delegated/vector_high/MATERIAL_SLICES.png)。

## 5. 当前结论与未完成项

指定旧坏终点的反射逃逸证据仍见[BASIN_ESCAPE_REPORT](BASIN_ESCAPE_REPORT_ZH.md)，不被外部负例抹除；它也不能替代外部物理标定、跨网格资格或独立确认。本轮清楚观察到搜索路径、当前prior/selection、噪声球与模型差异是不同环节。独立24的整组执行、封存和评分已完成，结果在独立确认报告中单列。全局费用仅按已核实目录范围列在[COST_LEDGER](COST_LEDGER.md)，不以低残差或发布完成替代科学验收。
