# R12：NN 接口保留与材料选择边界

本报告合并两组已完成开发证据：R8四公开对象的冻结 NN 接口比较，以及两个高对比度 true3D 几何加一个 cross-grid 配对。它们没有验证新增的独立24 secondary estimator。两组均先封存完整预测、再离线评分；本整理轮没有跑新实验。

## 1. 四个公开对象：DR保留了NN信息，但收益不普适

准确 ID 是 `r8_test_0001,0004,0022,0027`；0001为6源，其余12源。0027是已有暴露的 heldout-family 反例。它们不是 R11 同编号对象，也不是独立确认集。

直接重用原 R8 physics NN（epoch86、原10通道 features）；四个重推断结果与历史 public NN 初始化完全一致。没有新训练或按初始Phi拒收 NN。NN 本身携带原监督训练先验，因此 NN 与无NN方法的信息来源不同。

四臂为 uniform-direct、NN-direct、完整物理 NN 双场→DR、NN→原 CIE。统一原 scalar-TM backend、box[.02,20]、原 smooth/TV prior和最低原物理Phi选择。DR为scale6/omega1、3000步/audit50、grouped v2；CIE保留原 native Fourier/DPCS 与3000辅助步、`initialize_physical_current=False`。四臂共同最多3000迭代/120CPU秒原 corrector；proposal费用额外实测，不是等总支出比较。

|R8对象|初始NN amplitudeNRMSE|uniform-direct|NN-direct|NN-DR+corrector|NN-CIE+corrector|
|---|---:|---:|---:|---:|---:|
|0001|0.273914|0.668228|0.076345|0.077289|0.076404|
|0004|0.130735|0.368955|0.080315|0.080656|0.080315|
|0022|0.243888|0.289828|0.289830|0.066429|0.065756|
|0027|0.264002|0.460948|0.120932|0.120458|1.012975|
|四对象均值|0.228135|0.446990|0.141855|0.086208|0.308862|

DR相对NN-direct的主要材料改善来自0022。0001/0004分别小幅变差0.000944/0.000341；0027与NN-direct接近，而原CIE接口进入了材料明显更差的低Phi候选区域。不能把四例均值的差异解读为每例胜出。

原物理目标也保留：0022 NN-direct约0.1004618，而NN-DR/NN-CIE约5.35e-6；0027 NN-direct/NN-DR约4.2866e-4，NN-CIE约4.6267e-4。0027 CIE低数据目标仍对应材料误差1.012975，表明数据目标与图像恢复必须分列。

在共同corrector前，按原Phi选出的候选为：

|对象|初始NN NRMSE|DR raw NRMSE|CIE raw NRMSE|
|---|---:|---:|---:|
|0001|0.273914|0.107522|0.205175|
|0004|0.130735|0.074037|0.130735|
|0022|0.243888|0.055806|0.073886|
|0027|0.264002|0.166300|1.049422|

这支持有限开发观察：完整物理双场DR接口在四个指定对象上没有像该CIE接口在0027那样丢失材料优势。它不是独立确认，也不是反射的单因素归因——本组没有NN→AP，DR与CIE的current表示、辅助目标和proposal费用不同。没有人为错误prior的新臂，不能推出任意错误NN先验的鲁棒性。

成本：四对象累计 NN-direct **43.438CPU s**、NN-DR **78.432CPU s**、NN-CIE **104.459CPU s**；uniform-direct **39.273CPU s**。NN-DR/NN-direct逐例CPU比为1.239、1.644、4.549、1.580。实际共享预处理共0.763496CPU s（含其计时范围内的Torch导入、加载、features和推断），无新训练；初始NumPy及solve前module imports不在这些局部clock内。DR另有4个H Cholesky、192 setup RHS、126000 iteration RHS，不能仅报物理LU。

证据：[neural完整报告](delegated/neural/RESULTS_ZH.md)、[冻结协议](delegated/neural/FROZEN_PROTOCOL.json)、[预处理receipt](delegated/neural/PREPROCESS_RECEIPT.json)、[all16 seal](delegated/neural/ALL_PREDICTIONS_SEALED.json)、[scores](delegated/neural/OFFLINE_SCORES.json)、[cost summary](delegated/neural/COST_SUMMARY.json)、[图像](delegated/neural/materials.png)、[已付费Phi/CPU图](delegated/neural/objective_cpu.png)。

## 2. 高对比度true3D：raw DR好，原Phi corrector使其退化

seed2026111213在随机几何和结果前冻结。两个独立非球形三维几何分别含chi8和混合chi2/8；再将第一个几何用8³ forward生成、6³ inverse重建。两网格标签都是同一连续分区的解析体积平均。配对不是第三个独立样本。

每个逆解有216个实材料、648个复矢量分量、60个原R8电偶极照明和180电接收通道，初始chi=.25、box[0,8]。原 R8 CIE/Dual current/material/physical与quadratic prior保持；仅由本补充wrapper传入上界8。基线120cycles/inner12，AP/DR3000/audit50/grouped v2；共同corrector最多1000迭代/60CPU秒。保留 inherited projected-IFFT range SVD坐标；不称为逐字2019复现。

本组使用原vector物理目标

\[
\Phi=D+R,\quad D=\frac1P\sum_p\frac{\|F_p(m)-Y_p\|^2}{\|Y_p\|^2},
\]

\[
R=2\times10^{-4}\operatorname{mean}_{edges}(m_i-m_j)^2
+10^{-6}\operatorname{mean}_{cells}m_i^2.
\]

这里D没有1/2，且按源归一；不能与 scalar-TM 的 .5全球归一Phi直接比较数值。材料指标也采用本组relative L2，不与前节amplitudeNRMSE混算。

|对象|方法|raw材料relative L2|corrected材料relative L2|corrected Phi|完整CPU s|
|---|---|---:|---:|---:|---:|
|high8 matched|CIE|0.457843|0.129785|0.00072734924|4.268|
|high8 matched|Dual|0.454820|0.129783|0.00072734924|4.398|
|high8 matched|DR|0.061631|0.129783|0.00072734924|38.064|
|high8 matched|AP|0.103885|0.129786|0.00072734924|44.011|
|mixed2/8 matched|CIE|0.618636|0.227237|0.00037537091|7.606|
|mixed2/8 matched|Dual|0.608627|0.227232|0.00037537091|7.544|
|mixed2/8 matched|DR|0.078315|0.227239|0.00037537091|45.797|
|mixed2/8 matched|AP|0.216416|0.227260|0.00037537091|42.864|
|high8 8³→6³|CIE|0.397404|0.432218|0.00341592185|7.780|
|high8 8³→6³|Dual|0.395368|0.432220|0.00341592185|8.255|
|high8 8³→6³|DR|0.491405|0.434481|0.00341507045|37.880|
|high8 8³→6³|AP|0.424022|0.432217|0.00341592185|28.394|

共同corrector后，两个matched对象四臂材料误差与Phi近乎一致，DR没有最终材料收益且费用更高。保留raw结果能区分“没有找到好材料”和“找到后被同Phi选择/corrector偏好覆盖”。

## 3. 保存状态分解：何时是data/prior权衡

本分解没有重跑forward/adjoint或读取新truth。直接用保存的残差平方和廉价原R求Phi分解；AP/DR raw物理current已存，可用原S读取全球残差。baseline raw精确物理current未存，不伪造其精确残差；但已知的加权范数下界足以认证其球外。

|matched DR|raw D|corrected D|raw R|corrected R|材料relative L2|
|---|---:|---:|---:|---:|---|
|high8|1.0210215e-4|1.3576744e-4|7.1731734e-4|5.9158180e-4|0.061631→0.129783|
|mixed2/8|1.0034028e-4|1.3032556e-4|3.5855561e-4|2.4504535e-4|0.078315→0.227239|

这两例中，corrector降低Phi靠的是R的下降超过D的上升；不是数据和prior同时改善。物理CN99也从满足变为不满足：

|DR|raw物理残差/δCN99|corrected物理残差/δCN99|状态|
|---|---:|---:|---|
|high8 matched|0.990166|1.143895|球内→球外|
|mixed2/8 matched|0.983948|1.125702|球内→球外|
|high8 8³→6³|5.868818|4.569687|仍球外|

还有另一种负结果：mixed AP 的D从1.9413650e-4降至1.3032811e-4、R从2.6823625e-4降至2.4504280e-4，**两项均改善，材料误差仍从0.216416升至0.227260**。因此不能将所有材料恶化解释为单一prior权重问题，更不能由D/R同时下降推出材料恢复。

真实raw选择位置与费用也保留，不能按事后发现好候选的时间免费截断实际跑完的proposal：

|DR对象|原Phi选中的projection迭代|候选已付费CPU s|完整固定proposal实花CPU s|
|---|---:|---:|---:|
|high8 matched|2000|14.272|34.734|
|mixed2/8 matched|2400|24.908|40.745|
|high8 8³→6³|50|0.574|35.972|

raw-only 10CPU上限内的Phi所选材料误差分别为：high8 DR0.062917、CIE0.457843、Dual0.454820；mixed DR0.094606、CIE0.618636、Dual0.608627。CIE/Dual完整proposal实际仅花约2.60–3.44CPU s，并非同样花满10CPU。此前包含共同corrector的前缀另存；不能用raw-only比较冒充最终等支出胜利。

分解产物：[逐臂JSON](delegated/vector_high/RAW_CORRECTOR_DECOMPOSITION.json)、[完整中文表](delegated/vector_high/RAW_CORRECTOR_DECOMPOSITION_ZH.md)、[raw-only CPU前缀](delegated/vector_high/RAW_CPU_PREFIX_COMPARISON.json)、[零新Maxwell成本receipt](delegated/vector_high/DECOMPOSITION_RECEIPT.json)。实际新增仅外部接收算子装配一次及360列存量current读出，约0.131CPU s。

## 4. cross-grid与有限DDA资格

相同high8几何的8³/6³ clean response相对差异为 **0.261567**，远大于原1%随机噪声；以sigma定义的CN99球没有包括该离散模型误差。它表明该配对存在强模型分辨率差异，不证明全部合法coarse材料都不可能进球。

chi8对应epsilon_r9；内部最大cell相位约1.885（6³）/1.414（8³），是很粗的high-index离散挑战，不是连续Maxwell精度证明。mixed连续材料含2和8，但6³体积平均最大值7.80691，原样报告，没有据结果改几何。

本组累计276.860 solver CPU s、3913 physical LU、234780 source RHS、169260 adjoint RHS；投影另有6个H Cholesky、1080 setup RHS、1080000 iteration RHS。记录包括实际archive费用；nested timers/counters不重复相加。12个预测与前缀在truth评分前封存。

证据：[high3D报告](delegated/vector_high/RESULTS_ZH.md)、[冻结协议](delegated/vector_high/FROZEN_PROTOCOL.json)、[all12 seal](delegated/vector_high/ALL12_PREDICTIONS_SEALED.json)、[scores](delegated/vector_high/OFFLINE_SCORES.json)、[cost ledger](delegated/vector_high/COST_LEDGER.json)、[材料切片](delegated/vector_high/MATERIAL_SLICES.png)、[Phi/CPU](delegated/vector_high/PHI_CPU.png)。

## 5. 从开发负结果到预注册secondary规则

Primary保持原Phi，保证历史可比性。新增secondary只在已有、合法、完整物理、cap内、CN99球内候选中最小化**同一个R**，按Phi破同分；无合格点明确unavailable。它不生成新候选、不加forward、不向轨迹反馈，不是GT oracle。

注册文件 [SECONDARY_SELECTION_FREEZE.json](delegated/confirmation/SECONDARY_SELECTION_FREEZE.json) 如实记录：24观察已生成，确认inverse任务开始0，确认GT评分0时注册。其科学动机来自上述**已暴露开发结果**；不能写成在所有数据生成前就冻结，也不能把这些开发结果当作规则的独立验证。

该secondary同时用于独立24八个方法，并保持240CPU整child、144proposal/96corrector与原primary。附加R/残差检查、快照及写出仍计费。不同方法的审计池覆盖率、无候选频率、失败与atomic overrun资格必须一起核查。它是新增candidate-selection estimator，不能将未来收益全归于DR反射。

确认192运行的整组seal与离线评分现已完成，详见确认完整报告。secondary覆盖率与主规则分别报告；它未消除全部失败，也不证明DR普适恢复或原Phi数学错误。正面与负面机制的综合定位见COMPREHENSIVE_EVALUATION_ZH.md。
