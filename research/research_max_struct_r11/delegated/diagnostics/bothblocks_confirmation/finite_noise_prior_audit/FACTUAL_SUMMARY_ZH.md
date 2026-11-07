# 有限噪声球与既定prior：全部16×5离线审计

仅复用16个既有Y_clean与80个封存Y_prediction；使用原Regularizer((24,24))的同一R。原80候选、门限、成功判定、在线Phi规则及算法均未改。整个审计0 forward、0 adjoint、0 optimization；核心全80审计进程外部wait4 CPU1.261185秒，后续编辑排版CPU未外部采样。

16/16真实GT预测都在各自实际观测y的CN99球内：||Yclean−Y||/sigma范围23.8334–24.8372，共同门25.1654。80候选中65个也在对应球内，因此有65个实际观测球同时包含GT预测和候选预测；其中38个候选材料门FAIL。全部对应case/method在全表明确列出，没有挑选有利对象。其余15个候选数据不合格，尽管各自GT在球内；它们必须与有限噪声兼容的材料偏差分开。

65/80候选Phi低于GT；这65个同时具有更低的数据项和更低的R。没有出现“候选数据项反而更高、只靠prior下降才低于GT”的实例。70/80候选R低于GT，但其中5个仍被较差数据项压倒。这个有限比较不能把所有偏差唯一归因prior，也不能把低Phi当材料恢复证明。

|方法|GT+候选同CN99球|同球但材料FAIL|候选Phi<GT|候选R<GT|median R候选/R真值|
|---|---:|---:|---:|---:|---:|
|cie_single|12/16|8/16|12/16|13/16|0.760182|
|cie_dual|13/16|8/16|13/16|14/16|0.756880|
|global_dual|13/16|7/16|13/16|14/16|0.755050|
|fullG19|14/16|7/16|14/16|15/16|0.756203|
|fullGreenBothBlocks|13/16|8/16|13/16|14/16|0.755263|

全表同时记录||Ypred−Yclean||/sigma、材料L2/relative-L2/amplitude-NRMSE，以及low/high/background各区RMSE和均值。mixed低相真值为2.5；各方法跨8个mixed对象的low区估计均值中位数约2.426–2.467，但low区RMSE中位数约1.644–1.686。区域均值接近真值不能说明空间结构保持；也不支持“所有低相均被完全抹去”的说法。

这些是已实现候选在实际噪声实现、当前prior和固定门下的有限成对事实；不证明全局唯一性、所有噪声不可辨识、某个算法不能继续改进，也不证明一个唯一的优化路径或prior因果解释。核心解释由根线程完成。

可复现脚本：experiments/diagnose_bothblocks16_noise_prior.py。TRUTH_OBJECTIVES_ALL16给出GT数据项/R/Phi；CANDIDATE_DATA_PRIOR_ALL80给出全80分解、三个sigma距离和材料距离；OBJECT_PAIR_CLASSIFICATION_ALL16便于逐对象定位；AUDIT_RECEIPT与EXTERNAL_PROCESS_COST保留零FW与费用。所有candidate重构Phi与原封存Phi完全一致（最大绝对差0）。
