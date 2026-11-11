# R12 基准合同与证据范围

此表区分已完成的暴露开发诊断、公开实测预测、离散3D挑战与独立确认。新建反演cohort的材料评分在相应预测封存后进行。例外是原Pro源码重放：它继承公开toy的逐步GT诊断记录，未用GT选候选，明确为暴露开发重放。既有对象不能由一次新的封存变成盲测。独立24对象已192/192运行封存，已完成离线评分，见delegated/confirmation/CONFIRMATION_RESULTS_ZH.md。

|组别|实际对象/观测|物理与未知量|方法/预算|当前状态与可支持范围|
|---|---|---|---|---|
|原始源码重放|4个既有toy对象|原作者包内离散问题|direct/AP/DR及原200步续跑|已完成；只确认源码行为与开发反例，不作为新盲测。[回执](experiments/pro_reproduction/REPLAY_RECEIPT.json)|
|成熟基线接入|4个toy，scaled/clipped两套合同|继承scalar-TM与R8成熟读出|CIE/Dual/Global/A+/实际Green读出|已完成；原始/clipped分别保留，不混合物理合同。[费用与结果](delegated/baseline/COMPARISON_AND_COSTS.json)|
|R11坏终点干预|环CIE0001、壳G190004、桥G190022|24² scalar-TM，576实材料|同保存起点；direct、固定平滑multistart、AP、DR；付费240 body-CPU前缀|已完成的暴露开发路径证据；环/壳有近局部极小数值依据，桥未同等收敛；不是均匀初始化下优于所有成熟法。[报告](BASIN_ESCAPE_REPORT_ZH.md)|
|冻结NN接口|4个既有R8 ID，16预测|原十通道、epoch86权重、原物理与prior|uniform→direct、NN→direct/DR/CIE；共同corrector，proposal另计|已完成；无新训练；不是R11同编号对象或独立24确认。[报告](delegated/neural/RESULTS_ZH.md)|
|Fresnel公开实测|2001 single/twin，2GHz；2对象|24²/0.16m；27训练源×39实际接收行，711留出；空目标多极校准；实chi[0,20]|direct/Dual/AP/DR；180proposal+60corrector CPU上限；Dual另补完整500/500/500/1500合同|12臂外部组中8臂已完成，另2完整Dual完成；只有留出预测，无pixelGT。未知噪声协方差，.15球为模型差异诊断。[报告](EXTERNAL_AND_VECTOR_REPORT_ZH.md)|
|Austria原观测|既有R8 epsilon4/400MHz，**150²正演→32²反演**|20平面波、40接收器，1024实chi；目标chi=3；10%固定范数噪声|direct/Dual/AP/DR；另补原生完整Dual；同prior/box/初值/corrector；.12相对球|已完成，明确含网格模型差异；不是64²生成，也不是原论文全复现。粗GT残差.2252565在.12球外。[元数据修正](delegated/external/METADATA_AMENDMENT.json)|
|Austria匹配32|同32²目标/算子，新seed2026111232，1观测|32²→32²；同chi3、10%固定范数噪声；GT残差.0995131在.12内|direct/full native Dual/DR；同scale6、球、prior、预算|已完成；事后获准的OFFLINE_PHYSICS_DIAGNOSTIC，非新确认；分离原网格差异与仍存在的有限预算失败。[报告](delegated/external_matched32/REPORT.md)|
|moderate true3D|2个非球形几何：混合包覆立方体、连接双夹杂；加首对象8³→6³配对，共3观测|216共享实chi、648矢量电流分量、60偶极源/180接收分量；chi[0,4]；1%噪声|CIE、Dual、Global、A+、AP、DR；120原PR cycles或3000投影；共同120corrector后，统一再加1000iter/60CPU|18+18端点封存评分；较长共同corrector后各法趋同，短预算差异不能作为最终盆地优势。[报告](delegated/vector/REVIEW_READOUT.md)|
|vector_high离散高对比度|2个seeded几何：chi8、mixed2/8；首对象另8³→6³配对；3观测不是3独立对象|同216/648矢量离散模型；chi[0,8]；解析体积平均；chi8对应epsilon9|CIE/Dual120cycles、AP/DR3000步、共同1000iter/60CPU；固定原Phi选候选|12预测与CPU前缀均封存后评分；matched DR raw材料较好，但共同corrector后收益消失；粗网格不具连续Maxwell精度资格。[作者工作记录](delegated/vector_high/RESULTS_ZH.md)|
|独立DEV6与CONFIRM24|新6开发与24独立观测已生成，GT分离|24²/12源/48接收，继承scalar-TM|8法：CIE/Dual/Global/A+/实际G19/DR/directmultistart/冻结NN→direct；240 whole-child CPU|192确认与48开发均已封存；已完成离线评分，科学结论见确认报告和总报告。[执行交接](delegated/confirmation/CURRENT_EXECUTION_HANDOFF.md)|

## 共同比较限制

1. 原scalar-TM目标是`.5||Ypred−Y||²/||Y||²+R`；R8矢量组是按源归一化的`D+R`，无同一个1/2因子。不能把两者Phi数值或LU维数直接混为一个物理基准。
2. proposal步数、单步费用、状态归档与corrector费用不同。固定work、付费CPU前缀、整子进程CPU和wall均分别保留；没有统一“总费用相同”的普遍结论。
3. 合法源相关实测掩码、旧DPCS/Fourier、原prior、真实Maxwell校正器和失败端点均保留。没有用缺失接收器填零、用GT选early stop或把posterior/噪声可行性当作材料正确性。
4. 所有3D组都是小规模有限CM+RR dyadic-DDA挑战；高chi下粗cell相位尤其大。不能据此宣称连续Maxwell精度、完整2019原算法或跨模型通用优势。

费用范围见[COST_LEDGER](COST_LEDGER.md)，外部与矢量详细观察见[中文报告](EXTERNAL_AND_VECTOR_REPORT_ZH.md)。

确认调度继承了method-specific priority：五个baseline在Windows为BelowNormal，三种search继承controller优先级；单次中途快照为Normal，未记录全程逐task优先级。均固定独立P核、单BLAS线程并按process CPU截断；因此可报告固定CPU额度内的恢复，不能声称完全相同优先级下的wall-throughput对照。实测活跃核约1.5GHz，完整频率快照见confirmation记录，未在队列中途改电源设置。
