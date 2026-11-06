# 原3λ四立方体：实际执行状态

**后续更新：本页以下保留物理诊断前快照。随后预声明IGT source59的30/60差由CM.508374降至.359841，仍未达到10%噪声尺度；两RHS及费用保留、不扩大IGT生成/反演。父线程已批准并实际启动[CM九臂协议增补](CM_NINE_PROTOCOL_AMENDMENT_ZH.md)：cross-grid5臂+matched finite-grid4臂，统一packed16，九seal前不作材料GT评分。当前机读状态以`RESULTS/original_vector3d/STATUS.json`为准。**

**60³原尺寸数据已完成并冻结；五臂完整比较尚未完成，也未取消。父线程明确要求先处理53.5%的跨网格物理误差，GPU已完全释放，packed16实现待命。没有进行GT材料评分。**

数据保留原3λ范围、60源/180向量接收行、four-cubes χ=1/2/3/4、60³生成与30³全27000实材料未知量。生成端的严格零电流消元与44³包围盒卷积仅改变计算嵌入，h=.05、原绝对坐标、S/inc、自项和全部非零电流不变；一次独立完整/裁剪作用差为3.56e−16。反演端不使用已知支撑。物理始终为独立CM+RR点DDA，不能称为作者积分单元实现。

最终60源由52个已付Mac检查点、7个CUDA续算和1个已验证CUDA source59基准状态组成。全部复用状态在远端原矩阵上重新验证，原相对残差≤2e−6；平台RHS相对差约4e−17。新加入的圆对称复高斯噪声期望相对范数0.1，固定种子2026110619，不强制重标定样本范数；具体噪声实现与作者未公开细节明确区分。

最重要的现有物理证据为单次离线60³→30³真值响应对照：global L2=0.535254，mean-source RMS=0.535244，高于0.1噪声。该对照只做一次、保留60个合格完整RHS，并非反演输入、网格选择或GT调参。它限制了后续材料与论文复现解释；最终判断由父线程负责。

已完成的工程核验：同一source59完整RHS，i9/SciPy-MGS99.905wall秒、416步，CUDA驻留CGS2为5.890wall秒、368步；GPU状态代入原CPU矩阵的残差1.865178e−6。电流相对差1.364e−5、接收预测差6.591e−6。16.96倍包含正交化/检查节奏差异，不能声称纯GPU硬件提速。原矩阵、complex128精度、GMRES60/20cycle及原残差阈值保持不变。

一次实际完整current+material步骤回放，strided4=2.97547wall秒，连续[b,3,n,n,n]布局packed16=1.02073秒，收益2.915倍；material、gamma、current、field、alpha在该公开固定状态上差异全部为0。只比较此一候选，无扫描、GT读取或新完整求解。整个进程成本10.171875CPU秒，内部8.859秒计时已包含，不重复累加。

旧strided4 sourceCIE已保存MF6阶段113步的材料及完整物理预测；MF10另117步后为工程优化暂挂/终止，共230个完整辅助步、120个物理RHS（含60个air零RHS）。全部原始轨迹与616CPU秒在`RESULTS/original_vector3d/engineering_attempt1_strided4/`。暂停前wall622.14秒，总进程wall871.75秒包含暂停；旧程序没有导出完整Green计数，可能的末步剩余计数明确未知，不能填0。原始阶段未完成，未伪称收敛或计入五臂胜负。

准备好的五臂为sourceCIE/sourceDual/GlobalDual/AplusGlobal/LSIEsource，统一ρ=1、R=0、air0、box[0,4]、MF6→10、每阶段最多500、原文Eq16根号RMS<.001、无额外物理corrector。sourceCIE/LSIE的阶段内辅助目标单调性仅作现有轨迹接线诊断；Dual/A+不共享同一个阶段目标，不适用该要求。A+明确为R9向量变体，包含按完整固定J的A损失接受普通CIE种子提案，未冒称与R8冻结A+完全相同。

当前状态与去重成本见`RESULTS/original_vector3d/STATUS.json`及`ACCOUNTING_BEFORE_PHYSICS_DIAG.json/.csv`。计时数值/工程CPU合计3549.628秒，已知完整primal RHS245，其中60个为air零RHS；各进程wall不能相加当作总历时。旧共享预计算和GPU接线检查已包括在616CPU秒中，没有重复计费。所有物理失败/中断和已知计数缺口保留。

代码入口为`src/vector_original_extension.py`、`src/vector_original_cuda.py`与`experiments/vector_original_cuda.py`，配置`configs/original_vector3d.json`。远端独占任务目录为`D:/AI/MAX_STRUCT_R9_ORIGINAL3D_20261106`；新packed16代码包已上传但五臂没有重启。最终释放记录`GPU_RELEASE.json`实测compute-apps为空、GPU0%、显存0MiB。禁止在父线程物理诊断结论前自动重启。
