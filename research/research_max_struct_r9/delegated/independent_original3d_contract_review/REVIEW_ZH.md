# 原尺寸CM九臂：只读实现合同审计

本次限定静态审计未发现使普通基线变弱、或把A材料机制偷偷混入current更新的有效性错误。以下六点区分代码事实、已披露变体和解释边界；核心数学与最终科学判断由root保留。没有修改运行代码/预算，没有运行solver、打开GT/预测数组或读取在途评分；没有审阅small3D。

## 1. current/material的五臂区别与协议相符

实际入口`run_original_cm_nine.py:16`调用CUDA的`cuda_auxiliary_solve`。其`:305`–`:334`和CPU参考`vector_original_extension.py:315`–`:360`一致：

|臂|current残差缩放p|PR/线搜索聚合|material|
|---|---|---|---|
|sourceCIE|βc/(1+βcχ)，βc6→1|每源独立|普通CIE精确共享实χ读出，βm6→1|
|sourceDual|同普通CIE，βc6→3|每源独立|普通CIE读出，βm6→1|
|GlobalDual|同普通CIE，βc6→3|全部源/系数联合|普通CIE读出，βm6→1|
|AplusGlobal|与GlobalDual完全相同|与GlobalDual完全相同|局部A+材料读出，βm6→1|
|LSIEsource|p=1，β不参与|每源独立|ReΣconj(Q)J/Σ|Q|²并投影[0,4]|

所以GlobalDual与AplusGlobal没有A-current的混淆。sourceCIE与GlobalDual同时改变βcurrent与聚合；只有cross_grid内的sourceDual能进一步隔开这两项，不能将sourceCIE/Global之差全归因于global PR。matched_grid未包含sourceDual是冻结四臂合同的限制，不是实现偷偷漏跑。

材料仍全27000个实χ、共同box[0,4]、air0、R=0、ρ=1，没有inverse support mask。CPU的active_indices仅允许generation_only，实际CUDA逆构造为`OriginalVectorModel(30)`默认全域（`run_original_cm_nine.py:65`）。

## 2. A+是已披露的材料读出变体，不能过度归因

`vector_original_cuda.py:252`–`:294`使用不相交2³-cell/24-vector-component局部H、两轮坐标arc候选。它先提出普通CIE材料seed，仅在固定J的A-block损失降低时接纳，结束后重新算完整固定current的A-block损失再接受；不是用普通CIE损失放行A更新。CPU参考`LocalArcReadout`同样有这一结构。

该CIE-seed提案明确出现在代码`variant`及九臂协议中，不应称与R8原A+完全同一实现。GlobalDual/AplusGlobal比较隔离了“材料readout包”的变化，但同时包含A metric、有限sweep/candidate算法和seed接受，不能仅由此唯一分辨A metric或解析arc的贡献。未发现这一变体改变current，或给普通CIE附加较小MF/不同数据/额外停止限制。

## 3. MF6→10保留系数与联合状态，不是reset/白化

`CudaImplicitFourier.transport`（`:236`–`:241`）以原网格FFT-flat频率索引的嵌套位置，搬运每个分量、每个源的旧γ，其余新频率填0；检查嵌套性。apply/adjoint均用`norm='ortho'`及同一个Vr投影。每阶段只重置PR历史，不把γ重置，也不另做APCS SVD whitening。

J/Q/D通过线性增量保留，nested transport保持旧basis作用，因此没有在MF扩展时丢掉旧current。CUDA的`old_gamma=gamma`不是这里的alias bug：`:327`用out-of-place `gamma=gamma+...`，CPU版本则显式copy后in-place更新。这里只作代码合同检查，未运行数值恒等式测试。

## 4. Eq16使用根号RMS；它仍只是系数停止

CUDA`:335`–`:349`先按每源算`||Δγ_s||²/||γold_s||²`，再对60源平均、开平方，与配置阈值.001比较；CPU`:360`–`:374`相同。不是把平方量直接与.001比较。初始非零更新/零旧范数使停止量未定义并不触发停止；零更新/零旧范数按0贡献处理，这是实现的零范数约定，应保留说明。

每阶段至多500，未附加物理corrector；主要端点明确是stage2，不按GT或较低物理Phi挑stage。系数小变化不能自动称材料收敛、全局最优或完整论文物理复现。sourceCIE/LSIE的辅助目标不增只在已付轨迹作诊断；没有对Dual/A+不相同的交替目标施加额外拒绝门。

## 5. 全60源、共享算子与物理cache的范围自洽

一个resident进程共享G/E/S及一次SGram/EVD/Vr（入口`:65`–`:79`）；两个Y分别做Jd/Qd/Dd和数据权重。每次新arm会重新初始化m=air、γ=0、J/Q/D由该Y的receiver克隆，不继承上一arm或另一cohort的优化状态。

current均覆盖60列，packed16最后一块保留12列；physical validation也循环全部60源。这里是共享matrix-free算子/批量FFT，而不是“每个物理状态一次LU处理60RHS”的密集求解。GMRES逐源但统一采用真实未预条件残差核验。

入口`:81`–`:105`的cache只以完整材料bytes作key，预测与Y无关；命中后用当前Y重算loss，跨cohort复用相同m的确定性预测是合法的。air解析零预测实际RHS=0，首次物理求解、cache命中和来源都记录。此共享queue费用/摊销不等于每方法冷启动独立速度，协议和COST字段已作这一限制。

## 6. 噪声、两cohort与封存隔离无隐藏A优势；描述需保持精确

`public_data`（入口`:25`–`:51`）保持cross_grid Y逐元素原样；matched使用已付coarse clean加原noise乘`||coarse_clean||/||fine_clean||`，sigma同倍缩放。没有样本范数归一化，两个cohort为共享方向/seed的相关噪声控制，不是两个独立随机样本。每个cohort各自data_id、Y、Jd和wd，fixed receiver cutoff仍.1，不是本任务临时换成noise_ball99。

**隔离口径：**准备函数确实从`sealed_truth`命名容器读取已付`Y_clean`、`noise`和coarse `prediction`；源码没有索引材料χ、support或mask，也没有用它们更新材料。这是协议公开的观测生成步骤，不能笼统称整个run进程从未打开sealed目录。实际迭代接收的只是准备好的Y/其receiver状态；没有看到把clean/noise单独作为A+迭代oracle的通道。

入口`:145`最终端点取`endpoints[-1]`；`:165`仅九个成功端点齐全才写统一SEALED_ALL_NINE。数值失败会留下terminal failure，不伪装九个完成。源码中没有GT评分调用。没有读取实际在途输出，因此本审计不替代后续seal验证、远端加载代码核对、GMRES数值资格或root的最终成像判断。

## 复核范围

阅读三份指定源码、`configs/original_cm_nine.json`、基合同`original_vector3d.json`和`CM_NINE_PROTOCOL_AMENDMENT_ZH.md`。未读取任何large array、材料GT、在途预测、score或small3D结果；没有改代码/预算/raw。`READSET.json`保留本次源码快照大小/mtime，未做SHA检查。
