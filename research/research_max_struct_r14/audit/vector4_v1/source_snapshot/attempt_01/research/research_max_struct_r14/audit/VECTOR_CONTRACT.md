# Vector4 v1：首次运行前冻结合同

2026-11-15，主线程审阅已有 core equations 后授权实现并运行。24条：4个旧公开3D条件 × 4个metric组合和2个成熟sourcewise对照。来源为原R13公开NPZ；只复制公开文件，不打开/复制材料GT，不读取heldout，不在本轮评分。全部任务终止、预测封存后才可由主线程另行决定离线GT评分。

这是继承CM+RR finite vector-DDA的**几何对应实验**，不是R7 Arc复现、新Mie/PDE资格或一般continuum证明。三个独立旧几何，crossgrid8 frame重复第一个几何；全部已经曝光。

材料为每voxel一个实χ=εr−1，216实未知量，共享xyz/全部60源。J/Q为648×60，cell-major xyz；box[0,8]；原G=V Goff+sI，s=−1/3+i k³V/(6π)，原2×2×2空间块有24个currentcomponents。无GTsupport、无2Dstack、无polarizability替代、无新backend/GPU。

|ID|current|material|
|---|---|---|
|00|CIE_TSOM beta0|CIE_TSOM beta0|
|01|CIE_TSOM beta0|A_TSOM beta1|
|10|A_TSOM beta1|CIE_TSOM beta0|
|11|A_TSOM beta1|A_TSOM beta1|

四组合全部global实PR+、原exactquadraticlinesearch。前60cycles mf1/beta0=6/beta1=6；后60cycles mf2/beta0=3/beta1=1。原mf切换重置gam/旧gradient/direction；不移植2Dzero-paddingtransport。原projected3DIFFT加SVD-rangebasis、sourceSVD/DPCS、rho10、每源data/stateweights保留。各materialblock均用原`material`完整analyticgradient，L-BFGS-B相同maxiter12、ftol1e-12、gtol1e-8、maxls25、[0,8]bounds。不得调用硬编码[0,4]的`VectorPR3D.aplus`或2DArc。

两强对照直接调用原R13`run_vector_auxiliary`，方法`CIE_sourcewise`/`Dual_sourcewise`、outer120、inner12，保留其原不同beta schedule及sourcewisePR。对照仅加观测性的instance-local current/prepare hooks，捕获已完成state和计费；不改变其current/material算术或stage规则。00不冒称原Dual，因为其material后半beta也为3。

每条120整进程CPU秒，proposal绝对deadline72；随后同一原物理corrector至`min(120, corrector开始整进程CPU+48)`，保留严格48秒correctorallowance。setup/imports/SVD/分解/失败/序列化都计费。proposal原120cycles先结束时不添加cycles，也不把未用proposalallowance转给corrector。每2整进程CPU秒首次完成readout时审计，加initial/60/120/final原生审计。超限操作付费，其完成候选不得进入预算内选择。

原Phi=`mean_s ||Yhat_s−Y_s||²/||Y_s||² + 2e−4 mean_edges(χi−χj)² + 1e−6 mean_cells χ²`，无1/2，不使用2D全Frobenius归一化。全部合法、预算内、paidphysical候选以最小Phi选primary；secondary为原公开noiseball内最小R，无则primary显式回退。crossgridnoise球不覆盖模型误差。raw/corrected/primary/secondary、candidatepool、所有negative/failure/outbudget均保存。

固定fullstate在1/16/60/61/120及last保存inspectableNPZ，包含χ、J/Q、gam、旧precurrentgradient/direction、mf/stage/cycle、basis/DPCS、sourceweights、原prior/box/flags/counters；没有把m-only文件当resume。若预算/异常打断未完成readout，最后completedstate及in-flight状态分开标明。corrector使用同一个physicalmodel原value/gradient和box。

逐任务费用含physicalLU/sourceRHS/adjointRHS、basis/SVD/IFFT、densecurrentKB/SB、materialG、onehopforward/adjointRHS、setup/serialization/失败和CPU/wall。子进程记录OS可用的实际peakRSS；controller额外记录所有child lifetime的OS CPU合计，含最终serialization。旧历史counters不冒充新费用。

实际配置固定于`configs/vector4_v1.json`；冻结sourcecopies/provenance在`audit/vector4_v1/source_snapshot/`，inputs在`data_public/vector4/`。独立24任务允许2个本地AppleCPUchildren、单BLASthread。全部终止后controller写`predictions/vector4_v1/PREDICTION_SEAL.json`；seal只封存实际文件/任务状态，不代替科学review，不自动授权GT评分。
