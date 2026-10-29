# 接受、方向切向、驻点性与目标后果：R3 离线分析

日期：2026-10-30。只使用现有公式、已付费 development Jacobian/日志；**新增0 FE /0 RHS**。以下是有限维局部代数及其适用条件，不是新颖性声明；最终理论与物理解释由父线程 Astra 决定。四状态机制脚本仍为 **PDE NOT_RUN，等待父线程队列许可**。

## 1. 先分清四个不同问题

完整目标函数接受检查“这个有限步是否使实际目标下降”；采样 Jv 检查“这个方向的局部数据变化是否被 chart 表示”；驻点检查“是否还存在一阶可行下降方向”；目标准确性检查“恢复的目标是否接近真值”。前两项一般不能推出后两项。

本注中 J、残差 r 均为同一当前状态、σ=.01 白化、按公开参数尺度归一化的实坐标。精确和 chart 模型共享当前完整响应、prior、damping 与可行集。令 E=J−Jhat；此处符号与离线数据脚本的 DeltaJ=Jhat−J 相反，范数不受影响。固定二次 prior 的完整梯度 g=Jᵀr+p，近似梯度 ghat=Jhatᵀr+p。

## 2. 一个方向检查遗漏了什么

**局部恒等式。** 对无约束/共同自由变量中的 chart 步 δhat，记近似 QP 正规残差 η=Hhat δhat+ghat。则 F28 的完整残差满足

\[
e=H\widehat\delta+g
=\eta+E^T(J\widehat\delta+r)+\widehat J^TE\widehat\delta.
\]

证明：展开 `(H−Hhat)δhat+(g−ghat)+η`，使用 `H−Hhat=EᵀJ+JhatᵀE` 和 `g−ghat=Eᵀr`。

所以即使 `E δhat` 很小，仍未约束 `Eᵀ(J δhat+r)`；即使近似 QP 精确求解，也不能从方向 defect 单独推出完整正常残差 e 小。采样一个非零方向仅控制 E 在其一维张成空间上的作用；当可行材料维数更大时，对其余方向没有统一信息。步长因 trust region 缩小时，绝对 `||E δhat||` 还会随步长一起缩小。

这不是需要新数值玩具验证的事实：`E δhat` 和 `Eᵀr` 本来就是不同向量及不同作用，前者没有包括后者的代数项。有限次采样若未覆盖全空间且没有额外结构假设，也不能自动得到 `||E||` 上界。

**驻点后果。** 无约束时，`g−ghat=Eᵀr`。若 chart 梯度已小，完整梯度只有在遗漏的残差余切 `Eᵀr` 也小的条件下才会小。对固定凸可行集 C，采用投影梯度映射

\[
\mathcal G_\alpha(g)=\alpha^{-1}[z-\Pi_C(z-\alpha g)]
\]

时，投影非扩张性给出

\[
\|\mathcal G_\alpha(g)-\mathcal G_\alpha(\widehat g)\|
\leq\|E^Tr\|.
\]

因此有界问题应审计完整 KKT/投影梯度；不能把活跃坐标上非零的普通梯度直接判作未收敛，也不能因 chart 的 projected gradient 小而跳过完整检查。

## 3. 小模型误差为何仍可能改变目标更新

令 H 分块为 `(H_w,C;Cᵀ,H_x)`，`H_w≻0`，`S_X=H_x−CᵀH_w⁻¹C≻0`。在无约束或正确约化的共同活跃面内，F29 给出

\[
\widehat\delta_x-\delta_x
=S_X^{-1}(e_x-C^TH_w^{-1}e_w).
\]

这表明 relevant quantity 是经过墙消元后的残差，而非总 J 的某个盲截断指标；很小的 `λmin(S_X)` 可以放大这个残差。F30 是上界，分子没有统一小量或分母没有有效下界时，不构成小目标步误差证书。

活跃面可以变化时，F31 比共同活跃面等式更适用：若 H≽μ_H I，则在相同凸集合、精确 QP 最小点下，联合步差受 `(H−Hhat)δhat+g−ghat` 的范数除以 μ_H 控制。若 QP 未精确求解，还需加入相应优化/KKT 残差。若 μ_H 来自算法 damping，界仅属于这次阻尼更新；不能把它作为数据目标可辨识性或 posterior 收缩。

## 4. 目标相关的 Schur 一阶余切：独立复核

父线程提出的差分关系成立。固定 `Λ_W≻0`、固定 `Λ_X`，定义

\[
N=(A^TA+\Lambda_W)^{-1}A^TB,\quad
R=B-AN=Q_AB,\quad S=B^TQ_AB+\Lambda_X.
\]

对光滑 A、B 的微分：

\[
\boxed{dS=(dB-dA\,N)^TR+R^T(dB-dA\,N).}
\]

逐方向证明：`vᵀSv` 是
`min_h ||Ah+Bv||²+hᵀΛ_W h+vᵀΛ_X v`，最优 `h*=−Nv`；包络求导消去最优 h 的导数，得到 `2⟨Rv,(dB−dA N)v⟩`。对称双线性化即上述矩阵微分。此推导需要固定 prior；随状态变化的 prior 还需其导数项。

其可操作含义：应关注 `dB−dA N` 对 nuisance-profiled response R 的作用。它可以区分“很大的墙/目标切向误差但目标 Schur 后果小”和“总切向误差不大但击中目标有效响应”。这只是精确一阶关系；将有限 chart 差值代入微分，只得到 leading term，不能自动得到有限误差证书。

## 5. 完整目标下降仍不证明驻点或目标准确性

F32 保证被接受步的完整目标下降（按实现的接受条件），没有要求新点满足完整 KKT，也没有检查所有剩余可行方向。因此多次接受后仍可停在 trust-region stall，或近似模型导致的非驻点。

即使完整 KKT 成立，它仍是声明 objective 的局部性质。要从拟合/驻点进一步推出目标真值误差小，还需要正确模型、噪声控制、合适解分支及目标在允许墙类下的稳定可辨识性等条件。墙补偿、其他局部极小、模型误差和噪声拟合均未被接受条件排除。此轮离线导数分析不访问真值，因而不报告目标准确性结论。

## 6. 已付 dev 证据：每个差值都在同一状态计算

原始来源为 `results/online/dev_r0/{FIXED_CURRENT_CHART,ADAPTIVE_CURRENT_CHART}/endpoint_exact_jacobian.npz`，仅匹配 **z 与完整 data 完全相同** 的已保存 FINAL_JACOBIAN 中的 approximate J；fallback 文件内部的 full/approximate 对也始终共享 z、data。没有把两个方法不同端点的 Jacobian 相减称作 chart 误差。

派生输出：`ENDPOINT_CONSEQUENCE_RUN/SUMMARY.json` 与每条记录的 NPZ。总36条记录、33个有已保存近似矩阵的同状态对；包含 endpoint/FINAL 文件的重复表述，**不是36个独立物理世界**。以下只列各方法所选 endpoint，同一行中的差值才是比较对象。

| 同状态量 | Fixed endpoint | Adaptive endpoint |
|---|---:|---:|
| 完整 objective |207.231574|205.561236|
| 完整 projected stationarity（源码归一化）|.688186|.00978654|
| chart projected stationarity（使用相同完整对角归一化）|.0249107|.0173350|
| whole tangent 相对 Frobenius 误差|.625903|.0157950|
| wall tangent 相对误差|.381040|.00432043|
| target tangent 相对误差|1.004741|.0299666|
| 数据 Schur 相对 Frobenius 误差|.865160|.0308071|
| 完整数据 Schur 最小特征值|.00505027|.00517992|
| `||Eᵀr||`（总 / target）|18.8226 /9.29553|.452398 /.418387|
| 全19维1e-3阻尼下，目标无约束步相对差|.511716|.0761058|

数据 Schur 使用公开墙 prior、目标精度0、无算法 damping；步差为同状态局部 QP 的后果，**不声称是实际有界 solver 接受的步**。不得把约.005的局部最小特征值当成整个目标盒的全局稳定常数。

一阶目标有效公式采用 `DeltaJ=chart−full`，以 full A、B 为展开点：

| 有限 Schur 差与一阶 leading term | Fixed endpoint | Adaptive endpoint |
|---|---:|---:|
| `||DeltaB−DeltaA N||_F / ||R||_F`|1.06411|.0325773|
| 实际 `||Shat−S||_F`|440.505746|14.416109|
| leading term Frobenius 范数|365.120911|14.793901|
| `||actual−leading||_F / ||actual||_F`|1.30320|.0334799|

因此 fixed 端点的有限变化不能被这个一阶项定量替代；adaptive 端点该一阶项与实际差较接近。这是两个已付端点的事后解释，不是 rank、tolerance 或方法选择依据，也不是新实验成功标准。

另读 adaptive `joint_start_2` 已有 STEPS 日志：34次方向审计中25次所有频率通过、9次fallback；最后一次各频相对 defect 约.0244–.0320，通过既有.1门槛，但原始方向数据变化范数仅约3.30e-13–4.14e-13，最终 stop 为 TRUST_REGION_STALLED，完整 projected stationarity=.00978654。这支持“方向审计通过并不自动达到预声明.005驻点门槛”的有限运行观察。它不证明阈值应该改成别的值，也不提供目标真值准确性。

## 7. 成本、失败与权限边界

本次只做小型离线矩阵运算，0 FE /0 RHS。正式重算由 `run_suite.py` 包装：`results/jobs/dev_endpoint_consequence_offline_complete/RECEIPT.json` 记录约.469 CPU秒、.790墙钟秒、轮询峰值RSS32.24MB，退出0。

两次失败前置尝试均保留 wrapper receipts：`dev_endpoint_consequence_offline` 为 cwd 路径错误；`dev_endpoint_consequence_offline_repair` 为旧 fixed endpoint 缺少 measurement 字段。随后统一读取合法 `inputs/online/dev_r0/observation.npz`，成功运行。两次失败均未运行 PDE。没有读取 blind truth、评分文件或 generator，没有变更方法或阈值。

## 8. 科学图与可复用表

`CONSEQUENCE_FIGURE/same_state_consequences.png`、同名 PDF 和 `same_state_consequences.csv` 已生成，图中区分 fixed dev4条、adaptive dev3条、adaptive refresh-trigger26条，并显示两个所选端点的完整/chart驻点表。33条均为同状态 full/approximate 对；图注明确 endpoint/FINAL 的记录可以重叠。已目视检查 PNG；修复初版对数横轴次刻度文字拥挤后，最终版刻度、图例、表格和公式均可读，无裁切。PDF由同一绘图对象导出。

图的两个面板分别画实际目标 Schur 后果与一阶预测余项，不把相关散点图当成机制因果分离。绘图只读上述离线 summary；两次绘图 wrapper 合计约3.657 CPU秒，0 FE /0 RHS。最后版 wrapper receipt 为 `results/jobs/dev_consequence_figure_layout_repair/RECEIPT.json`。四状态机制实验仍等待独立运行许可。
