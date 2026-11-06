# 原尺寸 matched-grid 短原物理校正（development）

三个起点是 air、已公开 GlobalDual stage2、已公开 AplusGlobal stage2。原九臂 GT 已暴露，本组不是新的确认实验。单驻留 complex128 packed16 CUDA 算子；GMRES 仍逐源，未使用或宣称共享 LU。

共同原目标为 mean_source ||SJ−Y||²/||Y||²，无 1/2、无先验，χ∈[0,4]。L-BFGS-B 固定最多 8 个 accepted 步、24 次请求、12 次 line-search，ftol=1e−12、gtol=1e−8；所有拒绝 trial 都付正/伴随费用。按每臂所有已付状态的最低公开 Phi 选择，不按 GT。

|起点|初始→最终 Phi|初始→最终 PG∞|初始→最终 contrast NRMSE|accepted / 请求|新增正/伴随 RHS|CPU s / wall s|
|---|---:|---:|---:|---:|---:|---:|
|air|1 → 0.648519|0.0005831 → 0.0005832|1 → 0.975|8 / 11|600 / 600|118.906 / 120.328|
|GlobalDual_stage2|0.141197 → 0.0277534|0.0004009 → 0.0001406|0.4228 → 0.249|8 / 9|540 / 540|717.672 / 728.720|
|AplusGlobal_stage2|0.142384 → 0.0278043|0.0004056 → 0.0001284|0.4235 → 0.2495|8 / 9|540 / 540|714.875 / 725.892|

请求数与新求解数不同：initial 在 SciPy fun 前计算，fun 首请求命中相同材料缓存；air 初值已由梯度接线支付，解析 J=0 与伴随无 GMRES RHS。每次非零新 value+gradient 付 60 primal + 60 adjoint；回调读取已付状态，不再计算。上表为每臂此次新增实际归属费用，不能再加到整进程费用。两个 stage2 起点的既有 CIE 运行费用仍在旧九臂账；本表不是包含这些前缀的独立冷启动总成本，也不能把前缀重加为新消费。

梯度接线只在公开 air 点做预定单向 h=1e−5 检查；无必要 h/2，不复做旧恒等式。三个终点及每个 accepted 材料/预测全部封存后，才统一离线读材料 GT。contrast NRMSE 分母为 ||χ_truth||₂；paper Err 为 sqrt(mean((χ−χ_truth)²/(1+χ_truth)²))。PG∞ 是 ||χ−clip(χ−g,0,4)||∞，使用未重标定的原目标梯度。

本组只能检验有限 8 步物理方向是否继续改善这些公开端点，不能证明最终驻点、唯一性、全局最优，或依靠继续迭代拯救原 CIE。小梯度、残差下降与材料恢复分别报告。保留旧 CM+RR 高介电/原尺寸连续误差限制。

结果：RESULTS/original_vector3d/physical_refinement/{SCORES.json,TABLE.csv,GRADIENT_WIRING.json,RUN_PROCESS_COST.json,ACCOUNTING_FINAL.json}。图：FIGURES/original_vector3d/physical_refinement/{physical_curves.png,materials.png}。

材料图下排使用固定 z=-0.65λ（index8）切片，对应已知论文立方体 z=-0.7λ 附近；中心 z=0.05λ 原本位于两层目标间隙，首版图另存 materials_attempt1_central_plane.png。图仅离线展示，不选迭代/材料/评分；颜色范围固定 χ∈[0,4]。

结论边界：这两个已公开的辅助端点并未耗尽原物理改进方向；相同八步使二者的 contrast NRMSE 都降至约0.249，未显示 A+ 独有收益。air 八步仍为0.975。本组不能区分更长预算下的最终 basin/驻点，也不是原 CIE 加迭代的重新确认实验。以上为 matched-grid 有限 CM+RR 模型，不能转述为原尺寸 cross-grid 或连续 Maxwell 恢复。

费用：单驻留 GPU 主进程 1557.953125 CPU 秒 / 1583.523520 wall 秒，包括导数接线、共享 setup、全部三臂、拒绝trial与输出；每源 GMRES，无 LU。实际1740 primal+1680 adjoint RHS，其中接线60 primal，28个非零新VG共3360 RHS。一个解析air gradient无GMRES RHS，四次完整状态缓存调用不重复收费。perarm严格请求11/9/9，新增非零VG10/9/9；另有接线的解析air与value-only FD。

canonical ACCOUNTING_FINAL.json 对全GPU进程只加一次，局部图/运输/汇总独立相加，费用保持已知下界。第一次本机评分的Windows路径解析guard在GT读取前失败；0物理RHS、CPU未知、wall0.348秒，原始预测未改。绘图首版中心空隙切片与第二版色条裁切图及各自费用均保留。尚未计量的连接/AST预检CPU和进程输出尾部标null/下界，不发明零成本。
