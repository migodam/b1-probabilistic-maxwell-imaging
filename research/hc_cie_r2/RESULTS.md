# 逐项结果与原始数据

下表的分母分别是材料、端点或观测条件。一个材料的多源、多预算及多噪声条件始终配对。噪声球半径用各观测合同中的 $\delta$。

|问题|实际范围与结果|可以直接复查的证据|
|---|---|---|
|小模型能否精确消掉非线性|三格软数据，五格精确数据纤维，32 个同秩随机对照，多接收布置、峰值 1–100；五格特定几何可线性恢复，直接真实材料目标另有区间认证的驻点|[模型、优化与区间结论](TOY_AND_CERTIFICATES.md) · [五格全部起点](results/toy/hard5.json)|
|大网格有多少线性方程|21 观测、10 独立材料、570 个照明和截断条件。物理 Gram 秩分别为 575/576、1023/1024，精确余核仅为常数|[谱与误差界](SPECTRAL_RESULTS.md) · [完整谱的数组键](results/spectra/array_index.json) · [逐条件 CSV](results/spectra/spectrum.csv)|
|这些方程中有多少高于原测量噪声|五个代表观测、两种接收截断、三档公开电流幅度。R13 半环最多 45/6804，R13 全环与 SRAP 同侧线均为零|[逐条件有效秩](results/spectra/noise_rank.csv) · [实际协方差和谱](results/spectra/noise_rank_arrays.npz)|
|保留余项或实性是否有用|231 个截断投影条件；3 个真实球内材料另做复核；六个材料上实际部署 36 个前端和校正臂。实权没有提高四个半环校正后的恢复|[部分消元](results/spectra/partial.csv) · [球内反例](results/spectra/feasible_followup.csv) · [前端逐案](results/real_frontend/scores.csv)|
|可以排除多局部的修复|48 个父材料 × 五个历史端点的基础下界。八个 Dual-CIE 端点逐对枚举，两案普通秩界不足，三原子的几何界将其排除|[基础下界](results/support/archived_minimum.csv) · [所有位置对的最优距离](results/support/all_pairs.csv) · [两案保守数值界](results/support/all_pairs_verified.json)|
|有限跃迁能否恢复整图|十件暴露过的材料。修订版四类起点、七方法，保留原来的 7/40 个合格端点，未增加合格解|[端点与候选](results/imaging/score_v2/OUTCOMES.csv) · [前景和背景](results/imaging/score_v2/STRUCTURE.csv) · [实际墙钟](results/imaging/score_v2/WALL_PREFIXES.csv)|
|真实负曲率能否跨盆地|八件、32 臂。五个共败 CIE 端点都找到真实负曲率，有限盒内跳步未救回这些对象|[曲率及活动面](LANDSCAPE.md) · [方向数据](results/landscape/score_v1/DIRECTIONS.csv)|
|已知背景内的两格变化|24 个参考材料，每件两种变化。10% 差场时前两名 MUSIC 42/48；四候选单格拟合、交叉拟合、闭式读出各 47/48|[原始分数](results/interaction/offline.csv) · [闭式与单格见证](results/support/monitor_offline.csv) · [按独立父材料的配对不确定性](results/interaction/paired_object_statistics.json)|
|跨网格与孔径|六个 SRAP 父材料，两阵、两种变化、四测量条件。10% 跨网格差场的闭式读出 16/24，普通 MUSIC 同为 16/24。两者在同一个体上可能定位不同|[差分监测](results/cross_changes/offline.csv) · [粗化真材质与测量](results/model_gap/gaps.csv) · [原物理拟合](results/model_gap/scores.csv)|
|经典合成与公开实测|Austria 同形三观测 × NIE/Dual，FoamTwinDiel 一实物 × 两起点。实物的公平终点取实际五个留出源。没有实物像素真值|[训练与留出原始成绩](results/classical/score_v2/OUTCOMES.csv) · [实际复端口数组](results/classical/score_v2/HELDOUT_PREDICTIONS.npz)|
|三维同一改变能否被看见|四个已知三维参考、两种阵列、三个噪声条件。10% 差场的半球字典 4/4、面阵 0/4；干净面阵的两侧修正 4/4|[三维评分](results/vector_change/scores.csv) · [各极化和完整候选](results/vector_change/report.json) · [两侧子空间复算](results/vector_change/two_sided.json)|

同色标重建、物理 Pareto 和失败图可从相应章节点开。所有源均在同一批量里进出的任务给整批墙钟。使用旧端点的研究另外标明其既有费用。[汇总费用](COST.md) · [全部图像的文件目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/hc_cie_r2/figures)
