# 复用哪个算子，放到哪台机器

**先选与材料、端口和边界一致的 forward，再选执行设备。** 查询命令同时给出实际函数、设备、工作区源码和计费边界。完整清单有 26 个入口；表中的 CPU/GPU 表示已有实现适用的位置。[算子登记表](OPERATORS.csv)

```bash
python3 research/experiment_handbook/lookup.py --operator list
python3 research/experiment_handbook/lookup.py --operator paper_cie_fft
python3 research/experiment_handbook/lookup.py --operator srap_vector_ap
```

|物理或干预|已有入口和用法|比较容易出错的地方|
|---|---|---|
|24² scalar-TM 密集 VIE，原物理目标、梯度、corrector|CPU `tm_dense`、`physical_corrector`；已有同类 CUDA `tm_cuda` 的 `SharedOperators` 和 `ObjectiveLane`。混合 TRF/SVD 有 `trf_hybrid`|共享 G 与背景，一个独立材料对应自己的场和分解。分解、forward/adjoint RHS、候选的原物理审计，以及候选被拒绝的代价都要计|
|FFT-TSOM、DPCS/APCS 与 CIE|小网格 `tsom_native` 或已实现的 `tsom_cuda`；文献 Austria 用 `paper_2d_forward` 的 `acquisition/fft_forward/make_model` 及 `paper_nie`、`paper_cie_fft`|冻结数据确定性电流、Fourier 模式和增长、sourcewise/global PR、β、归一化和停止量。确定性电流被锁定与仅截 Fourier 空间是两种限制|
|A+、G19、R14 双度量| `aplus` 是四格不重叠局部 Green 块、实盒条件读出；`g19` 用 `LowRankDenominator`；`stateful_dualmetric` 保存 PR 历史、阶段和辅助状态|记录实际构造基的内部 G 或观测 S、是否计入材料相关度量导数。current/material 两块度量、MM/Arc、current 初值各自可能改变轨迹|
|R12 AP/DR、R13 事件| `dr_ap` 的 `DenseFieldCurrentProjectors` 与 `run_projection`；阈值和起跳用 `material_events`|DR/AP 用相同的物理双场和共享材料读出；预分解、投影 RHS、AP 准备、反射或加速、原物理验收都收费|
|SRAP 二维半环、直线与 Fresnel| `srap_2d` 的 `route_basis/BallProjector/run_cie_arm`；实测 `measured_ports` 的空场拟合和每源接收算子|同端口数还要核对坐标、波型、缺失通道、掩码和 adjoint。跨源共用或逐源噪声球跟随具体程序。60² 时区分建立投影的费用、实际投影次数和到达阶段|
|固定 Green 的多源材料峰值证书| `z_alg_sector_inertia` 用 `prove` 认证矩阵负特征值数和源空间谱界，并以 `certify_dark_obstruction` 核对严格暗向量|共享材料把标量不等式推广到照明组合。若二次型有 ν 个负方向，界只去掉至多 ν 个源方向；同网格与正定版本比较。旧环上限 5 时，逐源计算失去此界，12 源联合仍严格排除。[推导和复算](../z_alg_secant_tasks/multiple-illuminations.html)|
|固定 Green 的端口对称性与逐组发射| `z_alg_orbit_exclusion` 核对收发闭合并回放区间判据；`z_alg_exposure_prefix` 以实际发射器为单位回放停止规则|按保存阵列的坐标和相位确定源组；使用预设误差上限，同时保留超合同发射费用及事后固定前缀对照。[复算入口](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/z_alg_secant_tasks/highdim/reflection/README_ZH.md)|
|三维 Maxwell|R1 `tensor_dda` 检查小型各向异性、横向极化与 Mie。R8 `vector_paper` 用三维 Fourier PR；SRAP `srap_vector` 及 `srap_vector_ap` 用 dyadic FFT、CG 与端口投影|R1 张量自由度与 SRAP 标量 χ 乘三分量场不同；SRAP 三维观测的所有源共用一个 Frobenius 噪声球。分清发/收位置、极化通道与材料数，单列 GPU 显存、谱基、CG 残差及物理误差|
|NN 初值与 A17 投影交换|R8 `nn_features` 配 `nn_material`；A17 修复版 `a17_projected` 配同秩基和等时全梯度|NN 材料要按求解器建立合法状态。把同一初值同时给 direct、AP/DR 和 CIE 并记录转换损失。A17 须比较实际盒投影后可行收益和终点材料|

### 怎样在已有队列上运行

先查完整观测并恢复原数据。下面是对现有二维程序的调用示例，使用新指定的结果目录。这里的输出保存到已挂载外置工作区的本轮目录。阅读历史结果可直接使用只读复算器。[逐步实例](WALKTHROUGHS_ZH.md) · [SRAP 的原始复现和物理依赖](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/srap_cie_array_study_20261121/REPRODUCE.md)

```bash
python3 research/experiment_handbook/lookup.py srap_main24 --list-views
# 在具备原研究物理依赖的 Inv_SLAM 环境中执行。
SRAP_OUTPUT="$PWD/.artifacts/my_paired_run/predictions" \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
Gaussian/.venv_nn/bin/python research/srap_cie_array_study_20261121/run_study.py \
  --cohort confirmation --ids s2d_010__pw180_10x20 --arms dual,srap --workers 1
```

### Mac mini、i9 与 RTX 4060

|工作|设备选择及原因|发布时写什么|
|---|---|---|
|读冻结数组、逐对象评分、作图、小尺寸机制干预|Mac mini M4 可立即完成，也适合做单线程配对。已有计时争议先用保存状态重放|实际 BLAS 线程、wall/CPU，是否调用了新物理求解。按磁盘读取计费时和迭代 body 计时区分|
|多对象的独立 CPU 反演|XINAN i9-14900HX 有 24 核/32 线程。可同时利用 P 核、E 核及 SMT 提高总吞吐。等墙钟的同一个对象，各方法给相同核心类别、线程和可比负载；分批并轮转方法顺序|记录每对象完成的阶段、原物理审计次数、单臂 CPU/wall 与整批历时。避免把异构并发的超时分配当成物理优势|
|重复的 576 维 dense LU/adjoint、三维 dyadic FFT、已有批求解或训练|远端 RTX 4060 的已有 CUDA 实现优先。先检查占用、可用显存和所需精度。一次小作业若传输和建图时间占主导，可继续用 CPU|记录 device、实际 dtype、同步后的 GPU/wall、显存峰值、传输、预计算和摊销对象数；单列冷启动与共享热缓存|
|对设备表现有冲突|同一公开算子、同一对象、同一调度先看完成工作量。旧 U11 在两机甚至运行了不同秩和阶段；修正版本历史上，两机采样的单线程审计/FFT 用时仍差 2.18–2.39 倍|如需要解释 i9 的 1.6 GHz，读活跃核心有效频率、封装功率、供电/温度限制与实际耗时。全核页的平均时钟不能解释单线程排名；本手册没有更改系统功耗设置|

[旧 U11 机器与版本对照](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/srap_cie_array_study_20261121/engine_dev/HARDWARE.md) · [源代码与停机量勘误](BASELINES_AND_ERRATA_ZH.md)

遇到新物理类别时，先用现有参照验证最可能改变结论的环节。二维近场的源依赖 mask 不可交给任意矩形 Y 的“全端口”模型。三维可复用球的 Mie 或一次有必要的细化检查；已有鉴定只覆盖原网格、电尺寸、材料和误差预算。没有任务相关差异时，直接复用现有回执。
