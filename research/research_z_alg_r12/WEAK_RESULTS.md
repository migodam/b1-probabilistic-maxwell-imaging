# 576 维低对比材料的近似歧义与补测

## 数据和物理

使用 GN–NN R4 历史 TEST 的九个独立母体。原 TEST 三个分组和三个对比幅度 \(A\in\{0.1,0.2,0.4\}\) 的每个组合取一例。原始含噪数据及逐尺度 discrepancy GN 端点直接读自档案。本轮研究复用了已揭晓的 TEST。 [GN–NN R4 报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_gn_nn_preservation_r4/) · [原始数据与模型 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/gn-nn-preservation-r4-20261107)

物理模型为二维标量 TM 盘脉冲 VIE。24×24 个材料未知数和 48 个接收器，原 12 发射位置均匀分布于 \(r=1.6\lambda\) 圆周。\(k=2\pi\)，材料盘半径 \(0.018\lambda\)，原相对噪声幅度 \(0.001\)。回放使用原模型的块算子和 full-J，保留 \([0.02,20]\) 的材料盒约束。

补测使用相同接收器、波数和非线性 VIE。原半径的新增位置在原发射站位中间，整体错开 15°。近场对照把同样的新角度放到 \(r=0.72\lambda\)，源到最近材料盘边的最小距离为 \(0.05146\lambda\)。一发、四发分别在各自等频、等时长的单位单极源系数下比较，系数平方给出发射能量的相对标尺。旧位置重测使用独立的噪声。新增观测由同一原档材料通过同一物理模型生成，接收噪声为复高斯方向归一化到旧档推定的逐源绝对范数 \(0.001\|Y_{\rm old}\|_F/\sqrt{12}\)。不同重建方法在每个选源臂看到完全相同的新增数据。

## 有限材料差和来源分解

在原 GN 端点形成实际全 Jacobian，选取奇异值小于最大值 \(2\times10^{-4}\) 的材料子空间，在其中求低粗糙度方向。沿该方向正负扰动并投影回材料盒。试探的首个目标扰动幅度为材料 RMS \(0.05\)，必要时沿旧资料的强方向校正，并按脚本中的固定序列缩小试探。端点是否相容只看完整 VIE 对原观测的相对残差。生成这些端点和选源时使用公开旧观测与算子，真值用于在选源后模拟新数据和评分。

九例的两端都满足 \(0.0011\) 的旧资料阈值，实际残差介于 \(0.0009764\) 和 \(0.0010970\)。候选间材料 RMSE 为 \(0.04898\) 至 \(0.09993\)，中位 \(0.07395\)，相对原材料对比幅度的中位数为 \(0.495\)。旧数据中两个候选的响应差除旧噪声范数的中位数为 \(0.2663\)。

对每个新源，将离散入射场正交分解为 \(E_{\mathrm{new}}=E_{\mathrm{old}}c+b_\perp\)。在两个有限材料固定时，完整 VIE 关于入射场线性。因此用两端各自的正演计算 \(b_\perp\) 对应的有限响应差，再与旧源组合的响应差相加。三者以复向量相加，相对闭合误差在 \(3\times10^{-12}\) 内。表中的能量标尺为各自新观测噪声向量范数。

|选源条件|一发响应差中位数|四发响应差中位数|四发旧源空间外响应差中位数|
|---|---:|---:|---:|
|原位重测|0.2571|0.2508|约 0|
|原半径新增角度，按候选差选源|0.2957|0.2825|0.0773|
|近场，按候选差选源，同发射能量|2.2329|1.7607|1.7224|
|同组近场位置，将入射场范数降至远场水平|1.4131|1.1188|1.0932|

等入射强度消融逐角计算 \(a=\|E_{\mathrm{far}}\|_F/\|E_{\mathrm{near}}\|_F\)，近场系数为 \(0.6328\) 至 \(0.6541\)。接收噪声的绝对幅度固定。因正演对入射场线性，乘 \(a\) 同时给出真实全波响应差及其投影分量。此时近场发射能量为远场的 \(0.4005\) 至 \(0.4279\)。等场强与等发射能量是不同的物理比较。四发和一发使用各自的合并新噪声范数，表内无量纲数字可比较各预算内的不同几何。

四发同半径的随机、局部 Jacobian 和均匀选择响应差中位数为 \(0.2736,0.2765,0.2825\)。近场按离材料盘最近和按局部弱空间灵敏度的四源集合在九例中相同，响应差中位数为 \(1.7248\)。近场随机和均匀方位的中位数为 \(1.2801\) 和 \(0.8552\)。

以四发新数据联合旧观测，在 \(1.10\) 倍名义合并噪声范数内，原半径选源和原位重复均有 8/9 保留两个示例材料。近场候选选源有 5/9 排除两端、4/9 只留一端。单发近场相应为 1/9、3/9、5/9。响应差与噪声范数的比例衡量总能量，不足以推断统计分类的错误率或所有材料分支的可辨识性。

[有限候选材料图](figures/weak_candidates.png) · [补测和旧源空间分解图](figures/weak_source_geometry.png)

## 同观测材料恢复

四次近场补测后，共享这批数据的普通冷启动正则 GN 的中位 RMSE 从旧观测的 \(0.03159\) 降至 \(0.02652\)，Otsu-IoU 从 \(0.6054\) 到 \(0.6804\)，逐样本 RMSE 改善 6/9。普通冷启动正则 GN 用同半径候选数据的 RMSE 为 \(0.03615\)，用原位重复数据为 \(0.03783\)。一发近场候选数据的普通冷 GN 中位 RMSE 为 \(0.03551\)，改善 3/9。

归档逐尺度 GN 端点热启动的正则 GN 四发近场后中位 RMSE 为 \(0.01984\)，旧端点为 \(0.01847\)，IoU 中位均为 \(0.7029\)。从近场选出的有限分支出发做相同正则校正，9/9 通过 discrepancy，中位 RMSE/IoU 为 \(0.03584/0.5417\)。这里测到的新信息能排除旧相容材料，但有限分支作为反演初值没有转化为相对普通 GN 的恢复优势。

强基线复用 GN–NN R4 原 seed 7 Voxel 和 J-full Cascade 的冻结材料提案。两个网络按原 12 源的严格 GN 锚点推理，保持训练时的输入结构。用归档权重在 CPU 重放全部九例，材料与原提案的最大坐标差分别为 \(7.15\times10^{-7}\) 和 \(6.71\times10^{-7}\)。同观测对照采用归档提案，在每个源臂直接评估其联合 VIE 残差，并从该提案出发接受与有限分支相同的 14 步上限正则 full-J 和 \(1.10\) 倍合并噪声停止。另以 8→12→24 网格变化率、全材料 Jacobian 和零材料正则做纯物理 voxel 冷启动。

|方法|仅旧观测 RMSE / IoU|四发原位重测|四发同半径|四发近场|
|---|---:|---:|---:|---:|
|R4 Voxel 原提案|.02071 / .7606|.02071 / .7606|.02071 / .7606|.02071 / .7606|
|R4 J-full 原提案|.02664 / .6429|.02664 / .6429|.02664 / .6429|.02664 / .6429|
|R4 Voxel 同校正|.03443 / .6149|.02880 / .7029|.03310 / .6767|.02756 / .7068|
|R4 J-full 同校正|.03825 / .5723|.03258 / .6319|.03471 / .7000|.03109 / .6842|
|冷启动无正则 voxel full-J|.04345 / .4628|.04044 / .5284|.04329 / .5644|.03491 / .6017|

两套原提案在表中各组均有 0/9 通过 discrepancy。四发近场下，残差中位分别为噪声范数的 20.77 和 42.65 倍。两套校正和纯物理 voxel 在表中各组均有 9/9 通过。四发近场下，Voxel 和 J-full 的联合残差中位分别为噪声的 1.016 和 1.042 倍。冷启动无正则 voxel 在四发近场也降低了材料 RMSE。一发近场下，这三个校正或优化基线的 RMSE 中位分别为 .02986、.03566、.05289。

构造有限材料候选的中位成本为 0.43 秒、16 次 576 维 LU 和 960 个复 RHS。两有限端点预测全部 36 个补测条件的中位成本约 0.036 秒、2 次 LU 和 72 个复 RHS。原入射场空间分解的诊断成本单列。四发近场的 Voxel 和 J-full 同校正分别用 2 LU / 80 RHS 和 4 LU / 160 RHS，当前单线程 CPU 实测中位约 0.047 和 0.062 秒。原严格 GN 锚点在其历史环境的中位成本另为 45.16 秒、233 LU / 12,684 RHS，原网络在锚点之后的准备与推理时间另计。冷启动无正则 voxel 在四发近场的中位为 0.86 秒、61 LU / 1,744 RHS。

## 文件与复现

- [有限候选与补测主程序](src/weak_r4.py) · [同观测强基线](src/weak_baselines.py) · [绘图](src/weak_figures.py)
- [四发原始结果](results/weak_r4.json) · [材料数组、旧观测与全部模拟补测](results/weak_r4.npz) · [四发统计](results/weak_summary.csv)
- [一发原始结果](results/weak_r4_one.json) · [一发数组](results/weak_r4_one.npz) · [一发统计](results/weak_summary_1source.csv)
- [四发强基线详细结果](results/weak_baselines.json) · [四发基线数组](results/weak_baselines.npz) · [四发基线统计](results/weak_baselines.csv) · [一发强基线详细结果](results/weak_baselines_one.json) · [一发基线数组](results/weak_baselines_one.npz) · [一发基线统计](results/weak_baselines_one.csv)
- [独立复核输入](results/weak_inputs.npz)包括原几何和算子、全部新增源场、九份材料真值和 alpha 掩模、原观测与本轮新观测。新增源索引 0–11 对应同半径新增角度，12–23 对应近场角度，24–35 对应原位重测。

公开 R4 Release 的原始观测、原 GN 端点、真值和 checkpoint 在附件里。保留 R12 目录的层级，在 R12 目录执行以下命令。若已在本地恢复原 R4 档案，直接将 `--archive` 指向其根目录。

```sh
git clone --depth 1 --filter=blob:none --sparse --branch gn-nn-preservation-r4-20261107 \
  https://github.com/migodam/b1-probabilistic-maxwell-imaging.git r4-source
git -C r4-source sparse-checkout set research/research_gn_nn_preservation_r4/src

gh release download gn-nn-preservation-r4-20261107 \
  --repo migodam/b1-probabilistic-maxwell-imaging --dir r4-downloads \
  --pattern gn-nn-r4-data-public.part001.zip \
  --pattern gn-nn-r4-data-evaluator.part001.zip \
  --pattern gn-nn-r4-anchors.part001.zip \
  --pattern gn-nn-r4-frozen-python-dependencies.part001.zip \
  --pattern gn-nn-r4-models.part001.zip

for part in r4-downloads/*.zip; do
  unzip -oq "$part" -d r4-source
done

python src/weak_r4.py --archive r4-source
python src/weak_r4.py --archive r4-source --budget 1 --output results/weak_r4_one
python src/weak_baselines.py --archive r4-source --replay-checkpoints
python src/weak_baselines.py --archive r4-source --budget 1 --replay-checkpoints
python src/weak_figures.py
python src/weak_figures.py --results results/weak_r4_one --summary-only
```

需要 NumPy、SciPy、threadpoolctl、Matplotlib，以及重放旧 checkpoint 时用的 PyTorch。历史字节级预测可另外下载 Release 中的 `gn-nn-r4-predictions.part00*.zip` 解压到同一根目录。若这些预测文件存在，强基线使用原预测；`--replay-checkpoints` 同时运行冻结 checkpoint 并报告重放差异。主实验按旧观测形成有限候选并选择采集，再用新增数据做反演，复算时不会重训网络。
