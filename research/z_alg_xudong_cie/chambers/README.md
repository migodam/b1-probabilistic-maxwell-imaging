# 暗室观测输入和独立评分

开发使用 Fresnel Twin 8 GHz 原 `half180_tx`。确认保留另外两个 Fresnel 实物的 6 GHz 观测，以及 Manitoba E 型聚乙烯件的 3.5 GHz 观测。每个接收端口的坐标沿用旧模型，训练与测试的**发射源和发收对**不重叠。每个训练源使用原档案中该源的全部实测接收端口，因此这四个固定拆分没有独立的仅接收端留出。以下为一基端口号。 [训练布局回执](training_receipt.json)

| 观测 | 用途 | 网格 | 训练发射端 | 留出发射端 | 训练 / 留出复读数 |
|---|---|---:|---|---|---:|
| `fresnel_twin_8GHz_half180_tx` | B2 开发 | 60² | 1–9 | 10–18 | 2,169 / 2,169 |
| `fresnel_int_6GHz_half180_tx` | B3 确认 | 60² | 1–4 | 5–8 | 964 / 964 |
| `fresnel_ext_6GHz_half180_tx` | B3 确认 | 60² | 1–4 | 5–8 | 964 / 964 |
| `manitoba_e_3p5GHz_even_tx` | B3 确认 | 40² | 偶数 2–24 | 奇数 1–23 | 276 / 276 |

## 训练输入

`training.load_training(case_id)` 返回公开训练阵列 `E/B/E0`（同一阵列）、`Y`、全接收几何算子 `S`、逐源零基 `rx_ids`、`training_rx_mask`、`valid_mask`、`centers`、`grid_shape`、物理波数 `k`、等面积脉冲圆盘半径 `radius`、经验残差尺度 `delta` 和 `fft_geometry`。`Y` 为接收行 × 源列，`rx_ids` 为源 × 各自的接收行。`training.make_model(bundle)` 装配历史 R8 逐源接收的矩阵自由完整 VIE，成熟 Dual-CIE 直接使用这里的 `model`、`Y` 与 `fft_geometry`。Manitoba 在物理正演使用原保存 FFT 符号，Fresnel 使用旧 60² 等面积圆盘几何。 [训练 API](training.py)

Fresnel 重用旧训练与留出档案中的空场圆柱谐波拟合和原 `half180_tx` 掩码。校准不看目标散射。`delta=0.15 ‖Y_train‖_F` 沿用旧经验合同，未验证为硬误差上界。 [原训练拆分](../../limited_aperture_routed_dr_20261120/run_la_study.py) · [Fresnel 采集论文](https://doi.org/10.1088/0266-5611/21/6/S09)

Manitoba 使用第一次空场和 3.5 英寸金属标定柱。逐端口复增益为该柱的 PEC 解析散射除以测量散射。只在训练掩码复算校准，按原始增益式回算的 E 型训练读数与本地存档逐元素一致。`delta` 取第二次空场漂移范数 0.01178485 加留出 2 英寸钢柱残差范数 0.18548350，约占训练散射场范数的 26.95%。这是经验敏感性量。介质目标三维转移误差、真实标定误差和二维近似误差未得到统一上界。金属留出本身也没有参与标定增益拟合。 [存档标定与算子](../../z_alg_optical_inertia/validation/manitoba/physical.py) · [Manitoba 采集论文](https://doi.org/10.1109/MAP.2011.6138442)

打包程序在项目外置盘的数据目录读原始及既有校准输入，只写出偶数源的训练 NPZ。成熟引擎的 `load_training` 只读这些训练文件，不加载留出目标或名义材料。原始读数、私有完整标定包、本地训练包与残差数组沿用原数据授权，发布只带代码、统计及图。

## 独立评分

`score.score_material(case_id, material, freeze_record=…)` 接受冻结的原生材料，独立完成完整 VIE，再报告训练复场、未见发射端复场的绝对和相对 Frobenius 误差及每源误差。用原始掩码排除未采集通道，预测禁止增益重拟。计入实际正演费用和线性状态缺陷。确认评分入口需要记录已冻结算法的文件。独立确认的三件实测使用相同评分器。逐例训练及留出复场见本轮的[逐件结果](../results.html)。 [评分 API](score.py)

名义评价只使用离线的公开尺寸参考。Fresnel 沿用历史 R8 图解坐标、论文 Figure 5 的泡沫和塑料圆柱，名义对比度 0.45 与 2。Manitoba 沿用存档 Figure 6b 的 E 型截面，名义对比度 1.3 取自论文在 5 GHz 的期待值。两室均无与复场坐标注册的逐格材料真值。内部取固定名义支撑，背景取名义零对比度格，IoU 对预测与名义支撑应用同一个固定对比度阈值（Fresnel 0.225，Manitoba 0.65）。这些数用于同一参考下的几何诊断；实测横向比较以留出复场为准。 [Fresnel 名义构造](../../z_alg_optical_inertia/fresnel/physical.py) · [Manitoba 名义构造](../../z_alg_optical_inertia/validation/manitoba/physical.py)

## Twin 8 GHz 开发终点

冻结的成熟 Dual-CIE 第 533 次更新用全部 18 个源重新正演。训练相对误差为 **25.19%**，独立源误差为 **26.30%**。相应绝对误差分别为 0.322106 和 0.334111。训练误差超过 0.191771 的经验半径。与同一未注册名义参考比，内部相对误差 28.57%，名义空气 RMSE 0.08006，支撑 IoU 0.8240。终点正演花 0.36 秒，18 个源右端，状态相对缺陷 \(8.29\times10^{-10}\)。 [终点评分](dev_twin_terminal.json) · [阶段评分](dev_twin_stage_scores.json) · [稀疏残差与材料图 PDF](dev_twin_terminal_map.pdf)

## 重放

从项目根执行。需要 NumPy、SciPy、threadpoolctl，重画图另需 Matplotlib。Manitoba 原包和校准引用位于本地 `.artifacts`。Fresnel 的原始维护者和授权见项目原档案。以下命令只整理训练视图及重放开发评分。

```sh
python -B research/z_alg_xudong_cie/chambers/training.py
python -B research/z_alg_xudong_cie/chambers/development.py \
  --local-terminal-map .artifacts/z_alg_xudong_cie/chambers/dev_twin_terminal_local.npz
python -B research/z_alg_xudong_cie/chambers/figure.py
```
