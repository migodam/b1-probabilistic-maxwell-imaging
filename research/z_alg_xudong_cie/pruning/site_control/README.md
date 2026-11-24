# 单格位适用性对照

这是 12 件复杂盲测之外的独立、前瞻性适用性对照。四件新的强介质分别落在 16 个 24×24 格位中随机抽取的互异单格位，其余区域为空气。正对比度位于 \([2.5,4]\)，两轴实际支撑跨度须在 11–19 格。图形包含两种简单结构和两种非规则结构。冻结协议和代码保留各独立种子、缩放、抖动及盲标识。公开测量中没有材料类型或实际格位。 [冻结协议](FROZEN_PROTOCOL.json)

四件均沿用精确保存的原生 96²、400 MHz、半环和全环 10 发 20 收。每个阵列的 2% 和 10% 共用一个复高斯噪声方向，并沿用 \(10^{-6}\lVert Y\rVert_F\) 半径余量。合计 16 个拟合读数矩阵，单位是四件独立材料。先对实体非零格用原生 LU 生成读数，再用独立 LU、逆范数与 Arb 场残差复核真材料是否落入对应原生误差球。

**16/16** 个误差球通过区间复核，最小严格余量为 \(3.13187425\times10^{-6}\)。生成和复核合计 1.44 秒墙钟、1.08 秒 CPU，含四次生成 LU 和 60 个原生源右端。独立复核另含四次 LU、644 个逆矩阵右端、60 个原生源右端及 120 个接收伴随右端。 [公开范围内的生成与验证统计](ACQUISITION_RECEIPT.json)

默认位置 `.artifacts/z_alg_xudong_cie/site_control`。在线格位排除与 Born 排名只读 `public/manifest.json`、`public/geometry96.npz` 和 `public/observations/<case>.npz`。`offline/truth_key.json` 及真材数组由独立评分在全部在线预测封存之后读取。此协议的结果用于检验单格位前提成立时的证书与定位表现。

```sh
python -B research/z_alg_xudong_cie/pruning/site_control/generate.py --root /path/to/site_control
```
