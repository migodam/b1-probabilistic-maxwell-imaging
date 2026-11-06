# Right3：暴露对象上的posthoc结构诊断

**EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC。** 三条right/rank19/inner8已完整封存后按root授权评分。原168材料指标全部复用，不重打分、不改原结果；GT仅使用原24已解封缓存。新增评分Maxwell forward/LU均为0，没有读取conditional4新GT。

## 末stage、原Phi起点与终点

| 对象 | GP19原生stage3 NRMSE | 原Phi选中corrector起点 | GP19最终NRMSE / fixedIoU | 冻结PGP19最终NRMSE / fixedIoU | GP19主材料 / CN99 |
|---|---:|---|---|---|---|
| 0001 annulus | .139481 | stage3，.139481 | .101668 / .739130 | .098895 / .755556 | PASS / PASS |
| 0004 slotted shell | .550686 | stage1，.507647 | .515137 / .086957 | .515137 / .086957 | FAIL / FAIL |
| 0022 bridge | .403915 | stage1，.260731 | .436076 / .116788 | .337003 / .330357 | FAIL / FAIL |

GP19仍是1/3材料成功，与这三个暴露对象上的PGP19相同。0001保留明显救回；0004没有改善；0022的材料结果反而更差。GP19的最终Phi依次为5.24975e−5、.17393422、.01039669；0022虽低于PGP19的.01420888，GT误差却更大。**去掉左投影不是这两个失败对象的充分补救，也不能用较低Phi替代材料证据。** 这是posthoc三例结构比较，不是新确认结论或对任意GP/PGP路线的普遍判定。

0001与0022的right corrector达到3000迭代上限；其保存PG∞分别约4.83e−8和1.82e−6，不能写成完全converged。0004以投影梯度准则在560次迭代停止，保存PG∞约8.88e−9，但仍是高Phi坏数据拟合点。三个对象的全部四阶段、公开起点及final均保留，不用GT选择“最好stage”。

## 已付费用

| GP19对象 | 完整child GetProcessTimes CPU(s) | 原物理576-LU | 原物理复RHS | 19维小LU | 小系统复RHS |
|---|---:|---:|---:|---:|---:|
| 0001 | 281.640625 | 3110 | 74592 | 33761 | 810264 |
| 0004 | 137.765625 | 585 | 13980 | 27367 | 656808 |
| 0022 | 282.234375 | 3078 | 73812 | 33612 | 806688 |
| 合计 | 701.640625 | 6773 | 162384 | 94740 | 2273760 |

这些inverse已由root/producer运行，本scorer没有新增inverse。评分及整图渲染另付CPU1.170654秒/wall1.585744秒、0 forward。原CIE/A+/PGP19费用在对照表中明确标`REUSED_FROZEN_REFERENCE_NOT_NEW_SPENDING`，不因重复展示再收费。内部R10/R8/aux/corrector时间是nested视图；不同维LU不混计。端点不同停止次数，不能称同CPU比较。

交付：`FINAL_COMPARISON.csv`、`ALL_STAGE_START_FINAL_SCORES.csv`、`COST_LU_RHS.csv`、`ALL3_PUBLIC_SEAL.json`、`RECEIPT.json`；`MATERIAL_COMPARISON.png/.pdf`包含三个对象的GT/CIE/A+/PGP19/GP19，已目检。源码入口为`experiments/score_right_projection.py`。

本次仅完成right3，conditional4保持未读GT，等待root另行授权。
