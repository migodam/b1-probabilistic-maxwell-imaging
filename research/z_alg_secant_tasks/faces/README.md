# 规格盒两侧极值复算

每行对应一套保存的材料中心、每槽名义容差半宽和物理端口。72 行沿用 `array/results/` 的六个中心、四种采集和各自原 10% 噪声半径。`results.json` 保留两端距中心预测、端点间距及相容余量的区间。所有场距都是冻结离散算子的复端口向量二范数。

`face_parameters.npz` 中 `<case>__center` 为中心材料，`<case>__ROI` 为布尔任务掩码。盒的合法材料范围是 `[0.5,30]`，按中心余量剪裁名义半宽。纯 ROI 对照保持区外为中心。补偿端点的区外绝对材料值存于 `<case>__<width>__<arm>__LS_outside_±1`，补充 GN 存于相同结构的 `GN_outside_±1`。区内参数在复算时以记录中的二进制中心和半宽重新精确组装为 `c±r`，四舍五入后的打印小数不用于认证。

每个端点由保存的完整 `G/S/B` 用 128 位 Arb 复数球重解 VIE。`remaining_radius_about_F0` 是同一中心读数球的余量 `δ−max(||F₊−F₀||, ||F₋−F₀||)`。正下界保证以 `F₀` 为中心、半径为该下界的全部读数都保有完整先验任务跨度。`remaining_radius_about_midpoint` 给出中心可自由选择时这一对端点的余量。三组尝试都保留，依次为固定区外、一次区外有界线性拟合，以及线性结果失配时对两侧分别追加一次普通 GN。

从项目包根目录运行，需要 Python、NumPy、SciPy、python-flint 和 threadpoolctl。该复算读取公开冻结数组，生成 `faces/results.json` 和 `faces/face_parameters.npz`。

```sh
python faces/face_modulus.py
```
