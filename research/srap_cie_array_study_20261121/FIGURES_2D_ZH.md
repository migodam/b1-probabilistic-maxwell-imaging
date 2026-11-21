# 二维阵列研究图集

## 确认结果

| 图组 | 内容 |
|---|---|
| [24 对象浏览页](figures/confirmation_atlas_30s/index.html) · [整册 PDF](figures/confirmation_atlas_30s.pdf) | 每页一个对象、七种阵列、参考材料、Dual-CIE 与 SRAP-CIE 的 30 秒 raw/light |
| [有限柱面源截断配对图](figures/line_cap_examples_30s_light.pdf) | 预先固定 001、012、024，在同侧直线、对侧直线、同角弧上比较原版与奇异谱截断的 30 秒 light |
| [误差随时间正反例](figures/confirmation_trajectory_examples_30s.pdf) | 确认组 30 秒 raw 配对相对 L2 降幅最大和最小的案例，同时画实际物理残差及路由诊断 |
| [补充确认组 PDF](figures/followup_materials_30s.pdf) | 两种 120°、180° 及翻转照明四阵列在 8 个薄桥或连续材料上的 30 秒 raw/light |

## 历史与实测

| 图组 | 内容 |
|---|---|
| [历史材料对照](figures/legacy_materials_30s.pdf) | 12 个合成、6 个实测，Dual-CIE、SRAP-CIE 及旧 U11 路由 AP |
| [历史两机场景](figures/legacy_cross_host_30s_light.pdf) · [历史机制轨迹](figures/legacy_trajectory_examples_30s.pdf) | 三个固定场景的跨机器轻校正图，以及路由、直接投影、动量重启和图相容投影的历史对照 |
| [历史两对象孔径×密度](figures/development_materials_30s.pdf) | 180°/120° 及 10×20/7×14，含截断对照 |
| [Fresnel 匹配掩码](figures/fresnel_twin_materials_30s.pdf) · [0–1 色标](figures/fresnel_twin_materials_30s_low.pdf) | 两种 8 GHz 掩码每臂各 733 个训练通道，并共用完全相同的 TX/RX 留出目标 |
| [Mac 实测延迟组装](figures/fresnel_lazy_materials_30s_full.pdf) · [0–1 色标](figures/fresnel_lazy_materials_30s_low.pdf) | 八种 Fresnel 布置的 Dual-CIE、SRAP-CIE、SRAP-CIE lazy，均标注未截断峰值与可用的留出复场误差 |

## 布阵与算子

[七阵列](figures/planar_arrays.pdf) · [补组布阵](figures/followup_arrays.pdf) · [Fresnel 掩码](figures/fresnel_twin_masks.pdf) · [接收谱](figures/receiver_spectra.pdf)

材料图用同一 0–6 介电对比度色标。标为低对比的实测页用 0–1 色标，峰值取原始材料数组。L2 为相对材料误差，IoU 为参考峰值一半阈值下的区域交并比。实测无逐像素真值，图下 TX/RX 为相应留出侧的复场相对误差。轨迹用实际完成物理审计的时点，空心标记是算法按原物理目标选出的 raw 候选。

在本研究目录且评分已完成后运行。相同脚本也接受 `--budget 15`。图集脚本读取现存预测和评分。

```sh
python figures2d.py --section schematics legacy development measured_extra fresnel_lazy legacy_cross_host confirmation
SRAP_DATA="$PWD/data_followup" python figures2d.py --section followup
```
