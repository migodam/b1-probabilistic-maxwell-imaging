# Weak-scattering controls for the Born data test

A material with peak contrast 0.04 produced full-wave data outside the entire nonnegative first-Born cone at both 1% and 2% noise. The exactly Born-generated readings on that same shape and the full-wave readings on an ultraweak ellipse both admit explicit Born fits. All six readings use the original 96², ten-illumination, twenty-receiver Green and port maps. The comparison contains two independent shapes, with two noise levels per arm. The exact Born arm reuses the two-disk shape and has its own prespecified noise draw. [Frozen design](FROZEN_DESIGN.md)

| Generator | True peak | Clean generator vs. same-material Born relative error | Proved cone lower distance / \(r\), 1% / 2% | Explicit Born fit upper residual / \(r\), 1% / 2% | Decision |
|---|---:|---:|---:|---:|---|
| Ultraweak ellipse, full VIE | \(9.992\times10^{-7}\) | \(8.862\times10^{-7}\) | 0 / 0 | 0.999879 / 0.999924 | Explicit Born fit in both balls |
| Modest paired disks, full VIE | 0.040 | 0.060105 | 5.348457 / 2.748230 | 6.073441 / 3.154508 | Every nonnegative Born image excluded |
| Modest paired disks, exact Born | 0.040 | 0 by construction | 0 / 0 | 0.998794 / 0.998844 | Explicit Born fit in both balls |

Here \(r\) is the rounded-up realized noise norm plus \(10^{-6}\|Y\|_F\). The explicit fit scales the offline reference shape by one nonnegative real coefficient, chosen by closed-form least squares and checked by interval arithmetic. A fit in the ball proves compatibility. An upper residual outside the ball pertains only to its displayed shape and coefficient. For the modest full-wave disks the lower bound provides the exclusion. The direct \([0,4]\) search returned a zero lower bound on those readings. The cone result also excludes that box because its Born data are a subset of the cone.

For both VIE controls, the original material's full-wave prediction lies in both stated balls. The interval calculation bounded \(\|\operatorname{diag}(\chi)G\|_2\) by \(1.0551\times10^{-5}\) for the ellipse and 0.422341 for the disks, using the embedding Fourier symbol. It then propagated the checked current defect through the inverse and receiver map. This gave upper bounds \(3.49\times10^{-18}\) and \(1.305\times10^{-11}\) for the clean simulated data's errors, in absolute field units. The same disk material has a clean full-wave/Born field discrepancy of 6.01%. The controls show that the gate responds to a mismatch between the measured field, the data tolerance and the Born model. A small peak alone does not settle that mismatch. The results apply to these two materials and the stated acquisition.

## 中文摘要

超弱椭圆的原生全波观测在 1% 和 2% 的误差球内都有 Born 拟合。峰值仅 0.04 的双圆盘在同一原生离散算子上，干净全波与同材料 Born 场相差 6.01%。其带噪全波观测距整个非负 Born 数据锥的严格下界分别为半径的 5.348 和 2.748 倍。在同一双圆盘材料上改用 Born 生成为对照，两档读数亦都有可逐项核验的 Born 拟合。这个对照支持直接检验测量精度下的模型相容性。不能只按峰值判断 Born 是否够用，也没有在两件材料上估计总体通过率。

## Reproduction

From the project root, run the following with NumPy, SciPy, python-flint and threadpoolctl.

~~~sh
python -B research/z_alg_optical_inertia/physical_revision/born/controls/replay.py
~~~

The run reuses the unchanged Born test in [run.py](../run.py), the native FFT current replay in [analyze.py](../../modes/analyze.py), and that script's existing [curve operators](../../../submission/curve/native.py) and [96² geometry](../../../synthetic/public/fine_geometry.npz). The completed single-threaded run took 7.55 s wall / 7.28 s CPU including geometry, all generators, interval forward checks, twelve direct data bounds, explicit upper fits and array writes. An earlier run computed the bounds and arrays before hitting a JSON boolean serialization error. The completed run corrected serialization and records the implied box exclusion. [Results and phase timings](results.json) · [Saved arrays](arrays.npz)
