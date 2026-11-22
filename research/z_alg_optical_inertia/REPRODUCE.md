# 复现指南

从实验目录运行下列命令。使用同一套装有 NumPy、SciPy、threadpoolctl、matplotlib 和 python-flint 的 Python。保存的复数数组和 FFT 乘子按原二进制值读入。归档内保留了逐项验证所用的原始输入。

## 96² 细网主结果

交接 ZIP 中的合成观测、细网算子、离线合成材料和冻结候选足以复算本节。保留 `synthetic/public/`、`synthetic/evaluation/` 和 `results/` 的相对路径。复算前可另存原结果目录，因为脚本会写入同名回执。

```sh
python -B fine/verify_fourier.py
python -B fine/long_ratio.py
python -B fine/all_scalar.py
python -B fine/spatial_followup.py --verify
python -B synthetic/nonvacuity.py
python -B fine/plot.py
```

- `verify_fourier.py` 从 `results/theory/fine_scan.json` 选取已存候选，用 `synthetic/public/fine_geometry.npz` 和公开 Y 证明标量多源下界。物理变量仍是 9216 格。若要重跑标量数值搜索，先执行 `python -B fine/fourier_bound.py`。
- `long_ratio.py` 严格验证 96² 接收核内的暗电流；随后 `all_scalar.py` 构造大正虚部乘子下的拟合反例。两者合并覆盖所列全环案例中任意全局实标量虚部乘子、正实部乘子和单位范数单照明组合。
- `spatial_followup.py --verify` 使用 `results/theory/fine_spatial_followup.json` 封存的分块提议，复算第 03 件半环 10% 的细网严格排除。若要先重做两件分块提议的细网数值预筛，运行 `python -B fine/spatial_followup.py`，它读取 `results/core/gradient_certificates.json` 中的粗网提议。
- `nonvacuity.py` 离线读取 `synthetic/evaluation/optfresh_01.npz` 到 `optfresh_06.npz` 中的真材料。用原 96² 算子验证真材料进入膨胀球，并重新检查标量证书余量。细网空间证书自己核对同一生成误差补量。主图以 `fine_feasibility.json` 的 `inflated_radius` 为分母，重写 `figures/fine_peaks.png` 和 `figures/fine_peaks.pdf`。 [真材料可行性与误差账本](synthetic/FINE_FEASIBILITY_ZH.md)

## 32² 对比、局域试验与历史输入

- 新合成 32² 数值扫描可用 `python -B core/scan.py --cohort synthetic --output results/core/scalar_synthetic_replay.json`，直接读取包内 `synthetic/public/manifest.json`、对应 Y、(S,B) 和 `coarse_green.npz`。全局扫描、固定惯性分块梯度和独立区间验证的完整顺序见[粗网复现说明](core/README_ZH.md)。
- 交接 ZIP 按原路径带有五份旧观测和上一轮 180 点，解压后可直接进行配对复算。历史输入位于同级 `research/z_alg_secant_tasks/highdim/`。需要 `replay.py`、`public_vie_inputs.npz`、`observation_radii.json` 及 `results/<观测 ID>.npz` 和 `.json`。逐点配对另读 `sector_inertia/results.json`。五个 ID 为 `r10_twod_0001`、`austria_eps4_400MHz`、`s2d_007__pw180_10x20`、`s2d_007__pw360_10x20`、`s2d_010__pw180_10x20`。
- 若要从种子重新生成六件合成材料的观测，`synthetic/generate.py` 调用历史 `research/srap_cie_array_study_20261121/data2d.py` 和 R8 `src/paper_physics.py`。交接 ZIP 的 `research/` 树含所需 SRAP 与 R8 源码，包含原项目的 `research/delegated/max_struct_r1/r8_worktree/` 路径。再生成不属于上面的冻结输入证书复算。 [合成数据及生成合同](synthetic/README.md)
- 区域扫描还使用 `research/limited_aperture_routed_dr_20261120/.runs/public/` 的三份旧观测。历史端点和离线材料的逐案出处在 `results/roi/endpoint_metrics.json` 及端点评估脚本中。 [区域审计和先前端点复现](roi/REPORT_ZH.md) · [双频复现](multifrequency/REPORT_ZH.md)

## Fresnel 实测

附件根目录的 `DEPENDENCIES.md` 列出可离线运行的各项。交接 ZIP 保留敏感性表、图、证书摘要与复算代码。原始 Fresnel 测量和逐端口审计文件不转存。完整重跑需要从原采集附档或历史镜像取回 `FoamDielIntTM.exp`、`FoamDielExtTM.exp` 和 `FoamTwinDielTM.exp`，放入项目的 `research/delegated/a2_literature/data/2005_iop_21_6_S09/`，对应的历史 R8 `src/fresnel_adapter.py` 已随 ZIP 提供。 [实测输入、掩码和复算命令](fresnel/REPORT_ZH.md)
