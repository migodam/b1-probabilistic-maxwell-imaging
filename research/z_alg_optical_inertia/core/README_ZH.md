# 全局峰值与分块光学乘子的复现实验

实验检验在冻结的端口算子和联合 Frobenius 观测球下，虚部等式是否减少惯性负方向，是否提高可排除的全局对比度上限。输入为 5 份已暴露存档和 6 件新材料的 24 份配对观测。每件新材料有全环、半环以及 2% 和 10% 两种带噪视图。数值提议只接收公开算子和复场。选定乘子后的代理误差比较会读取离线误差账本。

## 算子和评估量

令 `U`、`Dmu`、`Dlambda` 分别为逐格上限、非负实部权、实数虚部等式权。代码用 `Z=(Dmu U+i Dlambda)/2`，计算 `K=S* S+Dmu−ZG−G*Z*` 和 `L=S*Y+ZB`。`G`、`S`、`B`、`Y` 均取原始二进制复数组。标量的 `lambda` 正号增加 `Im_H G`，`lambda_over_mu` 记为 `t`。原存档和新合成均匹配 `H1` 圆盘脉冲的虚部 `πka J1(ka) J0(kr)/2`。R10 的脉冲半径为 `.018`，400 MHz 网格使用 `h/√π`。阵列共享均匀脉冲尺寸，未额外乘网格面积。`dark_thresholds.json` 给出 `H1`、共轭 `H2` 的实际算子误差、暗子空间临界上限和接收奇异值。精确离散公式的虚部在浮点存档中有小负特征值。每个被采纳的严格证书独立处理 K 的符号。

统一尺度为残差 `||SJ−Y||_F`。新材料的严格结果以同一个 32² 有限模型和选定球为条件。噪声加本案代理差距使用预存的 32² 面积平均材料代理误差。没有全盒的统一跨网格误差上界。基于同一套公开端口的暗向量和拟合电流可证明部分单源松弛的不可排除性，单源条件包括任意单位范数的照明组合。`single_source_obstruction.npz` 保存核投影前的提议向量与二进制拟合电流，配套 JSON 记录区间余量。

## 数值搜索与严格检查

- `scan.py` 使用 22 个正、负、零 `t`，6 个 `μ`。先分解与接收无关的 `A`，再通过接收空间的低秩更新数惯性及求源空间余量。空间臂先测四种几何分割方向，再在最高数值通过点的下一档执行固定 4×4 分块坐标扫描。分块按物理位置从低 `y` 到高 `y`、每行从低 `x` 到高 `x`。
- `refine.py` 对半环 `.5`，以及全环、半环 `.75` 和 `1` 加密成 43×31 的标量对照。优化对象是同一公开复场和观测半径。图中的分块族曲线在已有分块搜索点纳入较强的标量解。
- `gradient_blocks.py` 从较强标量提议出发，在三个匹配步长上比较固定惯性的分块光学梯度与标量梯度。`gradient_controls.py` 复核数值新增排除，并进一步加密局部标量对照。
- `certificate.py` 用向外舍入的 Arb 重新组装冻结矩阵。由移位 Cholesky 反算缺陷和负子空间证实 K 的惯性及谱隙，再由源空间 Rayleigh 界证明截去最大 `ν` 个特征值后的正谱和。只有区间下界严格大于相应半径时才登记排除。`certificates*.json` 保留实际墙钟、CPU、特征分解与严格乘积次数。
- `dark_witness.py` 用高精度接收 Gram 求解的区间包络把提议投影到精确 `ker(S)`。`single_source.py` 同时对大光学比率构造公开数据拟合电流。

## 复现入口

从项目根目录使用具有 NumPy、SciPy、threadpoolctl、matplotlib 和 python-flint 的 Python 环境。历史存档在同一仓库的 `z_alg_secant_tasks/highdim/`。令 `D=research/z_alg_optical_inertia`，依次运行以下命令。

```sh
D=research/z_alg_optical_inertia
python -B "$D/core/scan.py" --quick --output "$D/results/core/quick.json"
python -B "$D/core/scan.py" --output "$D/results/core/scalar.json"
python -B "$D/core/scan.py" --spatial --cohort archive --scalar-input "$D/results/core/scalar.json" --output "$D/results/core/spatial_archive.json"
python -B "$D/core/scan.py" --spatial --cohort synthetic --scalar-input "$D/results/core/scalar.json" --output "$D/results/core/spatial_synthetic.json"
python -B "$D/core/dark.py"
python -B "$D/core/refine.py"
python -B "$D/core/refine.py" --caps 0.75 1 --array all --output "$D/results/core/focused_upper.json"
python -B "$D/core/verify.py" --cohort archive --output "$D/results/core/certificates_archive_scalar.json"
python -B "$D/core/verify.py" --spatial "$D/results/core/spatial_archive.json" --spatial "$D/results/core/spatial_synthetic.json" --focused "$D/results/core/focused_scalar.json" --prior "$D/results/core/certificates_archive_scalar.json"
python -B "$D/core/verify.py" --case austria_eps4_400MHz --case optfresh_01__pw360_10x20__n02 --case optfresh_05__pw360_10x20__n02 --focused "$D/results/core/focused_scalar.json" --focused "$D/results/core/focused_upper.json" --prior "$D/results/core/certificates.json" --output "$D/results/core/certificates_focused_upper.json"
python -B "$D/core/dark_witness.py"
python -B "$D/core/single_source.py"
python -B "$D/core/single_source.py" --case optfresh_04__pw360_10x20__n02 --output "$D/results/core/single_source_obstruction_block_case.json"
python -B "$D/core/gradient_blocks.py"
python -B "$D/core/gradient_controls.py" refine
python -B "$D/core/gradient_controls.py" polish
python -B "$D/core/gradient_controls.py" verify
python -B "$D/core/summarize.py"
```

[严格排除及费用索引](../results/core/summary.json) · [全局标量网格](../results/core/scalar.json) · [固定惯性梯度](../results/core/gradient_blocks.json) · [精确暗向量和大乘子阻断](../results/core/single_source_obstruction.json)

## 上轮 98 个源方向耗尽点

运行 `python -B "$D/core/legacy98.py"`，可把上一轮 180 个参数点按原观测、对比度上限和实部乘子逐一对应。零虚部乘子的惯性在 180 点全部一致。在原来 98 个 `ν≥s` 的点中，扫描 11 档正虚部乘子后有 8 点恢复到 `ν<s`，没有点得到正残差下界或超过相应旧噪声半径。Austria 和三个 SRAP 视图各恢复两点，原 R10 网格没有源方向耗尽点。三个 SRAP 视图恢复的点甚至可以达到 `ν=0`，源空间剩余正谱和仍为零。[逐点记录](../results/core/legacy98.json)
