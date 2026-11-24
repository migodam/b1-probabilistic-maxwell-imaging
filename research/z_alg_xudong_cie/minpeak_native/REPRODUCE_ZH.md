# 原生 96² 最小峰值夹逼复算

本目录的物理问题是在原数据误差球内，寻找非负、无损、各照明共享的对比度图，并给 `max(m)` 一个有严格上界的可行见证。六个物体相互独立，1% 和 2% 是各物体的配对条件。只有 `submission/curve/observations/` 中当前行的 Y 进入在线候选生成。原 2% `n02` 与 `pair2pct` 逐元素相同。1% 采用同噪声方向的单独已存中心。两档都使用原登记的扩大 Frobenius 半径。 [输入生成与冻结](../../z_alg_optical_inertia/submission/curve/README_ZH.md)

## 从冻结成果复算

依赖 NumPy、SciPy、python-flint、threadpoolctl。生成对照图另需 matplotlib。下面的命令从本目录运行，第一步使用公开的已存区间证据。第二步只使用原始公开算子、当前观测和冻结的候选材料，重新计算原网格物理一致性。只有第三步读取离线真值。 [归档下界的逐行回执](lower_bounds.json)

```sh
python -B lower_bounds.py
python -B verify.py --inputs runs gn gn_short continuation --output verified_replay.json
python -B score_offline.py --inputs runs gn gn_short continuation \
  --verified verified_replay.json --output bracket_replay.json --figure bracket_replay.png
```

运行原生正演验收时，活跃格包含所有材料值严格大于零的单元，原方程使其他单元的电流严格为零，因此消元保持 9,216 维原方程。直接求解活跃状态和转置端口，再用通过 Neumann 检验的逆界和 128 位区间 FFT 围住状态与接收转置缺陷。令状态缺陷 d=Mb−Ax，端口转置缺陷 D=S−WA，严格包络使用 `Sx+Wd`，其误差不超过 `‖D‖F ‖A⁻¹‖₂ ‖d‖F`。用扩大半径和预测残差上端比较。浮点旧正演既不提供半径，也不提供严格可行性。 [独立验证源码](verify.py) · [继承的逆界构造](../../z_alg_optical_inertia/synthetic/nonvacuity.py)

## 在线搜索路径

第一条路径是从非负 Born 初值开始的原生伴随 L-BFGS-B。24×24 和 48×48 只约束材料变量在格块中相等。每个实际残差和伴随都用原生 96² Green FFT。01 的 2% 首次路径用三个 70 次预算阶段。

```sh
python -B search.py --object optfresh_01 --noise 2 --out runs_replay
```

第二条路径是原生 Gauss–Newton。固定两个中心高斯幅度搜索分支，DPCS/Dual-CIE 在部分病例作为额外的已有初始化。每步用 10 列电流状态和 20 行伴随接收敏感度得到复端口雅可比。带盒的局部模型通过 144 步投影线性迭代求解，只有原状态残差下降的线搜索步可接受。由当前材料排名和邻域扩展选择细网格支撑。限制支撑预算便于独立验证；下界仍来自全体材料。

```sh
python -B cie_seed.py --object optfresh_01 --wall 65 --out cie_replay
python -B gn_recovery.py --object optfresh_01 --noise 2 --cie cie_replay --cap-path wide --out gn_replay
python -B gn_recovery.py --object optfresh_01 --noise 1 --cie cie_replay --smooth 2.4 3.2 \
  --starts-iterations 22 --refinement-iterations 26 --peak-calls 58 --out gn_replay
python -B gn_recovery.py --object optfresh_02 --noise 2 --cap-path wide --out gn_replay
python -B gn_recovery.py --object optfresh_02 --noise 1 --smooth 2.4 3.55 \
  --starts-iterations 17 --refinement-iterations 16 --peak-calls 45 --out gn_replay
python -B gn_recovery.py --object optfresh_05 --noise 2 --cap-path wide --out gn_replay
python -B gn_recovery.py --object optfresh_05 --noise 1 --smooth 2.4 3.15 --cap-path wide \
  --starts-iterations 20 --refinement-iterations 24 --peak-calls 55 --out gn_replay
python -B gn_recovery.py --object optfresh_03 --object optfresh_04 --object optfresh_06 \
  --smooth 2.4 3.55 --starts-iterations 20 --refinement-iterations 24 --peak-calls 55 \
  --cap-path wide --out gn_replay
python -B gn_recovery.py --object optfresh_03 --object optfresh_04 --object optfresh_06 --noise 2 \
  --smooth 2.4 --starts-iterations 12 --refinement-iterations 9 --peaks 0 --out gn_short_replay
```

第三条路径在已可行的当前材料上构造原生端口雅可比，并对更低峰值盒求切向步。先尝试降峰 0.15，失败后减半。通过同一个原生误差球后才接收。线性切向重投影实际用真正演评价。

```sh
python -B continuation.py --object optfresh_01 --object optfresh_05 --noise 2 \
  --input gn_replay --output continuation_replay
```

针对 02 的一次技术对照用原生 Green 投影到 16×16 状态块，逐步提高物理残差权重。其投影场拟合不能计作原生可行性。每份投影候选在原生算子再次评分。 [投影对照源码](lowdim.py)

```sh
python -B lowdim.py --object optfresh_02 --noise 2 --out wri_replay
```

归档结果以各观测的 JSON 收据和实际数组为准。并发与机器性能改变墙钟/CPU 截止附近的迭代数。验证从冻结候选出发，可脱离搜索重现。离线真材料只提供另列的可行上界和结构分数。发现在线候选的峰值大于真峰值时，仍可以采用离线真材料作为当前更紧的数学上界。研究过程是否找到了温和替代材料以在线候选与离线真值的实际峰值对比判定。
