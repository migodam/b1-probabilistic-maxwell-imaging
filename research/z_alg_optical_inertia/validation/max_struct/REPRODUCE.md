# 从公开复端口复算 R14 材料迁移

在本目录使用安装了 NumPy、SciPy、python-flint、threadpoolctl 和 matplotlib 的 Python。同一测试覆盖 7 件材料、9 个上限、每上限 11 个实部乘子和 13 个光学比值。[冻结方案](FREEZE.md) · [机器可读冻结](frozen.json) · [原观测与噪声合同](public_provenance.json)

```sh
python -B scan.py
python -B verify.py
python -B summarize.py
```

- `public_inputs.npz` 保存一份共用的公开 G、B、S 和全部七份原样复读数、几何与噪声参数。R14 原研究的每件 NPZ 重复保存几何，本包只去除了这种重复。 `scan.py --export-from PATH` 可从 R14 的 `data_public/independent24/` 重做这份公开打包。
- `scan.py` 对一组权重只分解一次 K，共享解七份数据的最多 84 个右端，再分别形成七份源矩阵。它仅生成数值提案，留下全部通过、未通过及耗尽照明方向的条件。
- `verify.py` 用冻结的前二规则选提案。每组入选权重只做一次 576 阶严格惯性检查，复用七件的任何已选读数。核查由 128 bit Arb 从公开数组重建 K、验证位移 Cholesky 缺陷及负子空间、独立界定 Schur 残差源矩阵谱。最后输出的是证书公式的上下界。上界不控制未知的物理最优残差。
- `offline_scalar_scores.json` 是本轮在授权历史工作区从原封存材料得到的逐例标量。公开复算线使用这些记录做图和覆盖汇总。若在有权读取原材料的工作区重评分，先完成 `verify.py`，再运行 `offline.py --reference PATH`，最后运行 `summarize.py`。离线步骤只写评分标量。公开包里的新数组均属于原合法输入。
- 完整网格和验证回执都保留。若要保留原机费用，在单独复本内复跑。若搜索从未得到某类正规候选，配对表的区间列留空；零证书值与未验证不能互换。

[原 MAX-STRUCT R14 生成与计分合同](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_max_struct_r14/BENCHMARKS.md) · [光学惯性论文](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_optical_inertia/submission/)
