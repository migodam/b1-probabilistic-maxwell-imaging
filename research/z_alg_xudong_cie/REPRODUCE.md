# 复算入口

论文和新的成像实验各有独立目录。Paper 1 的 LaTeX 包可独立编译。Track B 的 Python 实验沿用同一 B1 项目的历史算子、数据文件和研究目录。项目根执行以下命令，设置 NumPy/SciPy BLAS 为单线程。需要 NumPy、SciPy、python-flint、threadpoolctl、Matplotlib。在线反演入口只加载源、接收器、传播核、训练读数、噪声半径和固定参数。离线目录含真值与保留端口。分发包内的已保存场、材料、计时与证书允许直接查看原结果。

- [B1 原生 96² 峰值区间](minpeak_native/REPRODUCE_ZH.md) · [七件 MAX-STRUCT](minpeak_maxstruct/REPRODUCE.md) · [超阈值像素数](area/README.md)
- [四件开发病例原始电流轨迹](diagnosis/REPRODUCE.md) · [两个实验室的测量与校准](chambers/README.md) · [新合成数据与原生验收](holdout/README.md)
- [冻结的方法、参数与代价](benchmark/FREEZE.json) · [根据开发病例补充的双初值对照](benchmark/PORTFOLIO_FREEZE.json) · [有限支撑的物理推导](pruning/THEOREM.md)

## 新合成数据和两种实测

若使用随包提供的合成公开数据，直接从训练开始。重新生成时，在空的任务数据目录运行 generation.py；它将观测和离线参考分别存储。原 Fresnel 与 Manitoba 的实验室采集须使用各自来源的测量与校准文件。测试使用原有效通道；实物图上的几何和材料只作标称参考。

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
study=research/z_alg_xudong_cie
python "$study/holdout/generation.py" --freeze-record "$study/benchmark/FREEZE.json"
python "$study/benchmark/run.py" --workers 4
python "$study/benchmark/reuse.py"
python "$study/benchmark/run.py" --config "$study/benchmark/PORTFOLIO_FREEZE.json" --workers 2
```

每次候选使用原端口评分，只有在预算内完成的材料才能入选。合成 32² 候选按原像素边界复制为 96²，在原传播核上再解 Maxwell 方程。误差排名和实测留出均在选图封存后计算。若要重新计时，应使用新的输出目录，记录并行任务和线程设置。冻结数据与代码决定物理问题；按墙钟截断的迭代数会随硬件和并发负载改变。

```sh
python "$study/holdout/verify_truth.py" --freeze-record "$study/benchmark/FREEZE.json"
```

详细评分命令、逐种输出及可下载包见本轮交接。比较中在线搜索、正演选图、证书、失败尝试和复用几何的分摊时间分别列出。离线生成、证明生成材料在球内及揭盲评分另计。

完整的盲评分和表格在全部方法输出后运行。方法按物体分开读取训练数据和离线参考。直接阅读压缩包内的表图可以跳过重算。重新运行合成或实验室程序前，应将包放在新的目录，以保留原有付费回执。

```sh
python -B "$study/holdout/offline_score.py" --wait
python -B "$study/benchmark/compare.py"
python -B "$study/chambers/confirmation.py" --scores "$study/benchmark/chamber_scores.json" --plot "$study/benchmark/chamber_comparison"
```

定位试验中的 16 候选区全部使用同一批原始带噪观测。给出具体可行传递矩阵后，可用下列命令重放能量、互易性和材料之间的差距。它们读取实际保存的材料盒排除证明和传递见证。

```sh
python -B "$study/pruning/contraction_control.py"
python -B "$study/pruning/reciprocal_control.py"
python -B "$study/pruning/constitutive_margin.py"
```

[压缩包结构、依赖项及原始暗室文件入口](ARCHIVE_REPLAY.md) · [16 个候选位置的冻结实验](pruning/site_control/RESULTS_ZH.md)

追加四件的单变量比较、原生采集、同预算控制和赛后评分使用独立目录。[四件材料的复算顺序](mechanism/REPORT_ZH.md)

32² 表示误差与原生数据可行性的对照在全部在线结果封存后运行。[最优块平均、固定 0001 的可行见证及区间正演](holdout/projection/README.md)
