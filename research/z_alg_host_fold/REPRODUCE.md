# 复算这一个已知宿主实验

英文论文的定理、读数和时钟记录对应 2026-11-21 研究包中的 36 圆盘仿真。打开本目录的 `index.html` 阅读中文摘要，`manuscript.pdf` 阅读英文论证，`PRECEDENTS.md` 查对最接近的论文。模拟使用 30 个已知介质单元、6 个可独立变化的材料、3 个外部线源和 1 个复接收器。

## 使用发布时的两份观测

在解压目录（含 `src`、`results` 和 `figures`）打开终端。实验时使用 Apple M4、Python 3.12.13、NumPy 2.5.5、SciPy 1.18.1、python-flint 0.9.0、Matplotlib 3.11.2，以及支持 C++20 与 128 位整数的 Apple clang++。小矩阵运算在 CPU 上执行。

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install numpy==2.5.5 scipy==1.18.1 python-flint==0.9.0 matplotlib==3.11.2
python src/certify_fold.py
python src/point_fibers.py
python src/noisy_cover.py
python src/noisy_diagnostics.py
python src/normal_section.py
python src/host_noise_margin.py
python src/host_mechanism.py
python src/branch_figure.py
python src/noisy_plot.py
```

上述顺序沿用已存的 `results/host_operators.npz`，其中包含两份精确指定的总读数、噪声中心和浮点候选。正演系数的证明从源代码内的十进制几何、出射 Hankel 函数和完整宿主线性系统重新开始；新生成的复球系数写入 `host_operators_interval.json`。输出写回当前解压目录。要重建搜索和模拟原始读数，可在另一份解压目录先运行 `python src/host_model.py`。不同平台的 SVD 可以改变法向量的尾数，重新生成的读数定义的是相邻的一条线。

- `point_fibers.json` 给出两个端点的全部复根、材料分类和条件商上界。每个合法根又经过完整 36 状态正演回代。
- `noisy_cover.json` 与 `noisy_cover_cells.npz` 记录全域和局部范围内的计算费用、未决小盒、内外包与相同体积对照。`noisy_diagnostics.json` 单独保存两条可行材料路径以及通过全部 15 项低维检查的错误盒。
- `normal_section.json` 给出严格凹曲线及排除其他材料的整线证书，`host_noise_margin.json` 给出 222 个全域小盒的完整正演余量。`branch_summary.json` 给出最坏误差下界和行列式、宿主功率恒等式的核查。
- `host_model.json` 保留实际搜索的坐标、已尝试的偏移和原正演对照。`host_mechanism.json` 分别测试线性材料、仅有槽内反馈及完整槽间传播。

## 独立同伦对照

运行过的版本是 Julia 1.11.7 与 HomotopyContinuation.jl 2.24.1。下面读取上一步生成的两个中点多项式系统，跟踪两组候选，然后把候选送回同一套物理系数与区间证明。

```sh
julia --project=. -e 'using Pkg; Pkg.add(PackageSpec(name="HomotopyContinuation", version="2.24.1"))'
julia --project=. --threads=1 src/point_homotopy.jl
python src/point_fibers.py --hc
```

`hc_times.tsv` 分开记录系统构造、该系统首次跟踪和热运行。`hc_roots_minus.txt`、`hc_roots_plus.txt` 是成熟同伦的实际输出。Julia 进程含两组系统的构造、首次求解与重复求解，共用时 21.15 秒。首次软件安装和编译依赖是环境设置，不计入单个目标的科学比较；Python 表中另列了盒筛选器一次性编译和两方法共用的物理准备。耗时随负载改变。

英文排版源为 `manuscript.tex`，在图已生成后可用 `tectonic manuscript.tex` 输出 PDF。
