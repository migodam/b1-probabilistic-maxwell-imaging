# 复现入口

从本研究包根目录运行以下命令。首次查看可从各节链接直接读逐行结果和失败字段。全量重算建议在研究包副本中进行，因为部分入口会写回各自默认结果文件。历史大图的重放和候选重算命令另有输出目录。

原区间实验用 Python 3.13、NumPy 2.5.3、SciPy 1.18.1、python-flint 0.9.0 和 threadpoolctl 3.7.0。NumPy、SciPy 负责双精度全波求解和候选搜索，python-flint 的 Arb 球算术用于把舍入误差包入实数或复数区间。Matplotlib 只用于画图。各实验目录的回执有本机实际版本、计时范围及物理调用计数。公开重放可使用 Python 3.10 或更高版本和兼容的 NumPy、SciPy。用新环境执行时，可按下列命令安装依赖。

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install numpy scipy python-flint threadpoolctl matplotlib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

## 历史带噪大图

第一条命令只用包内冻结的物理矩阵、候选和带噪读数，重放五份观测的 45 个材料对和两个等预算对照的 28 个候选。生成的 JSON 给出两端距原观测和参考场的实际距离、端点场差、材料任务差及原球和辅助球半径。随后两条从保存的历史中心重新生成弱方向、全波端点和同预算 GN。无需重新生成历史带噪场。

```sh
python highdim/replay.py --output replay.json
HIGHDIM_USE_FROZEN=1 python highdim/run.py --output recomputed
HIGHDIM_USE_FROZEN=1 python highdim/compare.py --output recomputed
```

复核从两个 CSV 开始。`PAIRS.csv` 有中心、方向、幅度和各球余量。`BASELINES.csv` 有保存的历史端点、原观测残差和离线参考差。对照 GN 的“旧球内跨度”按各条历史半径重新计算。缺少球内中心、无对称材料余量或更大步幅移出数据球，均保留在逐案例的 JSON。原参考像素评分需历史工作树授权和材料原件，使用 `highdim/score.py` 重做。公开包提供原有离线评分。`HIGHDIM_WORKSPACE`、R10/R8/SRAP 历史发布入口和每个原件见复现说明。 [高维复现说明](highdim/REPLAY_ZH.md) · [逐候选 CSV](highdim/results/PAIRS.csv) · [历史基线 CSV](highdim/results/BASELINES.csv) · [GN 费用与跨度](highdim/results/GN_BUDGET_COMPARISON.json)

## 已知宿主、同读数材料对及增测

按序生成六个设计中心的完整状态和端口算子、搜索有限幅度材料对并用区间法验证、在增测端口重算原分支。`--enhanced` 在 16 槽中心重新按各臂自己的物理 Jacobian 寻找见证。第一次和第二次 `finite_chord.py` 的输出各有独立文件名。保留默认输出位置时后续程序能找到冻结输入。`plot_material_pairs.py` 依赖前面已生成的数组，重画材料分布。

```sh
python array/array_secant.py
python cert/finite_chord.py
python cert/finite_chord.py --enhanced
python cert/verify_augmented_ports.py
python cert/plot_material_pairs.py
```

从 [采集合同](array/DESIGN_ZH.md) 核对同能增测的源和端口。方向、幅度、端口差和公共场漂移在 [方向扫描 CSV](array/results/pairs.csv)。双侧和单侧搜索的成功、合法盒回溯失败及区间认证边界在 [材料对说明](cert/README.md) 和 [基线逐行回执](cert/results/finite_chord.json)。[增测本身的反例](cert/results/enhanced_counterexample.json) 与 [增测审计的旧分支](cert/results/certified_pair_augmented_ports.json) 回答不同的可见性问题。

## 72 个规格盒与仪器增益

使用 `array/results/physical_arrays.npz` 冻结的二进制算子。全盒界先全局筛查，面见证则把左区固定到任务两极，逐端用完整复场球认证。极面初筛和细分搜索可排除正侧极薄层。增益压力测试再比较公共或各自的复增益能改变多少。

```sh
python bounds/box_modulus.py
python faces/face_modulus.py
python bounds/task_caps.py
python bounds/refine_caps.py
python bounds/gain_pressure.py
```

面见证说明如何精确重组 ROI。逐盒结果还保留三组失败、每次补充搜索和自由选择读数中心的检验，后者未新增合格端点对。细分预算另列未决盒、严格上界和可行下界。增益压力回执区分各自增益与公共增益。 [规格盒面见证](faces/README.md) · [逐盒场距和余量](faces/results.json) · [规格盒界说明](bounds/REFINED_CAPS_ZH.md) · [细分预算回执](bounds/refined_caps.json) · [增益压力回执](bounds/gain_pressure.json)

## 分离组件、全端口和连续圆柱

前两条重算非对称六槽组件的主间距阶梯和 40,000λ、80,000λ 两组件补充阶梯。第三条将数值候选包进区间并验证根存在，另核对跨组端口差。30λ 两槽脚本有自己的种子、全耦合求解和严格证书。`full_port_radius.py` 读取两个区间证书，不重新求根，把全部端口下的最大与最小材料对场距、共同中心半径界重新计算。`cylinders.py` 独立使用连续圆柱解析展开和双精度耦合解，输出各角阶的局部与跨端口读数、材料差、边界连续性及成本。

```sh
python multiplicity/separated.py
python multiplicity/separated.py --sizes 2 --spacings 40000 80000 --output multiplicity/results/intermediate
python multiplicity/certify_separated.py
python multiplicity/conditioning/conditioning.py
python multiplicity/full_ports/full_port_radius.py
python multiplicity/continuum/cylinders.py
```

包内的 [历史六槽小输入](multiplicity/historical_inputs.json) 已含组件材料、几何、共同旧读数及初值。若需要从原档重新提取，`separated.py --historical <目录>` 可在 `--inputs` 的目标文件缺失时读取历史 `host_operators.npz` 和 `point_fibers.json`。同一原档可额外执行 `python history/check_pair.py <目录>` 复核原两根及奇偶公式。历史出处和报告见 [分离组件说明](multiplicity/README.md)。

- [主间距阶梯](multiplicity/results/numeric_ladder.json) 与 [补充阶梯](multiplicity/results/intermediate/numeric_ladder.json) 保留各距离的 Newton 停止和未收敛种子。[根盒证书](multiplicity/results/krawczyk.json) 逐一给出根和跨端口差的区间。
- [30λ 条件复核](multiplicity/conditioning/README.md) 的 [逐间距结果](multiplicity/conditioning/results.json) 包含近间距的失败。[全端口说明](multiplicity/full_ports/README.md) 与 [共同球上下界和所有逐对场距](multiplicity/full_ports/results.json) 明示参考范数和共同观测中心。
- [连续圆柱说明](multiplicity/continuum/README.md) 与 [角阶、逐种子和边界回执](multiplicity/continuum/results.json) 给出 M=0 至 16 的谱收敛数值证据，以及离散脉冲近似和连续模型各自的读数代价。连续展开目前无严格无限阶误差尾界。

## 连续界面、端口标定与镜像停测

单独复现固定孔径下的连续圆柱比较，以及用存档复数决定何时停止发射，可从公开包根目录执行下列命令。前两条重解连续界面并重画方位图；余下命令用归档的源、接收器、实际复数和联合误差上限重放 13 种布局、229 个视图。后者不需要历史材料标签。

```sh
python multiplicity/aperture/aperture.py
python multiplicity/aperture/plot_aperture.py
python highdim/reflection/replay_symmetry.py --stage frozen
python highdim/reflection/exposure_prefix.py
python highdim/reflection/port_gains.py
python highdim/reflection/plot_replay.py
```

其中高维试验对照同一顺序下读整批、固定半批、事后最有利的每布局固定前缀以及逐次读到足够即停。10 份实际噪声超额的视图和所有不具备阵列镜像的布局都保留。[几何与发射组的定义](highdim/reflection/README_ZH.md) · [连续圆柱的参数和失败记录](multiplicity/aperture/README.md)

## 材料峰值与接收器选择

单独复算论文中最有区分力的一项严格结果，只需对预先保存的环形观测重新认证。输入为归档复场、独立的误差半径表、Green 及收发算子，乘子直接指定。该入口已经独立重放，所有证书数值与原记录相同。

```sh
python highdim/sector_inertia/certify.py --case r10_twod_0001 --cap 5 --mu 0.001 --output replay_u5.json
```

从头执行同算子下的基线、原定五份观测、追加的一个环形材料、隐藏电流对照和新增接收位置的比较，可用下列命令。初始几何单独保存后才读取新位置的数据。五份原观测对应四份不同材料。[正定比较和依赖](highdim/sector/README_ZH.md) · [不定矩阵扫描及严格认证](highdim/sector_inertia/README_ZH.md) · [配对接收设计](highdim/sector_geometry/DESIGN_ZH.md)

```sh
python highdim/sector/sector.py
python highdim/sector_inertia/scan.py
python highdim/sector_inertia/scan.py --case r10_twod_0001 --caps 4.25 4.5 4.75 5 5.5 6 --output highdim/sector_inertia/ring_cap_extension.json
python highdim/sector_geometry/geometry.py
python highdim/sector_geometry/design.py
python figures/plot_peak_constraints.py
```

## 重画综合图

```sh
python figures/plot_research.py
```

程序仅从上述机器结果重画。输出包括 [分离组件 PNG](figures/multiplicity_scaling.png)、[分离组件 PDF](figures/multiplicity_scaling.pdf)、[规格盒 PNG](figures/specification_range.png) 和 [规格盒 PDF](figures/specification_range.pdf)。[实验结果](EXPERIMENTAL_RESULTS.md) 对照汇总各自的条件、实际失败与费用。
