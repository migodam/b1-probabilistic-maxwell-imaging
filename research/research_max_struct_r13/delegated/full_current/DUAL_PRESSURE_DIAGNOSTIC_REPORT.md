# 近似 AP 对的定向接收伴随修正诊断

范围：仅四个已经暴露的 public case；每例比较 AP128 与 DR128 后单次 A-shadow。没有读取 fresh/GT，没有远程操作，没有修改在线 C2 文件。数学证明与最终解释属于 parent；实现仅测量 `THEORY_DERIVATIONS.md` §4.1 的已指定代数量。

实现：`src/dual_pressure_diagnostic.py`。结果：`dual_pressure_diagnostic_v2/{case}/DUAL_PRESSURE.json` 与两份 `{route}_FIELDS.npz`。

## 计算约定

每例固定同一个原 `DenseFieldCurrentProjectors`，原 G/E/S/Y、原缩放 s=6、原 [.02,20] 材料盒、原 raw CN99 球，`grouped_ports=True`，W=I。两条路线都从同一个已付费物理初值出发；路线之间不传递辅助状态。

对每个所得 a∈A，计算 `c=PC(a)`、`d=a-c`、`a'=PA(c)`、`e=a'-a`。实际 λ 立即取自此次 `PA(c)` 的 `last_A_info['multiplier']`，同时读取 `pseudoinverse_limit`。源码的 T=sS、`1+mu*T H^-1 T*` 衰减对应半平方球约束的 λ=mu，无 factor-two 转换，没有猜测其它 key。

在 m_c 处调用原 `model.state(m_c)` 得到并保留一份 Maxwell LU，然后把 `S* r'` 与 `eX+sG*eQ` 合并为 2ns 个 RHS，调用一次同一 LU 的共轭转置三角求解：

`p=A_m^{-H} S* r'`，

`correction=-A_m^{-H}(eX+sG*eQ)`。

记录的是这个实际方向上的修正及放大，不是 resolvent 算子范数估计、全局误差上界、未来轨迹界或恢复证书。用于定位事件的 t_geom 始终从实际 c/d 和合法材料盒计算；只有 first event cell 是内点时才附带未修正的 pressure-only 比值，且明确标为 exact-pair 公式的反事实参考。

## 主要观测

| case / 路线 | ||e||/||d|| | ||correction||/||dX|| | ||correction||/||e|| | t_geom | first event cell 的 p-energy 降序名次 |
|---|---:|---:|---:|---:|---:|
| A r10_twod_0004 / AP128 | 0.00468 | 0.01082 | 2.287 | 5.691 | 2 |
| A r10_twod_0004 / DR→A | 0.81243 | 1.11346 | 1.287 | 1.183 | 10 |
| B r12_confirm24_0010 / AP128 | 0.01853 | 0.18031 | 8.439 | 7.807 | 19 |
| B r12_confirm24_0010 / DR→A | 0.60633 | 1.00548 | 1.407 | 1.232 | 17 |
| C r8_test_0001 / AP128 | 0.06170 | 0.64473 | 7.726 | 7.384 | 1 |
| C r8_test_0001 / DR→A | 0.22523 | 1.11984 | 4.492 | 63.051 | 435 |
| A r10_twod_0001 / AP128 | 0.00958 | 0.41735 | 32.402 | 32.158 | 283 |
| A r10_twod_0001 / DR→A | 0.19704 | 1.18303 | 4.265 | 93.192 | 526 |

A r10_twod_0001 的 AP scalar pair defect 约0.96%，但实际定向修正是 dX 的41.7%，其 `||correction||/||e||≈32.4`。这就是仅报告 scalar pair defect 会遗漏的可观测量。四个 DR→A 状态的定向修正均约等于或大于 dX；不能据此直接套用无修正的精确 AP 陷阱解释。

p-energy 最大位置与 first event cell 也并非普遍一致。AP 的 first cell 除 C 外均为材料边界事件；DR→A 的 C first cell 也是边界事件，能量排名本身不包含 Q_c 能量、边界相关项与定向修正。这些观测没有测量事件是否改善重建。

## 适用条件与数值限定

八个状态的实际乘子均有限，未进入 pseudoinverse limit，projector 的可达残差底低于正半径 δ，记录的 Slater 检查通过。JSON 同时保留原 `last_A_info`、ball residual/δ、互补量、材料幅值法向缺陷、`eQ-sG eX` 相容性和原物理状态残差。

作为适用性诊断，current KKT 相对残差约 `8e-8–7.2e-7`，重构 `dX+λs p` 与 correction 的相对差约 `1.4e-7–1.5e-6`。这些是实际投影/线性求解的数值残差，不能写成 machine-epsilon 精确性宣称，也不是本实验标题结果。未计算全局条件数、极点排除证书或 resolvent 范数上界；LU 成功及有限定向解的范围已明确标注。

## 成本与失败记录

每例共用一次 projector setup：1个 H Cholesky、1个端口 Gram eigensolve。每条路线128步；DR 路线另形成一次 A-shadow；两份诊断各增加一次 PC/PA。因此每例累计259次 A 与259次 C 投影。

每例3个原物理 LU（初值+两个 m_c）；每个诊断把两块 adjoint RHS 合并一次求解，无额外 LU。A/B/A1 每例 source36 + adjoint48 =84 complex RHS；C source18 + adjoint24 =42。成功批次合计12 LU、294 complex RHS。诊断 adjoint RHS 已包含在物理计数里，不能再相加。

每份 JSON 分开保留 projector 原始计数、诊断新增 G/S forward/adjoint 次数/RHS、物理模型计数，以及合并的 `total_dense_operator_cost`。合并表还显式补入原模型每次 FULL_FORWARD 中的 G@J/S@J、projector setup 中的 B*B 和端口 Gram；不把共轭转置矩阵的 materialization 计作一次算子应用。H/Maxwell LU/RHS 与稠密 G/S 乘法是不同成本类别。

成功批次逐例 CPU 约0.564、0.554、0.529、0.571秒，总计约2.22秒，BLAS1、串行；包含每例输入/模型/projector setup、两条路线、定向求解与 NPZ 序列化，排除模块导入和最终 JSON 序列化。

第一次 A r10_twod_0004 尝试在两条路线与两份 NPZ 已完成后，因 NumPy bool 的 JSON 序列化失败。修复仅为通用 NumPy scalar encoder；失败字段保留在 `dual_pressure_diagnostic_v1/`。该次3 LU/84 RHS 的已完成调用另记于 `FAILED_ATTEMPT_ACCOUNTING.json`，CPU 未在失败前捕获，明确为 unknown，不从 wall time 猜测。随后 A 重跑，B/C/A1 各首次在成功批次运行。故包括失败的已知物理总量是15 LU/378 RHS，而 **2.22秒不是含失败的完整 CPU 总计**。

未运行完整反演、恢复评分、文献新检索或机器精度 identity campaign。所有字段保存 c/d/e/p/correction、实际材料、两个 RHS、几何阈值、p/correction 分布及实际 PA 乘子，供 parent 独立审查。
