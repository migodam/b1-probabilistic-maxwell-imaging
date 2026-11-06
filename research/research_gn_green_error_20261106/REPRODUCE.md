# 复现

本轮是已公开GT下的离线误差诊断，复用原GN端点。没有重新优化、生成新GT或改变历史门槛。

从B1仓库研究目录运行，使用已有NumPy、SciPy、Matplotlib、threadpoolctl环境。原主机命令：

~~~sh
REPO="/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r2_publication/repository"
PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
cd "$REPO/research/research_gn_green_error_20261106"
"$PY" src/analyze.py
"$PY" src/summarize.py
"$PY" src/tsom_chart.py
"$PY" src/plot_results.py
~~~

输入直接使用同一仓库的research/research_z_alg_r4_diversity：

- predictions/endpoints.npz：导出的物理材料m，不是reciprocal变量。
- RESULTS_SUMMARY.json：按task_id加入shape、contrast、noise和旧PASS/FAIL。
- public_data/operators.npz：G/S/E0和几何，已有pulse积分系数。
- truth_after_freeze/<observation_id>.npz：读取m；alpha只是归一化形状。
- src/physics.py：保持原完整非线性物理。

src/analyze.py固定筛选CPU native_trf与reciprocal_trf，共144端点；保留全部12形状、3材料增量和两种噪声。重新运行analyze会支付180次完整LU，而只改图或查看阈值曲线时不应重跑它：

- 已保存error_arrays.npz、fixed_spectra.npz及cumulative_mode_energies.npz，可重画和改变统计口径。
- summarize.py只处理已存系数与表，无新物理。
- tsom_chart.py只从已存谱构造旧28维Two-fold正交chart，并分析已保存误差；0新LU，不运行优化。
- plot_results.py只读表和数组，无新物理。
- 主.1门与预声明.03/.3敏感性不变；额外阈值响应标作从保存系数得到的上下文，不是新拟合或GT调参。
- exact/numerical intersection、in_then_ob、ob_then_in、tsom_projected_span四类输出含义见MATHEMATICAL_DEFINITIONS.md。不要将后者的重叠百分比强行归一化成100%。

结果位于results/，主图位于figures/。单operator的output_fraction仅指同一个G作用后正交SVD分组的能量；twofold串联后的传播单元不保证正交，另有cross energy。

表格主聚合为sum(component squared error)/sum(total squared error)，另存逐案例百分比、普通范数比例、RMS、原数据误差、PASS/FAIL条件与real-material静态控制。不会把多噪声／多参数化当独立材料样本。

