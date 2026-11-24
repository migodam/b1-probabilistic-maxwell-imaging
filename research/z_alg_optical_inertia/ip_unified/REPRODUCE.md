# 计算复现与原始数组

将 GitHub Release 的 **Z_ALG_Unified_Evidence.zip** 解压到一个新目录，从解压根运行下列命令。可在不另下载 v3 数据的情况下复算本文新增的物理论证、读数、有损实验和历史矩松弛对照。路径与当前仓库相同。使用独立复本可保留归档内的原始用时和历史失败回执。

~~~sh
python3.13 -m venv .venv
. .venv/bin/activate
python -m pip install numpy scipy python-flint threadpoolctl matplotlib sympy cvxpy scs clarabel
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
study=research/z_alg_optical_inertia
n="$study/ip_unified/numerics"
~~~

## 先读结果，再重做最小的严格核验

统一文稿图可以从冻结数组和严格结果直接生成，解析反例与有损标量移位全族上界也可以独立重放。生成的图保存到同包的 ip_unified/figures。

~~~sh
python -B "$n/controls/exact_controls.py"
python -B "$n/controls/read_historical.py"
python -B "$n/native/required_currents.py"
python -B "$n/feasibility/hidden_current.py"
python -B "$n/loss/pure_shift.py"
python -B "$n/loss/pure_shift.py" --weak-loss
python -B "$n/loss/directional.py" finite
python -B "$n/plot_paper.py"
~~~

| 证据 | 原始数组和费用的位置 |
|---|---|
| 原生无损光学模型 | synthetic/public/fine_geometry.npz、physical_revision/array/observations/arrays.npz 和同目录 96 份阵列清单中的带噪读数。正文的原生全盒比较从六件物体各取 32 发 32 收、2% 一行。半径来自随包 manifest。 |
| 无损下界和局部补偿 | ip_unified/numerics/native 的原生负子空间、正空间解、中心电流和端口、投影阵列、全部成功与失败的 JSON，以及 COST_LEDGER.tsv。源矩阵证据、区间以及重放入口在 numerics/source。 |
| 全盒物理上端 | numerics/feasibility 的 search_runs、projected_runs、cie_runs 内含实际材料、原生电流、预测端口和费用。strict_endpoints 保存区间结果及单独低峰盘体观测。upper_endpoints.csv 和 JSON 汇合候选状态、实际用时及尚未判定的盒。 |
| 真实有损全波试验 | numerics/loss 的 finite_0_inputs.npz、finite_1_inputs.npz、native_inputs.npz 共存 39 个按材料配对的带噪复端口记录，并存材料及对偶候选。格点合同和种子见 loss/FREEZE.md。严格回执、搜索和区间费用与数组同目录；配对低峰见证的数组也随包提供。 |
| 历史四区矩与 \(q/Q\) 对照 | ip_submission_v2/energy_sdp/alt025.npz、submission/single/local_moment.npz、原两件 10×20 复读数和四区外向回执。原几何上的合法电流能量见证含能量脚本及保存的材料/相干输入。 |
| 历史小域环形对照 | validation/max_struct/public_inputs.npz 含共用精确表示算子与七件公开观测。offline_scalar_scores.json 保存原离线严格计分，包内使用该结果读取 r14_new_0018 的可行上端。 |
| 计算成本 | [统一费用索引](numerics/COSTS.md)、逐分支 JSON、COST_LEDGER.tsv 和数值运行日志。区间复算与浮点候选生成使用各自记录的计时边界。 |

六件原生物体的实部真值和受控冻结干净场原先已公开，保存在 synthetic/evaluation。可行性搜索只读已冻结的带噪读数和声明算子。在有损路径中，声明实部与损耗角用来生成新的前向数据；证书步骤使用读数、半径和材料合同。

## 分支复算

暗电流证明用冻结的 u100_best.npz、u100_best.json、u100_best_conditioned.npz、u100_best_strict.json 及原 32×32 观测。对每件材料，如果 \(u=1\) 盒内真物理电流能进入原观测球，未中心化的投影 \(\|\Pi_-(K)J\|_F\) 必须超过一个严格正下界。六个给定下界约为 7.65458 至 22.00851，均使用原坐标下的未加权 Frobenius 范数。该计算未判断是否存在这样的候选。[必要下界回执](numerics/native/required_currents.json)

物理对照重读已保存的原生电流和前向区间，使用同一个负空间投影。独立低峰盘体的真实响应进入它自己的冻结噪声球，且 \(\|\Pi_-(K)J\|_F\ge18.39955\)；另三件候选在材料盒内，但其读数残差严格在各自观测球外。[盘体与盒内候选](numerics/feasibility/hidden_current.md) · [区间回执](numerics/feasibility/hidden_current.json) · [联合幅度图](figures/current-modes.png)

### 原生无损证书与全盒材料上端

[无损补偿与完整原生复算命令](numerics/native/README.md#复算和成本)中的本地 Python 路径替换为上面环境的 python。完整复算使用 CVXPY 的 SCS 及 Clarabel。先读冻结结果时可直接完成源空间区间复核：

~~~sh
python -B "$n/source/profile.py" --case u100_best
python -B "$n/source/outward.py" --case u100_best --weight .001 --count 13 --objects 123456
python -B "$n/source/spectral_envelope.py"
python -B "$n/source/cap_resolution.py"
python -B "$n/source/cap_resolution.py" --transition
~~~

全盒原生上端可从已存候选重放，无须先再次优化。生成新候选的三个独立入口及预算见 [全盒上端说明](numerics/feasibility/README.md#复现)。按原试验路径，直接多尺度搜索和 CIE 初始化用于第 01、05 件，投影建议用于全部六件。投影或 CIE 辅助优化的残差仍须经原 96² 完整正演重验。

~~~sh
python -B "$n/feasibility/enclose_native.py"
python -B "$n/feasibility/summarize.py"
~~~

### 有损

[有损报告](numerics/loss/REPORT_ZH.md#复现)列出全部命令。按该处命令从 loss 目录执行，或从解压根将各脚本路径写成 "$n/loss/脚本名"。run.py finite 和 run.py native 分别重建独立有限网格以及沿用公开数组的原生复波形。verify.py 从存档 Y、宣告的合同与前向模型复算单侧界。native_states.py --replay --focus 重放正文用的两条原生候选。directional.py all 比较固定源基在全部相关有损条件下的开销和距离。

### 四区矩与历史对照

~~~sh
python -B "$study/ip_submission_v2/energy_sdp/verify.py" \
  --proposals "$study/ip_submission_v2/energy_sdp/alt025.npz" \
  --output "$n/controls/historical_moment_replay.json"
python -B "$study/ip_submission_v2/energy_bound/source/verify_physical_current.py" \
  --proposal "$study/ip_submission_v2/energy_bound/receipts/pixel_profile.npz" \
  --proposal "$study/ip_submission_v2/energy_bound/receipts/coarse_12x12.npz" \
  --output "$study/ip_submission_v2/energy_bound/receipts/physical_current_replay.json"
~~~

这些命令重读固定旧矩和物理候选。历史暗室读数按原发布页提供的采集方入口获取。小域环形 R14 的原封存材料仍由其原数据合同管理，公开包提供完整的测量与在线算子、历史离线区间标量及评分入口。

## 历史出处

本包复制了复算所需的旧版公开代码及数值路径。需要浏览历史独立报告或复查更早的实测、旧权重和比较族时，使用以下原发布项。

- v3 的原阵列、物理反演、Born 对照和早期校准说明。[v3 公共证据包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/optical-inertia-ip-v3/physical_revision_evidence.zip) · [v3 暗室取数及校准](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_optical_inertia/physical_revision/reproduce.html)
- 四区矩松弛、相干能量与原联合对照。[历史矩松弛证据包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/optical-inertia-ip-v2/Energy_SDP_Evidence.zip)
- 七件 24² 环形 R14 迁移数据，以及更早的公开验证汇总。[历史验证包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-validation-v1/Z_ALG_Optical_Inertia_Validation.zip)
- 原生有界全波搜索和 CIE 初始化器的历史方法。[原 CIE / 原生搜索证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-xudong-cie-20261124/Z_ALG_Xudong_CIE_Research.zip)

可从完整研究工作区重建证据 ZIP。打包入口使用公开历史白名单，保留同一目录结构。数组、附图和数值日志对应存档运行，计算机和软件不同可能改变数值提案或计时；严格核验从所选数组重新建立单侧端点。

~~~sh
python -B "$n/package_evidence.py" --output Z_ALG_Unified_Evidence.zip
~~~

按审读意见追加的负空间检验仍读取相同 96² 原始算子与误差球。下列程序分别保留弱正模、对所有可能的稀疏逐格权重计算上端；完整命令和逐项费用见 [正模态与稀疏行](numerics/native/deflated.md)。实际权重变化的外圈/内部对照见 source 中的 edge 文件。

原生变权对照通过端口/核反射误差与正负惯性上下端核算。紧凑有损源分组分别给出接收椭球上下端与匹配状态权重下的界。

~~~sh
python -B "$n/source/edge_core.py"
python -B "$n/source/edge_core.py" --interval
python -B "$n/source/edge_core.py" --interval --wide
python -B "$n/source/edge_symmetry.py" --source
python -B "$n/loss/source_partition.py"
~~~

匿名论文源文件与数值数据分别下载。将 IP_Anonymous_LaTeX.zip 解压后在该目录运行 `tectonic manuscript_ip.tex`，或用 XeLaTeX 编译两次；包内有文稿所需的全部矢量图、期刊类和参考文献。
