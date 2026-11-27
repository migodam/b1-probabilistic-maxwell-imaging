# 复现、输入和逐步运行

[源码与原始指标](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/hc_cie_r2) · [历史输入附件，90,741,531 字节](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/hc-cie-r2-20261128/HC_CIE_R2_Replay_Inputs.zip) · [逐候选和暗电流状态，124,526,330 字节](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/hc-cie-r2-20261128/HC_CIE_R2_States.zip)

**仅复查结果**可直接读仓库内的 CSV、JSON、NPZ。第二个附件含 340 个有限跃迁运行和 32 个动力学运行的全部材料候选，另有 21 份首次谱实验的扰动场及余项。解压时保留目录，候选出现在原 `results/imaging/`、`results/landscape/`，扰动场在 `results/spectra/field_states/`。[数组索引](results/spectra/array_index.json)

**重算**使用带历史研究代码的主仓新克隆。若仅从零运行某项，通常只须第一个输入附件。有些生成命令默认写回复算克隆的 `results/`，用新克隆可与已发表数组独立比较。按实际时间选根时，机器速度和系统调度会改变当时入池的候选。从已保存的前缀复算有相同的输入。

~~~bash
git clone https://github.com/migodam/b1-probabilistic-maxwell-imaging.git
cd b1-probabilistic-maxwell-imaging
unzip -n /path/to/HC_CIE_R2_Replay_Inputs.zip
unzip -n /path/to/HC_CIE_R2_States.zip
export HC_CIE_R2_INPUT_ROOT="$PWD" MAXWELL_KILLTEST_WORKSPACE="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
R=research/hc_cie_r2
~~~

本机环境为 macOS 26.6.2 arm64、Python 3.13.13、NumPy 2.5.3、SciPy 1.18.1、matplotlib 3.11.1、SymPy 1.14.0、mpmath 1.3.0、PyTorch 2.14.0，另使用 `threadpoolctl`。CPU 复双精度运行。网页重建另用 Python-Markdown。历史物理依赖及数据授权另列。[输入说明](INPUTS.md)

## 三项决定性结果

~~~bash
python "$R/src/toy_replay.py" --source "$R/replay_inputs/toy" --out hc_replay/toy
python "$R/src/toy_certificates.py" --source "$R/replay_inputs/toy" --out hc_replay/toy
python "$R/src/toy_material_interval.py" --source "$R/replay_inputs/toy" --out hc_replay/toy
python "$R/src/toy_feedback.py" --out hc_replay/toy
python "$R/src/spectra.py" --workspace "$PWD" --states-output hc_replay/fields
python "$R/src/spectra_feasible.py" --workspace "$PWD"
python "$R/src/verify_pair_bounds.py"
~~~

依次得到小模型坏极小与消融、真实算子的完整谱与球内材料复核、0017/0044 两件全部位置对的保守下界。数组索引可回查实际噪声、材料箱、父端点和实现版本。

## 相同历史端点与真实物理审计

~~~bash
python "$R/src/imaging_run.py" --config "$R/configs/imaging_v2.json" --out hc_replay/imaging --workers 2
python "$R/src/imaging_score.py" --sealed hc_replay/imaging --out hc_replay/score --figures hc_replay/figures
python "$R/src/landscape_run.py" --config "$R/configs/landscape_v1.json" --out hc_replay/landscape --workers 2
python "$R/src/landscape_score.py" --sealed hc_replay/landscape --out hc_replay/landscape_scores --figures hc_replay/figures
python "$R/src/landscape_probe.py" --sealed hc_replay/landscape --out hc_replay/probe
~~~

针对已发表的存档，额外结构评分与墙钟重放默认读取现有封存及标签，分别为 `structure_scores.py` 和 `wall_prefixes.py`，需要第二个状态附件。逐臂前景、背景、原端口、材料全图、活动盒、不同起点及截止规则见[有限跃迁](EVIDENCE_IMAGING.md)与[动力学](LANDSCAPE.md)。零/普通初值、BP、已知支撑和成熟端点只在相应实验内配对。二维前端可单独重算。

~~~bash
python "$R/src/real_weights.py" --workspace "$PWD"
python "$R/src/real_frontend.py" generate --workspace "$PWD"
python "$R/src/real_frontend.py" score --workspace "$PWD"
python "$R/src/radiation.py" "$PWD"
python "$R/src/noise_rank.py" --workspace "$PWD"
~~~

## 已知背景的条件读出及独立网格

~~~bash
python "$R/src/interaction.py" generate
python "$R/src/interaction.py" solve
python "$R/src/interaction.py" score
python "$R/src/support.py"
python "$R/src/model_gap.py" all --workspace "$PWD"
python "$R/src/cross_changes.py" "$PWD"
python "$R/src/paired_change_stats.py"
python "$R/src/change_figures.py"
python "$R/src/model_gap_figures.py"
~~~

`interaction.py` 把 24 个历史 CIE 材料当作公开参考，只对新改变评分。`support.py` 另做历史端点的支持数下界、全对排除和监测单格/闭式读出。`cross_changes.py` 使用前一步 `model_gap.py` 产生的细网参考端口与粗网公开数组。设计和种子在 `configs/`，已发表的新观测按对象与几何登记。[算法](ALGORITHM.md) · [按父材料配对](results/interaction/paired_object_statistics.json)

## 真三维与实测

~~~bash
python "$R/src/vector_change.py" --workspace "$PWD"
python "$R/src/vector_change_subspace.py" --workspace "$PWD"
python "$R/src/vector_change_figures.py"
python "$R/src/spectra_physical_limits.py" --workspace "$PWD"
~~~

监测沿用原 TensorDDA 的三分量和两个极化。独立搜索的输入、真支撑对照与后评分的入口见[三维实验](VECTOR_CHANGE.md)。外部 Fresnel 13 源分割需先从原作者授权渠道取得 `FoamTwinDielTM.exp`，设置 `FRESNEL_2005_DIR`，按[原件准备](INPUTS.md)生成读数。V2 是原尺度主合同，V1 保留作目标尺度对照。

~~~bash
python "$R/src/classical_run.py" --config "$R/configs/classical_v2.json" --out hc_replay/classical --workers 2
python "$R/src/classical_score.py" --sealed hc_replay/classical --out hc_replay/classical_scores --figures hc_replay/figures
~~~

实测公布允许共享的派生预测、掩码及评分。平铺输入附件完成后，在独立目录实算验证了历史观测、两种旧母模型、三维端口、数值小模型及区间证书。[已执行的平铺复算范围](INPUTS.md)
