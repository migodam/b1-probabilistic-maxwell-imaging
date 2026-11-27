# 从原始复读数复算

完整附件含本轮材料数组、25 口复数散射矩阵、两端误差、搜索失败轨迹、封存提案与带噪采集，以及四套原生物理算子和既有 CIE 最小依赖。请在空目录解压。Git 仓库存放阅读页、代码、图与数值回执。原始 NPZ 以研究附件分发。[完整计算档案](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/mimo-material-compensation-20261128/mimo-material-compensation-20261128.zip) · [旧依赖与恢复说明](replication/README_ZH.md)

## 核心结果的一次冷重放

这些命令从解压后的根目录运行。实际冷重放使用 macOS arm64、Python 3.13.13、NumPy、SciPy、python-flint、Matplotlib 与 threadpoolctl；没有使用 GPU。向外误差中的浮点计算要求 IEEE binary64，严算部分使用 Arb。

```sh
python3.13 -m venv .mimo-replay-venv
PY=.mimo-replay-venv/bin/python
STUDY=research/mimo_material_compensation
"$PY" -m pip install -r "$STUDY/replication/requirements.txt"
"$PY" -B "$STUDY/replication/run_replays.py"
"$PY" -B "$STUDY/replication/restore_prior.py" .
"$PY" -B "$STUDY/controls/same_data_path.py"
"$PY" -B "$STUDY/environment/decision_limits.py"
"$PY" -B "$STUDY/environment/report_decisions.py"
"$PY" -B "$STUDY/environment/fluence_loop.py" check 32 96
"$PY" -B "$STUDY/environment/fluence_modes/reattack/verify_dual.py"
"$PY" -B "$STUDY/environment/fluence_modes/reattack/verify_loaded.py"
```

冷重放在临时空目录内复算 32² 材料路径、96² 同幅度两类下界、七件各 6,435 个切分、以及新增读数下的两个外区预算对照。不需要读旧目标真值。重新做离线材料误差评分时才另行恢复历史合成真值。工作目录为独立副本，收敛过程和费用会随机器变化。[代表性重放的原始回执](replication/replay_receipts/replay_receipt.json) · [四个物理算子来源](replication/prior-manifest.json)

## 各项结论的源文件

| 要复算的物理结论 | 入口与依据 |
|---|---|
| 原 Pro 控制、0.01 与 1.985 及其证明区别 | `controls/compensation_small.py`、`compensation_exact_certificate.py`、`structural_checks.py` 应在 controls 目录运行，原始脚本把代数检查写到当前目录；`same_data_path.py` 读浮点节点评估严格校正管。[397 段回执](controls/same_data_path.json) |
| 已有原生带噪读数的连续路径及 B 重拟合 | [路径执行顺序](paths/REPRODUCE.md)。每个 `reduction.npz` 都带原生材料、支撑和可验证的传播算子。连续折线、校正布局、中点复算与单点见证各自保存。 |
| 6,435 切分、条件 Fisher、固定对最大相干波前、重新攻击 | [拓扑复现与十个比较臂](topology/REPORT_ZH.md)。`design_frozen.json` 与 `new_observed.npz` 分开。 |
| 任务方向、形状自由反例与局部两端 SOCP | [外区预算执行顺序](budgets/REPRODUCE.md)。方向的顺序为所有实部再所有虚部；加权距离显式使用格面积。 |
| 单个被动探头的一侧全负载界 | [负载圆盘与原生正演](environment/passive_replay.py)、[两状态矩阵核验](environment/bounded_real.py)。每个待审的端口对须有匹配的原生 `cover` 回执。 |
| 两类所有像素可独立变化的 460 MHz 下界 | [公开规格、开路转移和独立裸波](environment/frequency_classes_REPORT_ZH.md)、[96² 核心增幅 2](environment/matched_amplitude_REPORT_ZH.md)。 |
| 同一份材料隐藏两个频率 | [共享两频寻解](environment/joint_frequency.py) 与 [按相干波前再优化](environment/joint_coherent.py)。`common_*`、`joint_wave_*` 的相容性以对应 `cover` 和 `kyp` 文件为准。 |
| 已采样品的宽盒反例 | [八个固定核心分支的再拟合](environment/postmeasurement/REPORT_ZH.md)。复算先用 attack.py 读此前封存的旧新读数，再用 expand.py 增加外部自由度。报告分别列出两球，不以两者平均值验收。 |
| 新样品主动判读 | [动作筛选](environment/measurement.py)、[低功率测量](environment/energy_loop.py)、[原生 96² 同增幅闭环独立复现](environment/low_energy_96_amplitude2/REPORT_ZH.md)。复算已封存读数可直接运行判读步骤。重演新采集须另开干净工作副本，保留原提案。 |
| 从整类方向产生一发一读的真实波权 | [单发功率分配与盲测](environment/fluence_modes/REPORT_ZH.md)。coherent_modes.py 32 96 独立复算整类下界；fluence_loop.py 依次运行 choice、acquire、decide、replay，后接两个网格参数。对封存包运行只读的 check 32 96 即可重放实际读数。重做 choice/acquire/decide/replay 四步须使用干净输出目录。四端版本见 [四相硬件比较](environment/coherent_readout/REPORT_ZH.md)。 |
| 原生 96² 一发后的宽盒两解 | [130 维响应对偶搜索与完整重算](environment/fluence_modes/reattack/DUAL_RESPONSE_ZH.md)。先运行 dual_response.py unit_415 unit_862，延续低规格难例用 --resume unit_862，裸空间 25 端重放用 verify_dual.py，另一路独立的加载 24 端重放用 verify_loaded.py。background_budget.py 从已采单标量及原生 460 MHz 工作点复核所需外区改动范围。 |
| 同机开销 | 区分搜索、验证、投影和采集。每个实验证据文件的 `cost`、`wall_s`、`cpu_s` 记录本进程用量。[独立原生内存计量](environment/resource_replay_96.json) |

以 96²、核心增幅 2、400 MHz 为例，下列顺序重新求解两件材料的完整物理波形并验证全负载上界。若只需重验存储矩阵，可直接运行第二条命令。

```sh
"$PY" -B "$STUDY/environment/passive_replay.py" 96 2.0 wave_96_2.0_1.npz --bands 400 460 --steps 4096
"$PY" -B "$STUDY/environment/bounded_real.py" verified_responses_96_2.0_wave_96_2.0_1_400.npz
"$PY" -B "$STUDY/environment/matched_amplitude.py"
"$PY" -B "$STUDY/figures/make_composites.py"
```

精确等场、同一带噪读数内的整段、固定材料对的未来距离与覆盖两个材料类的未来距离在回执中分别列出。检查前者用原观测的 \(r\)，区分两份可能的新读数用 \(2r\)。数字界只指向保存的离散物理和声明的误差合同。[科学报告](RESEARCH_REPORT_ZH.md) · [局限与交接](NEXT_PRO_HANDOFF.md)

[各独立结论的机器可读汇总](results/main.json) · [从回执重新汇总](results/summary.py)
