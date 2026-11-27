# 从公开档案复算材料补偿

本目录携带旧研究的最小输入闭包。4 个实际使用的名义端口算子与原始发布附件中的字节一致。旧带噪观测、公开材料锚点和补偿候选来自上一版发布附件。候选压缩包保留本轮脚本读取的原数组键，数值与类型原样保留。冻结的 Fourier 物理算子、Arb 校验器、Dual-CIE 读出和 Fourier 子空间代码均来自公共提交 `5f806c9ac2cd568d28ea23bea02f6ccd6ee32c93`，21 份源代码和配置同时与旧研究附件逐文件一致。[旧版发布页](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-blind-uq-mimo-20261127) · [逐文件来源](prior-manifest.json)

## 直接复算

将这一轮研究附件按原目录展开，使当前目录内存在 `research/mimo_material_compensation`。在空的工作目录内复算，或使用下面的冷重放程序自动创建隔离目录。Python 实测版本 3.13.13，三个物理重放均限制为单 BLAS 线程。

```sh
python3.13 -m venv .mimo-replay-venv
PY=.mimo-replay-venv/bin/python
R=research/mimo_material_compensation/replication
"$PY" -m pip install -r "$R/requirements.txt"
"$PY" -B "$R/run_replays.py"
```

冷重放只复制本轮所需脚本与任务证据，然后在空目录恢复 57 个公共依赖文件。它分别重算材料路径的全段证书、原生 96² 同幅频率实验、旧材料的穷举排名，以及新数据下的严格预算见证。程序不会重新寻优旧材料，也不读取旧合成真值。复算产物与逐进程计时写入 `replication/replay_receipts/`。可用 `--keep-workspace /path/to/an/empty/dir` 保留冷重放材料数组和工作树。

需要在展开的工作树直接调用各实验脚本时，先恢复旧依赖。恢复脚本只创建缺失文件，已有一致文件保留原样，发现不同就停止。数据固定恢复到 `research/z_alg_optical_inertia/blind_uq_mimo_theory/`，与新研究目录为同一 `research/` 下的两条路径。

```sh
"$PY" -B "$R/restore_prior.py" .
```

默认包 `prior-core.zip` 为 55,634,576 字节。另一个 6,073 字节的 `prior-scoring-opt-in.zip` 只含已经公开的旧合成材料 `n32_c02` 和 `n32_c07` 的 `material`，仅用于重生成这两个对象的追加观测或独立评分。复用追加读数、设计和材料时无需恢复它们。需要重生成时，才使用 `restore_prior.py . --include-sealed-scoring`。在线设计入口读取各对象的带噪 `initial_A.npz`，独立模拟采集在提案落盘后从封存合成材料计算新读数。

## 直接使用旧发布附件

已有旧附件时可绕过两份精简压缩包。三个旧文件名、下载位置、包内路径以及本轮读取的键均由 `prior-manifest.json` 登记。包内的历史前缀是 `blind_uq_mimo_theory/`，恢复后增加工作树前缀 `research/z_alg_optical_inertia/`。以下命令只准备当前实验所需的成员。最后一个评分附件只用于上述选择性重生成。

```sh
NEW_WORKSPACE="$(pwd)"
OLD_ASSETS="$(mktemp -d)"
REPO=migodam/b1-probabilistic-maxwell-imaging
TAG=z-alg-blind-uq-mimo-20261127
gh release download "$TAG" -R "$REPO" -D "$OLD_ASSETS" \
  -p blind-uq-mimo-study.zip -p blind-uq-mimo-operators.zip
"$PY" -B "$R/restore_prior.py" "$NEW_WORKSPACE" --original-asset-dir "$OLD_ASSETS"
```

- `blind-uq-mimo-study.zip` 为 300,536,885 字节，包含公共代码、原带噪观测和历史候选见证。[旧研究附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-blind-uq-mimo-20261127/blind-uq-mimo-study.zip)
- `blind-uq-mimo-operators.zip` 为 517,338,723 字节。这里只用 `n32_f400_a08`、`n32_f460_a08_c048`、`n96_f400_a08`、`n96_f460_a08_c048`。[原物理系数附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-blind-uq-mimo-20261127/blind-uq-mimo-operators.zip)
- `blind-uq-mimo-scoring-truths.zip` 为 3,817,966 字节。选择恢复封存评分材料时再下载，并同时在恢复命令加 `--include-sealed-scoring`。[旧评分附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-blind-uq-mimo-20261127/blind-uq-mimo-scoring-truths.zip)

精简包里的归档 `aperture_control/geometry.npz` 保存 `E0_both` 和 `S_both`。恢复时从发布的 32² 400 MHz 算子拷入 `G`，其数值与旧几何文件的 `G` 逐元素完全相同，省去一份 18 MiB 的重复矩阵。其他用于对照的较大旧实验结果不需展开。

## 依赖清单

| 本轮入口 | 旧版外部输入 | 本轮档案输入 |
|---|---|---|
| `controls/` | 无。Python 和 SciPy 的小模型，严格算术用 python-flint。 | 同目录的配置和已保存见证。 |
| `paths/` | 32² 和 96² 的 400 MHz 算子，原带噪 `c01`，旧 GN/稳定零区与扩展斑块，旧 B 观测、B 种子和联合 A/B 种子。读取 `imaging/shared`、`verified_wave`、`native_replay`、`aperture_control`、`aperture_joint`。 | `*_reduction*.npz`、`*_paths*.npz/.json`、B 轨迹和后续联合重拟合产物。 |
| `budgets/` | 同一对 400 MHz 算子，原 `n32_c02`、`n96_c02` 的带噪观测、查询掩码和原实材料锚点。读取 B 投影工具、模型和严格验证器。 | `freeze.json`，对偶和场数组 `raw_*.npz`，bridge 切分原数组。 |
| `environment/` | 四个 400/460 MHz 算子和两件原 `c01` 的半径/观测 ID。使用旧裸空间 disk 脉冲 Green 和 Arb 验证器，新增第 25 个端口系数在本轮脚本里构建。 | `phase_*`、`wave_*`、`cover_*`、`verified_responses_*` 的材料和物理回执。主动模拟补测需本轮对应提案、噪声和频段读数。 |
| `topology/` | 两个 400 MHz 算子；7 个旧固定对、两件原生旧扩展对及各自原 A。Dual-CIE 加深方案读取旧 Fourier 电流子空间、实材料读出和正则化。只有重生成两个历史主动对象的新数据会读独立评分包。 | `campaign.json`、`loop/*/design_frozen.json`、原始初读与追加读数、逐臂材料及验证器输入。 |
| `inputs/pro/` | 小模型用 NumPy/SciPy；原参考文稿的可选重新排版用 lxml、Pandoc 和 XeTeX。 | 同目录的 TeX、模板、文献标注。 |

各脚本的直接导入、文件读取表达式见 [逐脚本索引](script_dependencies.json)。冷重放从精简包导入 `NativeFourierAPCS`、`Regularizer` 和 `_readout`，并记录提供文件的真实路径。整个依赖链无需原私有外置盘、旧 Git 工作树分支和未登记的编译插件。

## 冷重放记录

| 独立重放 | 墙钟秒 | CPU 秒 | 验收输出 |
|---|---:|---:|---|
| 32² 材料路径 | 1.515 | 1.490 | 3 条折线、41 段严格通过。最坏段上界 0.9893r；原生独立中点也在球内。 |
| 96² 同幅频率 | 7.021 | 6.773 | 400 MHz 配对距离上界 0.3173；460 MHz 在两端点盒半径 0.01 下的严格距离下界为 1.2759，单位为两倍噪声半径。 |
| 七对旧材料穷举 | 7.720 | 7.477 | 每对 6,435 切分，原始端口公平性和高斯投影回执通过。 |
| 新读数预算见证 | 1.461 | 1.434 | B 与补偿臂的六个 32² 全波材料重算，旧/新双球、常量/自由 ROI 和 0.35 外区预算全部核验。 |

这些用时来自新临时目录的回放，不是原实验求解或采集费用。频率重放运行本轮 `matched_amplitude.py` 的严格端点盒证明入口。材料路径复用本轮已经封存的降阶锚点并对每个全段、原生中点复查。拓扑预算回放对六份材料重新做原物理分解和严格残差审计。仍需重新寻优时，按各子实验复现入口执行，并单列原始搜索预算。[冷重放回执](replay_receipts/replay_receipt.json) · [频率证书](replay_receipts/native_frequency.json) · [历史排名](replay_receipts/rank_summary.json) · [原附件恢复测试](replay_receipts/restoration_modes.json)
