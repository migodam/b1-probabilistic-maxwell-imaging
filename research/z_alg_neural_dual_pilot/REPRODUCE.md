# Z-ALG 学习双乘子复现说明

本包保留可直接运行的研究路径。以解压目录为当前工作目录执行下列命令。打开 `research/z_alg_neural_dual_pilot/index.html` 可读研究结果，打开同目录的 `interpolation.html` 可读连续权重与接收孔径跟进。包中有 72 件训练或验证材料各自的五档配对噪声读数，以及主比较、低噪声补充和独立跟进的冻结复场。主比较用 32² 生成提议，以 96² 全维区间证书或已证明的 96² 限制上界作物理判定。

## 包内入口

| 要核对的步骤 | 入口与数据 |
| --- | --- |
| 样本与训练 | [训练方法](METHODS_TRAINING_ZH.md)、[无标签输入](data/synthetic_observations.npz)、[生成回执](data/generation_receipt.json)、[冻结方案](results/training_recipe_frozen.json) |
| 两个模型 | [测量特征权重](results/neural_dual_weights.pt)及[训练记录](results/neural_dual_training.json)、[固定特征权重](results/no_feature_weights.pt)及[训练记录](results/no_feature_training.json)、[独立验证输出变化](results/training_checkpoint_audit.json) |
| 主比较的 28 物体 | 论文一 6 件原生材料的 18 个配对全环读数及半径，另含 SRAP 22 件的原始 10% 全环读数及半径。主比较共 40 个观测、98 个观测和上限组合。Y 的行、列分别是 20 个接收点和 10 个平面波。数据由 [holdout.py](src/holdout.py) 读取。 |
| 单列低噪声补充 | 相同 22 个 SRAP 对象各新增 1% 和 2% 端口中心，在 \(u=0.75\) 评价。由 [holdout.py](src/holdout.py) 的 supplement 路径读取。它们是同一材料的配对条件。[冻结复读数](data/srap_supplement.npz) · [区间半径生成回执](data/srap_supplement.json) |
| 搜索及严格证书 | [提议入口](src/propose.py)、[梯度方法](COARSE_METHOD.md)、[原生证书说明](NATIVE_VERIFICATION.md)、[批量严格回放](src/score_batch.py)、[论文一固定对照](src/fixed_baselines.py) |
| 数值导数和限制上界 | [32² 解析导数回执](results/coarse_gradient_check.json)、[96² 独立方向检查回执](results/native_restriction_training_check.json) |
| 状态负方向与冷启动 | [光谱诊断](GEOMETRY_INERTIA_ZH.md)、[几何压缩回执](results/geometry_inertia.json)、[独立冷启动回执](results/cold_latency.json) |
| 连续权重和接收孔径 | [跟进结果](FOLLOWUP_RESULTS_ZH.md)、[独立 32 件冻结读数](data/followup_observations.npz)、[外部 8 件冻结读数](data/followup_external8.npz)、[权重线段](followup/NESTED_SEGMENT_COVER.md)、[接收孔径](followup/MASKED_PORTS_REPLAY.md)、[Galerkin](followup/GALERKIN_COVER_ZH.md) |
| 固定 64 区权重转移 | [独立 16 件冻结读数](data/transfer16.npz)、[生成回执](followup/transfer16_generation.json)、[冻结独立比较](followup/TRANSFERRED_WEIGHT_FROZEN.json)、[逐案严格结果](followup/results/transferred_weight_transfer16.json) |

论文一的依赖来自本地保管的公开验证包，随附三份公共算子文件、原生证书与固定对照的源码及回执、18 份配对读数以及 SRAP 22 件的读数与半径。配对清单保留本试验用到的 1%、2%、10% 原行。原 SRAP 生成器、离线前向校验及 11 个快照依赖也随包提供，依照实际 Python 导入路径放置。包内运行直接使用已有依赖和公开读数。材料阵列、干净场和用户原请求留在原工作区。跟进的离线评分 JSON 只保留预先生成的标量标签。成员路径和体积见解压根下的 `PACKAGE_CONTENTS.json`。 [论文一公开验证包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-optical-inertia-validation-v1)

## 环境

使用 Python 3.13，并安装 NumPy、SciPy、PyTorch、python-flint 和 threadpoolctl。训练环境的 NumPy 为 2.5.3、SciPy 为 1.18.1、PyTorch 为 2.14.0。区间核验使用 python-flint 0.9.0。重画图需要 Matplotlib，重编网页需要 markdown-it-py。历史轨迹保留原机付费墙钟和 CPU。优化器有时间上限，其他机器可取得不同的候选数和候选权重。

~~~bash
export PYTHONDONTWRITEBYTECODE=1
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PILOT=research/z_alg_neural_dual_pilot
~~~

## 使用已封存证据

- 提议文件位于 results/proposals，逐对象 JSONL 包含乘子、诊断与独立计时。迭代详情保存在 results/optimization_traces 的压缩 JSONL。严格区间、上界、未知状态和各项费用位于 results/native。补充队列在 results/supplement 下保持同一布局。
- 论文一固定权重回放在 [fixed_baselines.csv](results/fixed_baselines.csv) 与 [fixed_baselines.json](results/fixed_baselines.json)。固定对照沿用论文一的绝对权重和历史共享正定间隙。
- 严格过线的判据是原生证书的下端点大于冻结观测半径。96² 限制矩阵的严格上界只用于判断某一套乘子无法过线。32² 硬分数或未决结果不产生排除结论。

## 如需重算

用单独解压出的工作副本，可保留原包和原完整轨迹。默认训练和数学诊断会写入 results。重训可传入另一输出目录。若要用重训权重重跑盲评，应在工作副本中先保留发布权重，再将重训的两个选中权重放到 results 的同名路径。

~~~bash
python -B "$PILOT/src/generate.py" --workers 4
python -B "$PILOT/src/train.py" --mode neural --out "$PILOT/retrained"
python -B "$PILOT/src/train.py" --mode no_feature --out "$PILOT/retrained"
python -B "$PILOT/src/check_coarse_gradients.py"
python -B "$PILOT/src/native_check.py"
python -B "$PILOT/src/native_check.py" --training
python -B "$PILOT/src/check_native_restriction.py"
python -B "$PILOT/src/validation_swaps.py"
python -B "$PILOT/src/recover_unknown.py" --training
~~~

已完整的逐对象文件会被续跑入口跳过。若要从零重做搜索和严格证书，在工作副本中先将原轨迹和完成标记移到另一目录，保留 results/holdout_open.json。批量验证可在每个队列的提议生成完后运行。它也支持 --watch 模式追随提议的写入。以下命令分别运行主比较和补充比较，省略分片参数时处理本队列全部对象。

~~~bash
python -B "$PILOT/src/propose.py"
python -B "$PILOT/src/score_batch.py"
python -B "$PILOT/src/propose.py" --supplement
python -B "$PILOT/src/score_batch.py" --supplement
python -B "$PILOT/src/fixed_baselines.py"
python -B "$PILOT/src/fixed_baselines.py" --supplement
python -B "$PILOT/src/recover_unknown.py" --run
python -B "$PILOT/src/recover_unknown.py" --run --supplement
python -B "$PILOT/src/diagnose_inertia.py"
python -B "$PILOT/src/cold_latency.py"
~~~

冷启动入口读取原比较已封存的提议，每行另建原生验证器。状态诊断以独立固定子空间核查状态负方向。主比较图读取完整的 [逐例回执](results/cases.csv)，生成 [严格排除结果](figures/certified_outcomes.pdf)、[配对在线费用](figures/compute_reliability.pdf) 和 [粗网筛选与原生判定](figures/coarse_native_outcomes.pdf)。费用曲线只对双方均严格通过的同物体、同噪声、同上限配对先求比值，再按物体汇总。

~~~bash
python -B "$PILOT/src/summarize.py"
python -B "$PILOT/src/paired_tables.py"
python -B "$PILOT/src/figures.py"
~~~

## 权重线段、掩码和 Galerkin

以下命令读取冻结的独立 32 件、重新采集的外部 8 件、原验证 12 件及历史 28 件的读数。四组分别报告，配对噪声保留原半径。默认复用包中的 30×30 原生端口响应、九节点区间响应和 Galerkin 缓存。`jensen_port.py --prepare --rebuild` 重新支付 90 个原生驻点右端。先保留已有缓存，再移走 `segment_enclosures.json`，可用 `segment_cover.py prepare` 复算九点缓存，额外支付 180 个原生右端。`galerkin_cover.py freeze` 写入冻结规则，运行 `all` 或指定批次时会按需制备压缩表。若要对初次制备计价，先在工作副本中移走 `galerkin_enclosures.json.gz`。接收掩码的生产路径从九节点响应中复用几何。三个独立原生核查另计实际原生右端。

~~~bash
python -B "$PILOT/followup/jensen_port.py" --validation
python -B "$PILOT/followup/jensen_port.py" --observations "$PILOT/data/followup_observations.npz"
python -B "$PILOT/followup/jensen_port.py" --historical
for cohort in validation fresh historical external8; do
  python -B "$PILOT/followup/segment_cover.py" "$cohort"
done
python -B "$PILOT/followup/masked_ports.py"
python -B "$PILOT/followup/masked_ports.py" --continuous-only
python -B "$PILOT/followup/masked_segment.py"
python -B "$PILOT/followup/galerkin_cover.py" all
python -B "$PILOT/followup/galerkin_cover.py" historical
python -B "$PILOT/followup/galerkin_cover.py" paired_summary
python -B "$PILOT/followup/plot_results.py" --aperture --diagnostic
python -B "$PILOT/src/render_report.py"
~~~

已冻结读数足以重复上述评分。下列离线入口可在另一工作副本重建低噪声补充及两个跟进批次。它们产生读数、区间半径和离线评分记录。网络提议和证书入口只读取公开读数与半径。

~~~bash
python -B "$PILOT/src/generate_srap_bands.py"
python -B "$PILOT/followup/generate_fields.py"
python -B "$PILOT/followup/generate_fields.py" --archived8
~~~

## 固定提案跨材料转移

[转移方案](followup/TRANSFER16_DESIGN.json)把一份既有的 64 区物理权重与归档固定权重统一放在 96² 全环的同一误差球中比较。独立 16 件材料各有三个配对噪声视图。逐案的接受只用读数和半径，适用的状态证明及端口响应共享。生成器可在另一副本重建读数和离线标量。复现已经冻结的比较时，直接执行下组评分命令即可。

~~~bash
# 独立的离线数据生成入口
python -B "$PILOT/followup/generate_transfer16.py"

# 只用已冻结读数及噪声半径
python -B "$PILOT/followup/transferred_weight.py" transfer16
python -B "$PILOT/followup/transferred_weight.py" retrospective
~~~

`transferred_weight.py prepare` 单独准备共享原生响应。已有缓存时沿用共享证明。`all` 依次评分独立 16 件和各历史对照，分开记账。已有权重无需再训练。

## 逐份带噪数据重新优化的 16 件对照

四个单线程进程生成 48 份一秒梯度提议，全部完成后六个单线程进程核验。发布包已包含选中权重、全轨迹与严格区间。原机的墙钟截点会受硬件影响，跨机复核同一份候选应直接使用所附的冻结提议。下列命令适用于 Bash 或 Zsh。[比较口径和费用](followup/ONLINE_TRANSFER_ZH.md)

~~~bash
for i in 0 1 2 3; do python -B "$PILOT/followup/online_transfer_propose.py" --shard "$i/4" & done
wait
for i in 0 1 2 3 4 5; do python -B "$PILOT/followup/online_transfer_score.py" --shard "$i/6" --proposal-workers 4 & done
wait
python -B "$PILOT/followup/online_transfer_score.py" --assemble 6 --proposal-workers 4
python -B "$PILOT/followup/plot_online_transfer.py"
~~~

## 重打包

[重打包入口](package_release.py)从原工作树中的本地公开档案取依赖。已解压的交付包使用其中随附的依赖。先用 `--check-inputs` 检查最终文件清单，需要生成本地 ZIP 时再去掉该参数。

~~~bash
python -B "$PILOT/package_release.py" --check-inputs
python -B "$PILOT/package_release.py" --output Z_ALG_Neural_Dual_Pilot.zip
~~~
