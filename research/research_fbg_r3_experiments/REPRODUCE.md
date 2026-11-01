# 复算入口与证据版本

本目录相对于 B1 仓库根为 `research/research_fbg_r3_experiments`。复算同一结果应使用冻结公开数组、对应源码档案和原 PLAN；重新生成材料是另一项实验。本文仅提供命令，没有执行新实验。

## 环境

Apple Silicon CPU 使用已有环境，系统 `/usr/bin/python3` 没有本研究所需 NumPy。实际版本记录在 `audit/SOURCE_MANIFEST.json`：Python 3.13.13、NumPy 2.5.3、SciPy 1.18.1；另需已有 Matplotlib、threadpoolctl。

```sh
REPO="/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r2_publication/repository"
STUDY="$REPO/research/research_fbg_r3_experiments"
PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
cd "$STUDY"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

CPU 批次通常四个独立 worker、每 worker 单 BLAS 线程。保存同样的并发、硬件与预算以比较 wall；有限 wall 停止意味着不同机器不保证逐位相同路径。不要安装 CUDA 来运行小网格。

## 数据、源码与分类

| 阶段 | 固定输入／源码 | 分类与复算限制 |
|---|---|---|
| 原报告代数演示 | `input/FBG_R3_Reproducibility/`；`audit/pro_replay/` | 已有公式和两单元演示；不是新成像。已有 replay 回执足够，不需重复误差地板检查。 |
| 初始 E1／E2 | `data/e1_public/`、`data/e1_case_ids.json`；五个 `starts_4/6.npz` | 24 观测、12 材料合同，来自三个固定形状族；4×4／6×6 分别重新组装。E1 因评分器提前预加载旧 truth 的真实偏差，属于已披露受控开发，见 `BLIND_PROTOCOL.md`、`audit/TRUTH_ACCESS_DEVIATION.md`。 |
| 旧 E1 后续分支 | `results/e1_*/PLAN.json`；`audit/E1_OBSERVABLE_SOURCE.zip` | 包含初始容量失败、observable backbone、CSI/CIE 权重、多起点与固定总预算。结果版本分开，不以修复后代码覆盖历史失败。最早未单独保存源码的变体，以封存回执审计为准，不能把当前代码重跑冒称逐字复现旧版本。 |
| 已封存独立新鲜小网格确认 | `data/CONFIRMATION_FREEZE.json`、`data/confirmation_public/`、`audit/CONFIRMATION_SOURCE.zip` | 方法冻结在生成前，6 个新变形对象、两布局、clean/noisy，共24观测；8方法×3起点共576预测已封存。属于先前 affine physical-anchor N1 的确认，不是后面的 native_tsom 新算法。软件时间隔离不等于密码学保密。 |
| 新 native_tsom 分支 | `audit/NATIVE_TSOM_FREEZE.json`、`audit/NATIVE_TSOM_SOURCE_FREEZE.zip`、`results/e1_native_tsom/PLAN.json` | 仅6×6，固定 DPCS 与原生 mf2 APCS 仿射空间、192个实 current 变量；普通 Fourier 与任务选模／软释放／有限 seed 同池同容量。旧 E1 truth 已披露，属于受控开发。本文件不预判该新批次最终 N1 结论。 |
| E3 成熟基线诊断 | 历史 R2 `data_contract/stage4_public/`，本目录 `results/e3_*/` manifests／receipts | 旧12个已披露576材料对象；新跑成熟基线是历史诊断，不是 fresh holdout，也不是 N1/N2 的576成功。 |

评分后公开的无生成 seed 数组为 `data/e1_truth_after_freeze/` 与 `data/confirmation_truth_after_freeze/`，各自含 `mapping.json`。E3 使用 `research/research_fbg_r2_adaptive_feedback/truth_after_freeze/` 和本目录 `data/e3_eval_metadata.json`。在线 runner 不应导入 scorer 或读取这些材料／mask／评价映射。移动仓库时使用附带的 `data/e3_eval_metadata_portable.json`：仅将 truth_file 改为同一个NPZ的basename，以 --private-dir 指向原R2真值目录；case/object映射和数组不变，原审计metadata也保留。

`experiments/prepare_data.py` 是原始生成入口；`fresh_confirmation.py freeze → generate → run` 刻意分成三个命令。`generate` 使用独立生成熵并拒绝覆盖已存在目录。重放已发表问题**不要重新 generate**；已公开数组就是所需数据，不需要任何 secret seed、key 或私有生成状态。

## 历史源码重放与 E1 运行

不要把源码 zip 解压到当前研究目录。可在仓库 `research/` 下复制本目录到一个新的同级 replay 目录，再将相应 zip 的 `src/`、`experiments/` 覆盖到该副本；保留 `research/research_z_alg_r4_diversity` 等原始依赖。这样物理模块的 repository-relative imports 仍可解析。

精确任务以对应 `PLAN.json::tasks` 为准，而非 CLI 默认方法列表。对每条任务调用冻结 `experiments/e1_current_release.py::infer(task)`，仅将 `tag` 改为全新 replay tag、`public_dir` 改为同一公开数组的本机路径。运行前创建该 tag 的结果目录；使用四 worker 的 `ProcessPoolExecutor`、主入口 `if __name__=='__main__'`，像原 runner 一样记录返回 rows 和总批次 wall，全部任务完成后才写 `status='PREDICTIONS_CLOSED_TRUTH_NOT_READ'` 的 SUMMARY。不要混合不同源码 zip 的 rows，也不要覆盖已有输出目录。

当前 native_tsom 冻结配置的直接入口如下；`replay_native_tsom` 必须是新 tag，待其全部预测关闭后再评分：

```sh
"$PY" experiments/e1_current_release.py --tag replay_native_tsom \
  --methods cie_matched n1_fourier n1_energy n1_random n1_normal n1_task n1_soft \
    n1_finite n1_finite_fourier n1_finite_random n1_finite_krylov \
  --starts 0 1 2 3 4 --sides 6 --layouts uniform jittered \
  --backbone observable --affine-mode native_tsom --chart cie --joint-steps 80 \
  --wall 4 --max-lu 120 --max-rhs 2400 --workers 4
```

更强成熟 CIE 的原 mf1→2→3 日程保留在旧 `e1_baselines_initial`／其他已记录 PLAN 中；`cie_matched` 固定 mf2 是与当前 native APCS 容量匹配的额外对照，不能替代更强既有基线。

已封存 fresh confirmation 使用 `CONFIRMATION_SOURCE.zip` 与原 `results/e1_confirmation/PLAN.json` 重放。其生成／执行脚本故意不覆盖原目录；在 replay 副本用上述任务级新 tag 重放即可。读取已封存 SUMMARY 及重新评分不需要再次运行全部576个预测。

E0 新接线检查入口是 `experiments/e0_formula_checks.py`，现有回执为 `audit/E0_PRODUCTION_WIRING.json`；通常直接审计回执即可。`audit/N1_PURE_SPAN_OWNER_RECEIPT.json` 如存在，其 owner-reported 来源不等于原生执行日志，详见独立审阅。

## E2：精确有限核与强对照

从研究目录执行：

```sh
"$PY" experiments/e2_multibranch.py --tag replay_mechanism
"$PY" experiments/e2_multibranch.py --summarize-prefix results/e2_mechanism
```

第一条是新运行，第二条只读已有封存结果。四个匿名 case、五公共起点、k4/k6、三个公共支持、Q6/Q20／同核 block-TRF／同核六起点／随机／全材料多起点／匹配实际校正 seed 数的对照均在 `results/e2_mechanism_plan.json`。每次最多3个共同 full-material TRF seeds；实际 root、seed、LU、RHS 和小核费用分别保存。

入口 `src/algebraic_multiroot.py::run_multibranch(model,Y,m0,config)`；精确核为 `F0+O solve(I-diag(d)C,diag(d)E_A)`。所有源／接收用于排序，选定测量只用于求根。子集根、完整数据根、物理校正候选不同；不得把未达完整门的根改名为同数据物理解。原材料与全校正候选保存在 `_materials.npz`／`_candidate_materials.npz`，过滤及实际根间距离在压缩 ledger／summary。

## 封存后评分入口

评分器只读封存预测与已付费 residual，不新增 Maxwell 求解。先确认所有目标批次的 SUMMARY／E2 closure 完整，再使用明确 closure token。以下操作会写新的评分文件和 `results/COST_LEDGER.csv`，重放评分请用新 tag／replay 副本。

```sh
"$PY" experiments/score_and_plot.py --batch e2 \
  --results results/e2_mechanism_results.json \
  --private-dir data/e1_truth_after_freeze --authorize-closed PREDICTIONS_CLOSED \
  --tag replay_e2 --classification POST_CLOSURE_CONTROLLED_STUDY

"$PY" experiments/score_and_plot.py --batch e1 \
  --results results/e1_confirmation/SUMMARY.json \
  --private-dir data/confirmation_truth_after_freeze --public-dir data/confirmation_public \
  --authorize-closed PREDICTIONS_CLOSED --tag replay_confirmation \
  --classification PROSPECTIVE_FROZEN_SMALL_GRID_CONFIRMATION
```

旧 E1 或新 native_tsom 改用其已封存 SUMMARY、`data/e1_truth_after_freeze`、`data/e1_public`，classification 应含 `NOT_FRESH_BLIND`。未完成的新批次不能先评分。NRMSE≤.15、IoU≥.7、clean data≤.001／noisy≤.003 固定；online top1 按原物理分数，offline best@K 单列。对象、布局、噪声和初值不是独立样本；旧 Otsu-mask 门冲突只保留离线诊断，不改门营救。

## E3：已有 Windows CUDA resident queue

先读兄弟目录 `research_z_alg_cuda_pipeline/OPTIMIZED_OPERATOR_CATALOG.md`、`REPRODUCE.md` 及 R2 的 `REPRODUCE.md`。实际环境为已有 RTX4060 Laptop／Windows、`D:\python\python.exe`、PyTorch 2.6.0+cu124、complex128／float64；不安装驱动、CUDA、永久服务或另一个调度器。

本目录的 `experiments/e3_fullwave_imaging.py` 动态复用 R2 `FeedbackLane`、原 resident `SharedOperators/ObjectiveLane` 和兼容的 `hybrid_svd.py`；CIE port 在 `src/cie_tsom_cuda.py`。原生 source/receiver/current 算子共享一次上传，两条 lane 的有限 FIFO 完成一个即补一个，所有照明共享一份 LU。CPU 材料读出、等待、传输和最终物理关闭也计费。

实际两份清单在 `results/e3_limited_diagnostic_20261101/manifest_limited_digit2_attempt1.json` 与 `results/e3_mature_baseline_remaining_20261101/manifest_mature_baseline_remaining.json`。复制到**新文件**，只重定位已有 pipeline/R2/vendor/public 路径及全新 output_dir；保留 task、算法、45s/2000LU/60000RHS、精度与日程。原文件的 `resource_idle_confirmed=true` 仅是历史本次回执，不是未来 GPU 使用授权。

```powershell
# 在部署了本研究源码的现有目录；先检查资源，不中断其他 session。
D:\python\python.exe experiments\e3_fullwave_imaging.py replay_manifest.json --inspect-resources
D:\python\python.exe experiments\e3_fullwave_imaging.py replay_manifest.json
# 当前资源已协调、复制 manifest 的 idle 确认与新输出目录均有效后才执行：
D:\python\python.exe experiments\e3_fullwave_imaging.py replay_manifest.json --execute
```

不带 `--execute` 返回 NOT_RUN。两lane（或兼容四lane）不保证固定 SM 分配。共享准备、参考 LU、原始失败尝试与后续结果均应保留；不只合计成功任务。下载完整 finite-queue SUMMARY 和每任务 JSON 后，可用 `score_and_plot.py --batch e3 --results <两个SUMMARY路径> --metadata-file data/e3_eval_metadata_portable.json --private-dir ../research_fbg_r2_adaptive_feedback/truth_after_freeze --authorize-closed PREDICTIONS_CLOSED --tag replay_e3 --classification DISCLOSED_HISTORICAL_576_DIAGNOSTIC` 评分。

该 E3 是成熟基线差距诊断，**N1_E3／N2_E3 仍为 NOT_RUN**；历史原 8/12 与本轮历史数据重跑数分别报告。不得据此宣称新鲜576成功、N1/N2已扩大验证或科学升级门已通过。原 raster 的 row0 在上，绘图不改变材料数组。


## 最终native-affine新确认与事后组件

新确认独立于前一确认集：方法冻结06:45:53 UTC，之后生成六个新对象；648个预测于06:49:31 UTC封存，然后才读新目录GT。公共输入在 data/native_confirmation_public/，无seed真值在 data/native_confirmation_truth_after_freeze/。使用同一分布重新generate会产生不同对象，不是复算原结果。

当前代码与最终native源码一致时，可直接重放全部匿名输入，保持原计划预算，不读取truth：

~~~sh
"$PY" experiments/replay_plan.py \
  --plan results/e1_native_confirmation/PLAN.json \
  --public-dir data/native_confirmation_public --tag replay_native_confirmation --workers 4
"$PY" experiments/score_and_plot.py --batch e1 \
  --results results/replay_native_confirmation/SUMMARY.json \
  --private-dir data/native_confirmation_truth_after_freeze \
  --public-dir data/native_confirmation_public \
  --authorize-closed PREDICTIONS_CLOSED --tag replay_native_confirmation \
  --classification REPLAY_OF_SEALED_FRESH_NATIVE_AFFINE_CONFIRMATION
~~~

源码版本为 audit/NATIVE_CONFIRMATION_SOURCE.zip；较晚新增的replay_plan.py只是任务重放入口，没有改原求解函数。若对比更早变体，请按上述源码zip隔离重放，不能将新实现冒称原版。

新确认后只有一个针对已披露收益案例的24-task组件诊断：results/e1_posthoc_mode_controls/，入口 experiments/posthoc_mode_controls.py。它使用相同代码、预算和public starts，对 clean obs_7ad98676df61dbb6 与 noisy obs_b3b3447fe12c27c7 比较random/energy/normal/soft。明确分类 POST_DISCLOSURE_GAIN_CASE_COMPONENT_DIAGNOSTIC，不加进原648新鲜成功率，不外推统计结论。需要重放时同样使用replay_plan.py并给新tag。

全部任务入口使用旧公开观测和算子即可从未知材料出发计算候选；离线评分是另一步。最终主批次回执按execution_source和run_id去重，已引用的180个强基线native侧6网格记录不是新的物理运行。

