# Research W：复算与重跑边界

默认复现路径是读取已保存、已付费的数组，复算概率数值或重画结果。它不调用 Maxwell、LU/RHS、Book 或物理队列，不生成新观测，也不读取真值。数值读回通过只说明缓存数组与保存参考相符，不等于科学 GO、后验校准证明或连续 PDE 认证。

## 1. 当前目录与现有环境

工作树：

```text
/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1
```

本文命令从阶段目录执行：

```sh
cd "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1"
export WR_PYTHON="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
```

现有 Python 环境包含 NumPy、SciPy、Matplotlib。默认数组复算只需要 NumPy/SciPy 和本目录 `src`，验证器自行加入该路径。无需 `.deps`、skfem、历史 Maxwell 后端或安装新环境。

物理 launcher `tools/run_job.py:12–14` 则依赖以下路径：

```text
本阶段 ROOT/src
本阶段 ROOT 的父目录/research_w_continuation/src
本阶段 ROOT 的父目录/research_w_continuation/vendor
/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/Research_W/research-w-continuation/.deps
```

最后一项是现有 WR 环境的**绝对路径**，包含 skfem 等旧依赖，不是可迁移的相对路径。移动到另一台机器后，仅复制当前工作树或已有虚拟环境不足以保证物理重跑。必须明确恢复匹配的历史后端/vendor、依赖与源代码快照，并在新的授权工作副本中审阅路径适配。这里未修改正在运行的 launcher，也没有改写七个冻结物理/decoder 文件。缓存概率复算无需这些物理路径。

## 2. 默认：已付费 P1A 数组复算（0 新 RHS）

```sh
"$WR_PYTHON" tools/verify_w_probability_arrays.py
```

唯一数据读入清单：

- `results/P1A/rules/full_5.npz`、`full_9.npz`；
- 五个检查点对应的 `cdf_0_5.npz` 到 `cdf_4_5.npz`，以及对应 order9 文件；
- `results/P1A/data/OBSERVATIONS.npz` 中的 `y` 和 `parent_ids`；
- `results/P1A/reference/RESULTS.json`，用于逐条数值比较。

验证器不调用 `wr.prob_reference.main`，也不构造任何 Book。其积分使用纯代数 `wr.posterior_math.log_marginal_likelihood`：480 个实观测维度采用每实分量 sigma=.01 的 Gaussian 密度常数；在所有频率、发射源、接收器上共享一个 `gamma=a exp(i phi)`，仅积分一次；a、phi 的均匀先验求积权重归一化为 1，不使用 profile 最大值代替积分。GL11/21/41/65 按保存参考的全规则最大 log-likelihood 变化 2e-4 停止。

材料 full5/full9 权重总和应为 1；CDF 使用直接截断区间数组，保留原分支概率和原密度，分子权重不得重新归一为 1。验证器独立重建这些正权重与坐标，检查截断先验质量。每条观测复算 branch 事件概率、五个 order5/order9 CDF、归一化、位置均值/方差、对比度均值、绝对 log evidence、阶数差与保存资格标志。32 条观测、16 parent×2 噪声副本全部保留，不过滤不合格对象。

每次只写入新的 `AUDIT/array_replay/<UTC>/`，保留 `MANIFEST.json`、`RULE_NORMALIZATION.json`、`RESULTS.json`、`SUMMARY.json`；异常写 `FAILURE.json`。stdout 逐条显示读回结果和最终 CPU/wall、peak RSS、新 RHS=0。默认每批16个材料点，最多允许32，以控制临时数组；最终明确检查 peak RSS<1GB。比较容差为绝对 2e-9，归一化容差为 2e-14。原 reference、规则、协议、账本与物理数组一律不改。

程序结束码：0 为本次数值读回和内存检查通过；2 为读回/内存检查不通过；输入或执行异常返回非零并保留审计失败记录。这些结束码不是科学门判断。

## 3. 默认：从保存结果重画（0 新 RHS）

```sh
"$WR_PYTHON" tools/build_w_results.py --final --probability-figures
```

此命令仅生成/更新展示产物，例如 `index.html`、图表、CSV 和结果摘要；不会运行物理、补造缺失结果或读取新的真值。已有 Oracle/离线派生结果仍须按其原标签阅读。概率图只读保存的 reference/inference、point validation，以及已完成离线 scorer 的汇总表；不启动 scorer，不重新访问真值。若要保留原展示文件，请在证据副本执行重画。物理队列正在运行时，展示读回可能是 live snapshot；缺失项目必须保持 NOT_RUN。

`--final` 只用于父线程已经写入有效 `FINAL_OVERVIEW.json` 的情形，显示选项不改变科学结论。最终网页已按此入口构建，运行结果见RESULTS_REPORT.json；此显示选项不修改科学门或原账本。

## 4. 新授权后才能重跑物理

以下模块不是默认复算入口：`wr.prob_bank`、`wr.prob_data`、`wr.prob_reference` 的主流程、`wr.progressive_full_bank`、`wr.progressive_point_validation`、`wr.continuous_point_control`、`wr.green_reuse_control`、`wr.endpoint_validation`，以及科学队列与协调器。某些主流程含数组步骤，但也可能触发数据生成、缓存补点或求解，不能用于“免费重放”。

重跑 Maxwell/LU/RHS 必须获得明确的新授权，并由父线程在独立、可审计的新 campaign/工作副本规划来源、阶段上限、串行队列、源代码和配置快照、gate 与账本。当前 `configs/AUTHORIZATION.json` 和阶段/root ledger 仍是旧研究的真实费用依据。新授权不通过删除 `*_SEALED.json`、`CLOSED.json`、job receipts、FAILED/partial 文件或清空 ledger 实现；也不得复制旧名称制造一次“未支付”的重试。旧封账、失败、中断和已支付 RHS/CPU 永久保留。

现有 `tools/run_job.py` 会构造 Book、独占 job 目录并将整个子进程 CPU 统一写入 root ledger；物理操作另记 RHS/LU，嵌套 solver CPU 不能再次加进总费用。它不是默认数组复算的 wrapper。新物理工作仍要求单个 physics worker；缓存 tuple 复用必须满足精确坐标与原 literal source/public family 合同，不能把近似坐标或旧 LU 当成新候选的免费响应。

完整物理独立重跑还需要冻结的公共 Family 墙参数、原几何/网格约定、历史 backend/vendor 与依赖；已有数组读回不能证明这些输入在新机器上保持数值等价。任何真值生成或离线评分只能走其单独授权角色，在线代码不得读取真值或生成计划。本说明不提供绕过 gate/预算或恢复已封账 stage 的命令。

## 5. 已执行的一次数值读回

2026-10-23 UTC（本地 2026-10-24）默认数组验证结果：`ARRAY_NUMERIC_READBACK_PASS`，32/32 观察、12 个规则；与保存值最大绝对差 `2.2370993946196904e-13`。CPU `2.333737 s`，wall `2.4111774158664048 s`，peak RSS `63,733,760 bytes`，新增 RHS `0`。证据见 `AUDIT/array_replay/20261023T173334_761962Z/SUMMARY.json` 及同目录逐条结果。未读取 GT，未构造 Book，未执行物理或重画；本结论限于已保存数值读回。


## 6. 网页小包、复场Release与离线图库

`MaterialHypothesisPackage/EXAMPLE.json`保留原生格式；其引用路径在同目录镜像，下载整个小包即可获得实际650节点响应、joint权重、8份SMC轨迹及代表材料。完整Release另含全部观测、chain与材料，以及第一joint见证cal/joint×两世界×五频的20份完整复场。其余本地full fields不会假装已经上传；Release库存说明各项范围。

三个绘图入口同时放在本阶段`tools/`，兼容repository布局和Release并列stage/history布局，不依赖全仓库其他研究的工具。已存在图表建议在证据副本中重画；风险图可使用`--run-id <新的标签>`保留原图。controlled gallery只读取固定observation000的GT参数行，必须先检查offline/inference/point/control四个完成门，是明确OFFLINE_RENDER；不属于默认不读GT的概率数组复算。

记录环境见ENVIRONMENT.json。requirements-replay.txt只用于缓存统计/绘图，requirements-physics-recorded.txt另记录实际FEM依赖；记录版本不是跨平台数值等价证书。物理launcher仍保留实际执行时的路径合同，没有为便携性改写已执行实验。
