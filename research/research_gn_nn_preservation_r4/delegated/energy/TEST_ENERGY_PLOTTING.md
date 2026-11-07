# New independent TEST energy plotting handoff

状态：**绘图 CLI 已实现，仅静态语法检查；未运行、未渲染任何 TEST 图，未打开 TEST 结果、评分或 GT。** 实现依据仅为 `experiments/test_energy.py`、`src/energy.py` 及公开冻结方法清单 `audit/FINAL_RUNS.json` 的字段/命名。父线程在有限 TEST 预测封存、独立能量导出完成后运行，并负责真实图像检查与科学解释。

## 命令

仓库根目录，实际 `results/test_energy` 导出完成后：

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_gn_nn_preservation_r4/experiments/plot_test_energy.py
```

可用 `--input-dir` / `--output-dir` 指定实际导出和图片目录；默认分别为 `results/test_energy`、`figures/test_energy`。`--no-sensitivity` 仅关闭附录敏感性图。没有调度、等待、训练、重算、封存重验或成功阈值。主能量 CSV 不存在时立即 `NOT_READY` 退出，不生成虚构结果。没有运行任何输入数据测试或公式恒等式检查。

## 只读输入和实际字段

| 输入 | 用途与边界 |
|---|---|
| grouped_absolute_energy.csv | 读取 run/noise/quantity/stage/projection/cases，以及 before/after/delta/update/signed_cross 的已存 sum。按 cases 汇总全部 group/family/amplitude，不平均“组均值”而改变权重。 |
| mother_paired_bootstrap.csv | 只读 cohort=ALL 的已存 delta_mean、delta_ci_low/high、mothers。按 run/seed/noise 单独展示，不重新 bootstrap。 |
| current_realization_gap.csv | 只读已由 kind=NN 分支生成的 auxiliary/raw 配对标量。字段包括 state_defect_energy、current_realization_gap_energy、internal_realization_gap_energy、aux_receiver_update_energy、receiver_realization_gap_energy、signed_receiver_cross、physical_receiver_update_energy。 |
| data_residual_relation.csv | 用 obs_id 匹配已有 raw_1/final_top1 记录；读取 true_S_field_error_energy、observed_misfit_energy、noise_energy、signed_residual_noise_cross。它们与 receiver-update 字段是不同参考系，分栏显示。 |
| NOT_RUN.csv（可选） | 保留所有原因事件；可从现有数据记录确定 noise 的按 noise 计数，无法确定者单列 noise_unassigned，不猜测。producer 对空列表不写文件，因此文件缺席与“已确认零失败”不强行等同。 |
| RECEIPT.json（可选） | 只取导出状态、共同 GN 参考和已有空间元数据；不读取 seal、GT 文件或评分 JSON。 |

脚本不读取 per_case_energy.csv、GT 数组、operators、predictions/checkpoints 或 TEST_SCORES.json；不导入项目 physics/energy/network 模块。它仅在父线程授权的导出表可用后读取这些 GT 派生汇总，不回到 GT 进行任何操作。

## 固定的方法身份

主 seed7、两层级联固定为：

| 标签 | 实际 run |
|---|---|
| Voxel | v_physics_seed7_cascade2 |
| Jfull | j_full_seed7_cascade2 |
| Jnull | j_null_seed7_cascade2 |
| Gweak | g_weak_seed7_cascade2 |
| Gdark | g_dark_seed7_cascade2 |
| feedback1 | recurrent_feedback1_seed7_cascade2 |
| Gweak+const | g_weak_seed7_cascade2_const |

**命名歧义已解决：** Gweak+const 在表中 arm 仍可能是 `g_weak`，故严格按 run 区分，不把它混入原始 Gweak。其它 const/internal/finetune、hinge、refreshed-q 和 baseline 不混入这个预声明主图集，也不会因结果好坏自动替换方法。

三种子确认按冻结计划仅展示 Voxel/Jfull/Gweak 的 seed7/17/27。Jnull/Gdark 没有在此虚构另外两种子。主图与确认图引用同一 seed7 数据，method_counts.csv 对同一 run 不重复列成两个独立方法。

## 默认图片（十张，分开保证可读）

1. `strong_dark_absolute.png`：clean/noisy 分行，内部 strong19、接收 strong9、D509 分列；common GN→实际验证 raw_1→final 的绝对真实 Maxwell-current error 均值。
2. `strong_dark_paired_delta.png`：对应 signed raw/final−GN 均值变化及已有母体配对区间。
3. `weak_null_absolute.png`：内部 weak557、接收 weak39、结构 null528 的相同绝对口径。
4. `weak_null_paired_delta.png`：对应带符号配对变化。
5. `twofold_clean.png`、6. `twofold_noisy.png`：七方法逐行，真实数值交集＋余项47、I→O、O→I 三定义严格分列；每格保留 GN_R、GN_F、Δraw、Δfinal，负值不删。两个 GN 值分别属于 raw/final 可用母体集合，不能因覆盖不等而冒充同一个均值。
7. `auxiliary_to_physical_chain.png`：只连接匹配 NN 的 GN physical→auxiliary head→raw physical→final physical；列为 total、internal strong19、receiver strong9、D509。辅助头仅称 discrepancy，不称真正 Maxwell-solution current error。
8. `realization_state_receiver_accounting.png`：按 noise 分行，三列为原始 NN state/realization gap、相对共同 GN 的 receiver-update/gap/signed-cross、raw 真实 receiver-error/observed-misfit/noise/signed-cross。四个现成项并列显示，不相加做闭合验证。
9. `three_seed_confirmation.png`：三空间×两 noise，Voxel/Jfull/Gweak 的三个种子分别有自己的配对均值与区间；没有 pooled seed mean、seed-as-mother 计数或 pooled seed 标准误。
10. `core_threshold_sensitivity.png`：core5 的 τ=.03/.1/.3 已存 Δfinal 描述均值；不改变训练或选择 τ。敏感性空间秩不同，不当作同一空间中简单重标度。

所有标题明确 **NEW INDEPENDENT TEST、每个 noise 计划108母体、clean+noisy、公共 τ=.1 主口径**。图片上的实际 R/F 是 raw/final 可用数，标题的108永远是计划数，不自动代表完成数。没有把辅助矩阵/场误差代替真实物理误差；单空间和 Twofold 主图只读 quantity=`Maxwell_current`，忽略 material_complex_Green_embedding。

误差单位为离散电流平方范数（576格×12源），不是能量份额、物理功率或连续体积分。不同输出场/残差列使用其已存平方范数，均明确参考系。log/symmetric-log 仅为显示，正负和精确零保留，不设科学 PASS 线。

## Twofold 与统计口径

τ=.1 的 genuine numerical intersection ranks 为 `[0,0,0]/[0,20,509]`，余项47单列。顺序项为 ‖P_O P_I e‖² 与 ‖P_I P_O e‖²，不称相互正交的 joint intersection cells，也不强制归一化为100%。没有把 TSOM projected span 替换成真实交集。

导出的 bootstrap 只覆盖 primary CORE 单空间；**Twofold、余项和 τ=.03/.3 没有现成区间**。这些面板只显示已存配对均值，绝不编造 CI。若某个 primary CI 缺失，脚本可显示 grouped delta 均值，但明确标 interval NOT_RECORDED 并在 SUMMARY 中记录来源。所有 CI 始终限于同一 run、seed、noise 内的有效配对母体；失败母体数量继续显示，不当作成功点补零。

## Auxiliary→realization 匹配边界

`current_realization_gap.csv` 由 exporter 的 `kind=NN` 条件写入，脚本又仅限固定 NN registry，不把 CIE/current-baseline 的 raw incumbent 当作某个网络辅助状态的同一读出。Voxel 无 auxiliary head，标 N/A，不补零。

连接辅助/raw均值前，使用已有 gap 行、raw数据行的 obs_id 集合及两种能量组的 cases 确认同口径。**这是关联表的配对范围处理，不是物理恒等式检查。** 如果不能由这些表确认同一集合，则该方法的 chain 均值不连接，图上列出 unavailable/unmatched，SUMMARY 保留原始行数与匹配状态。如果 raw/aux 已匹配但 final 少了失败母体，保留其 raw/aux 链，final点留空，不拿另一个集合的均值续接。

当前导出表没有逐样本所有 projected auxiliary/raw/final energy，因此在集合不等时不能只靠 grouped sums 重构任意交集的均值。脚本不会补算或猜测；对应缺口留给父线程解释。

现成 realization gap 定义来自 exporter：state defect 使用 `q=E0+G j_aux` 的一致场；current gap 是实际 `J(raw)−j_aux`，internal gap 是 `G` 作用于该差的平方范数，并不直接等同于网络保存的任意 q_aux 头误差。receiver-update 项参照共同严格 GN，receiver-error/misfit 项分别参照 GT 输出和观测 Y；它们不能互换。

## 紧凑交付与失败保留

`figures/test_energy/method_counts.csv` 按13个实际 run×两 noise 列出计划108、raw/final覆盖缺口、可用数、辅助行、实际匹配行、NOT_RUN原因与无法确定 noise 的事件。`SUMMARY.json` 保存同一计数、六个主空间的真实误差/配对变化来源、未匹配链、缺失/非有限字段计数、空间元数据、图片列表和解释边界。源 CSV 中存在非有限或缺失数据的组不会被删去后重算成更好的均值；该数值标 N/R，coverage 与事件仍保留。

这些 availability/status 标签不是新验收门或科学 PASS/FAIL。脚本不读取分数来挑好结果，不替换失败，不平均种子制造样本量，不做公式验证、重复求解或新增测试。若真实数据揭示绘图布局问题，父线程在实际封存后运行并检查；本次仅交付静态实现与命令，尚不能宣称图已完成或目视通过。
