# RNN 与 TEST 协议独立静态审查

状态：**SCOPED_STATIC_REVIEW_COMPLETE；发现的2项均已由父线程修复并静态确认；未发现剩余实质实现问题。** 本次仅阅读指定理论、代码与公开协议身份文件。没有运行实验、反馈层、有限差分、公式数值检查、训练、上传或任何绘图；未读取 truth 数组、TEST score 或 TEST 私有数组。父线程保留最终理论有效性、原创性和实验解释的裁决。

## 实质发现及处置

1. **P2，理论等式缺少正交前提，已修复。** 初版 A17 dressing 把 `(I-R_U L)C` 直接写成 `C+R_U D_m G C`，但没有声明 `U* C=0`，一般还需 `-R_U C` 项。当前 `theory/A17_RECURRENT_FEEDBACK.md:15–21` 已显式声明正交新增模式，给出 `R_U C=0` 的理由并保留一般展开。这是理论叙述前提修正，**不是当前 RNN 代码的 Schur 实现错误**；当前代码本来就不求 A17 retained/Schur 块。

2. **P1，VAL 分数可误进入 TEST 图册并打开 TEST truth，已修复。** 初版 `experiments/plot_results.py:5–6` 接受任意 `--scores`，而 `src/plots.py` 总是选择预声明 TEST 图册；合法的 `--scores results/VAL_SCORES.json` 因而会进入 TEST truth 读取路径。当前 `src/plots.py:30–35` 在读取身份清单和 truth 之前要求 `METHOD_FREEZE`、非空且全部为 TEST 的 scored rows，VAL 选择图继续走独立的 score-only 入口。身份清单读取现在位于 `:41–43`，truth 读取位于 `:50–52`，都在该检查之后。父线程报告此前未将 VAL 分数传入这个 TEST 图册；本审查也未执行该调用，**不将静态入口缺口写成已发生泄漏**。

## 指定实现与理论的对应范围

- **固定材料内反馈与上游梯度不同。** `src/recurrent_feedback.py:172–189` 在 NN 电流增量后读出材料，`fixed_m=m` 后在内反馈循环保持同一个材料，最后再读出。`feedback_step` 的 `L*e=e-G*(m e)`、联合数据梯度和 `Lg=g-mGg` 位于 `:72–94`。`fixed_m` 没有被 detach：这正确地区分“本次显式搜索方向对电流取固定材料的偏导”和“训练对整个展开映射的上游材料路径反传”。没有发现将二者混淆的 stop-gradient。
- **hinge 是 majorizer，不是其精确线搜索。** `:74–87` 使用 `theta=max(0,1-delta/s)` 只缩放数据梯度；`:88–97` 的曲率分母保留完整 `||Sg||²/a`，没有再乘 theta。`:103–108` 分别记录 hinge/二次目标。`feedback_hinge` 的 rho=1 与公开 clean/noisy 门在 `:140–155` 固定；`src/predict_recurrent.py:84–94` 从公开 anchor 的 noise 显式传入。代码和理论对这两种步长的名称、适用范围相符。
- **共享实材料读出、固定归一化。** `recurrent_feedback.py:34–42` 对全部源共同作实 row-LS 并裁剪到 [.02,20]；`:68–71,163` 的观测与锚点电流归一化不随辅助电流扩大。反馈后的 field/residual 用同一步显式更新 `:98–101`，反馈结束后再做 row-LS `:188–189`。理论的有限剖面下降结论仍受其声明假设和精度范围约束，不能扩张到未经约束的 NN 提案。
- **没有把 A17 exact solve 当成 RNN。** 当前理论 `A17_RECURRENT_FEEDBACK.md:7,25–34` 区分矢量 DDA 原实验、当前 scalar-TM 约化代数及反馈架构启发。RNN 层只调用固定 G/G*/S/S* 乘法和 row-LS，没有 L inverse、LU 或 Schur 块求解；`recurrent_feedback.py:201–207` 明确记录 LU/RHS=0 及前向乘法成本、训练反传成本未估计。最终原始 Maxwell corrector 是另一条付费路径：`src/evaluate.py:1–13` 复用既有完整 VIE 评估器，不能将其科学验收归功于“RNN 已无逆地解完 Maxwell”。
- **两步训练与16步诊断分离。** `RecurrentFeedback` 在 `:144–165` 固定训练/正常预测为两次共享权重；`predict_recurrent.py:47,94,125–132` 正常导出仍标两次并要求外部真实 VIE 验收。只有 `experiments/recurrent_iterations.py:23–30` 显式改成预声明诊断深度。公开诊断合同 `audit/RECURRENT_ITERATION_DIAGNOSTIC.json:2–6,63–64` 是6个 VAL 观测、16步、至多96次离线物理解；`:35–52,57–61` 记录辅助/真实场缺口和公开残差，并明确不是不动点或 GT 收敛证明。未运行或检查该数值结果。

## TEST 分离、覆盖与图像身份

`audit/TEST_RESOURCE_PARTITION.json` 只包含公开母体/观测身份与全局顺序：local36母体/72观测、remote72母体/144观测；两端母体及观测无交集，总计108/216。这里只解析身份集合，不读取科学数值。`run_test_partition.py:24–53` 按完整公开 online manifest 核对母体配对、固定分配、FINAL_RUNS 与 METHOD_FREEZE；`:55–63,143–183` 保留平台/worker限制、既有baseline前置条件、有限FIFO与实际缺失案例。`run_evaluation.py:31–37,43–55,72–83` 仅处理所分配 TEST IDs、按原全局序号轮换方法、保存主机出处和中断状态。这个资源分区不构成算法速度或成像优越性的证明。

`src/score.py:25–30` 在实际 TEST seal 的数值评分前检查冻结；`:34–42` 保留 NOT_RUN；质量成功、预算内成功及总55秒账在 `:60–69` 分开。共同严格 GN 与独立冷启动 CIE 的参考身份在 `:90–111` 分开，避免把 CIE 的 uniform 起点误充共同 GN。

图册使用预声明 `audit/FIGURE_SELECTION.json` 的固定母体，实际 row identity 来自对应 TEST manifest；主图 `src/plots.py:68–82` 标注未训练图像/母体、类别、字体/规则、幅度和噪声。完整图册 `:83–100` 保留各方法与失败/缺失；条件于合格 top1 的数据/成本散点 `:101–112` 明示选择计数，不把它冒充全计划成功率。未读取或审美检查生成图，不能将静态身份路径审查写成图像已验收。

## 停止点

限定范围审查已完成。两个修复均仅作源文本确认，未追加实验或数值验证。没有凭本次审查确认训练完成度、RNN 收缩、16步稳定性、TEST成功率、A17数值证书迁移或新颖性。
