# 额外 direct-prior 强基线：准备完成，默认 NOT_RUN

这项检查只在根线程认为 CDEV64 的 direct Probit64/128 已竞争或支配当前候选时考虑。它从相同的 prior 点启动原始 Probit 优化，跳过 Native warm optimizer；用它检查阶段分配主张是否遗漏更简单的方案。**当前只有代码与模板，没有执行授权、物理查询、预测或真值评分。根线程可选择 NOT_RUN。**

## 冻结时机与样本

根线程必须先封存主 D 全 512 parent 的完整预测，再依据 CDEV 结果选择唯一 `cap=64/128` 与唯一 `mode=actual/q`，在任何 D 真值解锁或评分前创建本目录 `SUPPLEMENT_FREEZE.json`。模板的 `TEMPLATE_NOT_AUTHORIZATION` 永远不能执行。

补充只用 `parents_DLOW.json` 与 `parents_DHIGH.json` 各自**按原名单顺序的前 64 个** parent，共 128。该选择由 metadata index 决定，不看 GT、不筛样本、不生成新真值。它不是增加主 D arm，不能改变原 D 源码、比较族、主预测或 full baseline。

## 根线程需要填写的冻结

- `status=ROOT_SUPPLEMENT_FROZEN_PRE_GT`、`execution_authorized=true`，并明确 CDEV-only 选择与理由。
- 唯一 cap/mode，以及两条件 first64 的完整 parent ID 列表和原 manifest BLAKE2b identity。
- 原主 `D_FREEZE.json` 与 `ALL_512_PREDICTIONS_SEALED.json` identity；完成的 CDEV direct-prior `RESULTS.json` identity。
- `source_identities` 必须包含本 standalone `supplementary_prior_baseline.py` 的实际 identity。已冻结的 prior-only Core/Original 数学依赖还会按 `DEVELOPMENT_SOURCE_FREEZE.json` 逐项检查。
- 保持模板中显式的探索性统计合同。两条件分别作 20,000 次 parent-paired bootstrap，种子损失先平均；不合并噪声条件，不宣称补充 128 已获得原 512 的资格，也不宣称新旧比较族合计仍有未重新分配的全局 FWER 保证。

## 预算与物理保存

控制器只使用根线程中央 `shared.physics.Book('D',...)` 与 `Target(sigma=online.sigma)`。prior 点先支付 8 RHS，之后调用已冻结原始 Probit-only Core。64 cap 最坏 2,600 RHS/parent；128 cap 最坏 5,160 RHS/parent，128 parent 最坏 **660,480 RHS**。必须同时满足 C≤3.2m、D≤4.1m、总计≤8m 及当前实际授权；不能通过改标签、退款或降级父样本数救预算。

最终 q、chi、current、field、prediction、J、receiver transfer、残差与优化路径由原始 lossless mode/context 输出保留；中央 job 记录全部实际物理资源。没有存储额外 LU。每条新路径独立计时并包括 context persistence/cleanup；主 full/selected 的已有物理路径不重跑。中断、已有部分目录或未关闭 job 都需要根线程显式审计，不能自动重放。

## 将来执行顺序

下列命令目前均未执行。

```bash
../../.venv_nn/bin/python research/p3l_r6_task_optimal_maxwell/track_c/supplementary_prior_baseline.py anchors
../../.venv_nn/bin/python research/p3l_r6_task_optimal_maxwell/track_c/supplementary_prior_baseline.py predict
```

预测使用原始 frozen deterministic 三种子，q 模式按已知噪声比例 `.5/2` 调用现有 closure Action。所有 128 个 supplement 预测和上下文先全局 BLAKE2b 封存。只有补充 seal 早于主 `GT_UNLOCK.json`、且根线程已正确解锁主全局 seal 后，根线程才可调用 `score`。脚本在第一次真值打开前同时验证补充和原主全部预测封存；评分只读固定 first64×2 的真实材料。

输出分别比较 supplement 与主 KKT、原 full：材料 MSE/MAE、种子平均损失、单独 ensemble、1% NI 线性 contrast 与 20% wall contrast、配对区间、逐 parent NPZ。结果状态固定为 **`SUPPLEMENTARY_PRE_GT_EXPLORATORY`**；不触发自动 GO/TAP，不修改主 D 结果。

最终是否执行补充、如何解释强普通基线，以及是否另安排前瞻性数据，由根线程在任何相应 GT 解锁前决定。
