# R5-B 独立确认：冻结统计分析协议

状态：`FROZEN_BEFORE_FRESH_GT_ACCESS`。本协议只规定如何从逐 parent 损失得出统计量；没有读取新 truth、生成观测、运行物理求解或训练模型。对应机器可读文件为 `configs/STATISTICAL_PROTOCOL.json`。父线程负责最终科学结论，统计程序只输出证据与门槛状态。

## 1. 要估计什么

独立确认包括 A、B 两个互不重叠的 batch，各 256 个独立 parent，总计 512 个。五个冻结方法为 physics_probit、strong_cnn、k32、deterministic_matched、distilled_w128。每个学习方法使用已经冻结的三个训练 seed。物理 baseline 没有训练 seed，其逐 parent 损失在三个 slot 中重复，且程序检查完全一致。

先计算每个 parent 在三个 seed 下的损失，再对这三个损失作算术平均，最后对 parent 求平均。随机化/抽样单位始终是 parent；三个 seed 不当成三份独立测试样本，也不 bootstrap seed。因此区间描述这三个冻结模型的平均风险对新 parent 的不确定性，不外推为任意重新训练 seed 的方差。

主损失是完整 material 网格的 raw MSE，即逐格 `abs(prediction-truth)^2` 的平均；clipped 版本是预先规定的次要敏感性结果。不能因为 clipped 更好而替换 raw 主结论。预测先作 ensemble 再算 MSE 与先算各 seed MSE 再平均不是同一估计量；ensemble 只能在独立字段中描述，不进入主检验。

## 2. 四项固定比较与符号

令 S、K、D、C 分别为 strongCNN、K32、deterministic_matched、distilled_w128 的 parent seed-mean 损失，数值越低越好。

| 比较 | 逐 parent 效应 | 正值意义 | 冻结检验 |
|---|---|---|---|
| deterministic vs K32 | K−D | deterministic 更好 | 双侧，零效应 |
| distilled vs K32 | K−C | distilled 更好 | 双侧，零效应 |
| distilled vs strongCNN | S−C | distilled 更好 | 双侧，零效应 |
| 保留至少 80% 增量 | (S−C)−0.8(S−K) | distilled 保留 K32 相对 strongCNN 的至少 80% 增量 | 单侧，原假设效应≤0 |

前三项双侧区间完全高于零才支持命名候选更好；完全低于零明确报告候选更差。区间跨零为未解决，不能称等效。最后一项是固定的非劣效线性 contrast，而不是零比值检验。它避免弱分母造成不稳定推断。确认 retention 还要求 K32 增量分母 `S−K` 的点估计为正，且相同调整水平的双侧区间下界大于零。该附加前提只会收紧 retention 判定，不创造第五项可选择的确认性结论。

## 3. Parent 配对 bootstrap、置信水平与多重比较

固定使用至少 20,000 次 bootstrap，默认 RNG seed 为 2026102705。合并主分析在 A 内有放回抽取 256 个 parent、B 内独立抽取 256 个，再合并；每次保持两个 batch 大小不变。一次抽到的 parent 索引同时用于所有方法、三个 seed、raw/clipped 指标以及 retention 的分子/分母，保留配对相关性。不能分别 resample 方法或 retention 分子/分母。

唯一确认性 family 是 pooled TEST 的 raw seed-mean 四项比较。Bonferroni 控制目标 family-wise α=0.05，每项 α=0.0125。前三项使用双侧 98.75% percentile bootstrap CI，即 0.00625、0.99375 分位；retention 使用单侧 98.75% 下界，即 0.0125 分位。所有分位采用 NumPy `linear` 定义。四项使用上述不同侧数也仍然满足 Bonferroni 的每项 α 分配。有限样本 percentile bootstrap 是近似推断，不能宣称严格的有限样本覆盖率保证。

同时输出 retention contrast 的双侧 98.75% 区间方便阅读；它不替代预先规定的单侧判定。单个方法均值的同水平区间仅为描述，不能从它们是否重叠推导配对比较结论。

A、B 各自使用相同四项公式与相同调整置信水平，作为描述性 replication；clipped、ensemble、逐 seed 也都不是额外确认性 family。因为没有从 A、B、pooled 三者中挑一个显著结果，所以 pooled 四项不再按 12 项调整。若后续要把 A/B 各自结果提升为独立确认性主张，必须另行设计并在访问相应新数据前冻结 family；本次不允许事后升级。

## 4. Retention 比值何时可解释

比值为 `R = mean(S−C)/mean(S−K)`，两个均值使用同一组 parent。只有分母点估计为正、调整双侧 CI 下界也为正时，才输出比值点估计与 bootstrap 区间；否则标为 `NOT_ESTIMABLE_UNSTABLE_OR_NONPOSITIVE_DENOMINATOR`。比值推断不能救回未通过线性 contrast 门槛的结论。

通过稳定性前提后仍保留全部非零分母 bootstrap draws，包括尾部少量负分母 draws，不进行条件筛选或丢弃；同时报告非正分母 draws 数。若存在精确零分母 draw，比值区间不可计算，标为 `NOT_ESTIMABLE_EXACT_ZERO_BOOTSTRAP_DENOMINATOR`，不删样本、不截断比值。线性 contrast 仍正常报告。

## 5. p 值的解释与决策优先级

程序附带 null-centered 配对 parent bootstrap 近似 p 值。若观察效应为 t、bootstrap 效应为 t*，双侧使用 `(1 + count(|t*−t|≥|t|))/(B+1)`；单侧 greater 使用 `(1 + count(t*−t≥t))/(B+1)`，再报告 `min(1,4p)`。加一校正避免有限 Monte Carlo 下报告 p=0；B=20,000 时未调整 p 的最小分辨率为 1/20,001。

这不是 permutation/randomization 精确检验，也不是后验概率。centered-bootstrap p 与 percentile CI 不严格互为反演，在边界可能不同；冻结判定以规定的 percentile CI/下界为准，p 值只是伴随证据。不得择优选用 p 或 CI。

## 6. Batch 一致性与缺失/失败

报告 pooled 与 A/B 的全部点估计、区间、方向，并逐比较标记两个 batch 是否同为正值、是否方向相反。同向但宽区间不能叫两个 batch 都显著；一个 batch 显著不能替代 pooled 主检验；A/B 反向时必须限定 pooled 解释。A/B 只描述重复表现，不增加新的确认性结论。

缺方法、少 seed、非有限或负损失、重复 parent（包括跨 batch 重复）、确认 TEST 不是 256+256、物理 baseline seed slot 不一致、或 bootstrap 少于 20,000 均报错。不得删除失败 parent、重新抽样补齐、调整 seed、改 loss 或扩大样本挽救结果。统计文件必须包含全部指定 parent 的结果；前置推理失败由总运行协议报告并阻止未限定确认结论。

## 7. 逐 parent 损失接口与运行

输入 NPZ 不含 object/pickle，且只读取 allowlist：

| 字段 | 形状 | 内容 |
|---|---|---|
| parent_ids | [N] | 全局唯一 parent ID |
| batch_ids | [N] | 确认样本为 A/B |
| seed_ids | [3] | 冻结的三个训练 seed slot 标识 |
| methods | [5] | 上述五个 canonical 方法，输入顺序可不同 |
| raw_losses / clipped_losses | [N,3,5] | 逐 parent、逐 seed、逐方法完整网格损失 |
| split | [N] | 默认选 TEST；不混入 TRAIN/VAL |
| phase | scalar string | 确认阶段 independent_confirmation |
| ensemble_raw_losses / ensemble_clipped_losses | [N,5]，可选 | 另算的 ensemble 损失，仅描述 |

默认只选择 `split=TEST`、`phase=independent_confirmation`。其他 phase/split 需显式指定，并自动标为 `DESCRIPTIVE_OUTSIDE_PRIMARY_PHASE_SPLIT`，不具有确认性资格。标识、损失在进入本程序前由授权评估阶段生成；本程序不打开 truth/prediction 图像，不进行物理求解或训练。

```sh
python research/p3l_r5_b_independent_confirmation/src/statistics_confirm.py LOSSES.npz \
  --output research/p3l_r5_b_independent_confirmation/raw/STATISTICAL_RESULTS.json
python -m unittest discover -s research/p3l_r5_b_independent_confirmation/src -p 'test_statistics_confirm.py'
```

实现测试全部使用人工合成损失，验证配对抵消、固定 batch 分层、seed-mean/ensemble 分离、分母不稳定关闭、schema/phase/split 及样本数门槛。测试通过只证明统计实现符合本协议，不能证明任何新 material 成像结果。
