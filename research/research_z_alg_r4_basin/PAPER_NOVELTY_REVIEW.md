# R4 独立新颖性审查：先校正多个低谷，后按原目标排序

最新复审：已读取最终 V3 新鲜封存结果；GT 改善 `NO`，有限数据目标可达性有局部正证据，原创性仍 `HOLD`。首版与 V3 追加记录均保留，见末节。

状态：`HOLD / FIRST_STATIC_REVIEW / ONLINE_RESULTS_NOT_REVIEWED`。审查日期：2026-10-30。这是独立审查意见，不代替主代理的最终科学裁决。没有运行新实验，没有读取生成器真值，也没有修改理论、求解器或结果。

## 可以保留的主张

目前可成立的描述是：在固定有限维 scalar-TM VIE、全部 576 个独立实材料未知的公共后端上，利用精确方向 resolvent 的谱复用搜索完整合法材料区间；把多个**检出的**、经明确保留规则留下的低谷分别送入局部校正，再按共同原目标选择。它是一个待评估的有限预算求解器组合。

这不是新的 resolvent 理论，不是全局收敛定理，也不是尚已验证的 Maxwell 跨盆地优势。静态代码存在与上述理想机制的偏差，见下文；真实数值收益继续 `HOLD`。

## 贡献与既有机制分解

| 组成 | 独立判断 | 依据与边界 |
|---|---|---|
| 精确有限差、材料 Jacobian、二阶余项 | 标准隐式微分和 resolvent 恒等式 | 乘法顺序正确；公式本身不足以支撑论文的新理论主张 |
| 完整材料 ray 的有理响应与复 Schur 复用 | 标准矩阵有理函数数值结构，且有项目内 R10 直接前身 | 全 576 阶 setup 和每次三角 RHS 必须收费；方向是一维不等于材料更新低秩 |
| 多低谷先校正、后排序 | 区别于只取一个线 argmin；属于多起点/盆地跳跃机制家族 | 新的可评估点在于有理复用如何组织与保留种子，以及收益是否抵偿 setup 和分流校正成本 |
| 较高原始损失候选也先校正 | R3 已有明确概念前身 | 不可称首次提出此操作；须对照同方向、同 evaluator、同 corrector 的 raw-top1 消融 |
| 单材料可辨识但连续下降不可达反例 | 正确的解释性有限 VIE 例子 | 指定几何未实现；它证明逻辑分离，不证明真实 Hankel 576-cell 数据存在同类严格极小点 |
| 横向收缩 | 有条件的标准局部结论 | 不证明进入该邻域、沿全局纤维移动、GT 恢复或后验覆盖 |

普通回溯只缩步，扫描完整合法区间可跨越沿途上坡；普通全局线最小化仍只返回一个原始线最佳点。多低谷分别校正可与这两者不同。另一方面，“远种子 → 局部优化 → 原目标选择”本身没有脱离普通多起点/盆地跳跃框架。R4 如有贡献，须落在结构化种子、复用成本和实际有限预算结果上。

只沿 GN 方向不必然使多候选校正与 ordinary line search 完全相同：多个远端线低谷仍可分别校正。独立方向有助于避免 GN 射线的覆盖缺失；在边界投影把 GN 方向变为零的反例中尤其必要，但不是所有问题上不同机制的普遍必要条件。

## 文献证据边界

遵循 `scholarqa-research` 的证据分层，复用主代理已有检索；本审查没有追加搜索，未把限流误报为阴性检索结果。

`audit/literature_bundle.json` 记录多次 Semantic Scholar 429。唯一明显相关的 Habashy–Abubakar 2004 条目有题名、作者、年份及 DOI `10.2528/PIER03100702`，但没有摘要或全文支撑，当前只能记 `Tier C`。它不支持“该文已提出/未提出本机制”的实质判断。另两项低相关元数据也不建立最近邻覆盖。

主代理已核对的 primary anchors 可供其最终 claim ledger 使用：CSI（van den Berg–Kleinman 1997，[DOI](https://doi.org/10.1088/0266-5611/13/6/013)）、recursive linearization（Borges–Gillman–Greengard 2017，[DOI](https://doi.org/10.1137/16M1093562)，[arXiv](https://arxiv.org/abs/1608.06871)）、既有 R3-DIGITS varpro/direct-prior 来源（[arXiv v2](https://arxiv.org/abs/2403.15292v2)）。本独立审查未重新取得这些全文证据，不能为它们补造 Tier A 片段，也不据此宣布最近邻已排除。递归线性化、多频 continuation、辅助电流/varpro 与原数据 ray 扫描是需要区分的机制；它们不由相似关键词自动等价。

因此目前是 `LITERATURE_COVERAGE_INCOMPLETE`，既不能断言全球原创，也不能把不充分检索当成已证伪原创性。项目内直接前身已足以要求收窄公式原创性表述。

## 实现对论文机制的限制

静态快照中 `src/basin_solver.py` 的 `bridge` 在校正前以 `candidate_rms=0.25` 合并种子。这个阈值会删除材料不同的近邻低谷，与注释“只移除 materially identical”不符；近邻种子也可能位于不同校正吸引域。应由主代理选择：去掉校正前近似去重，或明确将贡献限定为“阈值去重后低谷”，保存被截断记录，并承认不能检验全部分离低谷的理想机制。

默认 ray corrector 是 `correct(..., trf=False)`，即原目标 L-BFGS-B 加物理伴随梯度。full-J 用于 GN 射线构造，TRF 则另有显式路径。故当前不能称“多个低谷各做 full-J 校正”。这不否定多低谷机制，但改变了校正动力学、成本与对照解释。

扫描和排序实际使用 normalized full data loss 加共同 smooth/TV；“原始数据 ranking”若被理解为纯残差则不准确。共同正则目标可合法比较，但报告须同时保留纯数据残差、正则值和总目标。纯数据零残差收缩定理不能直接解释此正则 corrector。

有限网格加 bounded refine 只支持 `DETECTED_VALLEYS`。不得写 all stationary points、global ray argmin 或 all basins。`ray_best_only` 是“检出候选中 raw-top1”消融，不是已认证的全局线最优基线。

## 解除 HOLD 所需的判别结果

先看同后端、同目标、同总预算的 multi-valley 对 raw-top1 消融，再看普通同 corrector 多起点。若只与 TRF 比较，局部优化器改变、种子改变和分配改变混在一起。成功路径应能追踪到一个原始分数较高但校正后改善的候选；如所有低谷回到同一材料簇，不能声称发现多个 basin。

收益要按全部 setup、LU、源/伴随/接收 RHS、Schur、三角 RHS、失败校正和墙钟计费。数据残差、材料 RMSE、IoU 和预先声明的联合成功门分别报告。历史四例已暴露 GT，只能是开发证据；后续冻结后匿名样本的适用范围仍是有限生成器与有限预算。

目前判断：`FORMULA_NOVELTY_NOT_SUPPORTED`；`ALGORITHM_COMBINATION_PLAUSIBLE`；`MECHANISM_BENEFIT_HOLD`；`GLOBAL_NOVELTY_UNRESOLVED`。

## V3 封存结果追加复审（2026-10-30）

追加只读取 `results/fresh_v3/SCORES.json`、`configs/FRESH_CONTRACT.json`、`audit/FRESH_PREDICTIONS_FROZEN.json` 和 `ALGORITHM.md`；没有新增实验、校验或检索。封存记录声明 18/18 预测在首次新鲜 GT 评分前关闭；关闭时间为 03:15:35 UTC，评分时间为 03:19:07 UTC。本审查接受记录所声明的边界，未独立审计生成进程或证明不可变性。字体替换被明示为首个生成失败且未产出 GT 后的生成器修正。

V3 已将校正前近似去重阈值从 .25 改为 `1e-7`，对超过预算的低谷以幅度顺序分散选最多四个，并记录截断。每 seed 同时分配 wall/LU/RHS；ordinary multistart 复用已付费 m1 anchor，新增 4/8/12/16。所有 L-BFGS 路径使用公开仿射缩放变量及链式梯度。因而首版对应问题不应再被称为 V3 的未修缺陷；仍须准确写“预算保留的检出低谷”，而不是“全部低谷”。校正器命名也已明确为完整物理 L-BFGS-B，而非 full-J。

六方法在三个新 digit 2/7/8 上联合成功均 **0/3**。所有数据门也未达 `residual<=.001`。以 online top1 看，uniform_bridge 在 7/8/2 的残差分别为 .015425/.111464/.045760，明显低于 detected-best-only 的 .762235/.780809/.911384；但材料 RMSE 分别为 10.8908/10.5933/11.3274，相比其 anchor 约 5.0503/6.2721/6.3376 全部恶化。oracle best-K RMSE 没有从新增候选得到比保留 anchor 更好的材料。这排除了“只是 GT 好候选被 online 排名选错”的本轮解释：有限 archive 内也未出现对应 GT 增量。

ordinary uniform_multistart 同样找到低得多的数据残差：7/8/2 为 .054634/.137025/.041534；它在 digit2 的残差优于 uniform_bridge。故有“多低谷保留可优于只校正 detected-top1”的有限数据目标证据，但没有跨三个目标一致优于普通重启、没有 GT 改善，也没有证明谱选择具有普遍额外收益。不同 arms 实际资源耗用不同；这是同上限下的冻结策略比较，不是所有算法均耗尽同一计算量的配对证明。

最终可报告的算法贡献是一个付费、可追踪的反例式经验观察：在这三项有限离散合成问题和冻结预算内，更多结构化 seed 的校正可到达原数据目标更低的材料，然而成像指标失败。它不能上升为新 Maxwell 理论、全局可达性、真正不同物理分支或概率解集。最近邻文献不足的限制保持不变。

最新裁决：`GT_IMPROVEMENT_NO`；`JOINT_RECOVERY_FAIL_0_OF_3_ALL_METHODS`；`FINITE_DATA_OBJECTIVE_ACCESSIBILITY_LIMITED_POSITIVE`；`ADVANTAGE_OVER_ORDINARY_MULTISTART_NOT_UNIFORM`；`NEW_THEORY_NOT_SUPPORTED`；`GLOBAL_NOVELTY_HOLD`。
