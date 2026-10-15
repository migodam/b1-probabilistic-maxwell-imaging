# 候选覆盖压力与保留歧义

本轮不把有限候选集合成功当成全材料空间覆盖。full625池的Best-of-K MAE=.0622242099，Coverage@K(MAE≤.05)=.40625；只有新增材料才可改变这两个值。

## 已执行删除与模型错配

同32parent、每parent两噪声，明确exposed followup；按删除后的池重新排名。full625、central_deleted125、annular_deleted125、combined_deleted225。固定2.4/3.2、raw、robust、sentinel与不测量都保留。

robust一block平均Selected-MAE：full约.103714，central_deleted约.12301，annular_deleted约.28291。对比完整池，推荐动作分别改变33/64、13/64、9/64个noise episodes。它们不是64独立对象。

central缺失29/30、combined缺失15/16相关episodes没有residual alarm。sentinel阈值来自χ²启发式，经历有限池选择与bounded gain，不具有声称的nominal coverage；无报警不保证世界齐全。固定3.2 sentinel不是新ranking胜利。

另有独立instrument drift负控制：future gain = old gain*(1+.03*(k−1.6))，原field与noise不变，拟合器故意仍假设shared gain。该表测试模型错配，不和in-model收益混合，不冒称材料色散实验。

未知结构family不在此有限径向prior支持内，真正的unseen nonconcentric family未运行。当前没有系统的prior-sensitivity新campaign，只有保留原prior和有限候选删除/扰动控制；不能据此称prior robust。所有原始删除、ranking、parent配对与drift记录在results/E2_COVERAGE_STRESS.json。

## 实际交付多个剩余世界

results/REMAINING_WORLDS.npz/.json保存整个625材料及固定3.2、robust、risk三策略的32parent×2noise末历史；每世界shared gain、χ²、材料参数、posthoc events/families和best+9兼容mask均保留。768个已执行trajectory step重放的χ²差为0；没有新噪声或新物理。

兼容K范围2–27、1–35、3–29；192条episode-policy中160条仍跨family。mask仅是有限目录的数据兼容启发式，不是posterior/confidence region。概率none，outside mass=null，整个未筛目录也保留，方便后续独立探索。

E3三个连续prior起点正确补入了一些接近truth的fullvalidated材料，使四例Best-of-K改善；仍不能保证发现全部结构。MAE改善与结构正确分开报告，不能把新ε变体个数称为新结构family。
