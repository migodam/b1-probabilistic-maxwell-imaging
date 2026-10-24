# 两批真正独立的新对象：冻结与独立评分

第一批512：seed20261025020，prefixR4L1；第二批256：seed20261025021，prefixR4L2。二者都不是R2已暴露的256+256。128新DEV单独用于方法开发。

FROZEN_LOCKED、FROZEN_CONFIRMATION、FROZEN_DENSITY_BOTH_COHORTS以及risk fresh_freeze_v3在新观测产生前绑定。模型、checkpoint、CFM NFE32、K128、canonical CPU batch64风险推断、失败规则和主要动作均冻结。第一批结果未用于重新选择第二批模型。

所有768 truth forwards、noisy observations、physics-only native96+probit64 contexts完整保存。在线进程仅接收观测和合法context；GT分别放inputs/offline_only。两批45风险与26密度实例均在线完成、物理PPC和pretruth优化诊断封存后，才有FIRST_HELDOUT_GT_SCORE_RELEASE与独立评分PID。

独立审计复核全部45 raw/clipped损失、fallback、源绑定、chronology与4,000parent-bootstrap：失败0，source mismatch0，损失差≤5.55e−17，区间差≤2.15e−16。三个训练seed的loss逐parent平均，不当作2304个独立对象，也不当ensemble图。

Spline主要方法由历史validation的.00712673对Gaussian.00713016的近似平局确定。新DEV曾更偏向Gaussian，但没有据此换主要方法。Spline/Gaussian最终没有显著风险差；不得在两批测试后挑赢家重命名。

主要科学结果两批分别报告；合并768仅作描述。其他45/26方法消融不宣称全族错误率控制。现在两批GT已EXPOSED，未来需要新的独立parent，不把它们再次当locked。

首64物理代表/Full-J优化诊断为预GT封存的固定子集，与全部768主风险比较分开，不从其GT结果重调网络或stopping。
