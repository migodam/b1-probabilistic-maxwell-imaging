# 经典Austria结果独立交叉审查

范围：只读 `CLASSIC_BENCHMARKS.md`、Windows最终与已付前缀紧凑评分、receipt/checkpoint metadata、`RESEARCH_REPORT_ZH.md` §10及 `PAPER_DRAFT.md` 摘要／合同表述。未读GT数组，未求解、重新评分或改动主报告；不涉及 original nine 或 physical refinement。最终科学判断由root负责。

## 仍需收窄的两处表述

1. **经典报告第9行的“残差约.0933→.0929”缺少80步纯数据项的独立metadata支持。** 保存的single CIE检查点给出含prior的Phi：`.004353397559→.004318402`；最终另外保存了纯数据相对残差`.0929265874`。不能把`sqrt(2Phi)`严格当成纯数据残差。建议与当前中文主报告§10一致，直接引用Phi下降。材料NRMSE `.1361716→.1851289` 的恶化已由离线紧凑评分支持；此修正不推翻该现象，也不需要重新求解。

2. **摘要末句“All candidate selection uses the same physical objective”及Methods首段“same box and prior”应明确适用合同。** 若仅指第一阶段匹配组，当前科学含义可成立；若覆盖本文全部扩展，则独立经典NIE/LSIE为R=0、另六臂smooth=TV=1e-4，且扩展之间归一化不同。建议限定为“within each declared matched task/prior group”。这是合同范围澄清，不是发现GT选择或评分作弊。

## 已核对事实

- §10的εr=3.5/4主表与最终16项评分的舍入一致；A+分别为幅度NRMSE`.18450365/.58300030`，source single CIE为`.18512892/.47475257`，Global-Dual为`.18962134/.57516534`。没有发现把contrast-relative L2、幅度NRMSE与逐像素相对epsilon RMS混用。它们分别以`||chi_true||`、`sqrt(N)*(epsilon_r-1)`、逐像素`|epsilon_true|`归一化，不能横向当成同一指标。
- A/A+为**global材料臂**，ordinary current与Global-Dual相同；source CIE/Dual及独立sourcewise论文参考另列。现有主表和§10正确区分，未把global对sourcewise差别全部归因于A材料metric。
- A→A+完整任务CPU：ε3.5为`2108.234375→1591.953125`（减`24.4888%`）；ε4为`4617.390625→2977.953125`（减`35.5057%`）。当前§10已修为“完整任务”而非纯材料内循环，准确。辅助receipt还包含current/readout/验证及不同轨迹，不能据此隔离材料内核的因果速度收益。A+相对source CIE并不更便宜。
- 全部16项corrector均`300 accepted / CAP`，optimizer未报告成功收敛；部分native stage也有MAXITER。没有整项数值失败不等于完全收敛。当前报告已保留该区别。
- 80检查点来自同一次连续corrector，未用GT改选线上停止；ε3.5多个CIE/A臂材料确实退化而原Phi略降。该现象不能唯一归因为噪声，也不证明GT是原目标最优点；约9.96%粗细模型误差仍在合同内。
- 300/600/900 CPU前缀为稀疏已存状态的公开Phi选择；存在最大约723秒覆盖缺口与输出尾部保守上移。报告已明确不是精确等CPU截断，不能据它给出严格等算量优胜或连续效能曲线。
- 两个公开材料条件不是16个独立新几何；不同prior的论文参考不是作者实现的逐值复现。ε3.5的A+与CIE差约`.00063`且IoU略低，ε4 A+更差。现有§10及摘要没有据此宣称A+普遍优胜、新全局凸化或新的独立成功机制。

总体：未发现主表数字或方法标签的实质错误；当前§10已修正完整CPU计时与Phi／纯数据项的两处范围问题。其余必须保留的限制是有限预算CAP、不同prior参考、稀疏成本前缀，以及公开开发对象而非盲确认。以上证据只支持冻结配方下的比较，不能推广为完全收敛算法排序或论文复现胜负。
