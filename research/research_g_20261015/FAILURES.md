# 失败与限制

54/72个新endpoint不相容，仍保留raw、corrected、full预测及费用；这些是成像失败，不是未执行solver。QP24/144未达优化器收敛flag，完整QP状态在COST_LEDGER，另索引FAILURE_LEDGER；full候选不能凭优化器flag安装。

E0初次无阻尼Newton数值饱和，完整初次输出保留，仅做backtracking数值修复。独立offline reader第一次将phase-relative路径当repo-relative，失败回执保留并计费，修正reader后只重跑零物理评分，物理搜索和online数组未重跑/改写。

唯一结构witness dS略高于.08，门槛加倍不计结构，且不是新孔。p250077已有真实洞仍未覆盖。两轮paired相同world不是独立复现，三parent和三步bandit不支持泛化/策略优越性。

所有方法Top4/8短缺，MAE最优仍base。G归一化coverage较宽不能替代原RMS目标。只有合法base的三个weak条件被纳入G，新任务不撤销原其它五条件失败。

G实验冻结时未接入T-guided，不能声称完成五臂比较；发布时T已公开交接，记录独立保留。有限Haar/局部mask并非全材料域完备支持；无新3D/GPU/NN/RL/A17消融，没有完整posterior概率或全空间coverage。
