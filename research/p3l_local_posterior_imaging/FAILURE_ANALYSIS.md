# 失败说明与不能排除的解释

本文件区分已经观察到的失败与尚未证明的原因。最终数值、费用和结束状态以 METRICS、参考资格审计及主线程裁决为准。

## 不能把混合失败说成物理没有不确定性

MALA、变化 decoder mobility、pCN、独立 HMC、非线性 spine-HMC 与原生 Elliptical Slice 已运行，但目前没有合格 FULL reference。较稳定的 likelihood 轨迹没有使材料函数自动混合；有些材料函数在有限链中甚至为常数，保留 UNASSESSED。128 粒子 SMC 到达 beta=1 后只剩一个祖先，同样不能当独立 posterior teacher。

这是 sampling/reference-limited 的证据，不是证明原 P1 只有一个物理世界，也不能证明其 posterior mean 与 MAP 相同。选中的深参考原计划有 16 个；预算优先用于已固定条件的失败归因和独立核验，其余未执行的条件不能冒充完成。

## Gaussian 参考失败不是 quadrature 精度不够

固定属性与联合属性的 Gaussian 材料均值均在暴露 development 数据上恶化。后者 GH9 到 GH17 的平均材料变化已很小，Full 重要性权重却仍退化。40 组强条件参考中 38 组 ESS<2；局部梯度修正很小的情形也存在这种失败。

目前的归因是局部 Gaussian 参考在有限范围上没有正确表达原 posterior。饱和 decoder 重新激活、属性耦合、非 GN 曲率与优化残留均可能参与；没有把它们统称成 Maxwell 高阶 feedback。新代理力配对试验正是为了隔离保留 decoder 的作用，结果必须按真实采样与成本判断。

## 训练执行了，但没有验证集学习增益

五个归一化神经密度各三个 seed，全部训练到 validation patience，并均选择 update0。实际 minibatch NLL 下降而验证 NLL 上升；独立源码/检查点审查没有找到已确认的缩放或恢复错误。它支持当前数据规模、特征与训练设置的泛化失败，不能证明任何神经 posterior 都无效。

简单正则化 Gaussian 的 validation NLL 好于这些初始化，但四个 exposed parent 的完整 MH 对照中，条件与无条件提案均 0/256 接受。正规化和较好弱坐标 NLL 没有修复 strong-reference 的真实 likelihood 失配。没有使用劣质链做 distillation，也没有训练 RL 来优化未经验证的概率奖励。

## 风险比较必须保留优化与坐标限制

native MAP 与 probit-coordinate mode 是不同点估计器。后者变好不等于材料 MMSE 变好；额外优化、属性 interior 正则化与坐标依赖都需要保留。development 的 32 个 anchor 中 27 个达到冻结 KKT 标准；128 个新对象的 native96 有 79 个达到该标准，不能宣称所有点已得到全局最优。

冻结评估不删除其余对象，也不根据 GT 延长某个反演。独立 parent 风险可以评价这些实际估计器，但不足以声称完全优化 MAP 与真正 Bayes action 的理论差异。

## 内存中断与费用

内存问题来自未受限的 transfer 对旧 AnchorState/LU 的强引用，以及多个 context worker 同时运行。已改为有界缓存与单一物理 worker，并保留所有已付费状态。109,488 RHS 的中断前费用、缺失尾部 wall 的保守上界及原失败回执没有删除。此后恢复并不使之前的实验变成零成本。

## 不能从近似求值对照推出 A17 增益

同 current space 的 A17 与普通 Galerkin 仅差浮点误差。16 个固定提案中 rank16 条件界提前判断 5 个决定，但它依赖浮点模型，未使用定向区间运算；setup、界计算与 fallback 均需计入。目前没有总 wall 优势，TM24 Full 本身很便宜。

## 文献与新颖性

ScholarQA 的两条发现查询都返回 429，已知主文献的定向核验不能代替全面 prior-art 搜索。LIS、Rao–Blackwell、Laplace、非线性传输、代理梯度 HMC 和 goal-ROM 均有成熟先例。本轮没有因为检索缺口就宣布原创，也不以失败回顾包装强 IEEE TAP 算法贡献。
