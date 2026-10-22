# P3-L 执行协议（2026-10-22冻结）

正式配置为configs/PROTOCOL.json。原P1 prior、noise、acquisition全部固定；本轮是用户独立授权12M complex RHS、48小时受控CPU、24 GPU-hours的新阶段，理论包的false执行标记不构成本轮授权。

先运行32个新development truth与强化物理反演；从仅依赖y的公开anchor特征中确定16个深reference。四条过分散链，1024 warmup、4096 production，所有拒绝保留。MALA对u=Phi(eta)使用精确uniform-prior transport，canonical q上的pCN/ESS作独立原理核验。未达到Rhat<1.01、材料函数ESS>=1000、mean-RMS MCSE<=.001且跨机制一致的reference标UNQUALIFIED，不能成为teacher。诊断通过也不证明全局混合。

所有成像比较以材料MSE/MAE为主，保留prior mean、充分优化TRF、MAP与全部失败。LOCAL-C采用固定Gaussian距离12及材料RMS .16，无noise-ball硬截断；FULL、LOCAL-C、SLICE、PROFILE严格分开。强条件Gaussian材料矩只对声明Gaussian参考准确。属性box完整保留。

最多512 joint-train、128 validation；若开发风险或可信UQ支持扩大，方法冻结后才生成并评价128 locked新parent。旧locked数据不访问。网络只用合法joint NLL或合格teacher；未合格teacher不蒸馏。几何/RL是计算机制对照，不改变目标。任何更新另存修订记录，不能改原P1来求正结果。

新增求解先由本轮独立ledger预约；失败、拒绝、LU、forward8 RHS、source-adjoint gradient8 RHS、receiver-J32 RHS均计费。每个新endpoint立即保存q、材料、prediction和backward residual；truth/anchor额外保存current/field/LU，链保留所有重复状态及RNG checkpoint。大规模采样的全部LU不持久保存，避免将可重建缓存当原始数据。48h为受控执行CPU wall，不含人工/agent思考及发布等待。

仅在可信概率与独立统计允许时主张>=5% MSE改善且parent paired95%区间排除0。没有改善也形成完整负结果；不把Best-of-K、候选数、短链或residual softmax用作posterior成绩。停止明确失败方法，预算是上限而非耗尽目标。
