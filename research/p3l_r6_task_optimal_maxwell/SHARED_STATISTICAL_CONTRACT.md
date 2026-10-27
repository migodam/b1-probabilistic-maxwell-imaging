# R6 共享统计与物理合同

原P1是固定24×24 scalar-TM/VIE、576单元、256 ROI、8源×32复接收。原始参考为Gaussian64×Uniform2；属性Gaussianization为u=Φ(ξ)，完整X66为标准Gaussian。P1原prior.py、原quintic与原background直接复用，未改Born、noise、source或receiver。

源序实化 `pack(y)=[Re32,Im32]` 逐源连接；f=√2 pack(F)，d=√2 pack(Y)，σ=0.00119211798514519。Phi=||F−Y||²_C/σ²=||f−d||²/(2σ²)，FULL posterior相对Gaussian latent的势为Phi+||X||²/2。Native坐标优化使用Gaussian64正则与Uniform属性box；Probit坐标点优化还含两属性Gaussian penalty，故两个MAP具有参数化区别。

材料平方损失是复材料的实Hilbert损失除以576。三个已训练seed先在每个parent内平均loss，再以独立parent做推断；ensemble图像risk另报。ROI风险分母为256。相区用共同完整物理anchor的场值定义，GT只用于独立离线评分。首次CDEV相区脚本误用truth分区已保留原输出并纠错；主风险/成本/选择逐项完全相同，详见PHASE_CORRECTION_RECEIPT。

输出角色：材料风险估计器、指定条件函数、低维概率积分、合法Maxwell具体世界分别评价。CNN/蒸馏和Gaussian provider均不提供本轮FULL posterior认证；noise-ball不是posterior硬支持。低维A活动X明确截断于[-4,4]、其余X=0，是新子模型，不是FULL marginal。

所有GT生成、在线反演/预测、离线评分分开。新的D条件为σ/2与2σ，材料prior/acquisition仍相同；旧模型是zero-shot、没有新噪声训练。D的512个独立parent清单在生成前冻结，两批均完成预测封存后才解锁GT；同一parent只统计一次。源/接收子采集不在已完成主结果内。

物理CPU1/complex128，实际RHS在LU之前预留。Full8源forward计8；receiver Jacobian额外32；失败、拒绝、line search和回退全部计费。保留原backward residual≤1e−10，它是冻结离散系统的数值检查，不是连续Maxwell认证。原Native96/Probit64优化器及容差不变；新低成本配置另有明确身份，不冒称原完整baseline。

新资源独立于已封账R4/R5/R5B：总8,000,000 complex RHS，受控CPU wall总和72h、MPS训练48 device-h，无云采购。历史缓存的增量成本和从零冷启动成本分开。默认只持久化无损坐标/复响应，最终端点另存实际材料、current/field/J，避免反复存dense LU。
