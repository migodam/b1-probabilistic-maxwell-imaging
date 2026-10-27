# KKT闭合 §6 与D冻结epoch的只读审查

新增物理RHS、模型推断、训练、D真值读取均0。仅复核源码、冻结/启动记录、已保存且已暴露的CDEV在线contexts。

§6的正号与尺度一致：目标为Phi+||z||²/2，r_w=white(y−F)，所以场梯度是z−J_w,zᵀr_w。无约束Gaussian场驻点给g=+z；属性box或Probit变换不改变这些场坐标。D的J用σ_new白化，而原风险feature residual用σ_train白化，因此g_feature=(σ_new/σ_train)g_new，驻点闭合应为r·z。代码把noise_ratio乘在z上，再作单位球投影，顺序正确。

单位球投影P(v)=v/max(||v||,1)是Euclidean度量投影，故对任意g,z都有||P(g)−P(z)||≤min(||g−z||,2)。固定当前q时，材料raw反投影对投影向量是线性映射L=δ diag(H′/τ)A，因此RMS界为||L||/24乘上前式。不需要真实驻点、df非零或属性内部条件。复材料按实Hilbert范数；δ的复乘只贡献|δ|。这只是raw材料特征界，不是归一化image特征、神经head输出或risk界。

噪声改变时需将公式明确写成e_feature=r·z−g_feature，比较B(q,g_feature)与B(q,r·z)。不能一边用训练feature单位的g、一边仍用未乘r的z。代码已正确乘r；报告“同一噪声约定”可进一步补这句以避免误读。

16个已有CDEV full/24+64在线context，分别作r=.5/1/2的固定点代数重尺度，共48次raw界检查全部通过，最小slack5.59e−5。它们不是新的噪声观测、不是D泛化证据、不是浮点定向舍入证书。原CDEV q-vs-actual风险区间跨零须继续保留；成本NI通过不推出q显著优于真实梯度。

D冻结实际六arm标签、Low24/64、global8 Bonferroni与alpha=.00625一致。KKT配置时间1793103226.975902、selected-weights记录1793103661.29012均早于Dfreeze1793104030.930937；DLOW/DHIGH生成与首anchor启动均在其后。冻结中的“此前D真值不存在”是当时root记录；本次只读审核可验证作业时间顺序，不伪称重新证明了过去的整个文件系统状态。初始冻结审查时source/checkpoint/development identities均一致；后续serializer真bug修复需另保留源版本与epoch，不覆盖这份初始核验事实。

JSON详情：`audit/KKT_EPOCH_PROVENANCE_QA.json`。这不替代主Codex最终科学判断。
