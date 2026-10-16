# 同proxy预算比较

| 方法 | event MAE | 材料CDF grid误差 | recall* | 4096次平均arm秒 |
|---|---:|---:|---:|---:|
| local | 0.012039 | 0.071142 | 1.0000 | 0.4998 |
| geometry | 0.026158 | 0.178926 | 0.9479 | 0.5615 |
| geometry_nuisance | 0.011720 | 0.075228 | 1.0000 | 0.5531 |
| global | 0.004364 | 0.028285 | 1.0000 | 0.5079 |
| hybrid | 0.007243 | 0.042557 | 1.0000 | 0.5231 |
| tempered | 0.016145 | 0.103910 | 0.9896 | 0.0127 |
| prior | 0.011206 | 0.074878 | 1.0000 | 0.4307 |

*固定参考事件质量≥.05、经验估计≥.005；这是本轮false-exclusion诊断阈值，不同于旧P2的.01 discovery。CDF仅保存grid最大差。SMC为64个加权粒子，MCMC为保留拒绝重复的链平均，ESS不可直接互换。


每个case/seed/arm共4096次likelihood（含初始化，geometry含8次差分）。1,024 checkpoint也在相同总调用位置。完整本轮proxy次数589824，含两次坏geometry与两次repair；不把失败隐藏为免费。

temperedSMC批量64比逐步MH便宜，不能只比较总sample数；其eventMAE约.0161更差。C1全局MH约.00436，hybrid约.00724；混合并未自动改善。warm arm秒数不含共享cold C1导入/模型加载，后者在每case PROPOSAL.json，整轮wall包含。当前所有MH臂都预生成q draws，local/prior存在额外非必要setup；这些实测时间属于本实现，不宣称优化后部署优势。训练与bank amortized费用追溯原P2，不能算从零免费。

E2共96states/896complexRHS，96次LU；原已支付cache矩阵算术另列。成功job wall1.6247s，首失败job纯wall未测，保守elapsed上界95.627s含主线程间隔，总physics上界97.252s≤3600。所有失败/重复均保留，剩余状态0。

E2 continuation有5个predictor+5个corrector端点，而每个MH arm只有5次proposal，不能将continuation较大位移说成同成本胜出。Jacobian的32个receiver-adjoint RHS/anchor由本次共享setup实际支付；普通local/prior若单独部署无需该几何setup。这里不是严格同wall的E2算法优越性试验。
