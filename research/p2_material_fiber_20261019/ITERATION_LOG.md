# 实际研究迭代

1. 假设：固定Gaussian .5截断阻碍恢复。四优化器按train/exposedB比较；coarse_TR取得4/5相容并被冻结，24fresh最终15/24相容。未通过者保留，不用GT替换。
2. 假设：材料度量与prior-aware谱有助有限探索。八臂真实.02/.04对照；dynamic/adaptive成功率较高。单一修订为刷新触发分离decoder与Maxwell secant；paired v2保留，held-out随后冻结执行。
3. 假设：非局部prior jump或tempering补充分支。global方向未稳定发现重要结构，短cold chain移动不足；继续保留负结果。
4. C数值错误：property=-6.94e-18。暂停、只读回放、保存原8state/320RHS；仅机器精度latent边界恢复。同stage ledger累计，未重置，失败parent的首次启动计算仍计费。
5. 学习假设：54个真实train更新可摊销方向提案。Gaussian/Flow三seed实际训练并物理验证，未优于非学习对照；不扩大网络。
6. 求值假设：A17有限exchange降低费用。同空间误差1.91e-15，但A17更慢/fallback更多，停止增量主张。
7. 独立审查：随机angular progress可能削弱普通Full-J。在B剩余额度补充sphere-only两强基线；新增84states/2176RHS，明确exposed。Full-J可发现替代，削弱软谱独占机制叙事。原主分析未改。

最后不是再调参数：整体REDIRECT，停止当前学习与A17优势路线。阶段余量封账。
