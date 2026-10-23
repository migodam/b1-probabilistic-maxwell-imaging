# 已保留的失败与不确定性

1. 原minimum_feature缺少严格统一定义；固定支撑资格不能替原宽prior全合法性。
2. measurement-only FWI未到d≤1。原R达到12000 RHS后封账；预算停止不是收敛。
3. 修复前conductivity上界一ULP导致一次LBFGS-B中止；旧receipt、部分端点与费用保留。
4. cal-only正确器保留材料距离却抹去任务差异；task-preserving大步cal=1.92998未合格。
5. seed1 joint=1.11267，未做通过叙事或额外救援；wide局部重优化3.00688未合格。
6. 首次P1B评分被未完成连续点控制的保护门槛拦截，0新增RHS、未读GT。补齐已预注册控制后使用新label完成，旧失败不删除。
7. SMC单次概率精度与reference资格不同；86个run至少一项branch/CDF超0.02，progressive observation21集成CDF=0.020374仍失败，全部进入评分。
8. proper branch评分未优于强prior，位置MAE无可靠改善，progressive风险差跨零。
9. response压缩的数组成本更高；高阶teacher成本不是免费，未满足20%端到端节省。
10. Green工程warm成本继承部分Full准备，setup嵌套，未实现部署fallback；不能以这些行直接证明standalone加速。
11. 最近文献检索两次429，原始失败包保存；不解释为没有先例，不作原创性认证。

逐次中止、拒绝及费用见FAILURE_LEDGER.jsonl、PHYSICS_LEDGER.jsonl及Release原始job receipts。B2最后的独立端点/argmax验证以最终正式门为准，不将进行中记录写成PASS。

Green账本完成事件漏写phase、而预约事件有phase，曾使通用审计器报告未配对。原FAIL快照保留。独立逐序核验后，只在stage/PID/scope/label的更早未消费预约上下文唯一时关联；歧义仍FAIL。没有修改原ledger或费用，也没有将live未完成预约当完成。见AUDIT/execution_integrity/GREEN_PHASE_SCHEMA_REVIEW.md。


## 最终尾项与封账元数据

圆腔cal=0.520776但dTsup24=1.379061<3，科学门槛失败；其最大项粗细审核通过仍不改变此判定。Joint最大项的独立网格审核亦完成。没有再搜索或修改P1prior。

最终封账helper先后暴露旧FAILED_PRE_EXEC receipt类型和早期缺失sampled RSS字段的兼容性错误，过程见AUDIT/closing。原账本从未修改；唯一旧startup收费得到独立匹配、缺失RSS明确列出，已写相同stage seals保留，不删除后重开。最终审计findings为空，但元数据关联/索引范围限制仍显式保留。
