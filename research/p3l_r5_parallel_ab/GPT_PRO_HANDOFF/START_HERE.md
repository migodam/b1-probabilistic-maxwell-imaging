# R5实验交接：两个不同的下一步

本轮实际跑完A/B、代码与失败都已保留。请先读[主报告](../RESEARCH_REPORT_ZH.md)、[机器指标](../RESULTS_AND_METRICS.json)和[裁决](../GO_NO_GO.md)，不要把开发结果改写为独立确认。

## 确实新增的实验现象

**概率：** 原P1的(0,64)二维子模型包含正质量rank0/2/1区域；上下平台各占.22086048先验。固定正秩响应图漏host层，正常root成功并不能挽救概率覆盖。正常coupled在sigma/4时，响应积分均值误差.000339 SD，而GN/Laplace约.1075 SD；另一方面高阶场正常例的Gaussian TV已经约.0043。不存在统一“新方法全面更好”的结果。

**物理：** 预定宿主接缝路径中，latent J为0，物理材料响应作用非零，真实Full-Maxwell响应步长斜率2.996758。这支持decoder三次退化，而非Maxwell物理无灵敏度。全局方向norm界又过松，六个R=3充分条件都不满足。

**风险：** Gaussian确定性统计可以取代K32特征并改善已有VAL风险，但直接积分在batch64更慢。一个214272参数raw-y统计蒸馏器保留92.96%fusion增量，额外provider时间2.54%，完整附加链54.57%。它压缩的是有用功能表示，不是合格posterior law。

**负结果：** 18个baseline训练构造的嵌套OOF中，四种线性RRR纠正仍无增量；Maxwell专属收益未建立。原CPU/MPS逐像素容差FAIL和初次内存失败保留。

## 本轮不能回答的事

- Full-66后验的全局质量与校准：NOT_RUN/未建立。
- Response atlas完整覆盖、严格TV/branch证书：未建立。
- 比全局归一化RTO、完整manifold lifting更省费用：没有证据，部分基线NOT_RUN。
- 新材料风险确认：NOT_RUN，旧VAL已用于模型选择；保留率探索区间[71.60%,135.63%]，不统计确认>=80%。
- Maxwell独立信息在强raw-y/decoder控制上产生稳定增量：未建立。

## 可立即继续的代码与冻结对象

A：actual-run原坐标/响应积分、分层体积、原始观察与节点、局部RTO失败和材料矩；见[A入口](../stage_a_response_coordinates/START_HERE.md)。

B：两个finalists各三seed的PROVIDER/HEAD、纯在线`Action.predict`、全部训练日志、嵌套fold名单、30次原始计时与实际VAL预测；见[B入口](../stage_b_material_risk/B_SUMMARY_ZH.md)。

物理与数值身份见[共享合同](../SHARED_STATISTICAL_CONTRACT.md)及[独立审计](../audit/FINAL_SOURCE_QA.json)。本轮已封账，未来工作必须新建授权，不刷新旧账本。

## 下一步顺序与所需资源

1. **先独立确认B，停止继续调模型。** 使用已冻结的确定性匹配和width128蒸馏模型、各三seed。新256+256的同配置R4强anchor历史费用估计约1.94M RHS，当前32768无法支付。未来建议另设约3.3M RHS的保守确认合同（按约6400 context RHS/parent、另计truth与失败），实测冻结对象级上限；这是待授权建议，不是本轮额度。
2. **A先完成混合秩reference，而非Full-66训练。** 优先把rank0/1层精确体积和transition数值误差闭合，比较合格RTO/分层方法。正常case的点wise精度不能补全未覆盖质量。
3. 只有以上各自有独立证据，才讨论response统计替代risk provider；本轮不融合。

更值得Pro推导的是[三个具体理论问题](THEORY_QUESTIONS.md)，而不是再次列举Flow/CNN/RL架构。
