# 失败与未解决事项

历史必要条件违反：

- results/PRIMARY_SUMMARY.json / distinguishable / H0 / Re(epsilon): ['variance_below_necessary_quantile_lower_bound']
- results/PRIMARY_SUMMARY.json / gain_mismatch / H0 / Re(epsilon): ['variance_below_necessary_quantile_lower_bound']
- results/PRIMARY_SUMMARY.json / gain_mismatch / H1 / Re(epsilon): ['variance_below_necessary_quantile_lower_bound']
- results/PRIMARY_SUMMARY.json / representation_mismatch / H0 / Re(epsilon): ['variance_below_necessary_quantile_lower_bound']
- results/PRIMARY_SUMMARY.json / representation_mismatch / H1 / Re(epsilon): ['variance_below_necessary_quantile_lower_bound']

原生四案例full-grid reference failure保留。历史bin-CDF与node covariance不同测度；原生/细化不是一个精确reference。四个旧treams oracle数组缺失。slogdet警告保留，Cholesky有限交叉核对通过。原生报告worker容量失败三次，父Codex完成，无模型/provider变化；成本未测。

E1/非同心/大型NN/locked/新发布NOT_RUN。完整失败日志results/FAILURES.jsonl与controls warnings；任何未解决数值条件不得通过放宽门槛抹去。
