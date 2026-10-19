# 最终资源封账

| 阶段 | 计费states / 上限 | complex RHS / 上限 | 受控wall秒 |
| --- | --- | --- | --- |
| A | 383/600 | 12728/24000 | 66.594 |
| B | 482/1000 | 8312/34000 | 96.625 |
| C | 847/1200 | 16472/36000 | 83.535 |
| D | 240/700 | 7200/20000 | 65.645 |
| E | 284/500 | 2016/18000 | 6.844 |

总计2236计费states、46728 complex RHS、319.243s受控wall、166.339s物理operation wall，GPU0。

2012次fresh full-state/LU与224次保守计费的ROM/factor-reuse尝试分开。small-core、pool和训练CPU有action或wall记录；cached audit次数可推回，dense surrogate callback未全面逐次计数，时间包含在受控总wall。历史26states/976RHS与其他研究均未改。端口setup1056 RHS不能从exact-T query-only列中隐藏。

RESOURCE_SUMMARY.json保留主分析时的2152states/44552RHS；RESOURCE_SUMMARY_FINAL.json包含后续84states/2176RHS强基线审查，不是重复授权。
