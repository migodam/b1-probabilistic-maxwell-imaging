# E0 history 子任务摘要

- 已完成：32 parent ×5 arms×12 steps，1920 rows，5折 parent CV，A/B/C/D + intercept；零新物理。
- 主结构标签：corrected physical p 与全部严格prefix full-valid ledger比较，δr≥.05R OR δh≥.025R，不依赖原material-distinct flag。
- 主结构 B−A log-loss=+0.003786，95%[-.031108,.040507]；B−C=+.003490，95%[-.033365,.040816]。没有稳定额外结构预测信号。
- 材料新颖/合法性有历史关联，但高维B的合法性log-loss .1891劣于intercept .1515；D history-only .0763。模型受32parent和线性规格限制；不能断言历史无用，也不能推出RL收益。
- static/adaptive相同384端点，原archive新增材料161/154、结构128/122；adaptive proxy84433 vsstatic22259，重复region131 vs96。
- A30 E2 4 parents新结构全零，未伪造分类性能；未混入Mie训练。
- 无propensity/反事实，不能识别新policy价值。
- 第一次数值Newton饱和输出保留，后加单调backtracking；根指令全prefix标签在新标签CV前冻结，原标签secondary保留。
- 两次分析累计 wall 0.802152s / CPU 0.771077s，GPU0 / RHS0 / states0；请root合并 E0_HISTORY_JOB.json 计费。

产物：code/e0_history_signal.py；configs/E0_HISTORY_SIGNAL_FROZEN.json；results/E0_HISTORY_SIGNAL.json；E0_HISTORY_GROUPED_PREDICTIONS.csv；E0_HISTORY_FEATURE_PROVENANCE.csv；E0_HISTORY_SHUFFLE_MAPPING.csv；E0_HISTORY_JOB.json；SEARCH_HISTORY_SIGNAL_DRAFT.md。

Root裁决：阈值为历史径向几何diagnostic，不能冒充新G的二维结构metric；历史弱预测不支持PPO/GRPO。当前结果足以支持保留轻量可证伪history arm，必须靠实际单位物理成本新结构数检验。
