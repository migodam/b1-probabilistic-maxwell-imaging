# 最终静态一致性复核

结论：未发现会改变主结论的数字或科学表述冲突。此为文件之间的静态一致性复核，不是重新实验、重新评分或独立科学裁决。未运行测试/solver/scorer，未修改在线数组、账本、主报告、网页或旧项目。

复核对象：RESEARCH_G_REPORT_ZH.md、index.html、GATE_DECISION.md；主要依据为results/E1_LEGACY/RESULTS.json、results/E1_DISCOVERY_SUMMARY.json、results/E1_OFFLINE/METRICS.json、RESOURCE_SUMMARY.json。只读补充核对results/STRUCTURAL_WITNESS_AUDIT.json、results/E0_HISTORY_SIGNAL.json和RESOURCE_SUMMARY.md，以核实witness细节及历史表述。未重算数值指标。

| 主要表述 | 实际保存证据 | 静态结果 |
| --- | --- | --- |
| 一个新增结构witness，paired两轮不计两个世界 | pooled_same_observation中三个case结构代表含base为1/2/1，新增0/1/0；仅p250073的bandit t1在first/quality各记录一次相同witness | 一致。报告明确两轮重复不是独立复现，也不跨观测相加 |
| 合法边界/支撑变化，无新孔洞拓扑 | 已保存witness审计dM=.0230400803、dS=.0817670181、XOR10；support42→52；两者holes_all/holes_ge6均0；deltaS=.16不分离 | 一致。报告保留近工程门槛、非全局拓扑/分支证明边界 |
| invalid→valid且retention约.138 | RESULTS中raw/correction1/correction2残差89.9456235/25.7212398/22.3179695；独立验证22.3179695；retention=.1380130655；solver/constitutive约1.4e-15 | 一致 |
| Best-of-2/4/8无改善，Best-of-4均值.04224802 | METRICS.summary各方法/两轮保存值均.04224802130330251；paired_differences中Best complex MAE差为0；historical/support coverage均1/3，G尺度coverage为1 | 一致。没有把G尺度coverage替代旧指标或宣称bandit成像胜利 |
| 四臂legacy比较，完整五臂未做 | RESULTS.methods使用multistart/static_scale/QD/bandit；T_guided为WAITING_FOR_T_NOT_IMPLEMENTED；报告和网页均明确T-guided NOT_RUN/未交付完整版本包 | 一致 |
| 本次quality替换未触发 | DISCOVERY_SUMMARY各arm的prefix_t=3行quality_replacements均0；paired离线差均0 | 一致。E0旧池369次替换被单独叙述，没有当作E1获益 |
| 首轮每方法36states/864RHS，端点2/2/2/3和新结构0/0/0/1 | RESULTS各case/method local_costs为12states/288RHS；DISCOVERY_SUMMARY保存的每case最终行与报告首轮表对应 | 一致。共享3base没有重复算进每方法 |
| 全局291states/6936RHS、GPU0 | RESOURCE_SUMMARY为291/384、6936/12288、GPU0；受控wall76.305075539、CPU72.441919 | 一致。报告76.305/72.442为正常舍入 |
| 历史没有稳定额外结构预测增益 | E0_HISTORY_SIGNAL corrected_structural_novelty B−A log-loss=.00378563894，区间[-.03110815497,.04050682010]；B−C=.00349007565，区间也跨0；causal_policy_gain=NOT_IDENTIFIED | 一致。报告同时承认合法性/材料重复的预测关联，没有把无结构收益扩大成所有历史信号均无效 |
| 不宣称持续优越性，不扩大到RL或非平凡3D | 报告、网页与GATE_DECISION均为NO_GO_SEARCH，明确单对象witness、暴露development、完整五臂缺失与机制归因限制 | 文本相互一致；本复核不替代root科学判断 |

范围说明：METRICS.resources_copied_not_recharged保留较早时点的资源快照（wall54.350957372、CPU51.14678），不是最终总账。报告采用RESOURCE_SUMMARY最终76.305/72.442，符合命名和费用范围，不属于冲突。后续新增数值分析若入账，最终总账与文档中的时间快照可能需要同步；本次静态复核自身不追加数值job。

无必须修正项。网页中“所有方法两轮平均复材料MAE”的精确含义由邻近Best-of-4表头和主报告明确限定，未读成全部候选平均误差；如将来独立摘录该句，宜保留Best-of-K限定。
