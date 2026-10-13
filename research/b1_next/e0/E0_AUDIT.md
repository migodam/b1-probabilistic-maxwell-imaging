# B1 E0 实施审计

六个案例、三条件事件及材料边际的25/41/65/97规则与独立二维分区自适应积分已执行，全部同一continuous-polynomial target经验条件通过。指定prior、noise、gain logdet、evidence、TopK余量、坐标变换、漏模式、metadata回放已核对。72历史条件坐标有5项必要条件不相容；旧摘要/旧source arrays保留。

新moments与direct CDF属于同一continuous density；原生节点moments与bin-CDF被标为不同近似，并非evidence归一化失败。全部216个历史quantile-cut CDF probes直接积分并保存。旧treams四原始oracle arrays仍缺，未重跑；未来source保存字段不追溯补齐。

隔离分支research/b1-e0-posterior-audit，main未修改。新增物理0/0，历史384/384；未借用C1。只有指定必要概率核对，没有旧全suite/smoke/哈希。原生subagent写controls；报告worker容量失败后由父Codex完成，不换model/provider。

新报告及裁决见RESEARCH_REPORT_ZH.md、GATE_DECISION.json；可复用资产见REUSABLE_ASSETS.md。同target积分参考不认证真实Maxwell/continuum/仪器概率；E1新增物理仍未授权；公开发布现已由用户授权。
