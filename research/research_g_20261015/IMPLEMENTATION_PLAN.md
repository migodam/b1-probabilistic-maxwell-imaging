# Research G 本轮执行范围

先完整阅读三个Pro交接文档，再做E0和实现统一搜索；本文件记录已选择的实施范围，不代替实验。

1. E0只读A径向1920条端点、A30二维轨迹与C1/A30源资产，区分raw、corrected、full-valid、archive。按parent检验prefix历史，保持两个物理target分开。
2. 冻结field-based archive：复杂材料距离、多尺度结构距离、success ledger、active archive、frontier及failure ledger。允许truth-free质量替换，禁止单链去重和重复发现奖励。
3. 使用用户明确允许的**完全相同、已核验历史后端**，复用A30 commit717d3与Research A commitf1的原adapter、原32实Haar GN函数。无需等T，也不重新实现T的forward/corrector/nucleation/A17消融。T-guided没有交付，状态WAITING_FOR_T。
4. 已暴露八条件按原残差/合法性筛选，只得到三个weak parent：250069、250073、250077。保留其余失败条件，不按truth挑对象。每个base在G重新付费full验证，再做小型四臂比较；这不是T新接口的完整五臂实验。
5. 两轮仅改变first-member与质量替换archive规则；同三parent、四方法、三轨迹、每轨迹raw+两次原corrector+独立full端点验证。计划291 states、至多6936 RHS；全局硬限制仍384/12288/10800秒，不能互换。
6. 在线材料与预测冻结落盘后，独立进程读取已暴露truth评分和绘图。TopK在线决定；Best-of-K不进入搜索。结果不理想也不增加parent、方法、步数或放宽阈值。

E0 archive-only修复重放是相同缓存池的选择诊断，不能当作新adaptive搜索反事实。E1使用真实新增full-state调用，不称缓存复用是物理加速。
