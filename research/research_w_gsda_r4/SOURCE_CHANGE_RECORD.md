# 实现修复与费用记录

仅在新R4目录工作，未覆盖旧R2/R3实验。历史出处见ENVIRONMENT.json、historical_evidence/和INPUT_EVIDENCE_AUDIT.json；本地R4源代码未自动提交，最终ZIP为交付快照，没有虚构新的commit。

- 开发前复用依赖目录的scikit-fem版本元数据读取失败，保留FAILURE_LEDGER第一条；随后使用现有A25虚拟环境。没有为此重建旧物理数据。
- 原生电流审计第一次把配置状态元数据传入CSIOptions而失败：3次FE/32RHS、2.754659秒已记录Book CPU。解析器修复后重试5FE/64RHS、4.080089秒。失败目录与重试目录均保留。
- 数字实验之前比较了等价的cached-combined材料装配；同一物理模型与积分不变。另一稀疏插值策略未显示足够整体收益，没有替换主实现。独立c000开发的TV.2/1以及空间续接结果均保留，未用数字真值再挑先验。
- Full-GN是另列的普通基线补充，从已封存、测量选中的B0终点开始。它不是重新发现盲目标；原B0输出不修改，成本累计B0+GN。
- 六个控制族第一次启动错误地解析虚拟环境解释器的符号链接，导致系统Python没有numpy。全部在物理模块导入前失败，为0FE/0RHS；实际外部墙钟和stdout保留，CPU明确未测。只修正解释器路径后重试，没有更换材料、数据、参数或优化器。
- 后处理发现吸收计数器只识别ABSORB而未识别真实事件EXACT_ABSORB。修正为规范化、按动作/起点/迭代去重后重新评分；B2/B3实际每case2次吸收+2次定时重开。原TRACE和封存结果完全不改。该修复新增0FE/0RHS。
- 后期明确了Book计时从实例/分区开始、未含Python启动导入的缺口。保留原计时并对最后四个FWI子进程另测RUSAGE_CHILDREN；不猜补旧成本、不重复相加。详见COST_ACCOUNTING_SCOPE.md。
- 物理封账后，BudgetBook.check拒绝旧研究额度中的新物理调用。新的独立部署须显式新预算。封账后的工作只有数组评分、人工视觉复核、报告及打包，不自动开始下一研究方向。

所有实际失败见FAILURE_LEDGER.jsonl；物理预约、完成、缓存与失败费用见PHYSICS_LEDGER.jsonl和receipts/。没有清空预算、删失败或进行SHA校验流程。
