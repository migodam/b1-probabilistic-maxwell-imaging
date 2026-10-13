# 真实资源费用与预算边界

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

独立阶段`B1-C1-JOINT-P2-20261013`由用户授权4000新增材料–频率attempt、7200受控数值wall秒，CPU优先，GPU大型训练0。旧B1 **384/384**与C1-P1账本保持只读。没有借用C1余额、改名历史ledger或新checkout重置费用。

物理总数来自唯一 [PHYSICS_LEDGER.jsonl](results/PHYSICS_LEDGER.jsonl) 的`BEGIN`预约，而不是只数成功响应。实际 **3459/4000**，剩余 **541**；其中1个pre-kernelsetup失败已retroappend补收，成功kernel3458，unfinished0。所有额外频率、selectedendpoint、teachertruth、SMC、correction、repeat与failure均计费。

| 类别 | cap | 实际attempt | remaining |
|---|---:|---:|---:|
| mechanism/qualification | 400 | 86 | 314 |
| response/reference/BayesianSMC | 2000 | 1847 | 153 |
| freshtruth/independentaudit/representatives | 800 | 758 | 42 |
| residualrefinement/MDNcorrectionreserve | 800 | 768 | 32 |
| 合计 | 4000 | 3459 | 541 |

主要用途是1250两频commonbanknodes、596genericSMCfullcalls、224独立field核验、256train/64validation/128freshparent频率评估、576两个demoeventresidualreference、86selectedrepresentative/MAPendpoints、192matchedMDNfullcorrection，以及mechanism86与retrospectivefailure1。cachehit548次保留paidprovenance，不收重复newphysics，也不能把coldsetup成本隐藏。

截至最终归档，kernel累计wall **12.5443s**、boundaryRHS **96488**。这些是Forward内部solve时间/RHS数，**不是端到端实验成本**。最终受控jobwall **2191.357s（36.52分钟）**，包括失败、重启、CDF积分、teacherrefinement、CPUtraining与evaluation；已包含gallery、全部失败与必要artifact作业，精确数值记入 [RESOURCE_SUMMARY.json](results/RESOURCE_SUMMARY.json) 及 [JOB_LEDGER.jsonl](results/JOB_LEDGER.jsonl)。两小时保护针对这些受控数值作业；重叠CPUjob conservatively相加，没有将并行wall免费抹掉。

15models×3000updates的trainingbodywall38.814s、CPU37.969s；冻结后的C1rawfresh评价bodywall42.004s、CPU41.440s。C1模块本身新增physics0，teacher/setup/freshtruth全部归B1唯一ledger，teacher不是免费监督。外部subprocessCPU在每个`CHILD_CPU.json`单列，不能把parent几乎零processCPU当完整费用。

originaljob start/end全局计数差在parallel零物理job中包含另一个job的预约，因此不能相加。为透明保留原row，另写 [JOB_ATTRIBUTION_CORRECTION.json](results/JOB_ATTRIBUTION_CORRECTION.json)：freshcontinuous第一次rawdelta86、batched192、contextresume192都属于其他concurrentphysicsowner，各自newphysics实际0。唯一`BEGIN`硬预算始终正确，future预约补owner tags；此修正没有减少任何总attempt或覆盖历史记录。

冷/暖proposal计时仅含modelreload/conditionexport/new64draws/logq，五次轮换顺序；不flushOScache、不包括fullcorrection。Full/diagonal/Gaussiancoldmedian1.446/1.685/1.357ms，warm.216/.218/.182ms。三种64-callcorrectionbodywall合计18.865s另列。SMC没有重复cold/warmphysicsbenchmark，记NOT_RUN，而不是补测后隐去费用。

GPUtraining0、oldlockedtruthreads0、DeepSeekcalls0、A17冻结稿修改0。剩余541只是授权余量，并不要求把它耗尽；当前参考和mechanism已足以作受限科研裁决，没有开启非同心大campaign或追加大型训练。所有原failedjobs与数值缺口见 [FAILURES.md](FAILURES.md)。


研究数值作业已关闭。完整控制wall不含对话/模型服务等待、人工写作与公开网络传输，不能当作整个chat的历时。最终阶段总账与逐job表见 [COST_TABLE.csv](results/COST_TABLE.csv)。五次cold/warm测量覆盖proposal acquisition；重复整套full-Mie inference的五轮cold/warm没有在剩余配额内运行，记NOT_RUN，不拿毫秒proposal时间冒充完整posterior速度。
