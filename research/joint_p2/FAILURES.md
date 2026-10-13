# 保留的失败、负控制和修正

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

失败不是事后隐藏的安装细节。它们影响原始摘要、积分可靠性、预算及可声明的科学结果。所有作业stdout/stderr与状态在 [jobs](results/jobs/) 和 [JOB_LEDGER.jsonl](results/JOB_LEDGER.jsonl)，物理预约与终止记录在 [PHYSICS_LEDGER.jsonl](results/PHYSICS_LEDGER.jsonl)。原B1 384/384和C1-P1结果不被重跑或重写。

| 实际问题 | 保留证据与处理 | 对研究结论的影响 |
|---|---|---|
| 最大合法半径逆变换的floatroundoff得到z<−3 | initialphysicsbank在548成功预约后中断；原failurejob保留；只对**精确rmax端点**回到精确z=−3，后续从paidcache恢复 | 没有扩展prior或裁剪任意样本；最初pre-kernelsetupfail通过追加记录补收1次attempt |
| 失败setup在早期adapter先decode后预约 | [ACCOUNTING_REPAIR.json](results/ACCOUNTING_REPAIR.json)与retroappend记录；实际kernel未启动，历史row未删改 | failure-inclusive总费用不能只数成功kernels；后续decodefail也先计费 |
| 环境没有pytest | 该0physicsjob失败receipt保留；改为stdlib/direct执行原必要六项density检查，不安装工具链 | 环境失败不代表概率算法失败，不增加physics，不把额外smoke当scientificacceptance |
| 密度fixture的float32/float64不一致 | 原失败保留，修正fixturedtype后六项概率/Jacobian检查通过 | 数学标准未改；没有通过降低精度/容差获取PASS |
| coefficient-coordinateadaptiveintegration未满足2e−6内部准则 | 原512cells及误差保留；改成physicalr/h坐标并保留原priorJacobian/atoms | 同target下数值方法修正；七例5/9独立积分收敛，不把failure改称evidence归一化错误 |
| trainteacherorder9的mean误差fraction .001433 > .001 | 原teacher与audit保留；用order13同target细化，前八例17阶comparison改善至7.816e−5 | 仅修积分，不重新生成truth，不降低CDF/SD标准，不把teacher当免费 |
| 原freshnode-resampledCDF混有node量化台阶 | 保留诊断输出；另做directtruncatedcontinuousintegral与模型draws比较 | 旧nodeCDF不能支持连续材料精度，特别不能仅凭eventmass达标宣称fullposterior |
| FreshdirectCDF的rawgrid存在.00265042非单调缺陷 | 原CDF数组与audit保留；最终改同密度positivepartition积分，单调缺陷0，不做isotonic或uniform-bin展宽 | 分区evidence/event差及新旧CDF差.01476339仍披露；当前comparison仅empiricalfixed-grid/子集refinement，不能认证所有阈值 |
| 逐观测CDF重算过慢 | 已完成计算与约9分钟受控wall保留；主动中断后改用相同Gaussian充分统计量批处理 | 只优化相同likelihood计算，不修改noise、prior或target；失败/中断wall也计费 |
| 独立geometrycontext使posterior更窄，初始adaptive2048cells不足 | 原contextjob失败与sensorstream保留；改context-aware partitions及8192cells同target同容差 | 不能因为加入sensor而忽略参考收敛；不挑选有利sensornoise重抽 |
| 环形/组合/prior-boundary世界补偿仍有多noise-scale残差 | 全部best-foundendpoint和边界参数保留 | 明确的负物理结果；没有调大noise或放宽prior制造ambiguity |
| generic64particleSMC事件质量seed变化大 | 两个合法full-targetruns均保留 | 高ESS与合法correction不认证finiteparticle的mass准确性 |
| 64callMDNimportancecorrection漏访小事件、MC区间宽 | branches、原draws、zero-countevents、bootstrap保留 | 不能把零访问当posterior排除；不宣称完整probabilitycalibration |
| full/diagonaldistillationeventMAE差95%区间跨零 | 64parentpairedbootstrap保留 | 停止fullcovariance在所有事件上更优的主张 |
| Jointtruth-pairtraining退化 | 三seedfull/diagonaljointcheckpoint和fresh结果保留 | 不删表现差的seed；distillation使用额外监督，不能混成同信息优势 |
| ScholarQAretrieval四项429、部分全文缺口 | 原providerreceipts与resolvedmetadata保留于文献账本 | 不是全面noveltysearch；没有首次/IEEE TAP原创贡献认证 |
| 汇总脚本读取不存在的`metrics`字段 | `p2-final-prior-refinement-and-comparison`的KeyError与约11.075s失败wall保留，新增physics0 | 修数据字段后使用实际seed-separated结果；不隐藏报表失败或凭空填数 |
| 并行零物理job的旧全局countdelta误归其他owner预约 | [JOB_ATTRIBUTION_CORRECTION.json](results/JOB_ATTRIBUTION_CORRECTION.json)补充三个job自身newphysics0，原row保持不变 | 总费用始终由uniqueBEGIN为3459；不能把jobglobaldeltas相加、降低总预算或重置ledger |
| Rare-event补充particles后PPC原实现按样本list均匀抽取 | 原`old_uniform_sample_list_diagnostic`保留在七个package；改按wholeposteriorweights抽取1024个样本 | 小事件补点不能获得与主要事件相同抽样权；修正不提供模型adequacy或calibration认证 |
| Radius-uniform细化后outerlogevidence仍含旧节点近似 | 以独立smooth-coordinateadaptive值更新，原值保留`old_log_evidence_coefficient_nodes` | 原prior没有变回coefficientprior；不同数值近似的target与evidence分别声明，不混入主方法比较 |
| 前八fresh观测通过未能发现全部reference误差 | 全128event/mean/SD13/17核对发现观测12/39/48失败；升级17并交叉21后通过，原13结果保留 | 主比较转acceptedreference，不能把旧子集审计称全panel认证；新增physics0 |
| Positivepartition脚本未匹配括号 | 约.05s的syntaxfailure及job记录保留，修语法后运行同target积分 | 报表/执行故障单列，不降低数学标准或重置预算 |
| Acceptedreference已保存CDF后，stringparentID引发空分组/NaN摘要 | 原失败保留；修ID类型后复用已完成数组，没有重跑积分或physics | 新比较与原指标分开归档，不能隐藏partialsuccess/报告失败 |
| C1acceptedmetadata使用`data`而非`reference`导致NameError | 原失败与约10.515sreceipt保留；修变量后成功，冻结draws及weights不变 | metadata错误不改写成模型失败/成功；重放wall计费，新增physics0 |

不进入的范围另列NOT_RUN：非同心3Dcampaign、一般多物体、旧lockedtruth、A17冻结稿修改、GPU大型训练、state/memoryNN扩展、无uniformbound的全域材料CDF认证。NOT_RUN不是已证明无效；预算内实际负结果与未授权/未实施内容分别保留。


发布整理还保留`p2-positive-partition-CDF-reference`语法括号失败、accepted summary中parent字符串导致NaN及C1 metadata字段名错误的0物理cost receipts。后续分别修复；accepted参考数组已保存，汇总恢复没有重复积分或物理。服务器丢失已完成会话句柄时，以持久化job receipt为准，未再次启动已完成作业。
