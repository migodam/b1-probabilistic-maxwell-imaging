# 独立预算与实际费用

本轮使用新的Research-A授权；旧B1 384/384、Joint P2 3459/4000和旧A30账本保持原样。未借用P2剩余541次。所有新调用、拒绝、重试和失败都保留。

|阶段|实际新物理|授权上限|实际求解RHS|执行时间范围|
|---|---:|---:|---:|---|
|E0|0|0|0|包内八偶极代数约.0143s；旧cache/interface回放与元数据读取另列|
|E1径向vector-Mie|2301材料–频率attempts|2304|64368个球层边界匹配线性小系统RHS|受控job body约171s，加已知早期评分进程约1.90s；最终数值JSON在results/RESOURCE_SUMMARY_FINAL.json|
|E2二维TM/VIE|328 full-state attempts|384|9214个576维复LU RHS，上限9216|物理/搜索两次执行累计约18.13s；保存数组后分析另.265s|
|GPU训练/新3D|0|0|0|NOT_RUN_NOT_AUTHORIZED|

E1与E2的RHS含义不同，不能将两个数字相加解释成同一种求解成本。E2只余2 RHS，不能再付一次八源求解；其56个state名额不能绕过RHS上限。E1剩3次，在研究完成后封口，不为了补图继续使用。

E1物理分配：实际后端audit9、开发轮1共164、轮2共68、轮3共36、新truth32、确认集1952（共享base32+五臂各384）、第二频率24、无补偿8、独立结构分层归因8，总计2301。确认每方法每对象12次，少于原包16次上限；增加第五自适应臂时没有扩大总预算。代理优化、失败、磁盘写入、加载与模型提案在对应job wall中，而不是按免费full-call记账。

E2第一次审计的1state/91RHS在账本中；使用可分辨的任意近似电流重做dual恒等式，容差未放宽。为给这次额外开销留空间，首两条件各少一个重复no-comp endpoint（共少16RHS），最终累计9214，不重置run ID或ledger。

历史共享资产：P2五阶张量bank共1250个已付费材料–频率节点，原15个C1 CPU模型和teacher、A30已付费背景与N3 LU都保留provenance。各方法均可用同一proxy；MDN历史额外监督单列，不能以本轮没有训练费用宣称它免费。A30载入旧LU的实际验证和receiver adjoint本轮重新计费。

冷/暖记录见 `results/COLD_WARM_SEARCH_TIMING.json` 和 `results/TIMING_SUMMARY.json`。固定前4parent、5次轮换，重跑搜索/选择，1200次完全相同材料端点cache请求，新物理0；warm没有更快。没有OS冷缓存清除，也没有把整个冷campaign重复5次。共同base、初始化、模型condition export不混成完整冷成像加速。

成本账本是层级式append-only：根 `COST_LEDGER.jsonl` 记录E1与E2子账本索引；E2所有state/RHS/action原始行在 `results/E2/COST_LEDGER.jsonl`。32个交接包的global cost snapshot重复引用同一账本，禁止相加。每个candidate的endpoint费用只归到真实对应预约。

口径限制：计时是明确包围的数值job body；人工阅读、推理、文档写作、发布网络传输不是CPU数值实验时间。解释器启动/外部命令初始化未逐项分摊，不以这些body数字声称完整端到端冷启动速度优势。初始四次评分进程的已知tool-wall约1.90s已单列补充；其他后处理、失败与重绘由job/专门回执保留。


最终封口：E1受控body累计 **172.650421s**，加单列初始评分tool wall约1.904269s；E2搜索/物理18.129854s、原保存数组后分析.264885s，独立进程离线重放约.958518s。最新正式包decoder元数据修复计入E1job；未增加任何物理/RHS。公开打包与上传属于发布I/O，另列，不冒充数值算法执行时间。
