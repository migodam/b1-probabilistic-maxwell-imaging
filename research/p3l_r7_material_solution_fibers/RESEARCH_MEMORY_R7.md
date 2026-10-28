# R7 研究记忆：历史证据、已排除路线与尚未成立的机制

2026-10-28；独立历史证据审计。此文只整理原报告、机器摘要、固定版本与已有文件身份，不作最终原创性裁决。新增 Maxwell RHS/LU/训练均为 0；没有读取 R7 GT、生成新观测或改写历史文件。

## 1. 阅读入口与身份

证据权威顺序是：原始固定版本报告及其限定、封账机器摘要、具体审计；handoff 和 R7T 是索引，不能覆盖较晚完成的原实验。详细逐条事实、原文摘录、源文件行号、完整提交见 [HISTORICAL_EVIDENCE_TABLE.tsv](audit/history/HISTORICAL_EVIDENCE_TABLE.tsv) 与 [JSON](audit/history/HISTORICAL_EVIDENCE_TABLE.json)。共 35 条历史事实、10 个固定工作树身份、17 份机器摘要读回。`FIXED_IDENTITIES.json` 记录报告 blob、报告最后修改提交及工作树 HEAD；10 份报告在当前树、原历史树、固定提交三者字节一致。没有运行 SHA256 检查。

|历史阶段|原工作树|本次读取固定 HEAD|
|---|---|---|
|F|b1-research-f|bb23534dd04983217668d301af6857f0db61885f|
|F3|b1-research-f3|bf27a17f5e3541fce8f1ec697966306476102d4e|
|P1|b1-p1-material-target|8e4519e687acfc5de7782c7b7ca0d055b013929f|
|P2|b1-p2-material-fiber|064d0aaed9772cadd573866eec84f159895b317b|
|P3-L|b1-p3l-local-posterior|c1171ab9c9006f5061f10293c7eadd5b977dd608|
|R2|b1-p3l-r2-risk-posterior|33471215c4e3005a9b810d5375e3ea2942d4d71b|
|R4|b1-p3l-r4-conditional-transport|092da74e320de38c086d07f24a1d98ca4892e34a|
|R5|b1-p3l-r5-parallel-ab|d813b4196ebb1340f051c0fd2804b7a49d256351|
|R5-B 确认|b1-p3l-r5-b-confirm|3e61f963f379ef4c6254d8b95f60297fe1bc19ab|
|R6|b1-p3l-r6-task-optimal|177a9ea37ada6dbad8ec67589bf3c6d6b99c0706|

当前 R7 分支起点也是 `177a9ea...`。上传 R7T 指定的科学版本是 `d1bdc8aceae4cd6a843f873b61b019d5eee3cd1d`，本地实际可以解析；不得把两者默认为同一版本。R7T 是 theory-only 输入，`inputs/PACKAGE_IDENTITY.json` 记录包身份，不将其中陈述升级成新实验结果。

R3 的完整原 TeX 和 PDF 实际保存在 `research/p3l_r4_conditional_transport/source_package/`。TeX 行 108–143 声明它是 2026-10-23 理论/实验设计，新增 Maxwell/训练为零，固定 R2 版本 `6eadb4690c962aaa95821fc4056b677f0f924a6b`。R7T `source/09_sources.tex:6` 所说“未取得 R3/R5T 全文”是 R7T 作者当时的范围限制；现在不能继续把本地 R3 写成缺失。此历史审计仅抽取 R3 的范围与成熟数学边界，不声称已逐行核验全部 4500 余行推导。

## 2. 共同模型与不能混称的对象

F/F3 是旧两 Gaussian P0 的低维条件实验。它们只涉及同一个暴露 development parent 在 ring/weak 两种背景，不是两个独立物体。F3 之后独立构造 P1，不能写成扩宽 P0/noise 救正例。

P1 至 R6 主要合同为：64 Gaussian 场坐标 × 2 Uniform 属性；Gaussianized 66 坐标通过 CDF 保留原测度；原 C² quintic saturation、24×24 scalar-TM/VIE、256 ROI、8 源 × 32 复接收、固定 gain。基础 sigma 为 0.00119211798514519；R6 的 sigma×0.5 和 ×2 是明确冻结的新噪声条件。Full 指原离散二维模型，不代表连续介质、3D vector Maxwell 或实测认证。

必须分开：同材料 decoder gauge、不同材料的有限噪声竞争端点、exact response fiber、posterior draws/质量、材料平方风险动作、优化驻点、计算速度。noise-ball 是资格诊断，不是 Gaussian posterior 的硬支持；同一 y 的确定性特征不增加 y 自身的信息。

## 3. 科学历史的正面与负面证据

下表的 ID 可在详细证据表中查到完整原文、源行号和提交；数字保持原报告口径，不重新解释最终 novelty。

|阶段 / 证据ID|保留的正面结果|失败、未完成与适用域|
|---|---|---|
|F / F01–F03|区分了参数坐标加速度与材料像真实法向曲率；两份 physics-only 逆解合法且 noise-compatible|严格驻点 FAIL；F2 四臂 NOT_RUN。旧 Maxwell 二阶向量缺失，仅 norm，不能恢复完整三向量归因。未确认非零 exact tangent kernel；learned Explorer NO-GO。|
|F3 / F301–F304|公共权重尺度商保留原材料 image，非均匀 rho prior 解析并核验；新冻结 quotient/KKT 门下两个 anchor PASS；weak 有任务尺度竞争端点|旧 F 八维 Gate2 FAIL 不改。边界 likelihood anchor 不是 posterior MAP。三个接近端点不是三个 modes；同56RHS普通 native 找到2个 weak 端点，高阶D为1，无增量。8步 MH 无混合/质量资格。|
|P1 / P101–P103|先材料构造，后独立授权物理资格；A residual22.2923<23.7854；26states/976RHS|初始材料报告的 physics NOT_RUN 必须与后续阶段区分。后续 B residual309.7873 FAIL；A非认证驻点/MAP。两个末点 Dg 数值 rank66。比较 A17/Full-J/ROM 等当时仍 NOT_RUN。|
|P2 / P201–P204|24对象中15个相容，archive有69个任务尺度不同替代材料，heldout29；保留真实 finite-noise 运动证据|heldout仅4/8合格。Best-of4 MAE 约1.35%是混合 archive oracle 诊断，不是匹配预算胜出。soft/dynamic未明确胜 hard/Full-J；Gaussian/Flow无收益；全局搜索短链不合格。E的误差判断为offline oracle；A17约11.1倍Full费用且17/48fallback。数据生成/反演同OS进程偏差保留。REDIRECT。|
|P3-L / P3L01–03|Probit坐标点相对native MSE降低11.71%；解析 supplied-Gaussian 材料矩可用；保留decoder降低代理响应误差|128对象近似平均 MSE .0615021 对点 .0110309，+457.5%；11组 FULL reference 全部不合格，SMC128单祖先。更准代理未改善材料MCSE且更慢。不能由错误分布平均推翻真正Bayes均值最优性。REDIRECT。|
|R2 / R201–04|首批新256对象 CNN风险降低25.27%，调整区间22.34–27.99%；另256对象确认static sensitivity小增量；1934/2048合法端点满足noise/ΔPhi竞争，实际步长.02/.04|endpoint不是parent/mode/质量。39个变差对象保留；仅138/256学习图过noise-ball。轨迹平均风险约+1.96%，policy未成功；FULL UNQUALIFIED。弱截断projector bracket非对合证据不是exact-kernel不积分定理。|
|R3|提出动态弱几何、完整条件输运、材料矩与刷新控制的理论和可证伪设计|theory-only，零新增实验。成熟Bayes/LIS、Gaussian更新、Frobenius等应用不能直接作原创声明；R4才实际检验部分机制。|
|R4 / R401–03|两批512/256新对象：相对R2风险改善13.10%/11.36%，相对强CNN5.36%/4.90%；Full66 properNLL与CFM energy score有真实相对改善；同8192步决定驱动重放71.397→24.153秒，逐步决定/路径一致|Gaussian global统计与Spline同样有效；绝对Phi失配、SNIS退化，posterior未合格。缺同谱decoder-only CFM臂，不能隔离传播增量。CFM密度稳定门FAIL，8求解NOT_RUN。111候选不是modes；四段路径仅2项通过完整支持门，直接法94。速度只属beta=.002短核，不是混合/MCSE加速。|
|R5 / R501–03|A正则区域局部材料矩更准确；B开发VAL保留92.96%增量，provider模块显著压缩|A全域覆盖/branch归一化/TV/等精度20%成本优势未建立；host RTO0/8，Full66等NOT_RUN。B新256+256确认受硬预算NOT_RUN；CPU/MPS3e-6 FAIL保留。A REDIRECT，B CONDITIONAL GO。|
|R5-B确认 / R5B01–03|512全新对象、冻结模型、封存预测后评分；student对StrongCNN MSE降低3.11%；Q64 .006854734优于K32 .006943309|>=80%增量保留率HOLD：点80.70%，调整区间60.28–113.12%，contrast下界−5.43036e−5。不能用ensemble救门。原先VAL92.96%不是独立结果。shared physics anchor仍主导完整费用。|
|R6 / R601–05|局部分层矩、真实Maxwell copula控制反例；两个新noise各256对象的原冻结风险/费用门通过|A仍无全域参考、task精度费用分离或Full66；B新目标未胜历史部署器；C/D原通过要同时附强Full-q与计时限定。GO-TAP REDIRECT；R5-B HOLD不改。|

## 4. R6 的限定必须和节时数字一起出现

原统计合同内，24+64 q候选对原完整系统：低噪MSE .005335107 对 .005463559，风险−2.351%，原窗口节时34.551%；高噪 .009555500 对 .009534329，风险+.222%，节时32.571%。两条件均通过原冻结1% NI与20%wall门，不能事后删掉 PASS。[R603；原报告:85–92]

强化同特征 Full-q 后，低噪完整96+64为 .005315940；候选比它高.3606%，1% NI contrast上界+1.26745e−5，强化比较资格不足。Direct-prior Probit128在独立冻结补充64+64对象上很有竞争力，不能把Native分配称为独有机制。[R604；原报告:100–104]

同窗口重放first16/噪声、三算法，共96路径/268544RHS，端点完全一致；Low节时27.65%/32.77%，Direct节时29.31%/36.34%，二者无确定胜者。低噪新窗口20%节时下界16.82%，不支持跨窗口至少20%的普遍说法。该小样本复测不能冒充全512新的确认。[R604]

封账为6,860,984/8,000,000 RHS、210,355states/LU、33训练阶段。unused1,139,016已封，不构成R7预算。R6-A截断2D/3D条件模型不是Full66 marginal；局部GH矩一致不等于全域漏质量被控。其copula例证明特定条件依赖不能独立化，并不独自解释旧Full66 Gaussian失败全部原因。[R601、R605；RESOURCE_SUMMARY.json]

## 5. 已停止的具体路线与保留价值

|不应按原主张继续包装|原始淘汰依据|仍保留的价值|
|---|---|---|
|P0 exact-local-fiber learned Explorer、高阶D优越性|F/F3 rank/尺度/同RHS增量不足|材料商、非均匀诱导先验、边界KKT和合法retraction基准|
|弱密度×单strong-Gaussian Full66成像|P3-L平均风险严重恶化，reference全不合格|解析Gaussian材料矩是输入分布上的正确工具|
|当前小Flow/大网络可自动修复posterior|P2、P3-L失败；R4相对学习进步但绝对响应失配|保留概率学习价值，不宣称理论不可能；逐类检验响应/覆盖|
|动态非完整路径天然给后验质量/风险收益|R2、R4可行运动有正面证据，路径平均和四段推进无收益|作为候选提案/可达性诊断，不带自动体积或概率权重|
|小dense TM24 A17/低端口秩独有加速|P2与P3-L总成本/回退负面；R4 ordinary GMRES也有效|标准current反馈和带任务的数值控制接口；不可从低latent维推低物理更新秩|
|R5蒸馏>=80%确认/端到端大加速|独立R5-B HOLD，共同anchor主导|有风险收益的功能模块压缩，需保留完整历史成本|
|R6 Native截断已是新Maxwell任务理论|Direct与Full-q强对照、窗口漂移削弱机制|固定合同内真实工程节省，后续作为强普通比较|

## 6. R7T：标准恒等式与待证机制分开

R7T `CLAIM_LEDGER.json` 24项完整保留。下列是输入包自己的分类，不是本审计新的最终 novelty 判断。

|R7T claim|分类|R7需要新增的证据|
|---|---|---|
|C01–02 材料pushforward posterior、给定材料的decoder条件核|成熟disintegration/Markov结构|实际decoder消元、atoms/overlap/tails处理，不能将同材料多latent重复计质量|
|C03–07 光滑商/原子否定Flow/S3普遍降维/exact-kernel不可积/非对合自动质量|若如此推论则错误；ledger明确否定|满秩、有限等材料校正、不同strata实测；弱谱分布与exact kernel分开|
|C08–10 旧端点、相对概率进步、材料风险|限定历史实验支持|保留正面事实，不升级成Full66校准|
|C11 current任意阶递推|已有物理数学的具体推导|递推本身不申报原创|
|C12 补偿后的首个响应阶与有限工作区间|Maxwell特化的新组合候选|非均匀真实材料、affine控制、强Full-J对照、系数误差与两个相邻固定距离；任务尺度档不能删掉|
|C13 共享补偿冲突分解|加权least-squares标准代数|同一个U、坐标、低阶材料系数；新采集冲突的预测和用途优于强普通指标|
|C14 原子branch质量阶与任务阶|奇异渐近的条件应用|正测度/支配/tail和真实成本，零测度弧不决定Full66概率体积|
|C15–18 协方差导数、条件期望风险、正prior下矩阵稳定性|成熟工具；零一阶协方差不等于无信息|Maxwell特有低成本近似及有限测量/高阶任务增量仍待证|
|C19–20 centered-task误差界与全局认证|有明确常数/覆盖条件的一般界；有限黑箱查询无一般全局保证|所有Z下界、漏块质量和任务包络必须可用；经验refinement不叫certificate|
|C21–24 R6截断、便宜任务积分、统一论文、同y特征信息|当前理论/成本/统一链NO-GO；同y确定性特征不增原Shannon信息|独立机制、全成本、条件任务价值，不机械拼装三线论文|

R7T `source/03_maxwell_jets.tex:23–29` 特别区分显式decoder链式法则与物理隐式反馈、多次散射阶数与材料Taylor阶数。`REVIEW_AND_SCOPE.txt` 的13条限定包括：实Hilbert补偿、共同物理基、零测度弧与概率体积、局部证书空不能阻止预定有限距离验证、两个积分接近不排共同遗漏。它们应保留在具体claim判定里。

## 7. 新 R7 parent 提案的独立性审计

父线程提供的八个候选为 `R7_LOCKED_000`…`007`；master2026102809；prior直接使用 `SeedSequence([2026102809,i,1])`，noise使用role2；E1 master2026102807另域。审计没有生成材料/观测/GT。

[机器审计](audit/history/R7_LOCKED_SEED_INDEPENDENCE_AUDIT.json) 在59份历史元数据文件中读到38507个parent标识、2545个scalar seed值：八ID无冲突，16条新prior/noise初始四个raw64词无历史scalar-seed初始流冲突，彼此无冲突。状态为 `PASS_METADATA_AND_INITIAL_STREAM_DISJOINT_WITH_SCOPE_LIMITS`。这只是可见元数据和初始化复用审计，不是随机独立性的数学证明，也未证明没有未收录manifest。必须直接把SeedSequence传入default_rng；将generate_state所得uint64再作为scalar seed会变成另一流，机器审计中uint64仅是身份marker。最终生产代码仍需记录NumPy/bit-generator版本、entropy及角色，并将新对象封存于在线/GT分隔流程。

DEV明确复用R6 CDEV前4个，已暴露，不计入新八对象。完整原raw注册表是 [DEV4_RAW_ANCHOR_REGISTRY.json](audit/history/DEV4_RAW_ANCHOR_REGISTRY.json)，只检查NPY headers，未读GT。历史原树下每个parent：

- `raw/online/{parent}.npz`：观测与sigma，源脚本 `src/run_c.py:17`。
- `track_c/contexts/CDEV/{parent}/native96_probit64/CONTEXT.npz`：y、最终x/q0、chi0、J512×66、预测/残差、H和native材料。
- 同目录 `original_optimizer/probit/RAW_MODE_CACHE.npz`：q、chi、current576×8、field576×8、prediction、J512×66、receiver_T32×576，四parent文件全部存在。

这些是原 `b1-p3l-r6-task-optimal` 工作树中的真实已付费raw资产；compact handoff不等于原资产。表中绝对路径指向原树，避免重复昂贵anchor。当前头的这些大数组通常未tracked，固定commit认证的是读取/保存源码与历史身份，不假称Git blob覆盖所有raw内容。下一步使用者可按实际需要做已有数组一致性核验，但此worker没有运行新物理。

## 8. 审计自身范围

本审计不运行新的ScholarQA/网络文献验证，不把R7T文献表当最终先例检索；final novelty由顶层Codex负责。没有读取R7新GT或旧科学大数组来选择方向；仅历史汇总JSON和DEV4头信息。正文中的百分比保留对应原父样本、seed-loss、baseline和区间口径，不跨cohort汇总。

仍需顶层处理：R7真正执行协议与预算、E1/E2方向和距离冻结、物理/decoder机制的最终解释、R7T条件命题逐项证明审查、与成熟算法和最新文献的原创性判断。历史各类成败不能替代这些工作。
