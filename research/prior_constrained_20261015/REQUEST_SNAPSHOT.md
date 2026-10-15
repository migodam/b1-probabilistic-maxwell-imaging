# Prior-Constrained Maxwell Multi-Solution Imaging — E0/E1/E2
# Execute a falsifiable pilot, not a large training campaign.

你在现有 B1 工作区执行一个独立、可复现的小实验。
输出目录：research/prior_constrained_20261015/
不要覆盖 T/G/A30 的历史产物、预算台账或结论。不要重新编号 B1/C1。
先完成 E0；仅在条件满足时自动进入 E1/E2。失败时交付负结果，不扩大预算寻找成功。

## 0. 问题与本轮裁决目标

检验原 A25/A30 材料分布是否支持至少两个结构不同、材料合理、对同一观测具有接近似然的竞争解释。
本轮不以强制生成 K 个候选为目标，不训练 RL、conditional diffusion、flow 或大型 NN。
N3 是强 Top-1 基线及初始化，不是 posterior，也不是所有候选必须贴近的硬先验。
必须区分：
(1) box 合法；(2) 数据 noise-ball 相容；(3) 点似然相近；
(4) 原材料先验支持；(5) 多个 posterior 高质量区域。
本轮可以建立 (3)+(4) 的竞争见证，不得把它直接升级成 (5)。

## 1. 读取与冻结来源

优先复用现有本地 checkout；必要时仅获取所需公开文本和已发布数组。
在 SOURCE_MANIFEST.json 中记录实际读取的 commit、文件 hash、数组来源、条件及坐标约定。
至少读取：
- migodam/b1-probabilistic-maxwell-imaging
  research/research_t_20261015/RESEARCH_REPORT_ZH.md
  research/research_t_20261015/ 中实际存在的接口、corrector、candidate arrays、预算文件
  research/research_g_20261015/GATE_DECISION.md
  research/research_g_20261015/ALGORITHM_SPEC.md
  research/research_g_20261015/ 中实际存在的 raw/corrected/archive 记录
- migodam/a30-maxwell-observable-material-manifold
  已知参考 commit: 717d3d26a290248bd2b9aeb55f3d03694b407b84
  research/delegated/a30_inputs/A30_INPUTS_AUDIT.md
- migodam/a26-om-feedback-prior-completion
  src/a25/dataset.py：generate_manifest、parent_material、chi 的组装
  src/a24/cases.py：backgrounds、noisy_observation
  实際数据 manifest / geometry / config / N3 checkpoint

不得用报告中四舍五入数字代替已有精确数组。
无法获得某项数据时写 MISSING，不猜路径、不补造结果。
旧 32 条件及 250065–250080 是 exposed development，不得重称 blind test。
不得打开原 locked test 的材料、几何参数、真值、clean response。

## 2. 冻结物理与真正的先验

保持原 24x24 scalar TM/VIE、8 sources、32 receivers、complex128、gain=1、同一 S/G、同一材料定义。
用原 convention 白化：每个复通道实/虚部乘 sqrt(2)/sigma，512 实通道。
ell(x)=0.5*||W(F(x)-y)||^2；tau=sqrt(chi2.ppf(.95,512))。
同一 y 的所有候选使用同一个已知 sigma；不得按候选预测幅度重估 sigma。
不得修改噪声、gain、支撑区域或 acquisition 来获得成功。

从实际 config 读取单位复相位 v、背景 chi_b、幅度 a。确认：
chi = chi_b + v*rho, rho real, nonnegative；
精确 P0 的 rho 来自原 parent_material 分布，包含其归一化和参数范围。
a1 是索引，必须读取其实际幅度，不能将字符串解释成 amplitude=1。
固定相位约束作用于 chi-chi_b，而不是绝对 chi。
rho≈0 是回到背景，不自动等于空气 chi=0。

始终分别标记：
- BOX：历史材料 box；
- PHASE：原增量相位、范围和背景约束；
- P0_EXACT：有合法原生成器参数、且端点确实由原生成器产生；
- P0_APPROX：仅通过统计/拟合近似先验；
- P1_RELAXED：显式扩展的先验；
- PRIOR_UNKNOWN。
不能将 PHASE 或低 denoising error 写成 P0_EXACT。
原 P0 若内部含四类形状混合，可以使用这一已存在的生成分布；
不得新增 hole-count 桶，不把真值 family/subtype/geometry 提供给在线算法。
本轮主要存在性结论必须由 P0_EXACT 见证支撑。
P1 只允许另表诊断，不允许偷换为原 A30 正结果。

## 3. 预算、访问隔离与实现约束

E0：0 新 full state，0 新 Maxwell RHS，0 GPU 训练。
可生成至多 2048 个 material-only 样本：只按 TRAIN 分布、用全新 seeds，
不得复用/反推出 locked 或测试 parent 的 seeds/geometry。
记录它们是新增材料样本，不是新增 paired Maxwell 训练数据。
只有 64 个独立训练 parent；512 condition rows 不等于 512 个独立几何。

E1：最多 64 个新 full-state attempts、2048 complex RHS。
E2：最多 192 个新 full-state attempts、8192 complex RHS。
全轮总上限：256 states、10240 complex RHS；任何更小的既有有效上限优先。
失败、拒绝、line search、full validation、reference construction 均计费。
导数/adjoint、factorization、wall time、材料生成、prior fitting、cache reads 分开计费。
到 cap 即停，不借用 T/G 余额，不自动追加实验。
这里是资源上限，不是要求必须耗尽。
数值结果没有严格 forward-error bound 时标记 FULL_VALIDATED_NUMERICAL；
线性残差小不等于输出误差证书，更不等于连续 Maxwell 验证。

设置严格接口：
online 仅可读 y、sigma、acquisition、公开 background、N3、TRAIN 材料先验及合法搜索历史。
offline evaluator 可读 development truth 或新 benchmark 隐藏见证。
新 benchmark 的隐藏两端、clean means、生成参数、truth identity 不得通过共同文件泄漏在线端。
所有候选冻结后才做 MAE/Best-of-K/见证覆盖计算。
人工构造 witness 的耗费单独列为 benchmark construction，不能算作免费在线候选。

## 4. E0：零物理审计

对全部可用 N3、T raw、T corrected、G raw/corrected、基础图和已有 pool：
1. 从精确预测缓存计算 residual、ell、对同一 y 的 pairwise log-likelihood difference。
2. 计算固定相位偏离、负 rho、背景下穿、幅度/峰值规则违反。
3. 将全域 DC/低频补偿、局部结构变化、core 外补偿分别计量。
4. 计算连续 rho 的多阈值/多尺度结构；不以 abs(chi) 或 abs(chi-chi_b) 的一张图识别空气。
5. 用训练材料、parent-disjoint folds 或 material-only 新 samples 校验先验可接受率。
6. 检查 prior-accepted raw proposal 经 corrector 后的 prior、结构、likelihood 保留情况。
7. 从已有同 acquisition 的合法先验材料及 clean-response 缓存寻找结构分离但预测接近的 pair。
   每个比较必须冻结同一个 sigma，不能用各自 SNR normalization 制造近等价。
8. 记录历史 truth 落在 95% noise ball 外的情况。这是预期有概率发生的事件，
   不称物理非法，不为消除它而改 tau。

必须报告对原 P0 membership 的证据强度。
无法拟合回原生成器，不等于证明不属于 P0；写 PRIOR_UNKNOWN/RETRACTION_NOT_FOUND。
若只有 learned prior 可用，先验验收不足以证明 P0_EXACT。

E0 输出：
E0_AUDIT.csv、PRIOR_SUPPORT_AUDIT.md、LIKELIHOOD_AUDIT.md、
CORRECTION_ATTRITION.md、CACHED_PAIR_SCREEN.csv、E0_DECISION.md。

## 5. E1：竞争材料见证，而不是一般重建加样本

先使用 E0 缓存见证；不足时进行有预算的 P0 内 pair search：
min ||W(F(chi_a)-F(chi_b))||
s.t. chi_a、chi_b 都为 P0_EXACT，结构距离超过预冻结阈值。
不能靠任意全材料补偿、极端先验尾部、扩大噪声或缩短 acquisition 制造 pair。

结构距离至少联合：
- 归一化 rho 的多尺度边界/占据差异；
- 连续 excursion-set 的持久结构差异；
- DC/总体幅度单独报告，不能靠其单独通过结构 novelty。
阈值由 TRAIN/material-only prior 的噪声稳健尺度预先冻结，不用测试 MAE 调参。
保持历史 normalized L2 与 changed-cell 指标，仅作对照，不单独作为新资格。

优先寻找 d_pred<=2 的固定 pair，并完整报告 d_pred，不将 2 当物理常数。
若无此类 pair，报告最小已找到的距离及结构/先验代价，不声称全局唯一。
对先固定的 pair，用真实模型生成 y=F(chi_a)+epsilon 或 F(chi_b)+epsilon；
不得把两者预测的中点冒充某个 truth 的正常观测。
冻结 32 个独立 noise seeds，利用已算 clean means 零 forward 重放；
报告两端同时满足 noise gate 及 |Delta ell|<=log(10) 的频率。
log(100) 仅作预声明敏感性表，不用于事后更换主门槛。
不得挑一个有利 noise realization 充当全部证据。

对最优见证，在剩余预算内做非零、预冻结尺度的 prior-local 扰动，
检查是否仅为孤立/脆弱碰撞；样本少时仅报告局部稳健性，不声称 posterior mass 已估准。
数值边界样本必须 tolerance refinement 或记录 BORDERLINE_NUMERICAL。

至少发现两组不共享 material parent 的合格 P0_EXACT pair 才进入 E2。
它们需结构分离、重复噪声下确有竞争、无明显极端先验/数值问题。
未满足则交付 NO_GO_FOR_GENERATOR_COMPARISON；
保留 ONE_SUPPORTED_CANDIDATE / AMBIGUITY_UNRESOLVED 状态，不追加训练。
没有找到 pair 不是不存在 pair 的证明。

## 6. 推荐算法：PC-SP（Prior-Constrained Structural Profiling）

实现一个小的数值搜索器，不训练新的条件生成 NN。

表示：
- 基本材料变量是固定复相位的实幅度场 rho。
- 多尺度 field/wavelet/DCT 是搜索与结构坐标；不能用固定 rank-32/64 作为永久材料空间。
- P0_EXACT 分支必须最终回到原生成器输出。
- 若原硬 rasterization 不可微，用 material-only prior secants / 可行有限参数变化构造提案；
  不把零几何导数解释成可辨识性，不悄悄平滑/更换原生成器。
- 连续 field surrogate 只能用于提案；端点投回 P0 后必须重新做结构检查和 full validation。
- 原生成器 retraction 失败时写 PRIOR_RETRACTION_FAILURE，不用 box GN 替代。

初始化：
保留原 N3；另构造 phase/prior-respecting、最小必要修正的基础候选。
不追求最低残差，不要求所有条件强行达到 8/8。
全局候选从 TRAIN/material-only prior 取得，不从隐藏 witness 或 development truth 取得。
每臂至少三分之一的 raw proposal slots 用于独立 global proposals，不能全围绕 N3。

局部方向：
在当前合法 anchor，用 A = whitened material Jacobian composed with prior-allowed directions。
用 Q 表示 shrinkage/多尺度 prior metric，必须正定，不能用 64 parents 的奇异 covariance
把未见方向全部硬设成不可能。
用 H 表示结构 descriptor 的局部变化。
优先方向解近似广义特征问题：
H^T H u = lambda (A^T A + alpha Q) u。
alpha/步长在训练材料与数值尺度上冻结；谱排序只是 proposal ranking，不是非线性可行证书。

沿所选结构方向逐步增加结构目标 t，
反复询问“达到 t 的先验合格材料能否保持近最优 likelihood”，而不是任意最大化差异。
有限步失败就缩小/终止该结构延拓，不无限放开 compensation。

共享 corrector C_P：
对所有生成器使用完全同一实现、阈值、budget 和 stopping。
在硬物理/先验、likelihood、补偿预算、结构保留约束内，
优先找可行点，再最小化相对 raw proposal 的改动。
固定结构保留率 eta=0.5：
沿 raw-seed 的结构方向，corrected-seed 的投影不得低于 raw 的一半；
另外检查最终离散/持久结构，不能仅靠投影冒充保留。
局部提案分别限制 edit 外补偿和 DC/低频补偿；
global 提案使用相对自身 prior seed 的修正预算，不施加相对 N3 的小半径。
这些限制由先验材料的自然变化尺度冻结，并对所有臂相同。
最多两次接受的 nonlinear correction updates；所有尝试端点均计费。
不能返回无限制全材料 GN 的结果并称先验修正。

likelihood gate：
先通过与 y 无关的 prior gate，再用同一 prior 下的公共 best-known likelihood。
原始 off-prior Full GN、隐藏 truth、隐藏 witness 不进入在线 reference。
online 使用共同公开 seed reference；最终额外使用所有臂合法输出的 union reference 做统一重判。
主门槛 Delta ell<=log(10)，敏感性 log(100)；同时记录 noise gate。
只能称 best-known-relative；除非有证书，不称已找到全局 likelihood optimum。
N3-anchor penalty 不得解释成独立 p(chi) 或用它重新计算 proper posterior。

## 7. E2：同一 corrector 的四臂最小对照

最多 4 个观察任务，来自 E1 的两组 pair，两端轮流作 truth。
冻结后在线运行：
M：prior-respecting multi-start；
T：有限低对比度/边界编辑 + C_P（标记 T-soft/adapted，不冒充历史空气孔 T）；
G：原 QD/archive 逻辑的 phase/prior-compatible 版本 + C_P；
P：PC-SP + C_P。

原强制 chi=0 的 T 单列历史诊断。若不在 P0，不把其自动失败当作新生成器优势。
四臂相同 prior、支持、N3、共享 prior reservoir、global quota、validator、noise。
每观察每臂最多 10 新 states；4*4*10=160；
其余 32 states 最多用于公共基础、必要的数值复核与小型 corrector attribution。
每臂最多 3 个 raw proposals，目标 K<=4 含合格基础图；不强制凑数。
导数准备与 scoring RHS 必须入各臂账；预算不够即记录，不免费借用别臂 Jacobian。
若 proposal / corrector 相同可缓存复用，另报告 standalone cold-start 成本。

在相同小组 raw proposals 上比较旧 corrector 与 C_P，独立归因 corrector 改善；
历史 corrector off-prior 的输出不进入正式 P0 archive。
PC-SP 在同一 C_P 下未优于 multi-start/QD 时，不能把 C_P 的收益归功于新生成器。
此规模只作排错与 Go/No-Go，不作显著优越性或跨分布泛化结论。

## 8. 评价与结果状态

保存 raw -> phase/prior retracted -> corrected -> full validated -> archived 全流程。
分别报告：
- N3 raw Top-1 MAE；N3 修正版 MAE/likelihood；
- 所有 raw 候选与正式有效集合的材料误差，不混合分母；
- valid-set Best-of-1/2/4；旧的包含未合格 N3 的 Best-of-K 另表；
- 两个隐藏 witness 邻域的结构覆盖，阈值预冻结；
- 各 candidate 的 ell、Delta ell、predictive distance、prior status/score；
- 最坏候选质量、结构差异、DC/域外补偿；
- correction retention、每种淘汰原因；
- 首个/第二个合法结构所需 states/RHS、成功发现概率与失败任务；
- parent/pair 级重复，不能把多 noise rows 当独立几何样本。

probability 字段：
posterior_mass=null。
可报告 relative_likelihood 和 prior compatibility；
不是从 archive 次数或 softmax(-ell) 生成“posterior probability”。
near likelihood、prior plausibility、posterior basin mass 必须分开。

允许输出：
SUPPORTED_ALTERNATIVES；
ONE_SUPPORTED_CANDIDATE_AMBIGUITY_UNRESOLVED；
PRIOR_CONFLICT；
PRIOR_RETRACTION_FAILURE；
BORDERLINE_NUMERICAL；
NO_SUPPORTED_CANDIDATE。
有限搜索无第二解不能报告 GLOBAL_UNIQUE。

## 9. A17 边界与小型数学检查

允许用精确有限差分恒等式作为预测/解释：
Delta j = (I - R0 H_h G)^(-1) R0 H_h E0；
Delta y = S Delta j，随后按原 W 白化。
所有 source 共享同一材料，H_h 为总材料变化，包括结构与补偿。
普通 full Jacobian 已含 resolvent，不能删 feedback 制造假的“无 A17”基线。
新增先验/active-face 约束后，不能直接套无约束 paired GN 保证。
没有 uniform nonlinear remainder bound 时，线性筛选失败只能记 unresolved。
本轮不强制加入昂贵 A17 学习/ROM 模块，不以此作为生成器成立前提。

## 10. 交付与停止

输出至少：
SOURCE_MANIFEST.json
PROTOCOL_FROZEN.json
E0_AUDIT.csv
PRIOR_SUPPORT_AUDIT.md
WITNESS_PAIRS.json / arrays（offline-only）
PAIR_NOISE_REPLAY.csv
CANDIDATE_LEDGER.jsonl
COMMON_CORRECTOR_SPEC.md
FOUR_ARM_RESULTS.csv
COST_LEDGER.jsonl
FAILURE_LEDGER.jsonl
E0_E1_E2_DECISION.md
RESEARCH_REPORT_ZH.md

报告第一页回答：
1. 原 A30 P0 下，找到了重要的竞争材料解释吗？证据到哪一层？
2. 旧 T/G 候选主要败在 prior、likelihood、结构还是 corrector？
3. N3 能否在小 MAE 损失下获得物理相容，是否仍是最有支持的单解？
4. 新方法收益属于 prior/corrector，还是生成策略本身？
5. 是继续多解生成、仅继续先验保护的 N3 修正，还是先验/benchmark 必须重做？

如果 E1 没有通过，不执行 E2，不训练更大的网络。
如果 E2 所有方法都只有一个合格结构，报告有限预算下未解决歧义；
不要降低先验门槛或增加随机 seed 直到出现正结果。
不得承诺在后台继续；本轮在现有预算内完成并封账。
