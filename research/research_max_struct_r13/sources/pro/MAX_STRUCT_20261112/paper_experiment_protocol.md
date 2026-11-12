# MAX-STRUCT：可直接交给 Codex 的论文实验协议

版本：2026-11-11 UTC。冻结历史仓库：`ec3a37da6f3ee90beb61e753067d4f639f4118f5`。

本协议基于本轮实际阅读的 `evidence_audit.md`、`theory_audit.md`、`literature_review.md`，以及 `normal_memory/*/{result,followup_result,integer_result}.json` 和 `offline_scores.json`。它规定后续实验与论文组织；**不把尚未执行的独立实验写成结果，不宣称事件算法已优于 G19/DR，也不要求达到某个任意成功率才继续研究。** 本次仅创建本地文档，没有改动仓库或远程设备。

## 1. 第一篇论文应回答什么

建议第一篇聚焦一个问题：**辅助电流表示与反射状态记忆，如何分别改变非线性 Maxwell 搜索能够访问的材料分支？** 主结果候选是有界共享材料本构投影的有限分支事件；表示 profile 与配对实验用于排除替代解释。暂不把一般观测可辨识性、极点学习和完整统计材料选择同时包装成另一套新算法。

工作标题可用：*Finite-Distance Constitutive Branch Changes in Nonlinear Maxwell Inverse Scattering*。如果事件干预没有稳定改善实际成本，论文应如实保留为机制和预测研究，并缩小算法主张；不把负结果藏在附录，也不因为它为负而停止表示机制、独立对象或真实物理研究。

### 1.1 当前可以承认的发现

- G19 重启开始使用 36 个 Fourier 列，后续所有11例都进入完整 Fourier 阶段。因此早期路径破坏与永久表示限制必须分开；最终仍存在固定 DPCS 约束。
- 在11个公开起点中，Fourier 成分通常占当前投影损失能量的大部分，而 retained-data 成分主导投影造成的接收变化。这说明电流能量、可观测性与材料影响是不同量；它尚不是成功预测器。
- AP/DR 相同投影与相同审计机会的五个候选池中，离线最优 AP proposal 也不能复制 DR 已产生的好材料。这个证据支持候选生成动力学差异，但没有隔离 G19 与 DR 的全部差别。
- 原 paired 实验给 G19 的物理候选审计主要在阶段边界，AP/DR 每50步一次。B 组 G19 被预算中断的最后阶段没有进入候选池，这是独立于表示空间的有限预算混杂因素。
- 两个 NN 例子在正确器后仍数据可行，却材料更差。最小 prior 的噪声球内选择也不能自动修复该问题。

这些结论的数字、源文件和条件以 `evidence_audit.md` 为准；不得把11例或新增三例追踪称为新盲测。

### 1.2 三个事件追踪到底支持多少

从相同已保存 AP3000 状态出发，本轮3000步续算如下。`r/δ` 指已付候选中原始物理 Φ 最小者；不是最后投影残差。

| 公开实例 | 首射线事件 t* | ordinary 首上下界翻转迭代 k（k=0为初始） | ordinary 最佳 r/δ | 连续事件 proposal 最佳 r/δ | 整数 proposal 最佳 r/δ |
|---|---:|---:|---:|---:|---:|
| A shell0004 | 5.73172 | 5 | 0.95770 | 0.97309 | 0.95514 |
| B0010 | 12.49514 | 12 | 30.07966 | 30.07966 | 30.07966 |
| C NN0001 | 14.07906 | 14 | 24.09591 | 0.92964 | 0.92975 |

离线材料评分中，shell ordinary/连续event NRMSE约为0.18708/0.17776；NN0001的ordinary、连续event、整数event分别为0.3348555、0.1331813、**0.128761149**；B0010都保留原始坏候选。既有published DR在该NN实例的raw/primary为0.107522/0.077988。这里raw是在原NN起点直接运行DR后的原始物理候选，primary还包含既有corrector/选择；本轮event是AP3000准备后的续算且没有corrector。因此只能说明尚未超过已有候选，不能当作同起点同成本的正式优劣检验。连续与整数 proposal 都可能改变 basin，因为 AP 对并不精确；整数 t 本身不会赋予“exact skip”语义。

续算计时只含各分支新增工作，不含共同 AP 准备和历史投影 setup；旧 `result.json` 的 CPU 又包含 AP3000与两条500步续算，**不能把整个数值冒充 AP-only setup**。后续竞争实验必须重新分项计时，而不是据这些数字宣称提速。

## 2. 数学主张与适用范围

理论以实、非负、有界、各向同性材料为主。标量 TM 的每体素矩阵为2×2L；真实三维各向同性 Maxwell 为2×6L，跨照明与三场分量共享同一个实材料。未知复损耗和张量材料需要新的投影几何，不能沿用同一个角度公式而不证明。

论文理论层次如下。

| 层次 | 可以写的命题 | 不可替代的验证/边界 |
|---|---|---|
| 基础引理 | 精确 resolvent 导数与有限差分；二次 profile 配方；DR法向漂移 | 归入已有工具，不宣称全新Maxwell恒等式 |
| 材料事件命题 | 对 c=PCa、a=PAc、d=a−c≠0，明确共享材料锥保持旧shadow的最大法向区间 | 区分实内点奇异值交叉与盒边界全局分支打平 |
| 扰动命题 | 近似AP对的一步射线误差≤(t+1)ε，及保持同分支时累计误差 | 逼近投影退化点时放大因子可发散；不能给所有真实轨迹精确步数 |
| Maxwell源—接收桥接 | 精确AP对的法向由接收残差伴随场产生；非零gap要求数据约束活跃 | 实材料、固定尺度、可逆resolvent、Slater与正确λ归一化；不是新的representer theorem |
| 精确跳步 | 在精确AP对和未变投影分支内，直接到下一个整数射线状态 | 只减少冗余投影；不改变普通DR选择的材料分支 |
| 近似事件proposal | 在真实近似状态沿同方向移动再继续DR | 这是可失败的有限扰动算法；必须测材料与实际成本 |
| profile机制 | 受限current目标=full-current目标+材料相关表示惩罚 | 旧G19双β/两块不对应单一共同目标；不得称为整个旧G19轨迹的精确Hessian |

特别保留 d=0 的错误数据可行材料作为负控制：纯可行性反射没有提供自动离开它的机制。噪声球内部数据约束局部不活跃，因此可行集切向角不能代替统计材料可辨识性。

## 3. 原创性与2026仿射DR外推近邻

必须正面引用 NFFT-SOM2018、Elser2007/2017、变量投影和 van Leeuwen–Yang2025，具体书目见 `bibliography.tsv`。解除 noisy DPCS、一般 discrepancy accumulation、SVD rank-one 投影、扩展场消元都不是首次发现。

本轮再次用**原作者 arXiv 记录的检索摘要**核查 Collard–Lindstrom，*Finite Recovery by Extrapolation for the Douglas–Rachford Method on Affine Subspaces*，2026-11-02，https://arxiv.org/abs/2611.02694 。它研究仿射子空间 DR 的 principal-angle 搜索维数，并用 LSM/least-squares extrapolation/GMRES 的仿射搜索空间在有限步恢复 DR 极限。PDF和HTML直接访问仍返回 DisabledError；本协议没有依赖二手AI解读或声称读过其全文。

**判断：**仅据原作者摘要，它不直接覆盖非凸共享材料锥从非零最近点对到首分支事件的阈值，也不等同于精确AP法向线段上的整数跳步。但它足以否定泛泛的“首次提出DR有限外推/精确跳步加速”主张。投稿前需获得全文核对；在此之前主新意只能落在材料盒约束事件、扰动范围、可预测性和实际算法效用，不能落在“外推”这个词。主来源为该段所列原作者 arXiv 记录。

## 4. 实验模块与唯一需要回答的因果问题

### M1：首分支事件的独特预测

对象是所有11个已公开paired起点，而不只三个正反案例。固定同一 G/E/S/Y、场电流尺度、材料盒和投影实现。现有三个状态直接复用；其余状态只按预先一致的准备方案生成一次共享AP序列，并保存100、300、1000、3000步的状态。准备成本单列，后续每个从该状态开始的独立部署总成本都要计入；离线机制分支之间共享准备只作为“新增算力”账，不能当冷启动成本。

每个快照在运行未来轨迹前写出 t*、首体素、活跃边界、ε/gap、射线投影裕度和预测适用标记。随后从同状态比较普通DR、AP、整数事件proposal、连续事件proposal。所有分支都按相同规则审计，不根据GT选择哪一状态/哪一事件值得运行。

评价分为三个不同结果：

1. 规定法向射线上的投影分支改变：这是公式结构验证，仅需少量点，不把机器精度恒等式当主要实验。
2. 实际DR轨迹上首大材料翻转的时间和位置：这是主要预测实验。记录所有微小移动，区分“任何材料改变”与“边界到边界改变”，不得事后改事件定义。
3. 到达数据可行候选和材料准确候选的成本：这是算法问题，不能由前两个结果代替。

对非整数 t*，从 z₀=a=c+d 出发的精确参考事件迭代是 floor(t*)；整数打平时报告 projector tie rule。实际误差与ε/gap、接近退化程度一起画图，不用GT拟合一个漂亮但无验证的阈值。

#### M1新增预注册诊断：源照明与接收数据如何形成法向gap

依据 `theory_audit.md` §12，定义 \(A_m=I-D_mG\)、\(r=W(sSa_X-Y)\)，数据球使用半平方约束 \(\tfrac12(\|r\|^2-\delta^2)\le0\)，其KKT乘子记作λ。对于**实各向同性材料、固定s、精确AP最近点对、A_m可逆且数据球满足Slater条件**，

\[
p=A_m^{-*}S^*W^*r,\qquad
d_X=-\lambda s p,\quad d_Q=\lambda D_m p.
\]

因此非零gap必须由正λ和活跃数据边界支持；内点体素的首事件可写成

\[
t_{*,i}=\frac{\|Q_{c,i}\|}{\lambda s\|p_i\|},
\quad \Re\langle Q_{c,i},p_i\rangle=0.
\]

盒边界使用原来的带交叉项公式，不能套内点比值。这里Q_c是辅助本构场，不一定等于物理总场Q(m)；r是辅助传播/数据状态的残差。它解释法向如何同时依赖照明、完整resolvent、接收覆盖、噪声度量与数据压力，不增加观测信息，也不预测GT质量。

在M1和既定孔径/照明/噪声压力实验中预先保存：数据球活跃裕度、正确归一化λ、源场行范数、接收伴随场行范数、交叉项、相应KKT残差与理论条件状态。检验改变接收/照明后**新AP状态**的实际法向和事件如何变化；不预设扩大δ必然增大/减小t*，因为状态、λ和p会一起变化。

实现不建立 \(H_m=WSA_m^{-1}D_{s^2+m^2}A_m^{-*}S^*W^*\) 或任何新巨大矩阵。优先复用投影已经得到的λ，使用一个现有G*和S*作用评估 \(A_m^*d_X+\lambda sS^*W^*r\) 的缺陷。只在预声明的小型诊断子集上付一次多照明adjoint block solve得到独立p；费用写入`event`/`diagnostic`，并记录RHS和Krylov工作。若直接取 \(p_{\rm proxy}=-d_X/(\lambda s)\)，它只是状态重表达，不能把两种t*一致当作新的独立实验证据。λ很小时保留 \(\lambda p\) 而不做不稳定除法。

真实近AP状态必须记录ε、投影误差和λ来自哪个投影调用。`PA(c)`返回的a′可能不同于保存的a；不能把a′的乘子和a的残差混用后宣称精确恒等式成立。δ=0时使用独立的数据等式乘子；Slater未验证/失效、resolvent奇异或数值近奇异的实例单独标注，不适用正标量λ结论。原式是熟悉的adjoint/extended-source结构在此几何中的桥接，不作为新的representer定理。本轮只补协议，不增加实验或把未计算诊断写成已验证结果。

### M2：反射记忆与大扰动的区别

在M1相同状态加入 `dr_reset` 与 `normal_random_control`。`dr_reset` 保留材料/场接口但清掉反射记忆；随机控制先构造同法向维度、同范数、同全局标量预算的扰动，并用冻结随机种子。另加不使用 t* 的同范数随机方向以及简单倍增射线搜索，分别排除“任何大步都行”和“只需要找一个足够大的t”。所有拒绝proposal和为选择t付出的投影都计费。

允许做小型状态/尺度敏感性：s取原值的1/2、1、2；材料上界取原上界及明确更宽的物理允许值。若只在贴近GT上界时有效，应解释为先验与边界效应。不同s改变投影几何，这是方法设计因素，不应解释为物理信息增加。

### M3：完整电流表示移植到相同G19动力学

不得只运行原始LSIE-NFFT，再与当前G19/DR比较后声称分离了表示。建立统一 current adapter，保留相同G19材料分母、完整导数、current块目标、PR更新、线搜索和β日程，仅改 admissible current family。

| arm_id | current表示 | 用途 |
|---|---|---|
| `g19_legacy` | J_D+(I−P_D)F_mf a，原3→6→10→12日程 | 精确历史基线 |
| `g19_full_apcs` | J_D+(I−P_D)a，从第一步完整APCS | 区分早期Fourier路径与永久DPCS约束 |
| `g19_nfft_full` | J_D+F_full a，等价完整J；照原块目标优化 | NFFT空间移植；直接挑战“解锁为新发明” |
| `g19_soft_low` | J_D+V_D b+(I−P_D)F_mf a，仍用原加权数据残差 | 只释放retained坐标，隔离显式噪声球的额外作用 |
| `g19_elastic_low` | J_D+V_D b+(I−P_D)F_mf a，显式全局噪声预算 | 测少量retained-mode校正是否足以获得收益 |
| `g19_preserve_anchor` | 使用R8已存在的physical-complement方案 | 排除重新发明已失败的warm-start修补 |

`g19_nfft_full` 保留J_D作坐标锚并不意味着锁定它。full F采用单位正交归一化，使坐标更换不偷偷改变Euclidean梯度/PR度量；为节省开销可直接优化δJ，但必须确认使用相同目标和正确链式梯度。与full APCS的未知量差是 retained coordinates，不是又改变G19材料式。`g19_elastic_low` 额外改变了约束，因此必须与`g19_soft_low`一起看，不能把其全部收益归因于坐标释放。

初始化另作一个小因子：原合法投影 vs 保留相同物理J（仅在该表示允许时）。无法完全保留的arm必须报告实际ΔJ而不是宣称“同完整状态”。“同材料起点”和“同场电流起点”分别作为两类配对。

### M4：profile偏置的有限材料后果

固定公开起点、G19/DR已保存好坏候选及它们之间预声明的材料插值。离线求指定二次current块的精确/高精度profile，分别报告full、full APCS、early Fourier三个值与差。只在这一小诊断集合付高精度消元费用，不能把它偷偷用于在线算法。

主要问题是是否出现**材料排序翻转**：相同m₁,m₂在full profile下的排序是否被表示惩罚改变，以及廉价指标能否预测这种变化。加入GMRES/Galerkin误差估计，避免把profile未充分求解当作偏置。真实物理目标Φ、固定current块profile和旧G19迭代使用的两个块分别命名。

若做物理切向补空间，先用一个预声明方向h而不是大量方向筛选；比较同秩随机、Fourier增加、Green奇异向量和材料相关方向。方向的求解、SVD与存储都计费，多照明补全可能是L列而不是“一列”。这部分是机制控制，尚不需要发展成另一篇复杂算法。

### M5：审计机会与corrector

先在原起点比较原阶段边界审计与统一审计。**即使不改G19，增加阶段内物理候选可见性也可能改变结果；必须先量化。** 对B组应保存进入mf12后的中间材料，解决原预算结束前最后阶段未审计的问题。

corrector单独分配预算，使用各方法已经选出的原始物理候选，不重新运行搜索。比较无corrector、原Φ corrector、显式噪声球约束的prior corrector、独立验证观测选择。它们回答的是不同估计器问题；不宣称最后一个在没有额外观测时免费可得。

记录 raw→corrected→selected 的材料距离、数据预测距离/δ、D和R变化、是否出入噪声球及离线材料损伤。保留C0004/C0022反例。原始候选永远保存在候选池；“不覆盖它”不等于有能力自动知道它更接近GT。

## 5. 成本与物理候选审计协议

### 5.1 两种比较估计量，不混成一种公平性

**机制比较A：相同候选审计机会。** 每条arm使用16个proposal工作窗口，每窗口分配相同新增优化CPU预算；在窗口末审计当前材料，加共同initial，共17个机会。固定阶段结束仍保存辅助状态，但不额外开放物理候选选择。若方法提前正常终止，保存终点并停止优化；后续重复材料可使用统一缓存，不强制重复昂贵正演。记录17个请求槽、实际唯一材料数与真实新增正演数。

窗口CPU不含物理审计，因此此设计匹配生成工作和审计机会，**不保证总CPU恰好相同**。昂贵材料的审计确实要付更多；最终以实际费用画性能曲线，并另做比较B，不能把相同窗口说成相同总成本。

**部署比较B：相同总CPU上限。** 总预算包括初始化、setup、优化、事件诊断、候选审计、选择及corrector。所有方法在同一预声明的累计CPU时间点请求候选；协作式停止在已完成的安全操作边界。预算将尽时不做“免费最终正演”，使用最后一个已付可审核候选。若一个不可中断的solve越过上限，记录超额与实际CPU，并按实际时间归属；不得丢掉失败或超额运行。

对两种比较均同时保存生产式选择（所有本方法已付可用状态）和共同槽位选择（只看匹配槽）；前者比较真实算法，后者隔离候选曝光。direct方法已经算过的物理状态可命中统一cache，但不可凭空多给其他方法免费正演。

### 5.2 记账口径

每个run报告独立部署成本与本次新增成本。共同G/S/SVD/投影factorization、NN推断、AP准备、候选物理审计、corrector都单列；若给amortized成本，必须同时给一次性总成本和摊销实例数。NN训练样本正演、训练GPU小时与推断分开，冻结权重不能写为训练成本0。

必要计数：G/G*作用的call与RHS数；S/S*同理；物理primal/adjoint请求与实际RHS；LU/Cholesky/eigh/SVD次数；Krylov迭代；projection A/C次数；噪声球乘子求值；被拒绝trial；peak RSS/GPU VRAM；CPU、wall与GPU同步后的kernel时间。缺失填null并说明，不能填0。

固定线程数、核心绑定和数值精度；以case为block随机方法执行顺序。预先抽取10%任务做重复计时，仅估计系统波动，不将最好时间用于结果。不同GPU kernel/布局/正交化差异必须写明，不能把算法和硬件优化混成一个速度比。

## 6. 独立2D确认集与成熟基线

### 6.1 数据划分

保留全部R4–R12与本轮状态为开发/机制集合。新生成96个材料对象，8个生成族每族12个；每对象的几何参数、连续材料函数和随机种子一次生成冻结。先生成manifest和hash，再运行算法，最终封存所有预测后由单独scorer读取GT。这个过程是新的独立确认；如果开发者提前看了结果，该批即标记exposed，保留数据并再建独立确认批，不能重新命名为blind。

8族覆盖：不对称多连通支撑；多尺度分离夹杂；细缝与近接触结构；弯曲细杆/桥；连续梯度材料；多层非同心介质；多孔随机介质；无明显分段边界的平滑随机场。与R8神经训练族重合的类别必须标注，不能仅靠换随机种子宣称结构泛化。

主设计使用32²反演，分别生成matched32²与cross-grid64²观测；与几何无关的固定12个照明和48个接收器。contrast按预声明分层覆盖低、中、高（例如χ最大值2–5、5–12、12–20），每族各档均衡；3%圆对称复噪声作为主条件，独立噪声seed按object配对。低contrast不删除，它是检查新机制是否破坏已正确路径的必要部分。

全部96对象使用cold与冻结NN两种初始材料。NN前向只做一次并复用，但每种部署成本都计该推断成本；同时报告NN-alone。主9算法为 CIE、Dual-CIE、GlobalDual、A+、G19、移植NFFT空间的G19、DR、事件proposal-DR、direct。共96×2网格×2初值×9=3456个run；不做所有因素的巨大笛卡尔积。

### 6.2 经典竞争者与压力设计

从manifest预先按每族3个选24对象，加入MR-CSI、CC-CSI与完整场ADMM/WRI类实现；与已冻结的主方法进行同成本比较。引用并尽量用作者实现，但统一单位、源接收几何、初值、材料盒、物理目标验算与成本账。不能用很差的自写裸CSI替代成熟正则化CSI。

另从manifest预先抽24对象进行平衡压力设计：1%/10%噪声、半孔径、照明数减半/增倍，以及材料损耗。每个对象只承担预声明的一部分组合，保证每条件包含各族，不用全排列增加重复计算。复损耗作为**算法外推与失败范围**，只有完成相应材料投影推导的算法才使用事件式；否则比较普通DR/G19/direct/CSI，并明确实材料定理未覆盖它。

在另外16个预声明对象上使用48²反演/96²生成检查分辨率稳定性；不能用放大像素图替代跨网格。再保留24个生成族/物理设置不同的额外分布外对象，用于最终冻结后的算法选择/预测指标确认。

### 6.3 统计与指标

材料指标同时报告旧amplitude-NRMSE、relative L2、支撑IoU和连续材料梯度误差，定义直接写入schema；对连续材料不强行二值IoU。以object为聚类单位做配对bootstrap置信区间，noise重复和多初值不能当独立对象增加样本量。

主要算法指标为相同实际成本下的数据可行率、材料误差分布、额外成功与被破坏成功的配对数、到噪声球的时间及失败类型。事件预测报告时间误差、首体素命中/邻域命中、误差区间覆盖率和预测失效频率；先定义邻域物理距离，不能看到误差后扩大。

ROC/AUC只用于开发阶段冻结的预测器；若只是探索散点相关，标注retrospective。对未在预算内到达噪声球的时间使用右删失语义，不能将其删除后只平均成功者时间。

## 7. 原始NIE/CIE基准：复现与新确认分开

NIE2016已成功访问共同作者上传的原始全文，已读推导与相关实验部分：https://www.researchgate.net/publication/296475773_A_New_Integral_Equation_Method_to_Solve_Highly_Nonlinear_Inverse_Scattering_Problems 。复现时沿用原2m域、20平面波/40接收、150²生成/64²反演、L=15及声明噪声。至少包括Austria的400MHz、εr=3.5与4，750MHz、εr=2；保留原日程和停止准则，再另外做成本匹配对照。原始FoamTwinDiel是独立实测复现项，不能变成合成替代数据。

CIE2019三维原配置按原文及R9已有合同/实现继续，不从零重复构造。先阅读：

- `research/research_max_struct_r9/delegated/original3d_feasibility/{REPRODUCE_ORIGINAL_CM,ORIGINAL_CM_NINE_REPORT_ZH}.md`
- `research/research_max_struct_r9/delegated/independent_original3d_contract_review/REVIEW_ZH.md`
- `research/research_max_struct_r9/delegated/3d_extension/physical_refinement/REPORT_ZH.md`

R9已记录CM+RR点DDA在原尺寸60³→30³的显著响应差，约0.535，超过声明0.1噪声；这不等于新算法失败，也不能用matched-grid恢复冒充连续物理验证。保留已有求解与费用，修正离散/自项/积分单元模型，同时继续明确标注的matched-grid机制实验。不要把这项历史物理问题变成停止整个研究的性能gate。

## 8. 真正三维矢量Maxwell与实测

### 8.1 新三维对象

计划24个新的三维介质体，覆盖非轴对称多连通壳、偏心夹杂、倾斜细桥、分层损耗与连续梯度；不是2D轮廓沿z挤出。主理论确认先用实各向同性材料，24个完整矢量照明/极化组合和三分量接收。统一24³或32³反演，与更细网格或独立求解器生成数据；具体网格由预先的解析体精度/资源测量决定，不能看算法材料结果再选择。

三维主比较保持精简：G19、G19-NFFT、DR、事件proposal-DR，加direct或field-ADMM在预声明12对象子集。先有相同物理模型的matched运行，再做cross-grid；每对象至少保留一个不同离散模型的外部预测验证。球/同心球用Mie，非球用独立边元Maxwell或可靠VIE实现。R9已有Mie验证与失败记录必须复用，禁止只重测已知通过项。

三维A投影不得建立巨大dense H。使用FFT dyadic Green作用、matrix-free SPD子问题、共同预条件器与热启动。每次投影记录残差和迭代次数；在可负担小3D对照上量化投影误差对事件预测的影响。近似投影改变理论条件，要给出数值误差账，不能用“CG收敛”四字代替。

### 8.2 Fresnel真实数据

先复用R8 `fresnel_adapter.py`、`fresnel_fetch.py`与已保存原始数据/hash及照明诊断；阅读 `research_max_struct_r8/delegated/fresnel/PRIOR_2005_DATA_AUDIT.md`。2001/2005公开对象已经广为人知，称为real-data validation而非新blind物体。2D至少包括单/双介质柱与FoamTwinDiel等不均匀介质。PEC对象不能简单塞进有限介电盒当同一个材料模型。

每频率保存原始场、单位、源/接收坐标、相位时间约定、背景扣除和校准参数来源。照明校准尽量使用独立参考数据；若只能在目标数据上估源复振幅/相位，将其列为共同 nuisance 模型和成本，并报告它可能吸收模型误差。不能逐方法调校准，也不能用已知目标材料挑频率。

预声明一组所有目标共有的频率/视角用于重建，另一组完整视角或频率用于外部预测；校准和选择不能窥视holdout。holdout提升提供额外观测证据，不能说在原相同信息量下已解决非唯一性。

三维Fresnel使用官方 https://www.fresnel.fr/3Ddatabase/ 与参数记录，先处理实际平面波照明、极化、几何遮挡和频率相关SNR。至少4个公开目标、每目标3个预声明频率块；报告数据预测、几何/材料误差的已知范围及校准不确定性。若材料标称值本身不确定，不用它生成看似高精度的唯一NRMSE。

## 9. 现实资源估计与扩展方式

以下是容量规划，不是已测速度或预算授权。先在固定公开校准实例测 kernel/forward/projector 成本，写入manifest；只更新算力估计，不据GT效果删方法。以现有paired全进程成本和本轮24²增量测量为量级参考，两者不可直接换算为32²/3D运行时间。

| 工作 | 建议资源 | 可复算的工作量/预估 |
|---|---|---|
| 11例机制与表示因子 | 16–32 CPU核，64–128GB RAM | 数百条24²/小型profile任务；预留30–100 CPU小时，所有准备/精确消元计入 |
| 96例主确认 | 同上，固定每进程1线程 | 3456run；若32²总上限600 CPU秒，每完整预算总计576 CPU小时；实际提前完成后如实少计 |
| 经典/压力/48²/额外分布外 | 32–64 CPU核；必要时一张GPU | 再预留约300–700 CPU小时；取决于高精度forward，不把上限乘积冒充实耗 |
| 原论文150²与30³/60³复现 | FFT/迭代求解；128GB主存，24–48GB GPU优先 | 复用已付数据/状态；单列新核验和缺失算法费用，不能重做整套历史账 |
| 新矢量3D | 128–256GB主存，24–48GB GPU；逐例批处理 | 约100–150条主/控制run；先预留100–500 GPU小时，使用下式按实测校准 |
| Fresnel2D/3D | 同一物理kernel；200GB以上结果空间 | 2D数十CPU小时至数百；3D约20–100 GPU小时，校准与重复求解另列 |

三维成本预测使用

\[
T_{\rm run}\approx T_{\rm setup}+N_A\bar k_H C_H(N,L)
+N_G C_G(N,L)+N_{\rm audits}\bar k_F C_F(N,L)
+T_{\rm selection/correction},
\]

其中每个C都由同硬件计时；禁止把一个GPU Green作用当一次完整forward。复杂度与实测时间同时报告。以complex128为例，30³矢量dense矩阵有81000²元素，**仅一份约105GB**；60³约6.7TB，未计factorization。因此3D主实现必须matrix-free，而不是增加内存硬扛。

样本规模是为了覆盖形态/物理因素与配对不确定性，非“通过某gate才能扩大”。若某机制失败，仍完成独立确认所需的负控制并调整理论解释；若算力不足，预先按整object block缩减或排队，保留随机化与所有失败，不按恢复好坏删case。

## 10. Codex实施任务与接口

建议新增独立研究目录 `research/research_max_struct_mechanisms/`，不覆盖历史R4–R12。以下是**待实现接口**，不是声称现有脚本已有这些参数。

1. `freeze_manifest.py`：固定repo SHA、物理模型、对象/噪声/初值种子、算法配置、预算、审计槽、硬件信息和truth隔离约定。
2. `current_adapter.py`：legacy Fourier、full APCS、NFFT full、elastic DPCS；只改变表示，保存每次初值投影误差。
3. `candidate_scheduler.py`：共同槽位、实际CPU、统一cache、final audit与中断语义；修复stage-end-only混杂。
4. `event_geometry.py`：实材料内点/端点阈值、tie rule、全局首事件、近AP误差与适用标签；只用当前状态。
5. `run_mechanism.py`：共享状态fork、普通/事件/重置/随机控制和精确分项成本。
6. `profile_audit.py`：固定二次块的离线full/restricted消元、求解误差和排序翻转；不提供GT给在线算法。
7. `run_confirm.py`：冻结确认集批处理、失败恢复和原子checkpoint，输出方法无关的已付候选接口。
8. `run_corrector.py` 与 `score_sealed.py`：分别执行估计器和离线GT评估，禁止通过分数再改本批配置。
9. `make_paper_figures.py`：只消费冻结输出，生成预声明主图及完整失败附表。

配置示例：

```yaml
schema_version: max_struct_mechanisms_v1
base_commit: ec3a37da6f3ee90beb61e753067d4f639f4118f5
truth_policy: scorer_only_after_predictions_sealed
physics:
  material: real_isotropic
  time_convention: repository_declared
  inverse_grid: [32, 32]
  forward_grids: [[32, 32], [64, 64]]
  n_illuminations: 12
  n_receivers: 48
  noise_relative_rms: 0.03
  discrepancy_definition: frozen_complex_noise_model
comparison:
  proposal_windows: 16
  include_initial_candidate: true
  cache_policy: identical_material_and_physics_hash
  deployment_total_cpu_seconds: 600
  corrector_cost_included: true
  no_free_terminal_forward: true
event:
  mode: approximate_proposal
  integer_parameter: floor_tstar_plus_one
  independent_cell_scaling: false
  retention: preserve_all_paid_physical_candidates
outputs:
  store_full_state: [initial, event_before, event_after, terminal]
  store_material_at_every_audit: true
  failed_run_policy: retain_with_reason
```

生产式event prototype必须与机制fork分开登记。可先固定少量检查点读取普通DR现有状态，构造a=PA PCz、c=PCa并作一次全局scalar proposal；新增投影和任何AP平滑都付费。并行保留一个不跳的分支会增加计算，不能计作单路径免费保险。主确认前冻结该规则；目前的“先AP3000再event”仅是机制准备协议，不能直接命名为低成本生产算法。

## 11. 精确输出schema

所有记录包含 `schema_version, experiment_id, case_id, run_id, arm_id, config_hash, base_commit`。UTF-8；数值不写字符串；不可知用null；复杂数组用NPZ或HDF5并在索引记录shape/dtype/轴含义/物理单位/sha256。

### 11.1 `case_manifest.jsonl`，每物理实例一行

```json
{
  "case_id": "confirm2d_family03_0007",
  "split": "independent_confirmation",
  "exposure_status": "unexposed_at_freeze",
  "geometry_family": "curved_thin_structures",
  "object_seed": 202611110307,
  "noise_seed": 202611120307,
  "physics_id": "tm32_cross64",
  "inverse_grid": [32,32], "forward_grid": [64,64],
  "material_parameter": "contrast_chi",
  "material_domain": "real_isotropic_box",
  "material_lower": 0.0, "material_upper": 20.0,
  "frequency_hz": null, "length_unit": "metre",
  "field_time_convention": "exp_minus_i_omega_t",
  "Y_units": "declared_by_forward_model",
  "whitening_array": "arrays/W.npz",
  "delta_absolute": null, "delta_relative_to_Ynorm": null,
  "delta_confidence": 0.99,
  "source_receiver_manifest": "arrays/acquisition.json",
  "truth_hash": "sha256", "observations_hash": "sha256",
  "forward_solver_id": "frozen_solver_and_discretization",
  "calibration_id": null
}
```

`frequency_hz`与δ在真实manifest中必须填数值或有明确不适用原因；模板null不是允许把缺失物理参数默认为0。使用仓库归一化噪声半径时必须同时存绝对值，防止本轮审计已经发现的单位混淆。

### 11.2 `runs.jsonl`，每run一行

字段：`initial_material_hash, initial_state_hash, initialization_id, current_representation, dynamics, beta_schedule, scale, material_box, audit_schedule_id, corrector_id, total_budget_cpu_s, terminal_reason, numerical_status, requested_audit_slots, completed_audit_slots, unique_physical_candidates, primary_candidate_id, raw_candidate_id, corrected_candidate_id, secondary_candidate_id, sealed_at_utc, truth_was_read`。在线run必须 `truth_was_read=false`。

`terminal_reason`枚举：`budget, stationary_rule, scheduled_steps, numerical_failure, user_interrupt, completed`。算法失败与基础设施失败分别记录，不能在统计时默默删掉。

### 11.3 `candidate_registry.csv` 与 `physical_audits.csv`

每候选主键`candidate_id`；字段：`run_id, parent_candidate_id, phase(initial/proposal/corrector), proposal_index, audit_slot, algorithm_iteration, stage_id, material_array, state_array, material_hash, eligible_common_slots, eligible_production, source_of_candidate, created_cpu_s, audit_finished_cpu_s, cache_hit, audit_cost_id`。

物理审计字段：`candidate_id, Phi, data_term_D, prior_R, residual_absolute, residual_over_delta, predicted_data_array, forward_solver_status, forward_relative_residual, primal_rhs_new, adjoint_rhs_new, memory_peak_bytes, audit_cpu_s, audit_wall_s`。辅助目标和物理Φ必须不同列。

### 11.4 `events.jsonl`

字段：`event_id, run_id, state_id, iteration, gap_norm, AP_pair_error, AP_pair_error_over_gap, t_star, predicted_cell, predicted_cell_coords, initial_material, initial_active_status, a_energy, b_energy, gamma, angular_interval, projection_margin, predicted_integer_event, tie_policy, prediction_type(exact_reference/approximate), numerical_projection_error, observed_first_any_material_change, observed_first_box_switch, observed_predicted_cell_switch, normal_ray_deviation_at_switch, event_compute_cpu_s, event_projection_count, jump_parameter, proposal_state_id`。

源—接收诊断追加字段：`lambda_value, lambda_constraint_normalization, lambda_projection_call_id, lambda_state_id, residual_state_id, data_ball_active_margin, slater_status, slater_witness_id, resolvent_regular_status, adjoint_field_mode(independent_solve/proxy/scaled_proxy/not_computed), adjoint_field_array, scaled_adjoint_array, source_field_row_norm_array, receiver_adjoint_row_norm_array, source_receiver_cross_term_array, source_receiver_kkt_defect, constitutive_amplitude_normal_defect, source_receiver_tstar_array, adjoint_solve_relative_residual, diagnostic_new_adjoint_rhs, diagnostic_Gadj_rhs, diagnostic_cpu_s, applicability_reason`。`slater_status`枚举`verified/unverified/failed/not_applicable_equality_data`；未验证不是默认为满足。数值求解成功也不能冒充一般解析可逆性证明。

阈值无限写 `t_star=null, t_star_status=no_normal_component`；不把Infinity写入非标准JSON。阈值与未来观察字段分开写入两个sealed阶段，以证明预测没有使用未来轨迹。

### 11.5 `cost_ledger.csv`

字段：`cost_id, run_id, parent_cost_id, shared_group_id, category(setup/init/proposal/event/diagnostic/audit/corrector/training/serialization), incremental_or_deployment, cpu_s, wall_s, gpu_s, primal_rhs, adjoint_rhs, G_calls, G_rhs, Gadj_calls, Gadj_rhs, S_calls, Sadj_calls, LU, Cholesky, SVD, eig, A_projections, C_projections, krylov_iterations, multiplier_evaluations, rejected_trials, peak_rss_bytes, peak_vram_bytes, completeness, missing_reason`。

`parent_cost_id`用于防止重复把子计时加到已含它的整进程成本。输出一个自动核算表显示process total与分项和差；该差是未分类成本，不强行归零。

### 11.6 `profile_diagnostics.csv` 与 `offline_scores.csv`

profile字段：`material_id, block_objective_definition, representation_id, full_profile, restricted_profile, representation_penalty, solve_residual, lower_bound, upper_bound, derivative_direction_id, directional_penalty_derivative, ordering_pair_id, ordering_flipped, diagnostic_cpu_s`。

offline评分字段：`candidate_id, metric_definition_version, amplitude_nrmse, relative_l2, support_iou, continuous_gradient_error, heldout_data_residual, calibration_uncertainty_label, material_success_definition, material_success, failure_class, truth_accessed_at_utc`。

`failure_class`允许多标签：`no_data_feasibility, feasible_material_error, prior_selection_damage, representation_bias_evidence, solver_access_failure_evidence, model_mismatch, numerical_failure, unresolved`。仅低残差不能自动赋`nonidentifiable`；要有实际候选/观测证据。对于实测缺少可靠GT的指标填null。

## 12. 决定性图表，而非平均NRMSE排行榜

1. **材料锥有限事件图。** 同一真实状态的投影能量随全局t变化，内点与端点两种竞争；显示阈值、观测首翻转、近AP误差，旁列B0010失败。增加源场/接收伴随场与活跃数据压力小面板，但只使用实际计算并注明proxy或独立solve的诊断。不要只放一张成功曲线。
2. **预报与实际的独立散点。** 横轴预测首事件、纵轴真实DR首大翻转；颜色为AP误差/事件裕度；标出未发生、错位置与退化点。与仅gap、Green谱、随机标量预测比较。
3. **表示×动力学×审计矩阵。** 同材料起点和同完整状态两种配对，标出每格可用候选数量与实际CPU。用于回答“是否只是fullcurrent/更多审计”。
4. **profile材料排序翻转图。** full与restricted的物理候选排序、DPCS与Fourier惩罚分开；给出能预测与不能预测的例子。
5. **已付候选云与corrector箭头。** 横轴r/δ、纵轴材料误差，候选生成/选择/正确器分色；显示C0022类在噪声球内变差，而不是只展示最终图像。
6. **成本—恢复与配对损益。** 全CPU、wall、forward/adjoint、内存四种成本；新增成功和破坏成功同时列。失败包含在分母里。
7. **独立物理验证。** 真3D三分量/极化、cross-grid与独立solver响应差、Fresnel训练观测和holdout观测；给校准误差带。

主表至少包括：算法/空间/目标/记忆/审计方式映射；新旧数据暴露状态；定理条件与反例；逐族配对结果；实际成本账。将完整每对象结果提供为机器可读附件。

## 13. 最强替代解释与对应反证

| 替代解释 | 必须如何检验 | 若成立应如何改写结论 |
|---|---|---|
| DR只是完整current空间更宽 | G19-NFFT/fullAPCS/elastic，保持G19材料块不变 | 收益主要属于表示；不归功反射 |
| G19好中间材料没被审计 | 统一阶段内paid审计和共同槽位选择 | 收益部分属于readout机会；修正算法pipeline |
| event只是大扰动/重置 | 同状态reset、同幅随机、简单倍增外推 | 缩小为有用的初始化启发式，不能称几何独有 |
| 事件只利用过紧材料盒 | 上界与尺度敏感性及连续材料 | 解释为先验驱动边界切换，不泛称Maxwell恢复定理 |
| AP准备本身支付了隐藏算力 | 全流程冷启动成本；在线现有状态检测 | 没有全成本优势就不宣称加速 |
| 低残差选错材料来自prior | 多候选同球、corrector分解、独立holdout | 区分材料选择和搜索；报告不可判定 |
| crossgrid误差比噪声大 | 独立forward/Mie与校准误差 | 限定matched离散模型结论，继续修正物理 |
| 为已见成功例设计机制 | 全部11例、冻结96+24新对象 | 若独立预测失效，保留反例并修正适用域 |

## 14. TAP / SIAM JIS / Inverse Problems 的现实定位

| 方向 | 第一篇可承载的核心 | 还必须补足的证据 | 当前不能声称 |
|---|---|---|---|
| 强TAP | 同实际成本改善高对比矢量Maxwell恢复；物理清楚的有限事件/表示控制 | 3D真矢量、损耗边界、原CIE基准、Fresnel可靠校准；成熟CSI/ADMM对照 | 三个公开2D续算已证明通用优越 |
| 强SIAM JIS | 有界共享材料事件定理+近似状态误差+因果机制实验+可实现干预 | 独立事件预测、表示与审计控制、实际成本收益或足够强的机制发现 | 一般DR角度定理/配方恒等式是新成像理论 |
| 强Inverse Problems | 与一般逆问题相关的新稳定性/不可辨识性/有限分支访问理论 | 超出单个有限矩阵模型的清楚条件、尖锐反例与界；可能需要新增测量设计或失配理论 | 当前事件阈值已证明材料唯一性或全局恢复 |

**优先建议：**围绕有限材料分支事件与表示归因先形成一篇SIAM JIS型数学—算法研究；如果最终主要价值是低成本强3D/实测表现，则改走TAP。Inverse Problems需要比当前三条代数/几何命题更深的可恢复性结论，不能只换叙事投稿。

论文结构：问题与反例 → 共享材料本构几何 → 精确事件与近似局限 → 最简事件proposal与表示控制 → 已付候选公平比较 → 独立2D → 3D/实测 → 不可辨识与失配。第一篇不必把所有R4–R12方法都解释成同一个理论；必要基线与历史证据如实保留。

若后来独立证明了连续Maxwell稳定性或有限数据材料商空间理论，再考虑第二篇理论论文。工程与理论两篇必须有不同主要结果和验证问题，不能把同一组恢复图拆成两篇。

## 15. 交付验收语义

最终交付包括冻结配置、完整run/candidate/cost/event/profile输出、全部失败清单、图表生成代码和逐主张证据表。验收的是可复现性、归因清晰度、理论条件与实际结果一致，不是某个任意成功率。

如果结果表明事件proposal没有稳定收益，但有限事件预测成立，应完成负结果分析、判断是否只是盒约束效应，并把投入转向经同空间/同审计控制后仍存在的问题；如果NFFT控制完全解释G19—DR差距，应承认并简化方法。如果两者都不能解决噪声球内错误材料，应将后续工作转向材料选择、额外观测或模型失配，而不是继续增加反射步数、Green rank或网络参数。
