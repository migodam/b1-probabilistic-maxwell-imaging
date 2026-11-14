## research/research_max_struct_r4/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT-R4：从 CIE 公式到空间残差传播
2: 
3: 本轮沿 R3 继续研究用户引用的 [CIE机制公式解析](sources/CIE_SESSION_REFERENCE.md)，明确验证原候选 A、B、C。**A 在一个历史困难弯条上带来真实材料恢复，并可缩成一次邻域交换；但六对象新鲜确认没有超过成熟 CIE，网格细化和分块平移平均也暴露了适用边界。B、C 的已测版本没有产生新成像收益。当前不能称已经得到成熟、普遍优于 CIE 的新算法。**
4: 
5: 有价值的结果是：识别了非局部残差权重改变材料读出的准确条件，完成了一个简单物理兼容算法、反相控制、新鲜确认、表示敏感性和失败后结构修正。没有用更低 residual、额外 prior、更多观测或 GT 选根代替材料恢复。R1/R2/R3 冻结证据不改写。
6: 
7: 直接阅读：[A/B/C逐项答复](ABC_VALIDATION_ZH.md) · [完整推导](THEORY_DERIVATIONS.md) · [理论状态](THEORY_STATUS.md) · [独立审查](audit/INDEPENDENT_REVIEW.md) · [复现](REPRODUCE.md) · [费用](COST_LEDGER.md) · [失败记录](FAILURES.md)。
8: 
9: ## 1. 真实模型、数据和评价
10: 
11: 复用现有二维 scalar-TM Hankel–Green disk-pulse VIE，未新造正演器。采用 `exp(−iωt)`、`λ0=1`、`k=2π`、背景ε=μ=1；实 contrast `m=εr−1∈[.02,20]`，相对对比源 `j=mE`，物理极化电流另乘 `−iωεb`。原始方程为：
12: 
13: ```text
14: j = X(Einc + G j),  X=diag(m)
15: F(m) = S(I−XG)⁻¹ X Einc
16: Phi(m) = 0.5 ||F(m)−Y||² / ||Y||² + R(m)
17: ```
18: 
19: 24×24 网格，576个实材料自由度，12个外部线源、48个外部接收器；ROI宽约.96λ0，网格间隔.04λ0、脉冲半径.018λ0。所有方法用同一完整多重散射 resolvent、公开 box 和历史 smooth/TV prior（各1e−4，width19.5，epsilon.01）。Fourier/APCS 实现来自成熟基线；当前小模型使用显式矩阵动作，不能把名称中的 FFT 当作本轮已实际发生的 FFT 加速。
20: 
21: 这是真实有限全波二维模型，不是三维矢量 Maxwell 或连续 PDE 证书。3D 新反演、R1 的 H(curl) 表示资格和非空单调性证书均没有被本轮二维结果替代。已有 R1 张量同场族和 Z-ALG-E1 连续材料图册的范围继续保留，不重复声称首次发现多解。
22: 
23: 开发对象是 R2 已公开弯条 `obs_b20d5a55677d4eca`，不能算新盲测。新鲜确认在算法冻结后一次生成六个独立材料（环、弯条、双岛各两个），opaque ID 不带家族，GT/clean Y/mask/seed 不进入反演或远端包。全部18个预测封存、按公共 Phi 固定候选选择后才评分。
24: 
25: 噪声为独立复高斯 CN，实/虚方差各σ²/2，σ为复幅度标准差。噪声判断用 `||F−Y||/σ≤25.1654222`（576复观测，99%半径）；历史归一化 Phi 与噪声 whitening 的尺度关系见 [数据合同](DATA_CONTRACT.md)。材料指标固定为 RMSE/12、预测自身 Otsu256 阈值的 IoU；沿用 .15/.7 门，仅作逐例描述，未重新调门槛或制造显著性。
26: 
27: ## 2. 从公式得到的 A：传播本构缺陷，而非近似正演
28: 
29: 令 `r=j−XQ`、`Q=Einc+Gj`，取实际 Green 矩阵中固定2×2空间单元块 Gin：
30: 
31: ```text
32: Kβ = I + X(βI−Gin)
33: Rexact = β Kβ⁻¹ r
34: ```
35: 
36: Gin=0 回到普通 CIE。只要变换可逆，其状态零集仍是完整原 Maxwell 方程；但离开物理解时，目标和路径会改变。因此有待验证的是有限缺陷的优化几何，而不是把等式重写本身当作成像贡献。材料梯度包括 Kβ 对 m 的导数，不能把它冻结后仍称精确梯度。
37: 
38: 分离自项，令 `D=I+X(βI−diag G)`、`L=D⁻¹XGoff`、`u=βD⁻¹r`，则：
39: 
40: ```text
41: Rexact = (I−L)⁻¹u
42: Rplus  = (I+L)u
43: ||Rexact−Rplus|| ≤ q²/(1−q) ||u||,  q=||L||<1
44: ```
45: 
46: 最后一式只在所用的同一个范数中成立。公开24网格整个合法材料盒的解析∞范数上界为 .16810（β=.2）和 .03784（β=1）。一次交换因此保持原状态零集，只需固定邻接矩阵动作和标量复除法，额外复杂度 O(Ns)，不需要小 LU。最终数据排序和 corrector 始终使用完整原物理 forward，未截断多次散射。
47: 
48: 算法保留成熟 mf3/6/10、各500步的双β调度（current β=1、后两阶段 material β=.2）、同 source metric 和 global PR。普通 global CIE 是直接机制对照；另给原 sourcewise CIE 追加 mf12、6500步，总8000步。三者最终都接同一原物理 L-BFGS-B corrector，1000步上限，物理目标统一缩放1e6以使数值停止尺度一致。
49: 
50: ## 3. 开发正例、反相控制和网格边界
51: 
52: ![开发对象与所有主要机制](figures/development_materials.png)
53: 
54: | 开发端点 | RMSE/12 | IoU | 原物理相对残差 |
55: |---|---:|---:|---:|
56: | 成熟 sourcewise CIE，8000步 | .557674 | .255814 | .0977066 |
57: | global CIE，共同 corrector | .389183 | .422680 | .107821 |
58: | 只吸收 diagonal self | .389418 | .422680 | .107821 |
59: | 精确四单元吸收 | .042654 | .928571 | .0009670 |
60: | 缓存版精确吸收 | .041691 | .928571 | .0009695 |
61: | 一次交换 Rplus | .041672 | .939759 | .0009679 |
62: | 同幅反相 Rminus | .558703 | .257813 | .0976930 |
63: 
64: 同物理数据、最终 prior 和 corrector 下，A 确实改善了这个对象的材料图。仅 diagonal self 和额外普通模态/迭代没有复制收益。它在 corrector 前已较好：一次交换 mf10 的 RMSE/12=.08327、IoU=.77451；corrector 后 .04167/.93976。global/self 的 corrector 则降低 residual 同时恶化材料。所有开发 corrector 都触及1000步上限，不据此认证吸引域或局部最优。
65: 
66: ![反相与必要细化](figures/validation_development_phase_grid.png)
67: 
68: π反相只令辅助 `Goff→−Goff`，真实 forward G 不变。两罚函数都保留 `||u||²`、`||Lu||²`、相同 metric 对角、图拓扑和逐次运算结构，只翻转 `2Re⟨u,Lu⟩`。反相失败支持符号/相位方向在该实例中有作用，不能外推为真实 Green 相位对所有对象最优。
69: 
70: 一次32×32诊断使用同192像素几何、同源/接收、同 Fourier 容量和 prior 数值参数。onehop/global **都恢复**：RMSE/12=.053068/.053102，IoU都 .935484，残差都约 .0009705。A 在24网格上的比较优势没有保持。这里同时改变有限脉冲离散、块物理宽度和噪声实例，不能把它称严格连续 PDE 收敛，也不能省略这个反例只展示24网格。
71: 
72: 跨平台复验使用原公开开发观测和冻结 onehop，在既有 Windows NumPy2.2.6/SciPy1.15.3 上重放；本地为 NumPy2.5.3/SciPy1.18.1。Windows 原物理残差 .00096835、Phi6.49305e−6，独立材料评分 RMSE/12=.041603、IoU=.951807，保留材料恢复正例；完整评分见增补。该次 solve wall约376.20秒，不能声称 i9 对此工作负载更快。第一次 wrapper 未加入已有 vendor 的 threadpoolctl 路径，在任何物理运算前失败，已保留费用并修正入口，未安装依赖或改算法。
73: 
74: ## 4. 六个新鲜对象：没有超越成熟 CIE
75: 
76: 全部预测先封存，随后首次读 GT。18个任务全部完成冻结辅助日程；没有失败任务被删掉，也没有短 wall cap 截断基线。
77: 
78: | 对象尾码 / 家族 | onehop RMSE/12 | global | mature native8000 | 三臂联合结果 |
79: |---|---:|---:|---:|---|
80: |98644322 / 环|.071621|.073369|.085988|均通过|
81: |4f224bb3 / 弯条|.761368|.741539|.743118|均失败|
82: |78134194 / 双岛|.073440|.072053|.069852|均通过|
83: |50fc483c / 弯条|.576876|.575755|.577784|均失败|
84: |fbac472d / 双岛|.073560|.072850|.071273|均通过|
85: |e5aa42ae / 环|.066619|.063476|.067396|均通过|
86: 
87: 三方法都是 **4/6**，相同四对象成功。onehop 相对 global 的材料 NRMSE 仅1/6更低，5/6更高；IoU是3高、2低、1相同。样本小，只描述逐例，不能作统计显著性或人群优越性宣称。
88: 
89: ![全部六对象与十八个冻结端点](figures/confirmation_six_materials_postseal.png)
90: 
91: 预先声明的三臂候选池按同 Phi 选最小，成本为三臂总和，仍只有4/6。两个弯条恰好都选中了三臂中材料误差最差者。42cb 的 A 将 Phi从 global 的 .004105 降到 .002193，但 NRMSE从 .74154 升到 .76137。这个结果明确区分了物理拟合与材料正确；没有离线按 GT 改选。
92: 
93: 两个失败对象的合法 GT 自身都满足声明噪声球，白化残差为25.1189、24.0092，Phi约6.54e−6、6.92e−6。因此当前高残差端点不能作为信息不可辨识的证据；它们没有达到已知可行的数据拟合区域。另一方面全部 corrector 都触限，末段目标仍下降，不能未经诊断就宣称已进入认证局部极小。
94: 
95: 详表、候选池择劣和每个任务成本见 [确认评分](results/CONFIRMATION_SCORES.md)、[独立确认审计](audit/confirmation/CONFIRMATION_COMPARISON_AUDIT.md)、[停止诊断](audit/confirmation/FAILURE_STOP_DIAGNOSTIC.md)。
96: 
97: ## 5. 失败后继续推导与结构修正
98: 
99: 材料读出给出一个准确的结构判别。实化线性化写为 `b−C dm`，C满列秩，M为正定残差 metric，Π为 range(C) 的正交投影。在无prior、内点、精确根的局部问题中：
100: 
101: ```text
102: D_M = (CᵀMC)⁻¹ CᵀM
103: D_M − C† = (CᵀMC)⁻¹ CᵀM(I−Π)
104: ```
105: 
106: 所有右端的读出都不变，当且仅当 range(C) 对 M 不变，即 `[M,Π]=0`。逐cell对所有源共享的正标量权重有这种抵消；跨cell权重一般没有。实际某一步改变还需 `CᵀM(I−Π)b≠0`。这解释了非局部 metric 可以把无法由局部实材料解释的 defect 分量送入材料更新，也解释了它可能有害。一个复源就可能产生这种余量；每cell只有一个实方程、C为方阵的玩具模型则会丢掉该机制。
107: 
108: 这是已有加权 LS 和 R3/D1 抵消关系的结构精化，不是新发明的投影定理。实际算法存在 prior、活动 box、非零状态 defect、不精确内层和双β，因此不能把这个局部等式直接当完整迭代的收敛定理。
109: 
110: 由确认失败进一步提出一次无调参修正：平均四种2×2分块平移的传播，得到水平/垂直邻边权重1/2、对角邻边1/4；边界不周期 wrap。权重来自分块包含概率，保留相同解析可逆界，但去掉任意 parity 起点。它平均的是 residual transport，不是四个平方目标。
111: 
112: 该版本只作两个**已暴露开发诊断**：旧正例和按 public 数据不合格、opaque ID最小选出的42cb。两个结果都失败：旧正例残差 .107660、材料NRMSE .38972；42cb残差 .143491，NRMSE从原onehop的 .76137下降到 .64140。这是保留的局部材料指标改善，但仍未成像且数据残差恶化，不能用GT优选它。没有调权重、β或停止规则救回结果。它不能改变原6对象/18端点的确认分母。这否定了“只需去除分块起点即可修复迁移失败”的具体解释，也表明原正例依赖更细的优化路径。
113: 
114: ![分块起点修正的全部两例](figures/partition_refinement.png)
115: 
116: 六次新原物理梯度诊断只评估已有失败端点，新增6LU/144RHS，没有优化或GT访问。box投影梯度范围约2.28e−7至8.94e−6，高于原1e−8停止目标；随后为第一个公共失败对象的 A/native 各给同额额外4000步原物理续跑额度。native实际3361步达到原投影梯度目标（7.46e−9），残差仍 .088672，表明这个基线的失败不能只归因原1000步太短；这仍不是局部极小或唯一性证书。A在251步因函数相对下降停止，投影梯度3.87e−6仍未达标，残差 .065584，属于数值停滞而非认证收敛。独立评分中，native材料NRMSE从 .743118恶化到 .750890；A为 .760862，也未恢复。新端点见[续跑评分](results/DEVELOPMENT_CONTINUATION_SCORES_V1.md)，原确认18预测和4/6不改。
117: 
118: ## 6. B 和 C 得到什么
119: 
120: **B：纯 port-only 消元会回到原目标。** 对自由全 current 的二次端口罚函数精确消元，只得到原数据残差的正比例缩放和材料无关常数。加入固定 `eta>0` 后才得到材料相关 metric：
121: 
122: ```text
123: C = S(I−XG)⁻¹
124: M = (eta I + lambda CC*) [eta I + (1+lambda)CC*]⁻¹
125: ```
126: 
127: 本轮实现含 resolvent 变化的精确 envelope 梯度；eta仅由公开 host固定。lambda=0 的 WRI/state控制和lambda=1各250辅助步，再接同 corrector。lambda=1残差 .093255，但材料NRMSE .533540，IoU .244094；lambda=0也失败。有限eta确实改了几何，但没有成像收益；不继续扫eta找正例。改变辅助 metric 与未缩放 prior 的相对作用也被明确保留，不能只凭 prior 名称相同排除这种影响。
128: 
129: **C：完整有限响应没有给出可用逃逸方向。** R3端点自由面33维最小曲率为正，而全空间有负曲率；box切锥投影破坏了协同符号抵消，使两方向曲率转正。四条固定完整有限射线无改善；进一步两次固定500步 box二次模型修正仍无有效射线。保留了原完整 Maxwell 有限响应，失败来自方向/可行性及预算范围，而不是人为 Born 截断。没有全 box/global 证书，也没有用GT差方向构造伪成功。
130: 
131: 另保留 DCTMC式非局部 T完成、照明加权局域化和 reciprocity/光学定理诊断。T可高度拟合通道却不能局域化成正确材料；其更低 residual 未改善材料。已测 reciprocity 修正落入极弱可见 current 方向，代价巨大，没有据此新增硬约束。它们均在 [研究循环](RESEARCH_LOOP.md) 和失败记录中说明如何分析、修正、转向。
132: 
133: ## 7. 计算成本与可复现性
134: 
135: 六对象确认复用现有 i9 四个 CPU worker 的有限 FIFO 队列，无新服务、无CUDA安装，实际GPU任务为0。计入共享生成一次，总 **18653次576维物理LU、445932个complex RHS**。前台远端实测 elapsed **1179.03秒**；记录CPU约3804.34秒。worker wall累计3900.01秒不是并行历时。
136: 
137: | 六任务工作 | 物理LU | complex RHS | worker wall累计秒 |
138: |---|---:|---:|---:|
139: |onehop|6202|148332|1147.09|
140: |global|6216|148680|577.95|
141: |native8000|6229|148848|2174.97|
142: 
143: onehop比普通global消耗更多内层计算，没有等算力优势；native得到更多普通模态/迭代，但仍不能把它视为所有CIE日程的充分收敛。小4×4 LU、稀疏邻接动作、标量除法、源/接收RHS和满维物理LU单独记账，不能混加成一个“forward次数”。32网格1024维、微型接线16维也与576维分列。拒绝候选、数值负结果、缓存回放、入口失败和重复上传都保留。
144: 
145: 全轮已记录物理总量为38757次LU、980383个model complex RHS：其中576维36698/931161，1024维2054/49116，16维接线5/106；另有非局部局域化30次LU/17280RHS，不混作物理LU。全轮开发、验证、确认、后续诊断及后处理的最终总账见 [COST_LEDGER](COST_LEDGER.md) 和 `audit/cost_ledger/`。金额、LLM/API收费、部分早期bootstrap CPU没有可核定数据，写null而不是零。未核验SHA；普通git基版本、冻结配置与执行源码快照用于来源记录。
146: 
147: 代码入口、公共数据重放、原始数组和环境版本均交付。既有Windows和本地软件栈实际版本见 `audit/ENVIRONMENT.json`。新的子任务明确配置 Sol medium（整理/单篇读取/绘图）或 Sol high（数值/评分/独立审查），父线程保留数学和最终判断，见 `audit/MODEL_ASSIGNMENTS.md`。
148: 
149: ## 8. 论文价值和下一步
150: 
151: 不能宣称发明了 Green 分解、材料相关 CSI 预条件、T重整化、Neumann展开、加权LS或完整物理Hessian。2000年的 [preconditioned CSI](https://doi.org/10.1109/36.851966)、[NIE/CIE与重写综述](https://doi.org/10.1109/OJAP.2021.3063248)、[CIE-I/SOM](https://doi.org/10.3390/jimaging5020027)、[DCTMC](https://doi.org/10.1103/PhysRevE.94.043318) 都有直接相邻先例；原文可得和未决部分在文献账本中区分。
152: 
153: **已成立的成像收益**限于旧24网格开发对象，经过反相、self、成熟扩模基线以及跨平台复验。**未成立的主张**是新鲜分布上的成像成功率或等成本优势。B/C的已测方向也未突破。最有研究价值的是“同物理零集不等于同材料更新”、本构切空间不变性判别，以及正例如何随分块/网格改变而消失的反例；目前更适合机制与失败边界研究，尚不足以独立支撑普适优越算法论文。
154: 
155: 下一步应直接分析 `CᵀM(I−Π)b` 在成功/失败路径中的物理方向性，区分可解释的材料改变量与错误反馈，不继续挑分块、扫β或叠加网络。已保存的状态、完整数据和新投影梯度让这个问题可以从小量有内容的诊断开始。不能用两个失败弯条证明数据不可辨识，也不能把新增照明的可能收益混入本轮同数据比较。

## research/research_max_struct_r4/CHAT_HANDOFF.md

1: # 给Chat / 后续研究Agent的R4入口
2: 
3: 先读 [中文报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/RESEARCH_REPORT_ZH.md)，再读 [A/B/C答复](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/ABC_VALIDATION_ZH.md) 和 [完整推导](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/THEORY_DERIVATIONS.md)。网页与全部图在 [这里](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r4/)。
4: 
5: ## 必须继承的结论
6: 
7: 1. A是实际四cell Green的本构缺陷传播，不是近似物理forward。精确block可缩成一次交换 `R=(I+L)u`，整个公开box有解析可逆界。旧开发弯条RMSE/12约.04167、IoU.93976；Windows独立栈复验约.04160/.95181。
8: 2. self-only、同幅π反相和成熟native8000在这个旧例失败，支持具体局部机制；但32网格的普通global也恢复，优势不是跨网格普遍事实。
9: 3. 六个新鲜对象、18个冻结端点、同prior/同corrector：onehop/global/native8000全是4/6，相同两个新弯条都失败。onehop相对global NRMSE仅1胜5负；三候选最低Phi选择仍4/6，两个失败弯条都选到最差材料臂。
10: 4. 失败后用四种分块平移平均构造8邻接传播，未调权重：两个已暴露开发例均失败，旧positive丢失；42cb NRMSE稍改善但数据残差更差。这不改原18和4/6。
11: 5. 一例成熟native再给4000步物理额度，实际3361步达到原PG精度，仍没有恢复，NRMSE还变差。A在251步FTOL停滞，PG未合格。原1000步确认和追加费用分开，不称局部极小证书。
12: 6. B纯port-only全current消元是原目标缩放；有限eta版本确实改变metric，但lambda1更低residual更差材料。C完整有限协同射线未逃逸，原因涉及box可行性和方向，而不是公式不精确。
13: 7. 有用的读出判别：`D_M−C†=(CᵀMC)⁻¹CᵀM(I−Π)`。局部权重抵消取决于材料切空间不变性；非局部耦合一般改变读出，但改善没有保证。它是已有加权LS/R3取消关系的精化，非首创新基础定理。
14: 
15: ## 证据入口
16: 
17: - [六对象原始分数](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/results/CONFIRMATION_SCORES.md)
18: - [独立确认与成本审计](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/audit/confirmation/CONFIRMATION_COMPARISON_AUDIT.md)
19: - [最后科学审查](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/audit/FINAL_SCIENTIFIC_REVIEW.md)
20: - [费用与维度口径](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/COST_LEDGER.md)
21: - [失败后如何推导/转向](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/RESEARCH_LOOP.md)
22: - [全部代码和数组目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_max_struct_r4)
23: 
24: 不要把单个公开开发正例、最低residual、额外cost、更多观测或GT事后指标改写成新盲测优势。preconditioned CSI、NIE/CIE、T重整化、Neumann和variable projection都有直接先例；具体算法的普遍收益和论文原创性没有成立。3D新反演NOT_RUN，R1单调性证书缺口未被解决。
25: 
26: 所有已评分GT现在都是公开开发资料；不得再次称独立blind。新的算法版本允许分析这些失败，但必须另立开发版本，并用新独立数据确认。

## research/research_max_struct_r5/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT-R5：材料读出中的残差传播有效在哪里，失效在哪里
2: 
3: ## 最终判断
4: 
5: **R5把旧弯条成功的作用位置定位到了非局部材料读出，并发现它依赖给定阶段路径；没有得到稳定超越成熟 CIE 的新成像算法。** 普通电流更新加上one-hop材料读出，已能复制旧正例；但两条新弯条仍失败。只在早期或后期使用这一读出也未复制成功。最有价值的产物是这组因果干预，以及解释其适用边界的有限缺陷导数和块谱分析。
6: 
7: 本轮不是九数字新基准或三维成像论文。完成了四个公开开发对象的8条轨迹重放、174个状态诊断、120个同状态反事实读出、7次带共同原物理corrector的有限干预，以及一个真实三维vector-DDA两体素谱反例。未拟合成功分类器、未增加照明、未换prior、未用GT在线选候选。168个检查点不是168个独立材料样本。
8: 
9: ## 四个必须回答的问题
10: 
11: ### 1. 为什么R4在历史弯条上有效？
12: 
13: 能够确定的是**作用位置和路径条件**，还不是几何充分条件。把one-hop分别放入current和material两块，得到：
14: 
15: | 电流更新 | 材料读出 | RMSE/12 | IoU | 完整物理相对残差 |
16: |---|---|---:|---:|---:|
17: | 普通global | 普通 | .389183 | .422680 | .107821 |
18: | one-hop | one-hop | .041672 | .939759 | .0009679 |
19: | one-hop | 普通 | .558514 | .257813 | .0976930 |
20: | 普通global | one-hop | **.041480** | **.939759** | **.0009692** |
21: 
22: 前两行复用冻结历史端点，后两行为新R5干预；本地重放旧global/one-hop的各阶段状态与历史一致。四臂使用同数据、576个独立实材料、同box/prior/初始化、同mf3/6/10及每阶段500步，最终共同1000步original-Maxwell corrector。辅助运算量不同，不能称等wall或等全部算量优势。
23: 
24: 材料读出的额外非局部项确实能改变路径，而current传播不是复制这个正例的必要部件。普通self项和反相控制此前均失败，排除了“只要改自项或加等强耦合就行”。这不是更多物理信息，也不是更高材料容量。
25: 
26: ![旧弯条位置消融](FIGURES/03_old_bar_factorial.png)
27: 
28: 在第一个自然阶段结束后，完整搬运m、J和Fourier系数，再按原程序扩模/重置PR。普通global→material-only最终NRMSE=.388613；material-only→global为.558738；全程material-only为.041480。**早期状态单独不够，后期更换读出也不够**，仅对这一固定切换点成立。不能据此宣称证明了吸引盆地、必须每一步传播、或其他切换点都会失败。两切换corrector均耗尽1000步，未获渐近收敛证书。
29: 
30: ![自然阶段切换](FIGURES/06_natural_stage_handoffs.png)
31: 
32: ### 2. 为什么两个新弯条失败？
33: 
34: R5没有找到能仅凭形状提前判定成败的规则，但排除了几个简单解释。失败对象具有与旧成功相同的正cross、同相残差占优和明显的metric导数项；这些量的存在本身不能解释成功。删除current传播的简化也未修复它们：
35: 
36: | 已公开对象 | R4 one-hop NRMSE | R5 material-only NRMSE | R5 IoU | R5原物理残差 |
37: |---|---:|---:|---:|---:|
38: | 弯条42cb | .761368 | .766323 | .106280 | .066004 |
39: | 弯条63ab | .576876 | .404823 | .403846 | .111383 |
40: | 正常环100b | .071621 | .072525 | .917355 | .001010 |
41: 
42: R4确认端点来自Windows，R5来自本地CPU；表中微小差异不能解释为算法优势。63ab的材料误差下降较大，但残差更差，仍未恢复；不得据GT在这些候选间选根。42cb几乎维持同类错误结构，正常环保持成功。
43: 
44: ![公开开发扩展](FIGURES/04_extension_materials.png)
45: 
46: 在one-hop最后辅助阶段的检查点，两个失败对象分别有37、23个材料上界活动cell，旧弯条和环为0；这些上界点的原物理负梯度也朝外，并非单纯“辅助目标顶住了一个物理梯度希望释放的box”。失败对象沿实际材料步的精确物理曲率为负，而两个成功对象为正；但这些检查点有非零投影梯度，不能叫驻点或证明局部极小。Gram中约1e−18的负小特征值是浮点现象，不是负物理曲率或精确零空间证据。
47: 
48: 因此当前支持的是：**相同局部残差滤波在不同全波状态上会把材料路径引向不同区域；块内相位结构没有对接收数据方向或材料正确性的普遍保证。** 两失败端点还有显著原物理失配，不能用“不可辨识的等数据解”替它们开脱。预算不足、非凸路径和离散表示影响未完全分离；R4的延长corrector、32网格及分块平移结果继续约束结论。尚不能给出几何层面的完整失败定理。
49: 
50: ### 3. 得到了什么可验证的新认识或简单改进？
51: 
52: **第一，根处加权LS不是离根实际读出的完整解释。** 固定辅助J，r=J−XQ，R=P(m)r，则
53: 
54: ```text
55: dR[h] = −P C_mat h + (dP[h])r
56: g_state = −C_matᵀPᵀP r + a_r
57: <a_r,h> = <Pr,(dP[h])r>.
58: ```
59: 
60: 被根处公式省去的a_r可与甚至大于总state梯度，因为不同项会抵消。普通CIE的这一分量a_r非正，因此其负梯度方向指向提高contrast；它也存在于成功轨迹，不能当失败分类器。实际读出还含prior、box、非线性内层和current历史。
61: 
62: **第二，正cross并不反对one-hop有效。** 168个冻结检查点的base-plus全域cross全部为正，范围.000299至.107936；空间上仍有负cell。反相控制自身的cross符号应取反，不能把base-plus诊断误称反相目标值。加权LS transport与真实材料步的cosine范围−.495至.450，范数大不等于推动了实际步。样本中34次实际步出现Phi下降但材料误差上升，另19次相反；这只是相关轨迹的描述，不是独立统计显著性。
63: 
64: **第三，滤波的谱结构可明确写出。** 对均匀材料2×2正方块，四个几何固定Walsh模态的Green特征值为
65: 
66: ```text
67: λ = (2ge+gd, −gd, −gd, −2ge+gd)
68: one-hop residual weights = |1+aλ|².
69: ```
70: 
71: 实际scalar-TM算子在β=.2、公开box上界m=20时，同相模态权重1.23674，两个交替模态.94677，棋盘格.90752。四对象三阶段末态的同相残差能量占.7095–.8344，贡献正cross，其余模态合计贡献负cross。它解释了总体符号，也解释反相为何改变优化；成功与失败都满足，所以没有得到新的可靠成败预测器。非均匀材料会混合这些模态，分析保留了真实混合项。
72: 
73: **第四，简单改进只获得局部资格。** material-only少了current块中的残差变换，保留旧例恢复，值得作为机制控制保存；它不是已成熟的普适新算法，也没有等实际成本胜出的证据。
74: 
75: ![采样轨迹](FIGURES/01_sampled_trajectories.png)
76: 
77: ### 4. 是否已具备二维、三维和公开基准上的稳定成本优势？
78: 
79: **没有。** 二维只有公开开发对象上的作用定位与有限干预，未开展新独立确认；R4六对象4/6共同结果和32网格优势消失仍有效。历史九数字成熟CIE已全恢复，本轮没有把它们重跑成新创新证据，完整历史映射/状态见[基准合同](BENCHMARKS.md)。
80: 
81: 三维复用了R1真实dyadic Green，以两个相邻体素、ε=3、两种合法横向入射极化作6维全波解。它给出一个有区分力的反例：β=.2时，沿两体素连线的同相电场模态权重1.34284，而横向同相模态为.85619。这里“纵向”相对于体素连线，不是非法的入射纵向极化。一个同相模态被加强、另一个被减弱；二维scalar解释不能无条件搬到三维。仅1 LU/2 RHS，**不是三维成像或连续PDE精度认证**。
82: 
83: ![三维极化机制反例](FIGURES/07_vector_pair_modes.png)
84: 
85: Fresnel/eispy2d完成官方来源、物理约定、文件格式和标定缺口调查，未完成实测文件下载/标定/反演。没有伪造实测材料NRMSE，也未以合成数据冒充实测结果。
86: 
87: ## 证据、计算与执行边界
88: 
89: 所有online候选由同一原Phi选择；辅助目标不代替Phi。GT在预测/诊断封存后由独立评分脚本访问；modal分析在评分后提出，明确为post-score development且自身不读GT。[原封存理论前六节](audit/THEORY_DERIVATIONS_AT_FIRST_SEAL.md)和[文档演化记录](audit/POST_SEAL_DOCUMENT_EVOLUTION.json)分别保留；后续理论附节不是预注册内容。具体记录见[复现](REPRODUCE.md)和[账本](COST_LEDGER.md)。
90: 
91: 新增物理共 **7717 LU、183686 complex RHS**，其中7716次为576维二维系统，1次为6维矢量系统；维度不可混同为等成本。二维各任务wall合计290.396秒，另加三维约.274秒；这是任务计时之和，不是并行历时，更不是包含文献/绘图/Agent的完整总时长。全部本地CPU；无需为这些小任务占用i9/4060。LLM及货币费用未获计量，保持null。失败尝试已收费，历史预计算只引用不重计。
92: 
93: ## 论文价值与下一轮选择
94: 
95: 标准恒等式、Walsh谱、dyadic纵横分解均不主张原创。较有价值的是**可复现的作用位置消融、自然阶段路径依赖、有限缺陷导数与受约束物理曲率的区分**。目前适合作为方法学机制研究/负结果证据，尚不足以支持IEEE TAP级“稳定优于成熟CIE”主张。
96: 
97: 下一轮优先的问题不是再调β，而是：在给定物理预算下，如何构造能与接收数据方向关联、并尊重box切锥的材料更新？公式起点应是完整有限响应与残差加权曲率，而非单独cross符号：
98: 
99: ```text
100: dF[h] = S(I−XG)⁻¹ diag(h) E
101: d²Phi[h,h] = ||W dF[h]||²
102:              + Re< W(F−Y), W d²F[h,h] > + d²R[h,h].
103: ```
104: 
105: R4候选C已表明，把无约束负曲率方向投影到box可能毁掉它的协同作用，不能照搬重跑。最小下一步应先在失败端点的合法自由面/切锥，检验是否还存在可产生有限Phi下降的耦合方向，并同时对成功控制检查；只有出现新证据才接更多corrector与新盲测。若这一步仍无效，应保留本轮“局部残差滤波不充分”的反例，转向物理可见材料子空间，而不是扩大图像或超参数数量。
106: 
107: ## 直接证据入口
108: 
109: - [轨迹与反事实完整分析](TRAJECTORY_ANALYSIS.md)
110: - [推导及成立条件](THEORY_DERIVATIONS.md)
111: - [研究循环](RESEARCH_LOOP.md)、[完整负结果](FAILURES.md)
112: - [R4原Chat候选A/B/C结论](../research_max_struct_r4/ABC_VALIDATION_ZH.md)：A局部正例/泛化失败；B测试profile无成像优势；C受约束射线未逃离。R5只深化A，不重命名或隐藏B/C。
113: - [文献与原创性账本](literature/PRIOR_ART.md)
114: 

## research/research_max_struct_r5/CHAT_HANDOFF.md

1: # MAX-STRUCT-R5 给 Chat 的入口
2: 
3: 先读 [中文报告](RESEARCH_REPORT_ZH.md)，再读 [数学推导](THEORY_DERIVATIONS.md)、[轨迹分析](TRAJECTORY_ANALYSIS.md) 和 [研究循环](RESEARCH_LOOP.md)。所有R5二维对象已公开，非新盲测。
4: 
5: ## 最值得继承的结论
6: 
7: 1. 旧弯条成功可由**普通global current + one-hop material readout**复制：NRMSE=.041480、IoU=.939759；仅改current失败（.558514）。本例收益定位于非局部材料读出，不需current传播。
8: 2. 自然mf3阶段后，两种完整m/J/A交接均失败：global→material-only .388613；material-only→global .558738。指定协议中不是只需一个早期好初始化；也未证明全局盆地/所有切换点。
9: 3. 简化仍未恢复42cb/63ab：NRMSE=.766323/.404823，环=.072525。63ab材料有所改善但物理残差更差，禁止GT选根。
10: 4. 168个base-plus全域cross均正；几何固定块谱说明同相残差被加重、交替残差被放松。成功与失败共享此结构，因此没有可靠预测器。
11: 5. 非根处 dR=−P C dm+(dP)r，metric导数不可省略。根处加权LS transport不是实际prior/box/非线性readout的完整模型。
12: 6. 两体素真实三维vector-DDA给出极化反例：β=.2同相纵向weight1.34284，横向.85619。是离散6维全波机制例，不是三维成像。
13: 
14: ## 证据与限制
15: 
16: 8重放/168检查点+6历史态/120反事实；7次带同原物理corrector干预；1个矢量小模型。只有4个二维材料对象。原始数组、source snapshots、失败和各run receipts全保留。新增7717LU/183686RHS（不同系统维度分开），全部本地CPU。没有新照明收益、额外prior收益、NN或参数扫描。
17: 
18: R4的六对象4/6共同结果、32网格比较优势消失和进一步corrector负结果全部继承。R5九数字/新盲测/真实测量/3D成像NOT_RUN。历史九数字完整single/dual表见[附录](RESULTS/HISTORICAL_DIGITS_REFERENCE.md)，不是R5运行。
19: 
20: 原Chat候选A/B/C：R4的[A/B/C完整答复](../research_max_struct_r4/ABC_VALIDATION_ZH.md)仍有效。R5深化A并定位作用位置；未把B/C重命名成新成功。
21: 
22: ## 下一步最有价值的问题
23: 
24: 不要再用cross正负、单独||transport||或相同零集给成像成功作保证。优先研究完整原物理residual-feedback在**合法材料切锥**中的协同方向，结合接收数据方向；先解释R4候选C的box投影破坏，再做一个失败例加控制。只在有限原Phi下降及独立确认有新证据后扩大。现有结果是机制定位/反例的论文素材，不是成熟“稳定胜过CIE”的算法。
25: 
26: 标准乘积法则、Walsh对角化、dyadic纵横分解均不主张原创；原创性仍需围绕具体可检验机制及有效算法单独建立。
27: 

## research/research_max_struct_r6/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT-R6：A路径扩大对照与FFT-TSOM接入
2: 
3: ## 最终判断
4: 
5: **A-material 已成为可运行的 CIE/FFT-TSOM 材料更新插件，能补偿两个历史困难弯条上的失败；目前证据不足以把它作为稳定优于 Dual-CIE 的默认替代算法。** 本轮实际完成21个材料×4方法＝84个主任务，另完成2个预指定材料×2方法＝4个真正FFT工程任务。数字1–9全部纳入；六个新对象由独立端生成，GT保持隔离时冻结求解日程，预测全部封存后才打开GT。实际生成回执早于最终config冻结约5分钟，不声称先冻结config再生成数据；失败与额外成本均保留。
6: 
7: | 分层联合成功 | CIE | Dual-CIE | Global-Dual | A-material |
8: |---|---:|---:|---:|---:|
9: | 公开历史数字1–9 | 8/9 | **9/9** | 8/9 | 8/9 |
10: | 公开历史困难对象 | 1/6 | 1/6 | 1/6 | **3/6** |
11: | 六个独立新对象 | 4/6 | 4/6 | 4/6 | 4/6 |
12: 
13: 联合标志沿用材料NRMSE≤.15、预测自身Otsu IoU≥.7及相应数据/合法性合同；历史clean门与新CN99噪声门分开，不混成同分布成功率。连续指标和所有对象见[完整逐例表](RESULTS/analysis/ALL_CASES.md)。A在数字5的NRMSE=.1142、数据拟合合格，但IoU=.6667，不能把8/9误读为它完全没有成像出数字5；Dual-CIE同例IoU=.7091。普通CIE在数字7则是明显的材料和数据失败。
14: 
15: **可以build upon FFT-TSOM，而且本轮已经实现并跑完。** A的材料读出不依赖显式Fourier矩阵。真实FFT版本替换了U/G及伴随作用，保留非周期自由空间Green、自项、接收和完整Maxwell验收。两例A的总CPU下降约16%和23%，但仍比同FFT的Dual-CIE贵；长轨迹终点有浮点敏感性，不能说仅是完全相同算法轨迹的免费加速。
16: 
17: ## 1. 实际改了什么，保留了什么
18: 
19: 固定公式是
20: 
21: ```text
22: r = J − mQ,    Q = Einc + GJ,
23: d = 1 + (β−Gii)m,   u = βr/d,
24: RA = [I + diag(m/d)Goff]u.
25: ```
26: 
27: Goff是实际Green在预定2×2 cell邻域内的非对角部分。A只用RA作材料块残差；current仍用普通global PR CIE，精确材料导数包括d随m变化的项。R5缓存仅减少重复求值。最终候选全都由完整物理F(m)和共同公开prior评分，不以辅助损失或GT选根。
28: 
29: 四臂同24×24网格、576独立实材料未知量、12源48收、box [.02,20]、host .1、原smooth/TV prior；mf3/6/10各500步，再自然续接mf12共1500步。共同原物理corrector最多3000接受迭代，原停止精度不改、没有45秒墙限。CIE为每阶段current/material同β；Dual-CIE为分块β，二者都保留成熟sourcewise PR；Global-Dual用来拆分global PR本身的影响。没有GT字体/support降维或训练。
30: 
31: mf12仍受DPCS affine约束，并非无约束current；同迭代数也不等实际费用。原Phi=.5||F(m)−Y||²/||Y||²+R(m)沿用基线，不冒称σ白化后的Gaussian likelihood；σ用于独立discrepancy。详见[数学推导](THEORY_DERIVATIONS.md)、[基线审查](BASELINE_AUDIT.md)及[冻结协议](BENCHMARK_PROTOCOL.md)。
32: 
33: ## 2. 历史困难对象：正面结果确实扩大，但范围有限
34: 
35: 以下NRMSE均用公开contrast幅度归一化，数字越小越好。
36: 
37: | 对象 | CIE | Dual-CIE | Global-Dual | A |
38: |---|---:|---:|---:|---:|
39: | FBG曲条7b8d | .57618 | .39313 | .39313 | **.04686** |
40: | R4旧弯条b20d | .55310 | .51072 | .38955 | **.04128** |
41: | 已公开失败弯条42cb | .75341 | .73918 | .75081 | .70840 |
42: | 已公开失败弯条63ab | .57463 | .58186 | .57730 | .40475 |
43: | 环e499 | .42483 | .42478 | .42480 | .42484 |
44: | 双岛437d | .07595 | .07649 | .07639 | .07641 |
45: 
46: 两个成功弯条上，A的原物理相对残差分别为.000154与.000966；Dual-CIE约.112与.101。它们不是仅辅助目标下降，也不是仅材料误差稍小。Global-Dual同样失败，因而global PR本身不能解释这两个正例。**新增的经验认识是：材料块传播的补偿效果可以在第二个公开弯条复现。** 这些都是历史开发对象，不能宣布两次独立盲测突破。
47: 
48: 42cb和63ab仍失败；63ab尤其说明指标不可混淆：A的NRMSE下降到.405，但物理残差.1114反而高于Dual-CIE的.0984。环的四种终点几乎一致。故不存在“弯条都受益”或“空腔都失败”的简单形状规则。
49: 
50: ![全部历史困难对象](RESULTS/scoring_primary/DEV6_DIFFICULT_all_objects.png)
51: 
52: ## 3. 六个新对象：没有新增确认优势
53: 
54: | 新对象 | CIE NRMSE | Dual-CIE | Global-Dual | A | 四臂结果 |
55: |---|---:|---:|---:|---:|---|
56: | 非等三岛 | .09500 | .09376 | .09471 | .09474 | 均合格 |
57: | 双层空腔 | .08110 | .08272 | .08258 | .08282 | 均合格 |
58: | 细弯条 | .03627 | .03475 | .03437 | .03532 | 均合格 |
59: | 偏心薄环 | .61473 | .61841 | .61436 | .61507 | 均失败 |
60: | 近邻双岛 | .06773 | .06583 | .06434 | .06438 | 均合格 |
61: | 开缝壳 | .61904 | .63312 | .63310 | .62596 | 均失败 |
62: 
63: 新六例中，A对Dual-CIE的NRMSE有3低3高，对Global-Dual只有1低5高；多数差值很小，不能拿全21例13个较低制造统计显著性。新细弯条、双层腔和近邻耦合对象都成功，也反驳仅靠“有曲率/有孔/有强邻域耦合”就预测A优势的解释。
64: 
65: ![全部独立新对象](RESULTS/scoring_primary/NEW6_CONFIRMATION_all_objects.png)
66: 
67: 六例为同离散合成数据，噪声为.1%复RMS，非独立网格或实测验证。现在GT随报告公开，后续使用只能标公开回归。新材料、噪声、先验和评分合同见[数据合同](DATA_CONTRACT.md)。
68: 
69: ## 4. 失败来自哪里：可行目标见证与停止状态
70: 
71: 仅在封存后，用生成端已保存的合法真材料和Yclean计算原Phi，不新增正演、不改候选。偏心薄环的可行Phi为8.70e−6，四臂却停在约4.85e−3；开缝壳的可行Phi为8.55e−6，四臂约3.81e−3至4.07e−3。算法残差分别约为真材料噪声残差的102–103倍、91–94倍。这说明它们**没有到达已存在的低目标可行区**，不能把本次失败直接解释为数据不可辨识；这同样不证明反演唯一。
72: 
73: 两者还必须区分：薄环四臂都达到3000 corrector上限，最终投影梯度约2.8e−6至5.9e−6，预算/慢收敛影响仍未排除；开缝壳四臂按投影梯度条件结束，范数约9e−9，是更清楚的差局部驻点证据。全84任务共76个梯度停止、6个迭代上限、2个相对目标变化停止；另外两个上限是历史42cb的Dual/Global，梯度仍为3.84e−6/1.66e−7。不能把所有失败都写成已收敛，也不能把所有失败都归咎预算。
74: 
75: 在四个成功新对象上，重建Phi反而略低于真材料Phi：公开平滑prior与噪声使真材料不必是正则化目标极小点。更小Phi不自动意味着更准确GT。原始见证见[离线可行目标诊断](RESULTS/scoring_primary/GT_FEASIBLE_WITNESS.json)。
76: 
77: ![固定六新对象的原Phi与物理LU](FIGURES/analysis/fixed_six_original_Phi_vs_LU.png)
78: 
79: 曲线只连接已付费的host/阶段真实物理验收及corrector点，辅助迭代间没有偷补forward或把代理目标当Phi。CPU曲线与原始CSV也随包保存。
80: 
81: ## 5. 真正FFT接入：已运行，成本收益有界
82: 
83: 已有NIE/CIE/FFT-TSOM先例不能记为本轮发明。新的工程入口用
84: 
85: ```text
86: UA = Qd IFFT(mask·A),
87: U*z = mask* FFT(Qd z),
88: GJ = crop[ padded free-space Toeplitz convolution(J) ].
89: ```
90: 
91: 材料读出保持原A公式，因此确实可以在FFT-TSOM上build upon。接收S和最终物理LU仍dense；没有把完整反演标成O(N log N)，也没有把二维FFT堆成三维。接线与浮点边界详见[FFT实现说明](FFT_TSOM_INTEGRATION.md)。
92: 
93: | 预指定对象/方法 | 显式Fourier CPU秒 | FFT CPU秒 | 显式NRMSE→FFT NRMSE | 完整联合结果 |
94: |---|---:|---:|---:|---|
95: | 数字1 Dual-CIE | 36.51 | 34.14 | .09656→.09592 | 均合格 |
96: | 数字1 A | 59.02 | 49.76 | .09598→.09644 | 均合格 |
97: | 旧b20d Dual-CIE | 39.48 | 38.02 | .51072→.38957 | 均失败 |
98: | 旧b20d A | 70.58 | 54.44 | .04128→.04049 | 均合格 |
99: 
100: A总CPU分别减少15.7%、22.9%，其辅助阶段分别41.76→35.51秒、53.53→43.97秒。两例FFT的物理LU反而增加，说明总时间节省来自辅助实现等实际差异，不能说所有计算项都减少。Dual-CIE在b20d上的Phi变高而NRMSE变低，再次显示残差与材料质量不同；浮点作用差异可能改变长轨迹，不能用这张图宣称FFT本身提高了信息量。
101: 
102: 在同FFT后端中，A仍比Dual-CIE多约46%和43%总CPU。这里只有两例工程比较，没有做大网格复杂度实测或六新对象FFT再确认。NIE原2016全文已取得，但独立成熟NIE重现未完成；故本轮诚实使用已审计的CIE/Dual-CIE主基线，详见[NIE审计](NIE_BASELINE_AUDIT.md)。
103: 
104: ## 6. 成本、九数字与完整负结果
105: 
106: | 主臂：21对象 | 总CPU秒 | 辅助CPU秒 | corrector CPU秒 |
107: |---|---:|---:|---:|
108: | CIE | 897.22 | 394.98 | 492.12 |
109: | Dual-CIE | 932.81 | 401.97 | 520.70 |
110: | Global-Dual | 956.57 | 466.37 | 480.51 |
111: | A | **1504.32** | **1049.64** | 444.91 |
112: 
113: A物理corrector总费用稍低，但材料内循环开销使整组CPU约为Dual-CIE的1.61倍；没有同实际CPU预算的普遍优势主张。主84任务157418 LU、3773976 complex RHS；四FFT任务6607 LU、158376 RHS。主有限四worker队列两段实际运行历时合计1140.37秒，另有修复间隔41.26秒；FFT队列55.36秒，不能用任务wall之和代替这些并行历时。生成、接线、失败、评分、画图另列[成本账本](COST_LEDGER.md)。
114: 
115: ![全部方法的材料误差与实际CPU](FIGURES/analysis/primary_material_nrmse_vs_actual_cpu.png)
116: 
117: 九数字全部图如下，保留CIE的数字7失败以及A/Global的数字5 IoU不足，不只展示成功样本。
118: 
119: ![数字1至9全部结果](RESULTS/scoring_primary/DEV9_DIGITS_all_objects.png)
120: 
121: ## 7. 科学判断与下一步
122: 
123: **有价值的经验结果**：只改变材料块的局部离流形几何，在相同数据、prior和global current控制下，额外恢复两个历史弯条；因此A可以作为有成本的补偿候选。现有数据不支持把它直接设为通用默认，更不支持论文级“稳定击败成熟Dual-CIE”的标题。FFT组合已落地，有两例有限工程收益；这属于实现贡献，不能占用已有FFT-TSOM/NIE的原创性。
124: 
125: **真正值得追的公式问题**仍是材料可解释方向与非局部传播的关系：RA=P(m)r 的精确导数为P dr+(dP)r。R5已排除仅看交叉项符号或H/L²尾就预测成功的简单解释；R6进一步显示相近局部几何可落到同一错误全局结构。局部四cell传播能改搜索路径，却未必能改变跨空腔的全局场模式，这是下一项待验证假说，当前不是定理。
126: 
127: 下一轮宜选成功7b8d与失败开缝壳这一对已公开对象，在相同已存状态上比较“材料切向正交补的传播”与“跨腔Green作用”是否产生不同有效材料步，再预注册一项简单修改并用全新材料确认。薄环另作预算/曲率诊断，不能与已近驻点的开缝壳混为同类。先回答这个有区分力的问题，不默认扫β、加NN或再堆更多同类图。
128: 
129: 本轮未新增3D成像、Fresnel实测或独立NIE复现；不以空模板伪装完成。R4 B/C的既有否证与边界保留在[R4 A/B/C报告](../research_max_struct_r4/ABC_VALIDATION_ZH.md)，R6只扩展用户指定的A。完整复现见[REPRODUCE](REPRODUCE.md)，独立审查见[audit](audit/INDEPENDENT_METHOD_REVIEW.md)。

## research/research_max_struct_r6/CHAT_HANDOFF.md

1: # MAX-STRUCT-R6 给Chat的入口
2: 
3: ## 首先读取
4: 
5: - [中文完整报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r6/RESEARCH_REPORT_ZH.md)
6: - [网页和全部图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r6/)
7: - [数学结构](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r6/THEORY_DERIVATIONS.md)
8: - [逐例数据与费用](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r6/RESULTS/analysis/ALL_CASES.md)
9: 
10: ## 已实际完成
11: 
12: 21个材料×CIE/Dual-CIE/Global-Dual/A-material＝84个主任务；另两个预指定对象×Dual/A＝4个真FFT任务。A是普通global current+仅材料块one-hop，继承R5精确重复求值缓存。固定mf3/6/10各500＋mf12共1500辅助步，共同原Maxwell corrector上限3000；576独立实未知量，无GT在线选择。
13: 
14: | 分层联合成功 | CIE | Dual | Global | A |
15: |---|---:|---:|---:|---:|
16: | 历史数字9 |8|9|8|8|
17: | 历史困难6 |1|1|1|3|
18: | 新确认6 |4|4|4|4|
19: 
20: A额外恢复两个已公开曲条：7b8d NRMSE .04686，b20d .04128；不是新盲测突破。42cb/63ab及历史环仍失败。数字5 A因IoU .6667未过联合合同，NRMSE .1142和数据拟合合格，不能只说“数字5完全失败”。
21: 
22: 新共同失败：82ac偏心薄环与40ac开缝壳。真材料可行Phi分别8.70e−6/8.55e−6，算法约.00485/.00381，残差约GT的91–103倍。说明没到达已知可行低目标区，不是本轮已证明的数据不可辨识。薄环四臂达到预算上限，开缝壳四臂梯度停止；不可混淆。
23: 
24: ## FFT与NIE
25: 
26: 可真正build upon FFT-TSOM：native Fourier action/adjoint＋非周期零填充Green，A读出不变；接收和最终物理LU仍dense。数字1、旧b20d的A CPU分别59.02→49.76秒、70.58→54.44秒，均保留恢复；同FFT的Dual仍更便宜。四个工程终点非逐位相同，不宣称通用纯速度等价。
27: 
28: NIE2016作者全文已取得，原文已有FFT-TSOM实际方法；本轮没有完成独立成熟NIE代码复现，因此主比较是可验证CIE/Dual，不用self-only冒名NIE。
29: 
30: ## 成本与判断
31: 
32: 主A总CPU1504.32秒，Dual932.81秒（1.61倍）；主157418LU/3773976RHS，FFT6607LU/158376RHS，setup/失败另列。全部新增已知物理164051LU/3932556RHS（其中12LU是16维接线小模型）。主四worker队列运行合计1140.37秒、另修复间隔41.26秒，FFT55.36秒；不要把重叠task wall之和叫历时。
33: 
34: 结论：A是有局部互补价值、已实现FFT的候选插件，还不是稳定优于Dual的默认方案；没有3D或实测新成像结论。当前最有价值的后续问题，是P dr+(dP)r如何把材料切向正交补的残差投到有效更新方向，以及四cell传播为何未改变跨腔错误结构。先对成功7b8d/失败开缝壳做状态固定的公式诊断，再新材料确认；不靠扫β或增加网络。
35: 
36: R4 B/C没有重跑，原否证保留。R1–R5报告/数组未改写。新六对象的GT现在随交付公开，之后不得重复标盲测。所有路径、失败与复现入口见同目录README/REPRODUCE/COST_LEDGER/FAILURES及audit。不要把本交接中的待验证假说当结论。

## research/research_max_struct_r7/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT-R7：低成本A+、有理材料动作与论文资格
2: 
3: ## 最终判断
4: 
5: 本轮已得到一个有实际计算收益的简单算法：**A+用Möbius状态项的二次驻点候选，加原先验的MM接受，替代昂贵的材料非线性读出。** 新30个材料上，原A和A+成功对象完全相同（21/30），A+总CPU降低约31.4%；相对成熟Dual-CIE仍多约11.2%CPU。A/A+在两个新对象上显著改善Dual/Global失败，平均NRMSE有改善证据；成功率的普遍提升尚未成立。
6: 
7: **精确有理搜索没有在所测真实失败锚点提供独有的跨盆地能力。** 实测与小型3D反演已实际完成，均没有建立A的可重复优势；3D连续物理精度和Fresnel天线校准仍有明确缺口。论文应围绕低成本条件读出、适用范围和反例，不应包装成“极点学习解决全局Maxwell反演”的统一框架。
8: 
9: 追加的20个选定家族新材料中，A+为12/20、Dual/Global各9/20、single8/20；曲条有重复收益，开缝壳仍失败居多，且壳体成功可由single解释。
10: 
11: 所有新代码、原始结果、失败、图、成本、证明、论文初稿与复现入口独立保存；R1–R6未改写。全部算法用同一公开数据和原物理Phi选候选，GT只在预测封存后评分。输入ZIP的算法建议只作为待审查材料，没有当作已运行证据或额外执行授权。
12: 
13: ## 1. 实际完成与分母
14: 
15: 30个独立材料×6方法=180个主反演，追加后选家族20个新材料×4方法=80个确认反演；另有6个相同对象×3方法的初始值对照18run、九数字A+回归9run、四披露对象×5readout的完整开发20run。30是独立材料簇，不把照明、迭代态或重复初始值当新材料。
16: 
17: 24个新对象同网格生成/反演；6个36网格生成、24网格反演，仍使用同一disk-pulse物理族。六形状各五个独立几何实现，五条件捆绑contrast/噪声/源数/孔径/网格；没有单因素因果设计。数据/先验/噪声约定见[BENCHMARKS](BENCHMARKS.md)。基线、原A先冻结；A+只用四个旧对象开发并在新对象A+反演前冻结，生成早于最终A+选择，未声称全算法冻结后才生成数据。
18: 
19: 三个公开失败锚点各做rank2/rank8搜索、原点局部对照和相同corrector续跑，另做一次rank32结构修正；两个真实Fresnel对象×3方法=6个实际实测反演；真实3D完成125/216材料DOF、375/648矢量系统的16个初始反演和16个长校正结果。它们不混入30个二维材料的统计分母。
20: 
21: ## 2. 新30对象：完整强基线比较
22: 
23: |方法|材料成功/30|同网格材料+CN99 /24|中位NRMSE|总CPU秒|LU|复RHS|
24: |---|---:|---:|---:|---:|---:|---:|
25: |成熟CIE|21/30|15/24|0.08722|1543.60|54873|1154196|
26: |成熟Dual-CIE|20/30|14/24|0.08584|1539.15|53408|1137060|
27: |Global-Dual|20/30|14/24|0.08653|1512.97|53923|1137684|
28: |原A-material|21/30|14/24|0.08177|2495.14|54384|1148820|
29: |A+|21/30|14/24|0.08233|1711.34|52156|1103040|
30: |原物理L-BFGS-B|10/30|9/24|0.31489|838.18|47754|1009764|
31: 
32: 全体方法都输出576个真实材料自由度；原current Fourier/DPCS表示未按形状缩减。
33: 
34: “材料成功”沿用幅度归一化NRMSE≤.15、预测自身Otsu IoU≥.7、合法box；数据CN99另列。成功数不是统计显著性的代名词。跨网格六例的sigma不包含模型误差，不能称完整误差证书。全部连续逐例指标见[评分CSV](RESULTS/scoring/primary/SCORES.csv)，全部180张材料图分为[1](FIGURES/primary/all30_group_01.png)、[2](FIGURES/primary/all30_group_02.png)、[3](FIGURES/primary/all30_group_03.png)、[4](FIGURES/primary/all30_group_04.png)、[5](FIGURES/primary/all30_group_05.png)、[6](FIGURES/primary/all30_group_06.png)，没有只挑正例。
35: 
36: ![全部对象的连续指标](FIGURES/primary/continuous_material_distributions.png)
37: 
38: 与同global-current的Global-Dual相比，A平均NRMSE差−.04382，对象配对bootstrap95%区间[−.10088,−.00568]；中位配对差的区间跨0，成功率差+1/30的区间[−.0667,+.1333]。平均误差改善主要由少数大幅恢复贡献，不能据此写成各几何稳定胜出。无多重比较调整，区间仅针对本冻结混合条件下的材料抽样，不是普遍Maxwell保证。
39: 
40: A+−A平均NRMSE差+.000604，区间[−.000923,+.002872]；同样21个成功对象不等于未来等价证明，也未预注册非劣性margin。A+每对象CPU少26.13秒，区间[−33.82,−18.62]。与Dual相比，A+平均NRMSE差−.03778[−.09292,−.00309]，每对象CPU多5.74秒[1.92,9.54]。完整[统计审查](RESULTS/scoring/analysis/ANALYSIS.md)及[独立次要A+对Dual比较](RESULTS/scoring/analysis/SECONDARY_APLUS_VS_DUAL.md)保留。
41: 
42: 两个具有实质幅度的正例：
43: 
44: |新对象|Global-Dual NRMSE|Dual NRMSE|A NRMSE|A+ NRMSE|A+数据残差|解释|
45: |---|---:|---:|---:|---:|---:|---|
46: |12曲条，contrast12/六源|.73511|.73500|.03965|.03914|.000942|普通single CIE也失败，支持本材料路径增量|
47: |20开缝壳，contrast18/半孔径|.29655|.29636|.06981|.06872|.001057|single CIE已约.07022成功，不能称A独占|
48: 
49: 相反，09细弯结构A+ IoU=.69663、Global=.71264，NRMSE都约.095–.098；这是阈值边界loss，不是灾难性材料失败。17的A+也未过IoU，而single通过；24跨网格曲条single失败但Dual/Global/A/A+通过。单β、双β和global的短板不同，不能只挑最弱一臂。
50: 
51: ![含开缝壳正例及其余四例的完整组](FIGURES/primary/all30_group_05.png)
52: 
53: ## 3. 相同实际成本改变了什么判断
54: 
55: 完整endpoint上A约为Dual的1.62倍CPU，A+约1.11倍。单看相同3000辅助/3000corrector上限不公平；direct6000没有辅助工作，也常很早停在差点。30个direct均报告数值收敛（28梯度、2相对函数），但材料只10/30通过；低CPU不等于成功反演。
56: 
57: 逐对象取双方已观测重叠CPU/LU/RHS窗口，按不超预算的最后已支付Phi incumbent评分，不按GT挑迭代。A对Global在成对共同最大CPU时为17/30对20/30，平均NRMSE差区间跨0，额外辅助工作抵消了一部分endpoint收益。A+对A在共同CPU为21/30对17/30。A+对Dual的共同最大CPU/LU/RHS为22/30对20/30，平均NRMSE差区间仍为负；但某些中间检查点的IoU通过会在最后更低Phi后丢失，不能把22改成最终成功数。
58: 
59: 预算材料为float32、每10次物理记录与稀疏辅助stage端点，存在检查点分辨率和unused budget。等LU没有给昂贵材料求根免费，等CPU也不是孤立硬件wall基准。完整对象级曲线和次要比较见原始统计文件。
60: 
61: ![材料与实际费用frontier](FIGURES/primary/pairwise_actual_cost_frontiers.png)
62: 
63: ## 4. A+数学与开发：保留收益，压缩代数成本
64: 
65: A只在材料读出使用RA=β(I+DzGoff)(J−DzH)，其中z=m/[1+(β−Gii)m]。零对角Goff使其对单个z严格仿射，多个坐标仍双线性。实材料的z在合法圆弧上。原prior用相切MM上界保留，四颜色独立批处理。
66: 
67: 精确正则化条件最小化需要五次驻点求根；最终A+只取**状态**二次驻点、box端点及旧值，用完整原prior-MM接受，因旧点在集合中仍有固定电流loss非增。它不是正则化条件全局最小、不是整体Phi单调，也没有消除Maxwell非线性。详细[推导](THEORY_DERIVATIONS.md)和[方法消融](METHODS_AND_ABLATIONS.md)列出适用条件。
68: 
69: 四个披露对象原A CPU284.21秒，最终A+197.78秒，辅助CPU193.72→104.94秒；7b8d单例仍略慢4.8%。所有20开发run保留，60,000次readout留存loss上升0；最终quintic根为0。projected8虽然局部下降，却使b20d从Phi约6.47e−6退到.005832，说明便宜内下降与外路径保持不同。
70: 
71: 封存后GT评分：最终A+和原A在四个开发对象都3/4材料通过，两旧弯条恢复保留；projected只2/4。九数字A+为8/9，旧Dual9/9保留为强历史参考；不把新FFT与旧dense差别当纯readout因果。六重复初值m=1时A/A+均5/6、Dual4/6，但这些不能新增分母，也不能GT选更好初值替换主表。
72: 
73: ![全部开发版本材料恢复](FIGURES/primary/disclosed_four_object_aplus_development20.png)
74: 
75: ## 5. 精确有理搜索的能力与边界
76: 
77: 背景resolvent给出的ΔF=OP(I−DKP)^−1DEP是已有标准精确有限响应。R7实现两cell四端口Gram和一般q小矩阵的精确端口二次型，保留完整背景多次散射；不是Taylor升阶。动作根据数据残差、源/收可见性、材料状态和合法box选择；两块各四cell是rank8，32cell两幅度是rank32，费用都按真实秩记录。
78: 
79: 在开缝壳、薄环和42cb上，所测rank2/8/32动作均未产生独有材料恢复：原12个完整校正端点（含原anchor与local控制）及另3个rank32端点共15个，材料通过0/15；它们只涉及3个已披露失败锚点，不是15个独立材料或15个新proposal。开缝壳的有限改善只有约1e−11量级；自由face的最小计算Hessian值约4.34e−9不包含活动约束方向，不是局部极小证书。薄环/42cb续跑原anchor已达到Phi .00492946/.00218249，有限动作之后没有更好结果。预算不足/慢收敛与真正差驻点因此不能合并。
80: 
81: ![有理动作与原anchor续跑](FIGURES/mechanisms/rational_vs_continuation.png)
82: 
83: 小响应精确不制造有用动作。神经网络若只拟合这些同一动作的响应，不会创造不存在的更好端点；学习更丰富proposal仍是未测问题。rational串行预算含前序对照，setup子计时也不是全部准备时间，不能据此宣称matchedCPU全局化优势。端点验证、所有拒绝动作和源/收列求解均记账。
84: 
85: ## 6. 实测和真正3D：实际做了，优势尚未成立
86: 
87: Fresnel四个真实介质TM文件已从官方IOP取得；两个对象在2GHz各跑Dual/global/A，正确处理total−incident、时间共轭、原几何、缺失通道及target-absent标定，留出1530个实际通道预测。单/双圆柱训练残差.14193/.10062，留出.26755/.23444，三臂各自材料完全相同，且触box、双圆柱corrector达到上限。起始corrector目标等于共同host，与拒绝辅助候选一致，但aux候选数组未保存，不能严格认证argmin。原horn照明近似与模型失配未解决，不报告虚假逐像素GT误差。详见[Fresnel审计](FRESNEL_AUDIT.md)。
88: 
89: ![真实测量的相同终点与模型失败](FIGURES/fresnel_material_endpoints.png)
90: 
91: 3D使用真正dyadic Green与3×3自项、六个严格横向平面波、外部矢量接收，125/216个独立实材料未知量；不是二维堆叠。16个初始反演和16个共同长校正显示主要排序由预算解释，充分校正后各方法接近：A在125-grid球/双球/非球括号NRMSE约.17844/.19965/.13300，没有持续恢复优势。三维采用其独立的quadratic+l2 prior、噪声/归一化和相对材料L2 NRMSE；“同prior/同Phi”只指各benchmark内部，不能将三维数值与二维幅度NRMSE直接合并。此处是vector-adapted joint/dual-like/A，不是成熟3D CIE，也未把二维A+强行当张量解析读出。
92: 
93: 独立Mie比较暴露binary几何alias：误差6.58%→27.59%；一次固定体积分数修正变成8.40%→8.20%，仍不支持连续3D成像优势。数据、初始失败、修正及成本均在[3D合同](BENCHMARK_3D.md)。
94: 
95: ![三维长校正后的实际材料](FIGURES/vector_continuation/vector_continued_slices.png)
96: 
97: ## 7. 是否达到论文目标
98: 
99: Q1：A属于材料依赖的局部本构残差metric改变。固定全局current对照说明它确实可改变非根材料路径；新30有平均误差和少数大幅恢复收益，但成功率泛化区间、single基线、跨网格、实测和3D边界都限制普遍主张。它不是新增观测或普遍全局求解器。
100: 
101: Q2：A+已把昂贵quintic/通用内优化简化为圆弧状态候选与原prior接受。30新对象保持所观测成功集合，CPU降低31.4%，这是本轮最明确的算法正面成果；不能无margin称理论非劣性。
102: 
103: Q3：有限有理响应提供便宜准确的动作评价，但本轮未带来独有逃逸或额外材料成功。动作空间与预算的失败边界比极点拟合阶数更关键。
104: 
105: Q4：本轮解析读出/evaluator不需要学习；NN/MLP/LNO训练NOT_RUN，不意味着学习proposal普遍无用。没有凭同数据变换解除不可辨识性。
106: 
107: Q5：当前最合理稿件是**A-material低成本条件读出与适用边界的二维算法研究**，附真实测量/3D负证据；不是“极点搜索+学习”的全局统一论文。已提供完整[PAPER_DRAFT](PAPER_DRAFT.md)，但IEEE TAP/Inverse Problems层面的广泛物理恢复优势尚未建立。下一步应先解决真实照明校准和独立物理误差、确认可重复作用区域，而不是增加网络或Taylor阶数。
108: 
109: ## 8. 成本、审查、复现与后续确认
110: 
111: 主180任务共316,498个576维LU、6,690,564个复RHS；其CPU合计9,640.37秒，task wall之和不是并行历时。开发、生成、验证、有理分支、实测、初始化、九数字及绘图/发布费用另外列，不混算不同维度LU。完整[COST_LEDGER](COST_LEDGER.md)保留未测RSS/金额/能耗和审计成本缺口，未知不填零。
112: 
113: 内部独立实现/评分与claim审查分别记录；实现A+的worker在claim review中披露作者角色，不能冒充外部同行评审。引用工作和证据等级见[PRIOR_ART](PRIOR_ART_AND_NOVELTY.md)。没有把标准MM、Woodbury、Möbius几何、GENOME、LSL或FFT当首创。
114: 
115: 基于主30的两个正例，已在完全冻结算法下另生成曲条/开缝壳各10个新对象，比较A+、single-CIE、Dual-CIE、Global-Dual共80任务，作为**事后选定家族**的独立确认；另列全部结果与费用，不并入原30分母，不重新选择已有对象。该扩展回答局部收益是否重复，不能把某个样本数当作科学结束gate。
116: 
117: ### 追加20家族确认：收益有重复，普通基线仍重要
118: 
119: 这20对象是在原30评分后选定家族，再按固定分布独立生成；算法、预算、prior和规则不改。原先3臂后依据原30结果、在targeted GT未开之前补single-CIE，全部80任务完成，没有筛样本或替换原30分母。
120: 
121: |方法|材料成功/20|材料+CN99/20|平均NRMSE|CPU秒|LU|复RHS|
122: |---|---:|---:|---:|---:|---:|---:|
123: |A+|12|9|.28525|1089.69|43594|849516|
124: |Dual-CIE|9|8|.35356|1061.57|46876|909312|
125: |Global-Dual|9|8|.35604|1027.61|45126|870528|
126: |single-CIE|8|6|.37375|855.19|43013|836184|
127: 
128: 曲条10例A+/Dual/Global/single为9/8/8/5，开缝壳为3/1/1/3。壳体中A+与single的成功集完全相同；不能把相对Dual的两个壳体增益说成A独有。曲条target03仅A+成功，是新的有条件算法路径增量，不是整个家族保证。
129: 
130: A+对Dual平均NRMSE差−.06831，配对20对象bootstrap95%区间[−.14213,−.01242]；成功率差+.15，区间[0,.30]，下界触零。CPU均值差每对象+1.41秒，区间[−1.65,+3.86]，不能以区间跨零宣称等价。对Global相同成本范围也保留平均误差改善；对single在相同CPU时仅10/20对8/20，平均差区间[−.16501,+.00548]跨0。完整endpoint对single更好，不能冒充等CPU优势。
131: 
132: A+的8个失败（1曲条、7壳体）全部保留，均未过CN99；2壳体报告数值收敛，6例达预算。每例都存在离线更低的GT可行Phi，不能直接称为信息不可辨识。详见[条件确认统计](RESULTS/targeted_analysis/ANALYSIS.md)和[成本及失败](RESULTS/targeted_analysis/COST_AND_FAILURES.md)。全部20×4材料图：[1](FIGURES/targeted/all20_group_01.png)、[2](FIGURES/targeted/all20_group_02.png)、[3](FIGURES/targeted/all20_group_03.png)、[4](FIGURES/targeted/all20_group_04.png)。
133: 
134: ![追加家族确认中的全部最后五例，含失败](FIGURES/targeted/all20_group_04.png)
135: 
136: 这次扩展支持“低成本A路径的部分曲条收益可以重复”，也否定“开缝壳普遍更好”。它仍不是独立物理/实测高contrast优势证据。原30与后选20分别报告，不计算一个漂亮的50对象混合胜率。
137: 
138: 
139: 代码入口与获取大型原始轨迹的方法见[REPRODUCE](REPRODUCE.md)；所有新GT评分后公开，之后不得重复标新盲测。失败链见[FAILURES](FAILURES.md)。

## research/research_max_struct_r7/CHAT_HANDOFF.md

1: # MAX-STRUCT-R7 给Chat入口（2026-11-04）
2: 
3: 先读[中文完整报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r7/RESEARCH_REPORT_ZH.md)、[网页与所有图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r7/)、[推导](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r7/THEORY_DERIVATIONS.md)、[论文初稿](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r7/PAPER_DRAFT.md)。
4: 
5: ## 已实际执行
6: 
7: 原30个新材料×6方法=180；后选家族20新材料×4=80，分母分开。另有初始化18、旧九数字9、公开开发20；三个旧失败锚点15个校正端点（含控制）；真正矢量3D初始16+长校正16；真实Fresnel2对象×3=6。没有NN/GT在线选根，所有新GT评分在预测全封存后。R1–R6未改写。
8: 
9: |原30材料成功|CIEsingle|Dual|Global|A|A+|direct物理|
10: |---|---:|---:|---:|---:|---:|---:|
11: |成功/30|21|20|20|21|21|10|
12: |CPU秒|1543.60|1539.15|1512.97|2495.14|1711.34|838.18|
13: 
14: A+是RA=β(I+DzGoff)(J−DzH)的条件Möbius状态二次驻点+原prior-MM接受；不是精确正则化条件最小、不换prior，保留完整loss fallback。普通global-current不改。相比A，原30成功集相同，CPU少31.4%；相比成熟Dual仍贵11.2%。A-global平均NRMSE改善CI负，但成功率差CI跨0，不能说普遍胜出。
15: 
16: 原30两个大阳性obj12曲条/20壳体；obj20 single也成功，obj12 single失败。obj09是IoU边界loss（.6966vsGlobal.7126），obj17是相对single边界。等成本按公开Phi incumbent检查点，不按GT选迭代；较低Phi可以丢失材料IoU。保留预算分辨率/float32限制。
17: 
18: ## 最后的家族确认
19: 
20: 在原30后选家族、再生成20独立材料，所有算法固定；不要和原30合并成事前50对象总体胜率。A+/Dual/global/single为12/9/9/8；曲条9/8/8/5，壳体3/1/1/3。壳体A+和single成功集合相同；曲条target03仅A+成功。A+对Dual平均NRMSE差−.06831[-.14213,-.01242]，成功率差CI[0,.30]。总CPU1089.69/1061.57/1027.61/855.19；对single等CPU mean误差CI跨0。八个A+失败完整保留，均有更好离线GT-Phi见证。
21: 
22: ## 其他路线和资格
23: 
24: - 精确rank2/8/32 finite response端点正确，未获得优于同corrector原anchor续跑的独有材料恢复；15端点0成功，只涉及3个旧锚点。低秩动作不够不是公式错误。标准Woodbury/Gram/MM/Möbius/GENOME不是新发现。
25: - projected8丢失b20d，说明内块下降不保外路径；最终A+保留两旧bar材料恢复。九数字A+8/9，旧Dual9/9历史引用，不重新标盲测。
26: - 真3D矢量DDA375/648维、125/216材料DOF已经实际反演；共同长corrector后方法接近，无持续A优势。binary Mie refinement失败，一次体积分数修正后误差仍约8%；不能称连续3D资格/成像优势。3D用独立prior和相对L2 NRMSE，不是二维A+测试。
27: - 官方Fresnel原始文件已取得，六次真实反演+留出预测完成；三臂同终点、触box、天线校准不足，无可信像素GT或A优势。原数据/全文不再分发，公开获取脚本及derived结果。
28: - NIE/LUCIE未独立成熟复现；MLP/LNO/LSL未运行。不把未做算法伪装baseline，也不证明学习普遍无用。
29: 
30: 最适合写A-material低成本条件读出与适用边界论文，不能称全局极点学习框架或泛3D/实测胜利。下一步优先独立物理/实测照明校准，保留有条件曲条收益并解释残余错误结构；别加无依据gate、扫β或用大网络掩盖。
31: 
32: 最终配置configs/aplus_final.json；完整复现REPRODUCE.md；费用COST_LEDGER.md、全部否证FAILURES.md、内部审查audit/。大型原始轨迹在[本轮Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r7-20261104)，所有分卷解压到仓库根。原30和后20的GT现在公开，此后都只是开发数据。

## research/research_max_struct_r8/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT-R8：神经初值如何进入 CIE/A+，以及物理模型如何限制成像
2: 
3: 本轮实际执行、评分和研究分析已完成；所有新增工作位于R8，R1–R7保持不变。750项原比较与126项机制扩展先分别完整封存，再统一开测试真值评分。已运行、探索性尝试和未完成的原尺度任务分别标记。
4: 
5: ## 先给判断
6: 
7: **本轮已经完成轻量NN训练、接入成熟FFT-TSOM/CIE及sourcewise/global A+、原论文Austria合同实验、Fresnel真实测量反演和真正的3D矢量材料反演。现有跨合同证据不足以支持“A+普遍优于成熟CIE”。** 更有价值的认识是：材料初值携带的信息会被硬DPCS、早期Fourier截断以及校正前的候选选择改变。一个好的NN材料初值不一定能完整进入辅助current迭代；更精确地传递current，也不保证到达更好的材料终点。
8: 
9: 此外，Fresnel发现了一个实质性物理约束问题：沿用合成任务的正contrast下界，会排除实测ROI中的空气。在同一已支付锚点上仅将下界从.02改为0，8GHz留出误差从.266降到.136，但内部高峰并未消失。这是可行材料集合的修正，不能记为A+或NN的收益。
10: 
11: **独立确认的主要正结果来自NN接直接原物理优化，而非NN接A+。** 同分布/未训练家族/跨网格的平均幅度NRMSE分别为.0965/.1624/.0843；NN接global A+为.2961/.5939/.3293，sourcewise A+为.3218/.5916/.3281。固定APCS偏置在部分对象改善，但没有稳定保住NN直接物理支路的收益。原A+在历史弯条上的作用仍存在，不能推广成普遍优势。
12: 
13: NN直接物理相对均匀初始化的同一物理算法改善了材料质量，但平均花费更多时间，不能称该比较中加速；相对本轮CIE/A+辅助管线，它在每个硬件/cohort层的最终平均CPU与NRMSE都更低。相同已支付预算的比较仍有反例和宽区间，不主张逐对象普胜或通用SOTA。NN训练分布是额外先验，费用和物理适用域必须保留。
14: 
15: ## 1. 实际完成了什么
16: 
17: |证据组|实际工作|不能据此宣称什么|
18: |---|---|---|
19: |轻量初始化器|1536训练、192验证；24737参数，120epoch；另有同容量Born特征网络及论文域适配网络|通用神经逆散射发明，或训练先验不算额外信息|
20: |独立确认|42新材料，24同分布、12未训练家族、6跨网格；750原任务+126机制扩展|876个独立样本，或开发对象重新变盲测|
21: |开发诊断|8验证材料，48原始轨迹+32仿射偏置对照|靠开发胜例证明泛化|
22: |历史回归|4公开困难材料、21新端点，复用19旧端点|新的独立确认成功率|
23: |原论文对象|59个实际端点，包含150²生成/64²反演Austria、较长corrector及域偏移修复|59个独立材料，或独立NIE实现就是作者成熟软件|
24: |Fresnel|真实2001单/双圆柱、2005 FoamTwinDiel；校准、逐源接收、五方法、NN迁移、NIE、多频及air0诊断|逐像素GT精度，或同算子synthetic冒充实测|
25: |真3D|TensorDDA全矢量，52个小尺度端点；独立Mie与网格、另3λ原电尺寸实际尝试|2D堆叠、原论文60³/30³已完整复现，或3D NN已训练|
26: 
27: 明细分别见[数据合同](BENCHMARKS.md)、[原论文实验](PAPER_BENCHMARKS.md)、[后续分析](PAPER_FOLLOWUPS.md)、[Fresnel](FRESNEL_REPORT.md)、[三维](VECTOR3D_REPORT.md)。同对象的算法、噪声、初始化和网格变体不增加统计对象数。
28: 
29: ## 2. 从公式得到的实际机制
30: 
31: ### 2.1 NN材料与辅助电流不是同一个初值
32: 
33: 原全波模型为`J=X(Einc+GJ)`，观测`F=S(I−XG)⁻¹XEinc`。FFT-TSOM的电流族为`J=Jd+UA`，其中`U=(I−Pd)F_k`，`Pd`是保留接收奇异子空间的正交投影。NN只给材料m0时，原辅助初始化仍是Jd；直接把物理电流jp投影到这个截断族，也会损失一部分信息。
34: 
35: 投影位移精确分成`δjD=Jd−Pd jp`和`δjN=Pu jp−(I−Pd)jp`。第一项来自硬DPCS与NN端口预测的不一致；第二项来自Fourier截断。两项在电流欧氏范数中正交，但Maxwell本构缺陷是`(I−X0G)(δjD+δjN)`，传播算子之后不再保证正交。因此“电流更接近物理解”不推出“本构残差更小”，更不推出材料恢复更好。
36: 
37: 本轮实现了简单接法：固定`Ja=Jd+(I−Pd)jp`，再优化`J=Ja+UA`、`A0=0`。它保留初始物理电流的整个APCS分量，优化未知量不增加，最终完整APCS可行族也不变；早期截断仿射族发生平移。jp复用已有初值forward，新增投影与存储照实收费。该接法同样适用于普通Dual-CIE，不能把其全部收益归给A+。
38: 
39: 这不是新的Maxwell恒等式。正交投影、学习current再接TSOM等均有先例；本轮可能成立的增量是明确定位投影损失、展示简单修复的作用及反例，并检验总成本后的恢复效果。[完整推导](THEORY_DERIVATIONS.md)保留条件和反例。
40: 
41: ### 2.2 校正前选根并不保持校正后的排序
42: 
43: 设C为相同original-Maxwell corrector。即使C对每个输入都下降，`Phi(a)<Phi(b)`也不意味着`Phi(C(a))<Phi(C(b))`。辅助阶段可能选到一个当时目标更低、但后续比NN直接物理优化更差的起点。
44: 
45: 因此既运行NN+CIE/A+，也运行NN+直接物理优化，并在开GT前冻结两支路终点按原Phi选择的组合。两条支路全部收费。这是标准的保留候选控制，不是免费增加一次corrector，也不是新的全局收敛定理。
46: 
47: 开发对象0003上，原NN global A+的幅度NRMSE=.4249，截断联合初始化=.7420，固定偏置后的Dual/global A+/sourcewise A+=.0315/.0326/.0326；然而NN直接物理优化=.0297且更便宜。0002、0005的偏置结果改善来自最终选择stage−1，实质上保留NN初值再校正，不能算辅助current解的贡献。
48: 
49: 反例0006更重要：原联合初始化=.0527，偏置Dual/global/sourcewise却=.6846/.5844/.6036。全部8个初始本构缺陷都改善，材料终点却不统一改善。这否定了“初始缺陷降低可单独预测恢复”的推断。没有用这些开发结果修改确认集评分或挑选测试对象。
50: 
51: ### 2.3 A+与sourcewise/global的作用位置
52: 
53: A+保留普通CIE current更新，仅对材料读出使用R7的Möbius条件结构和MM接受。`z_i=m_i/[1+(β−Gii)m_i]`下，零自作用传播使单坐标状态残差对z_i条件仿射。其接受只控制固定电流材料子问题，不认证原Phi或GT误差。
54: 
55: 本轮实际比较sourcewise/global × ordinary/A+，每一格同时有均匀与NN材料初始化。主确认还保留原A-material及直接物理优化。用户提出的NN接入sourcewise/global A+已经作为完整42对象比较执行，并非只在好看的对象上展示。
56: 
57: ## 3. 神经初始化与独立确认
58: 
59: 以下为**完整预算端点的平均幅度NRMSE；括号为继承复合材料标签的成功数**。标签要求NRMSE≤.15、Otsu IoU≥.7及box合法，只用于评分。不同总体不合并；相同端点预算上限不等于相同实际计算费用。
60: 
61: |方法|同分布24|未训练家族12|跨网格6|
62: |---|---:|---:|---:|
63: |均匀→原物理|.2497（2/24）|.4047（1/12）|.1959（2/6）|
64: |NN→原物理|**.0965（9/24）**|**.1624（6/12）**|**.0843（2/6）**|
65: |均匀→single CIE|.3443（8/24）|.6101（3/12）|.3275（2/6）|
66: |均匀→Dual-CIE|.3684（7/24）|.6383（1/12）|.3282（2/6）|
67: |均匀→Global-Dual|.3569（8/24）|.6254（2/12）|.3287（2/6）|
68: |均匀→原A-material|.3389（8/24）|.6359（2/12）|.3288（2/6）|
69: |均匀→global A+|.3406（8/24）|.6375（1/12）|.3283（2/6）|
70: |NN→Dual-CIE|.3282（7/24）|.5678（2/12）|.3283（2/6）|
71: |NN→Global-Dual|.2892（8/24）|.6048（1/12）|.3283（2/6）|
72: |NN→global A+|.2961（7/24）|.5939（2/12）|.3293（2/6）|
73: |NN→sourcewise A+|.3218（7/24）|.5916（2/12）|.3281（2/6）|
74: |NN固定偏置→Dual|.2901（7/24）|.5230（4/12）|.2599（2/6）|
75: |NN固定偏置→global A+|.3140（5/24）|.5183（4/12）|.2280（2/6）|
76: |NN固定偏置→sourcewise A+|.3129（6/24）|.5016（4/12）|.2283（2/6）|
77: 
78: 其余完整臂及全部逐对象值保留在[630主比较表](delegated/warm_analysis/PRIMARY_630_object_rows.csv)、[126偏置结果](delegated/warm_anchor_analysis/CONFIRMATION_SCORES.json)和[全部42材料图](FIGURES/warm_material_all_42.pdf)。NN原始初值的平均NRMSE为.1469/.2190/.1205，后续直接物理优化改善到上表值；不是只把网络输出当作最终重建。ID中15个NNdirect未达复合标签的对象，有12个仅IoU未过，因此不能把标签失败全称为高连续材料误差。
79: 
80: 固定偏置的结果与开发预测相符的是“有条件修复”，不支持“保留更多物理current必然更好”：未训练家族成功从原NN A+的2/12升到4/12，但仍低于NNdirect的6/12；ID global A+反而从7/24降到5/24。普通Dual也得到相近修复，不能归为A+专属结构。
81: 
82: ### 3.1 一个具体的路径反例
83: 
84: 预先选定展示对象0027中，NN初值Phi=.30965，global A+给出的校正前候选Phi=.01565，看起来更好。但NN直接物理优化在1814 accepted步按PGTOL正常停止，最终Phi=.00042866、NRMSE=.12082、IoU=.84967；global A+校正3000步到上限，最终Phi=.00046350、NRMSE=1.01370。sourcewise A+同样失败；固定偏置Dual在2789步达到数值PGTOL仍NRMSE=1.06415。
85: 
86: 该例中NN直接支路实际少于3000步，不能简单用它的6000步上限解释优势；它支持“校正前更低Phi可以丢失更好后续轨迹”。但没有证明辅助支路永远不能逃出某盆地。全组的严格符号差计数可能包括浮点微差，不能全部称为实质失败。
87: 
88: ![预定多孔对象：NN原始初值与真实下游端点](FIGURES/neural_initial_to_final_fixed_r8_test_0027.png)
89: 
90: ### 3.2 同预算与完整成本
91: 
92: 按对象取两方法最终实付成本的较小者作为共同上限，并只用该上限内已保存的原Phi incumbent，不按GT挑迭代。CPU、LU、complex RHS各自比较；稀疏checkpoint造成的未用预算与间隔均公开。NNdirect相对均匀direct的五个硬件/cohort层，同CPU上限下平均NRMSE差均为负，但它的最终CPU是均匀direct的1.06–2.74倍：均匀支路经常较早停在较差结果，不能把它的失败快停当作完成相同质量任务。
93: 
94: 相对辅助算法，NNdirect在五个硬件/cohort层的最终平均CPU约为：ID Mac32.03、ID i969.77、未训练家族Mac39.12、i9128.52、跨网格Mac30.82秒；对应均匀Dual为62.12、144.84、69.89、159.74、62.65秒。该比较包含部署特征/推理，未把训练当免费。后续补充的十个直接对照属于**看到结果后的描述性成本分析**，不是新预注册显著性检验。同CPU上限下也有小反例：ID/i9 NNdirect对NN Global-Dual平均NRMSE差+.00281，95%描述区间[−.06069,+.06401]；不能称普遍优势。不同硬件层包含不同对象，不能把层间差异归因于硬件。
95: 
96: 完整区间及费用见[数字备忘录](delegated/warm_analysis/FINAL_NUMERIC_MEMO.md)与[直接物理成本对照](delegated/warm_direct_budget_analysis/REPORT.md)。许多小分层区间跨零，未经多重比较修正；结论限定在本次冻结对象和模型。
97: 
98: ### 3.3 保留候选也不能用低Phi代替材料选择
99: 
100: 两条完整NN支路按最终原Phi选择的global A+组合，在ID达到NRMSE=.07754、13/24，比NNdirect的.09647、9/24改善；但未训练家族变为.34498，跨网格变为.22882，均比单NNdirect差。费用是单NNdirect的2.38–3.17倍。Global-Dual组合有非常接近的结果，不能把组合收益归为A+唯一作用。
101: 
102: 这不是GT挑根或隐藏了免费corrector：两个候选在开GT前已按Phi封存选择，两条支路全部收费。负结果表明，末端较低Phi也不能保证正确材料；不能为了保住算法主张改成离线best@K。
103: 
104: 额外120项只覆盖预定8对象（ID4、未训练家族2、跨网格2），实际完成均匀1、Born、DPCS、随机和Born-only网络消融。Born-only NN接global A+在ID4上4/4，physics NN为3/4，但两者在这2个未训练家族对象上都0/2；不据这个小子集回选正式权重或宣称42对象优势。全部负例与训练费用保留。
105: 
106: 网络输入是公开测量的Born/DPCS反投影、场能量和采集参数，不含GT、支撑、形状名或真实contrast。1536训练与192验证按独立材料种子划分，12个未训练家族对象另外保留。NN提供了训练分布先验，即便最后显式prior相同，也不能称没有额外先验。
107: 
108: 原physics网络训练CPU387.69秒，含自身训练/验证生成的边际训练成本402.56秒；Born及域适配开发另计。推理成本只是总成本的一小部分，真正要看完整辅助与物理迭代是否节省，以及节省能否摊销训练。训练耗时不能只取最佳epoch之前的一小段。
109: 
110: ## 4. 历史A收益保留了吗
111: 
112: 在7b8d、b20d两条已公开困难弯条上，NN global A+的幅度NRMSE=.04559/.04096，NN sourcewise A+=.04551/.04179，都保留原有恢复。NN Dual-CIE两例未恢复；NN直接物理优化仅恢复b20d。42cb与开缝壳40ac仍共同失败。
113: 
114: 这说明在这些历史对象上，NN没有自动替代A+改变轨迹的作用；它们存在互补的机制证据。但4对象不能证明普适优势；旧R6与新R8实现有FFT数值细节差别，严格初始化比较以同版本新增臂为准。[全部四对象及失败图](HISTORICAL_WARM_STARTS.md)不做成功筛选。
115: 
116: ## 5. 原始NIE/CIE对象与跨物理泛化
117: 
118: Austria不是仅画一个近似图案：实际采用论文2m区域、20平面波、40接收、10%复RMS噪声，并运行150²生成/64²反演。εr=4、400MHz下，CIE/Dual/global A+/sourcewise A+/独立NIE的幅度NRMSE约.441/.427/.501/.440/.441。数据残差接近.10，但内部环伪影和接近contrast20的峰值仍存在，真contrast只有3。
119: 
120: 原合同4096个独立实材料变量只有800复观测。Jacobian实核维数至少2496；恒秩、内点等条件下可讨论局部等数据族，但不能据此证明该有空气活动边界的特定Austria真值存在同box多解。跨网格误差又与噪声共同存在。这里越来越低的数据残差不能独立决定材料选择。
121: 
122: 原NN在论文采集上域偏移明显，因此实际补训256个论文域随机材料和64验证对象，包含不同家族但不含Austria真值。25025参数网络训练80epoch，与原训练数据平衡混合。它改善若干早期A+结果，但较长共同校正又可能降低residual、恶化材料；没有因此宣称域适配解决了问题。原确认权重未被替换。
123: 
124: 独立NIE按原式和停止规则实现，但不称作者认证的成熟代码。原论文预算、坐标、prior与适配差异均在对照报告中列明。
125: 
126: ## 6. Fresnel为什么失败，修正到哪一步
127: 
128: 本轮使用真实原始测量，而非合成替代。逐源mask和adjoint保留全部合法通道；用target-absent数据拟合满足背景方程的入射场，替换R7单一柱面波近似。单/双圆柱空场留出误差约.025/.024，旧近似约.391/.384；Foam约.115，旧近似约.644。独立圆柱/多体散射检查量化了网格误差，不能把校准拟合本身叫物理验证。
129: 
130: Foam原五方法的8GHz留出误差都约.266，材料都有不可信的环/尖峰。A+把A的材料内评价11135次降到2710次，辅助CPU76.93秒降到57.74秒；这只是该合同上的局部费用改进，不能掩盖共同的材料问题和总费用。真实数据无可信逐体素GT，因此不报告伪造的NRMSE或IoU。
131: 
132: 固定NN迁移到2001单圆柱后几乎回到相同终点，没有建立总成本收益。双圆柱存在不同材料终点而测量预测很相近，但这不是严格连续多解证书。独立NIE使用air0、无额外prior、无共同corrector的另一个合同，较平滑但留出误差约.348，不能直接称材料更准。
133: 
134: 从相同锚点加入6GHz数据，两频留出RMS从.3033降到.2324，6GHz改善、8GHz略恶化；峰值虽降低仍偏高。这是新增信息的收益与权衡，不能记作同数据A+提升。
135: 
136: 随后发现原[.02,20]下界排除了真实空气。保持同一锚点、prior、校准、通道与上界，只有下界改为0，300步物理校正将8GHz train .2255降到.0816、heldout .2660降到.1361；6GHz事后预测也改善。外缘区域98.27%可回到0，说明先前人为背景具有实质影响。但峰值从7.639升到7.964，内部尖峰仍未解决。此项为事后模型诊断，不能声称材料已恢复。
137: 
138: 进一步从共同初值.1重新运行air0 Dual-CIE和A+，保持原3000辅助/1000校正预算、prior、校准、通道和各自PR规则。两支8GHz留出分别.137237/.137128，CPU741.18/750.06秒，均达到1000步上限；材料相对差2.35%、实际通道预测差.1215%，峰值仍约6.93/6.91。修正后的材料图接近，没有可信A+恢复增量。两支都改善原合同的拟合，说明首先需要修正物理可行集合。现25项实测运行完整保留，总已测CPU10603.74秒，失败成本另计。
139: 
140: ![修正空气下界后的同初值对照，统一色阶且无像素GT](FIGURES/fresnel/air0_corrected_comparison.png)
141: 
142: ## 7. 三维不是二维堆叠，但还没有稳定优势
143: 
144: 三维采用真实dyadic Green、三分量current/field、60个外部电偶极及矢量接收，包含实材料四立方体、空腔/棒、有损包覆体。实scalar体素在sI3自项下可以推导A+条件结构；复介电不能直接使用实圆弧，本轮采用合法复参数读出，未假装复A+已完成。
145: 
146: 原投影IDFT系数PR下，四立方体LSIE/CIE/Dual/Global/A/A+的relative-L2材料误差约.4632/.4620/.4638/.4630/.4624/.4608；A+ CPU40.84秒，其余约17.47–22.09秒。误差差异很小，费用更高，未形成竞争力。8³反演/10³生成后误差约.249，各方法仍相近；空心立方体仍失败。
147: 
148: Mie球参考中，极化率保持的部分单元修正把5³/6³/8³误差降到6.38%/4.88%/3.88%。这是标准有效介质结构的应用，不能把这些误差贴到不同尺寸的全部成像对象。3λ论文原电尺寸也进行了真实16³生成/12³反演尝试，因内波长解析不足及费用保留截断结果；未完成60³/30³成熟复现。
149: 
150: 重复callback导致的多余LU经过缓存修正，同结果的8³ LSIE CPU从227.39秒降到101.28秒，LU653降到350。这是真实工程改进，所有最终臂共享，失败和旧费用保留。
151: 
152: ## 8. 成本、公平性与论文资格
153: 
154: 所有新确认对象的全部方法在同一机器运行；Mac/i9绝对CPU分层，同对象配对。材料、源、网格不同的LU不能直接当等成本。辅助FFT、材料内循环、feature/SVD、推理、训练、拒绝候选、两支路和失败均计费。并行线程wall之和不是研究历时；缺失计时为null或已知下界，不能填零。
155: 
156: 完整预测和候选选择先封存，再开测试GT。新集算法冻结后不按评分回调网络、β或预算。126项固定偏置扩展是在公开验证上提出、测试GT仍未读取时追加，保留时间及原750项，不伪装最初一次性预注册。部分公开测试Phi已在执行期间可见，材料评分始终隔离。
157: 
158: 目前最合理的论文定位是**学习初值的辅助状态传递、候选排序与模型适配边界的机制研究**，结合低成本NN直接物理优化的冻结基准证据和A+条件更新的费用分析；不是已经成立的通用更强Maxwell反演器。CS-Net、learned warm start、NeuralBIM、CIE-SOM-NET等已有直接先例，泛称“NN+CIE”不具有足够新颖性。
159: 
160: 实验规模、原数据和反例已足够形成可审查的完整初稿，但尚不支持以“广泛击败成熟CIE的A+算法”作为IEEE TAP/Inverse Problems论文主张。还缺一个原创且可预测的状态传递/选择机制，在新确认中既保留NNdirect收益又提供其没有的材料增量，并在合格实测/3D合同下重复。完整初稿不等于已满足期刊新颖性和接受要求，见[新颖性审计](PRIOR_ART_AND_NOVELTY.md)和[英文初稿](PAPER_DRAFT.md)。
161: 
162: ## 9. 下一步应继续哪条公式
163: 
164: 优先研究`(I−XG)δj`与材料切向读出的联系，以及`C(argmin Phi)`和`argmin Phi(C)`在实际预算中的差别。现有反例提示：仅降低current投影误差不够；需要保留有价值的物理分支，并找出无需GT、无需昂贵双路全跑的选择依据。此处尚无成熟预测器或全局保证，不以更高阶Taylor或更深网络替代这个问题。
165: 
166: 实测下一步应围绕合法air背景、照明/接收校准与离散误差，三维则围绕内波长/界面解析和模型误差下的材料选择。没有证据支持继续盲扫β、只加corrector或扩大网络深度。
167: 
168: 复现、全部数组、失败、费用和独立对象图分别见[REPRODUCE](REPRODUCE.md)、[RESULTS](RESULTS/)、[FAILURES](FAILURES.md)、[COST_LEDGER](COST_LEDGER.md)、[FIGURES](FIGURES/)。
169: 
170: 已发布并匿名核验：[研究网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r8/)、[原始Chat交接](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r8/CHAT_HANDOFF.md)、[52个完整数组归档包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r8-20261106)。归档共5,308,028,340字节；公开访问回执见[audit/PUBLICATION_RECEIPT.json](audit/PUBLICATION_RECEIPT.json)。

## research/research_max_struct_r8/CHAT_HANDOFF.md

1: # MAX-STRUCT-R8：给下一轮Chat的直接交接
2: 
3: 公开入口：[研究网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r8/) · [完整52包Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r8-20261106)。网页、图、raw报告、Release及匿名小包下载已经核验；全部52个资产的名称和字节数与本地归档一致。
4: 
5: **先读结论：NN已实际训练并接入sourcewise/global A+，但本轮最强独立二维结果是NN＋直接原Maxwell优化。A+没有在原论文、实测或三维上证明稳定普遍优势。** 完整证据、负结果和英文初稿均已保存；42材料的876预测先封存再评分，不是876独立样本。
6: 
7: ## 核心数字
8: 
9: |方法；平均幅度NRMSE|家族内24|家族外12|36→24跨网格6|
10: |---|---:|---:|---:|
11: |均匀＋直接物理|.24966|.40474|.19593|
12: |NN＋直接物理|**.09647**|**.16245**|**.08431**|
13: |均匀＋Dual-CIE|.36843|.63825|.32821|
14: |NN＋Dual-CIE|.32822|.56779|.32835|
15: |NN＋global A+|.29609|.59393|.32933|
16: |NN＋sourcewise A+|.32182|.59164|.32811|
17: |NN固定电流偏置＋global A+|.31395|.51831|.22800|
18: 
19: 主NN 24,737参数、1536训练/192验证对象、120epoch。GT只用于合法训练和封存后评分；没有按测试GT选根。辅助上限3000＋3000 corrector，直接物理6000，另提供同CPU/LU/RHS已支付端点。NN直接比均匀直接最终CPU更高1.06–2.74倍，是质量增益；相对本次辅助管线平均质量/费用较好，仍无普遍等成本优势。完整配对分层和置信区间见[NN报告](NEURAL_WARM_STARTS.md)。
20: 
21: ## 这轮真正发现了什么
22: 
23: 1. **有用初值会在接口中丢失。** `jphys(mNN)`进入`Jd+(I−Pd)F a`时，硬DPCS误差与Fourier截断分别可测；Maxwell缺陷由`(I−XG)δj`控制，current L2不能代替。固定`Jd+(I−Pd)jphys(mNN)`偏置部分修复，但所有开发例初始缺陷下降仍不保证材料改善。原理是既有投影结构的具体接线解释，不是新Maxwell恒等式。
24: 2. **更低起点Phi可导向更差材料。** 预定0027的NN直接1814步达到NRMSE.12082；先经global A+虽起点Phi更低，3000校正步后NRMSE1.01370。不是单凭预算上限可以解释，也未证明永远无法逃逸的basin。
25: 3. **最终Phi选更多候选仍会伤害泛化。** NN直接＋global A+双支路按最终Phi选，均值.07754/.34498/.22882，家族内改善、家族外和跨网格恶化，CPU为单直接的2.38–3.17倍。Global-Dual有类似现象，不能归为A+独有增益。
26: 4. **Fresnel先修材料合同。** 原χ≥.02把真实空气排除；修为[0,20]后从共同初值实际重跑Dual/A+，8GHz留出.137237/.137128，之前约.266。仍有尖峰且两法近同，不能称材料GT恢复或A+胜出。
27: 
28: ## 原论文、真实测量、三维都实际做了
29: 
30: - **Austria**：原150²生成/64²反演，2m ROI、20平面波40接收、400MHz、10%复噪声；εr=4五臂误差约.427–.501，没有稳定A+优势。独立NIE依原式实现，但不称作者认证的成熟复现。32²探索、长校正和无Austria标签的采集域NN修复全部保留。
31: - **Fresnel25端点**：真实2001/2005文件、空场多极照明、每源mask/adjoint、独立圆柱参考、FoamTwinDiel、独立NIE、NN回归和6+8GHz。新频率收益是新信息；无逐像素可信GT，不报伪NRMSE。
32: - **真矢量3D52端点**：TensorDDA/dyadic、三分量场、外部电偶极源、Fourier PR、四立方体/壳体/有损及跨网格；A+没有成本优势。Mie球部分体素误差细化到3.88%，只是指定前向控制。原3λ场景实际尝试16³/12³并明确截断，不能称60³/30³论文规模复现。3D NN及3D实测未运行。
33: - **历史四对象**：NN＋global/sourcewise A+保留7b8d与b20d，42cb和开缝壳仍失败；都是公开开发证据。
34: 
35: ## 论文判断和下一轮
36: 
37: 现在可写的是**学习初值的传递与候选选择失效机制，以及跨物理合同的可复核实验研究**，不是“普遍优于CIE的新A+”。CS-Net→SOM/TSOM/CSI、learned warm starts和CIE-SOM-NET都有近先例，不能把小CNN初始化当新发明。投影/仿射偏置/Möbius/MM也不能单独支撑高水平原创性。
38: 
39: 下一步最值得研究的公式是有噪声硬DPCS投影如何扰动本构状态，以及数据可见方向和训练先验方向如何在物理迭代中被保留。新简单机制必须以NN直接原物理作为强基线，不能只超过原均匀CIE。应寻找可预测反例/适用条件，避免继续堆候选却仅用同一Phi选根。实测应继续量化天线/离散/正则误差；原尺度3D还缺独立精度与充分预算。
40: 
41: ## 阅读与复现入口
42: 
43: - [完整中文结论](RESEARCH_REPORT_ZH.md)、[NN完整结果](NEURAL_WARM_STARTS.md)、[数学条件与推导](THEORY_DERIVATIONS.md)
44: - [原论文基准](PAPER_BENCHMARKS.md)、[Fresnel](FRESNEL_REPORT.md)、[3D](VECTOR3D_REPORT.md)
45: - [英文论文初稿](PAPER_DRAFT.md)、[新颖性](PRIOR_ART_AND_NOVELTY.md)、[失败与修正](FAILURES.md)
46: - [全部42材料图](FIGURES/warm_material_all_42.pdf)、[NN初始至最终图](FIGURES/neural_initial_to_final_all42.pdf)
47: - [成本](COST_LEDGER.md)、[复现命令](REPRODUCE.md)、[完整数组发布索引](delegated/cost_package/ARCHIVE_INDEX.md)
48: 
49: 科学文件与紧凑结果在仓库，完整数组/轨迹/训练数据及评分后合成真值在同仓库Release分包；第三方实测原件和全文不再分发。所有R1–R7冻结证据保持不变。未知失败CPU、启动、能源/货币成本仍明确未知，累计CPU小时不是并行历时。

## research/research_max_struct_r9/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT-R9：A/A+的增量、噪声硬锁与非凸搜索
2: 
3: ## 最终判断
4: 
5: **A/A+有可信的局部高对比正例，但目前不能说已经做到“普通问题普遍打平、困难问题稳定击败成熟CIE”。** A改变的是材料读出的空间本构残差度量；A+主要降低这个读出的计算成本。Möbius条件更新没有消除全局材料反演的非凸性，“扩大GT吸引域”也尚非已证明结论。
6: 
7: 第一阶段完成113次反演运行：47次公开开发/机制控制、66次新六材料比较；另有2次公开端点曲率诊断。结果给出了一个可复现的基线修复：按公开噪声选择DPCS秩，让新六例中的单β CIE、Dual、Global、A+都从2/6恢复到4/6。这一收益首先属于噪声与电流可行集的修复，**不是A+独有的恢复增量**。A+修正后平均NRMSE略低，但没有多恢复对象，平均CPU比修正后的单β CIE高约19%。
8: 
9: 另一条已知WRI/CSI型完整current松弛路线恢复了此前共同失败的40ac开缝壳，证明不同辅助路径确能获得原方法未取得的材料解。进一步简化为“低罚权初始化＋原物理校正”仍成功且更便宜。但这两条路线在新六例中都仅2/6，两分支组合也未超过4/6。**不能把一个开发正例包装成新的普遍凸化算法或论文级成像突破。**
10: 
11: 扩展已经实际完成：原NIE Austria的16项Windows比较、小型全矢量三维29项、原3λ四立方体九臂，以及三个原尺寸起点的短原物理校正。普通CIE在Austria ε3.5明显优于直接物理优化，在原尺寸三维明显优于LSIE，说明基线确有难以被普通局部物理优化替代的能力。**A+相对原A节省部分完整反演成本，但没有稳定超越匹配current的普通CIE；同样八步物理校正把三维Global和A+的材料误差都由约.423降至.249。**
12: 
13: 本轮最有解释力的是两个相反现象：三维小系数步会提前停在仍有明显原物理改进空间的端点；Austria继续降低原物理目标却损害材料质量。前者涉及本构缺陷经过完整resolvent向端口传播，后者受噪声、离散与prior影响，不能用“统一多迭代”同时解决。全部科学队列结束；没有将R9包装成新的实测或NN训练证据。最终代码、紧凑证据与大数组分别见仓库和[完整Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r9-20261106)。
14: 
15: ## 1. 继承与对照资格
16: 
17: R4–R8、Formula Atlas、FBG的实际代码和正负证据保留不动，本轮在隔离工作树、独立R9目录执行。复用R8 native FFT-TSOM、R7 A+、FBG-R3原Maxwell forward/adjoint及原物理corrector；没有另造弱CIE。审计读取630份R8完整回执，并保留R7的高对比正证据。
18: 
19: R8均匀初始化高contrast26例的材料成功数：单β CIE 6、Dual 3、Global 5、A 5、global A+ 4、sourcewise A+ 5。A+−Global的平均NRMSE约−.01505，但中位仅−.0000713，且没有A+独有成功；不能用均值遮蔽失效。R7高contrast18例A/A+各11/18、Dual/Global各9/18，属于不同采集/噪声/对象合同，正结果同样保留而不混分母。
20: 
21: 比较A读出时，global A/A+对应普通Global-Dual；sourcewise A+对应sourcewise Dual。只用global A+击败Dual不能排除current更新方式的作用。A+相对A的内层费用下降已在R8观察到，但不是相对所有CIE方法的免费升级。完整审计：[AUDIT_ZH.md](delegated/baseline_audit/AUDIT_ZH.md)。
22: 
23: ## 2. 物理、数据和评分合同
24: 
25: 第一阶段为**24×24、576个独立实contrast未知量的二维scalar-TM有限VIE**，exp(−iωt)、出射Hankel1、固定μ。复用原G/E/S及非周期FFT作用；源和接收器均在ROI外。材料m=εr−1，box[.02,20]，均匀初值.1；全部方法同一smooth+TV先验，最终以完整原Maxwell观测的Phi选已付费候选。
26: 
27: 新六对象为独立随机几何的ellipse、disks、annulus、comb、star、perforated，幅度2.5/12/18，复RMS噪声.1%/1%/3%，12源48接收、全/半发射孔径。一次生成、不重抽、不按结果过滤；六个对象才是分母，66次求解不是66独立样本。属于小型同离散确认，**不是跨网格、实测或3D证据**。已有R8的相关局限继续有效。
28: 
29: 主要算法/参数在新数据生成前冻结；补充single-CIE同噪声规则、简化弱阶段臂在生成后、所有确认评分前加入，依据是旧DEV证据，未读取确认GT调节。该时间区别保留在配置，不把补充臂称作完全相同的预生成注册。所有66个预测和两组portfolio的原Phi选择先封存，再首次读取六个GT。
30: 
31: NRMSE=逐像素RMSE/声明contrast幅度；支持分割使用prediction-only Otsu。沿用旧材料成功定义NRMSE≤.15且IoU≥.7及box合法，不为本轮改门。CN99噪声拟合、Phi、求解器停止、材料指标分开。辅助3000步＋共同corrector最多3000步；direct最多6000步。松弛/延拓预算单列，不伪称相同迭代即相同计算。
32: 
33: ## 3. 找到了什么问题，实际改了什么
34: 
35: ### 3.1 含噪DPCS硬锁
36: 
37: 旧阈值rtol=.001在四个R8采集算子中均保留17个receiver奇异方向，最弱保留方向使截断逆的条件数约394。对公开独立圆复噪声，硬读出电流的预期噪声量是
38: 
39: `E||ΔJD||² = Ns σ² Σ(i≤r) 1/si²`。
40: 
41: R8的3%噪声组，这一预期RMS量相对已观测JD范数的中位约91%；不是实际GT电流误差，但足以说明固定强模态阈值不能忽视噪声。R3已有无噪声/释放控制失败，因此它不是全部错误的唯一原因。
42: 
43: **实现修正：**公开σ按E|n|²=σ²解释，取完整数据99%噪声球，选最小连续奇异前缀r使丢弃数据能量落入该球。不同r共享同一个完整噪声事件，不从GT选秩。只改变硬锁/可优化current子空间；后续数据项仍使用完整Y，原物理验收也保留所有通道。额外SVD收费。减少r会释放一些原来锁死的current方向，不能称独立current自由度完全不变；材料仍576维，所有普通CIE也得到相同规则。
44: 
45: 这属于已有截断逆/差异原则的具体修复，不是新的Maxwell理论。代码：[baseline_adapter.py](src/baseline_adapter.py)。
46: 
47: ### 3.2 真正不同的辅助路径
48: 
49: 测试完整current松弛
50: 
51: `min(m,J) ||SJ−Y||²/(2yd) + λ||J−m(E+GJ)||²/(2ed) + R(m)`。
52: 
53: current在固定m时用完整Cholesky精确消元，使用包含原prior的包络梯度；区别于R2已有固定CIE权重profile。λ按初值的公开接收resolvent尺度归一，依次.01/1/100，各100次，再接原物理corrector。对照包括固定λ共300次、相同λ三次重启、数据目标同伦、普通direct。
54: 
55: 此处是已有WRI/CSI/variable-projection思想的机制控制。固定current子问题凸，不代表材料约化目标凸；Gram逆度量仍依赖m。参考[罚法先例](https://doi.org/10.1088/0266-5611/32/1/015007)及[约束松弛分析](https://doi.org/10.1088/1361-6420/adab89)。
56: 
57: 低罚权首阶段在40ac已进入好的物理区域，因此另试“首阶段直接交原物理校正”，复用已付费状态；单独方法费用仍补上全部初始化，实验总账不重复收费。没有增加NN、先验或观测。代码：[landscape_solver.py](src/landscape_solver.py)。
58: 
59: ## 4. 开发对象的决定性结果
60: 
61: |公开DEV|普通Dual / Global|A+|噪声修正|完整current路径|
62: |---|---|---|---|---|
63: |b20d弯条|NRMSE .390 / .488，失败|global .0406、sourcewise .0404，恢复|A+仍恢复，普通Dual/Global仍失败|三种最初延拓/固定松弛均失败|
64: |40ac开缝壳|.634 / .619，失败|.626，失败|秩仍17，不能修复|松弛延拓.0581、IoU .9618，恢复|
65: |0000正常控制|均恢复|恢复|各臂均恢复|各臂均恢复|
66: |0027高噪声困难例|固定秩约1.01，失败|固定秩约1.02，失败|四臂约.118–.120，全部恢复|direct及三种路径仍失败|
67: 
68: 40ac固定罚权重启仍NRMSE .4579；仅重启不能复制成功。弱阶段交接NRMSE .05835、IoU .96923，含已付初始化前缀的standalone归属费用19.6337 CPU秒，其中本次新增执行16.3525秒、已付前缀3.2812秒；对比完整延拓28.10秒。它不是另一次冷启动的独立计时，实验总账也没有重复收费。后两个dense profile阶段不是此开发正例的必要条件。
69: 
70: 原物理Hessian诊断也没有被当作魔法：40ac A+终点投影梯度9.43e−9，自由55坐标Hessian最小特征值+4.27e−9；b20d Global终点有一个−1.62e−8方向，但有限扰动＋同预算续跑与不扰动续跑都约NRMSE .492。活跃约束打开的临界锥未完全验证，不能宣布严格局部极小或全局证书。
71: 
72: ![四个公开开发对象及失败](FIGURES/development_materials.png)
73: 
74: 完整逐例表：[DEVELOPMENT_COMPARISON.md](delegated/scoring/DEVELOPMENT_COMPARISON.md)。
75: 
76: ## 5. 新六材料：修复有效，A+独有优势没有出现
77: 
78: |方法|材料成功|平均NRMSE|平均完整CPU秒|
79: |---|---:|---:|---:|
80: |原单β CIE|2/6|.55754|48.43|
81: |原Dual-CIE|2/6|.54021|46.06|
82: |原Global-Dual|2/6|.54211|47.96|
83: |原A+|2/6|.54352|52.67|
84: |噪声适配单β CIE|4/6|.18666|38.84|
85: |噪声适配Dual|4/6|.18882|40.81|
86: |噪声适配Global|4/6|.18766|42.98|
87: |噪声适配A+|4/6|.18552|46.37|
88: |直接物理优化|2/6|.40359|31.19|
89: |完整current松弛延拓|2/6|.34173|15.45|
90: |弱阶段＋原物理校正，含初始化|2/6|.34262|14.42|
91: |噪声A+＋松弛，两次完整付费|4/6|.18488|61.82|
92: |噪声Global＋松弛，两次完整付费|4/6|.18735|58.43|
93: 
94: 噪声修正额外恢复的是comb0003和perforated0005，各普通CIE同样取得；annulus0002与star0004仍共同失败。新普通控制0000、0001都保住，但六例不足以证明总体非劣性或成功率显著提升。A+均值略好来自少数连续误差差别，没有新增成功；不能称其恢复增量已稳定复现。
95: 
96: ![六例逐对象材料、拟合与成本](FIGURES/confirmation/paired_cost_matrix.png)
97: 
98: [全部六例GT＋十一臂材料图](FIGURES/confirmation/all_six_materials.png)包含所有失败。简化弱阶段平均更便宜，但0000与0003反而稍贵，不是逐例加速。两portfolio在DEV分别4/4与3/4，在确认均4/6；额外花费没有多恢复对象。保留这个泛化负例，不再根据确认GT修改组合。
99: 
100: 已付费轨迹的**事后**30/60 CPU秒切片中，噪声修正各CIE仍4/6，原各CIE2/6；例如30秒时噪声A+平均NRMSE .16482，继续到最终反而.18552。这提示Phi下降仍可伴随材料恶化，不能据此按GT选择30秒停止。切片只取实际已付费原Phi incumbent，不插值、不GT挑迭代；不是新的预注册性能点。见[PAID_CPU_SLICES.csv](RESULTS/summary/PAID_CPU_SLICES.csv)。稠密profile的current Gram/Cholesky全部计入CPU，不拿较少physical LU冒充等总算量。
101: 
102: ## 6. 对“凸化”和剩余失败的回答
103: 
104: **PROVED_UNDER_ASSUMPTIONS：**严格递增变换保留局部极小；局部微分同胚材料换元在驻点处通过合同变换保持Hessian惯性。因此仅用Möbius改写或对原Φ取单调函数，不能删掉原错误极小。交叉熵采样搜索、KL数据loss、熵正则运输、熵prior是不同机制；对概率或运输计划凸，不等于对m↦F(m)凸。
105: 
106: 完整current松弛确实改变中间目标，40ac正例证明这种路径有价值；其b20d和新确认负例不支持普适恢复保证；有限预算失败本身不作为目标非凸性的证明。已有特殊有限端口Galerkin凸化定理需要整个状态空间由匹配source/receiver Riesz基张成等条件，当前576内部自由度和48×12外端口不满足，不能强行套用。
107: 
108: 没有把共同失败直接归为不可辨识。匹配有限模型下可用保存的GT全波响应作**离线可行点见证**；见[audit/CONFIRMATION_FAILURE_DIAGNOSTIC.json](audit/CONFIRMATION_FAILURE_DIAGNOSTIC.json)。更低的已知可行Phi只能说明当前求解未达到该值，不能证明唯一性、全局最优或整个连续PDE盆地结构。后续重点应是如何绕开共同错误材料区域，而不是继续给A+添加更精细的一维读出步骤。
109: 
110: ## 7. 第一阶段定位与后续问题
111: 
112: 1. **A解决的是局部材料读出的空间残差耦合；A+解决的是该读出的部分费用。**它们尚未解决硬DPCS噪声、所有Fourier表示限制、sourcewise PR问题和全局非凸搜索。
113: 2. **本轮可采用的具体修正是noise-aware DPCS，且必须同样更新强CIE基线。**在新六例里普通修正CIE已解释主要收益；A+不应被宣传为普遍优于它。
114: 3. **保留低罚权完整current初始化作为便宜的互补候选。**40ac实际被恢复，简化有成本收益；但确认没有额外成功，暂不能作为普遍更强的默认算法。
115: 4. **目前最有价值的是机制和失效研究，而非新算法广泛优胜论文。**有明确噪声硬锁反例、合法路径互补、错误终点曲率、公开原始数组和成本；原创性仍需与既有CSI/WRI、差异原则及CIE近邻严格区分。该六例研究同离散、小样本；已完成的文献对象及真正三维按下述独立合同报告，也不替代未开展的新实测验证。
116: 
117: 下一轮应围绕已有公式中的两点设计少量可反驳研究：一是用噪声相容的软约束替代二分硬锁，同时保持已修好的普通CIE基线；二是研究低罚权阶段为何在40ac产生有效可观测材料方向、却把新annulus/star带到共同错误终点。不能因为知道GT形状就给它们定制prior，也不能把更低Phi包装为恢复。NN、全局凸化及更广泛统计优势在第一阶段均未被证明；跨维度结果由新增全矢量实验单独评价。
118: 
119: 完整推导见[THEORY_DERIVATIONS.md](THEORY_DERIVATIONS.md)，费用与缺口见[COST_LEDGER.md](COST_LEDGER.md)，运行与公开数组见[REPRODUCE.md](REPRODUCE.md)。
120: 
121: ## 8. 已完成的三维扩展与物理修正诊断
122: 
123: 新增29次小型全矢量三维反演已经完成：8³=512个实各向同性材料、1536个complex场分量，采用真正dyadic Green和60外部电偶极照明。与二维scalar-TM分表；不是将二维堆叠成三维，也不是原paper3λ电尺寸。
124: 
125: |三维公开对象|Global-Dual|A|A+|Direct|
126: |---|---:|---:|---:|---:|
127: |缩小四立方体|.25022|.25398|.25039|.25486|
128: |薄空心体＋杆|.29281|.29562|.29355|.30006|
129: |ε9双立方体|.47752|.47854|.47728|.47829|
130: |可分辨厚壳控制|.18066|.18063|.18063|.18063|
131: 
132: 表为χ相对L2误差，**与第一阶段幅度归一NRMSE分母不同**。在共同300步原物理校正后，没有看到清晰的A/A+终点优势。前三例同Mac平台的Direct CPU约79–86秒、A+约163–165秒；厚壳同Windows平台Direct253.72秒、A+511.38秒。辅助240＋校正300与direct300不是等迭代/等算量；完整已付费曲线、所有方法及停止信息在[VECTOR3D_EXTENSION.md](VECTOR3D_EXTENSION.md)。
133: 
134: 物理核查没有被低线性残差替代：ε9球的16³ Mie预定散射平面误差仍43.1%。一次标准源单元积分修正为43.4%，未解决主要误差。根据界面条件继续做了两个已知球离线正演，scalar arithmetic和局部法向tensor降到27.8%与25.4%；但其截面误差分别4.7%与25.6%，并非一项表示所有物理量都更好。该反例说明单一积分指标可以掩盖复场分布误差。
135: 
136: 这些是已知几何的前向资格诊断，没有把GT界面法向交给未知材料反演，没有把静态平均、IGT或tensor公式记为新发明。原3λ四立方体的全60源CM+RR 30/60观测差约53.5%；固定source59单元积分后由50.8%降至36.0%，仍很大。因此目前不能称已独立复现作者连续高对比结果。详见[VECTOR_DISCRETIZATION_DIAGNOSTIC.md](VECTOR_DISCRETIZATION_DIAGNOSTIC.md)。
137: 
138: 原Austria16项已在远端i9四个P核全部完成。原3λ的cross-grid五臂＋matched-grid四臂已在4060全部完成，见下两节；不能把换了数据生成离散后的改善记为同数据算法收益。
139: 
140: ## 9. 原3λ论文对象：CIE优于LSIE，但A+没有进一步优势
141: 
142: 新九臂为真正3D矢量，27000个实材料未知、60外部偶极源、180个接收分量，MF6→10，air初值、box[0,4]、R=0。全部18阶段达到原Eq16系数停止，无CAP或数值失败；主终点固定stage2，没有按GT选择迭代。它是CM+RR离散下的原几何/电尺寸算法适配，不能称作者代码或已取得连续物理精度的复现。
143: 
144: |方法|cross-grid χ相对L2|matched-grid χ相对L2|matched完整物理Φ|
145: |---|---:|---:|---:|
146: |source CIE|.41460|.42727|.13478|
147: |source Dual|.41088|未列入冻结四臂|未运行|
148: |Global-Dual|.41072|.42275|.14120|
149: |A+ global|.41182|.42347|.14238|
150: |LSIE|.91006|.92953|.49414|
151: 
152: 这里Φ是逐源相对数据残差平方的平均，没有1/2，与第一阶段数值不能横比。A+与Global-Dual在两组数据上材料、物理拟合均接近且略差。matched的实际共享队列wall为Global298.35秒、A+380.33秒；补入setup与复用归属后仍是299.22/381.21秒。A+的host CPU较低并不代表运行更快，GPU算术必须看wall与设备计数。
153: 
154: CIE相较LSIE显著减少χ误差，能分辨四块；但并非精确恢复。论文逐像素相对ε误差会受到背景占比和高ε归一化影响：air为.18966，cross LSIE为.17994、CIE却为.21627；仅看这个指标会掩盖CIE在高contrast块上的实际改善。故保留完整色轴材料图、χ相对L2和原物理拟合，不挑有利指标。
155: 
156: ![原3λ同离散控制的完整材料对照](FIGURES/original_CM_matched_grid_materials.png)
157: 
158: matched数据下已知同一GT的可行Φ=.010203，CIE系终点仍为.135–.142；其辅助data RMS却只有.017–.027。因此粗细网格差异不能解释全部缺口，Eq16小系数步也不是原物理收敛证书。复用`r_phys=d−S L(m)^−1 e`及已付范数，当前残差方向的传播增益估计区间为CIE[12.11,13.27]、Global[6.97,8.06]、A+[7.00,8.09]。不同β的状态范数不同，不能按这些数字给算法排序；Global/A+同规范的结果显示局部材料读出没有控制这一端口缺陷。
159: 
160: 该区间上端不是operator-norm上界，数值也没有严格线性误差包络；它是已完成端点的机制诊断，不是新定理或证书。推导、适用条件及完整九行数据见[THEORY §9](THEORY_DERIVATIONS.md)与[诊断原始表](RESULTS/original_residual_transport/DIAGNOSTIC.csv)。
161: 
162: 为区分尚有原物理下降空间与更难的错误区域，另冻结matched的air、Global-Dual、A+三个初值，各做相同8 accepted step上限的原物理L-BFGS-B；三臂已全部完成、封存后统一评分。这是GT已经公开后的开发补充，旧九臂预测、终点和费用保持不变。完整原尺寸合同见[原九臂报告](delegated/original3d_feasibility/ORIGINAL_CM_NINE_REPORT_ZH.md)。
163: 
164: |起点|χ相对L2，开始→八步后|原物理Φ，开始→八步后|本次新增wall秒|
165: |---|---:|---:|---:|
166: |air|1→.97501|1→.64852|120.33|
167: |Global-Dual|.42275→.24899|.14120→.02775|728.72|
168: |A+ global|.42347→.24950|.14238→.02780|725.89|
169: 
170: 两个辅助端点都还存在有效的原物理改进方向，A+没有独有收益。三臂均八步预算停止，不能称最终驻点。air的八步实际便宜得多，且Global/A+此前的辅助运行并非免费，故该表不构成等总成本的warm-start胜出证明。新增整进程1557.95 CPU秒／1583.52 wall秒，1740正向＋1680伴随RHS；完整独立增量账包含绘图、运输和工程失败，已知CPU下界1560.18秒。没有重复计入旧九臂费用。
171: 
172: ![三起点八步后的材料；完整共同色轴](FIGURES/original_vector3d/physical_refinement/materials.png)
173: 
174: 详细梯度、请求/缓存区别、全部已付状态和成本见[短物理校正报告](delegated/3d_extension/physical_refinement/REPORT_ZH.md)。这只是匹配有限CM+RR模型，不能转述为跨网格或连续Maxwell的准确恢复。
175: 
176: ## 10. 原NIE Austria：普通CIE的优势确实存在，A+没有稳定增量
177: 
178: 两种公开材料条件的16项Windows运行已完成，并在全部预测封存后统一评分。64²=4096个实材料未知，150²生成数据，400MHz、20源40接收、10%原噪声；保留原MF/β日程、完整native预算及同一连续300步原物理校正。A/A+的current与Global-Dual完全相同；独立NIE/LSIE参考为R=0，与其余六个匹配prior的臂分开解释，不冒称作者代码。
179: 
180: |方法|εr=3.5幅度NRMSE|εr=4幅度NRMSE|两条件完整CPU秒|
181: |---|---:|---:|---:|
182: |source CIE|.18513|.47475|1552 / 2603|
183: |source Dual|.19040|.48472|1616 / 2531|
184: |Global-Dual|.18962|.57517|1577 / 2799|
185: |A global|.18922|.58654|2108 / 4617|
186: |A+ global|.18450|.58300|1592 / 2978|
187: |直接物理L-BFGS-B|.60683|.58057|382 / 477|
188: |独立NIE，R=0|.18644|.47560|1514 / 2418|
189: |独立LSIE，R=0|.61526|.59125|586 / 687|
190: 
191: 因此“直接物理优化总比论文CIE强”不成立：ε=3.5时CIE系保留Austria结构，direct/LSIE大幅失败；但这也不是A+的独有收益。A+相对原A的完整任务CPU减少约24.5%/35.5%，支持该实现的整体工程收益；这些比例不是纯材料内循环计时，也包含不同轨迹的成本。相对ordinary Global，它没有稳定材料增量；ε=4时还比source CIE明显差。global与sourcewise的差异不能全部归因于A材料metric。
192: 
193: ε=3.5时A+比single CIE的NRMSE只小约.00063，Otsu IoU反而.87395<.87678，不能据此宣布优胜。ε=4多种CIE/A结果出现εr=21的假峰，材料幅度失败没有被约.62–.69的Otsu重合率遮蔽。完整幅度图和全部方法见[经典报告](CLASSIC_BENCHMARKS.md)。
194: 
195: 另一个决定性反例是**更长物理优化会损害材料**：ε=3.5多个CIE/A方法在第80步的NRMSE约.136，继续到300步后升到.185–.190；以single CIE为例，保存的原物理目标Φ只从约.0043534降到.0043184。该Φ含共同prior，不能将sqrt(2Φ)直接称作纯数据残差。ε=4也多有退化及box假峰。这不证明GT就是优化目标的最优点，也不唯一归因于噪声；该任务还存在约9.96%的粗细模型差异。它说明错误材料不总是“还没优化够”，也不能用更低Φ作恢复证据。
196: 
197: 所有16个corrector最终都是300步CAP，不能写成完全收敛；部分native stage也达到500步上限，表内全部保留。80点是预先保存的同一连续运行检查点，**没有事后把80改为线上最优停止**。300/600/900 CPU前缀仅使用已付保存材料、先按公开Φ选再评分；快照间隙可能很大，不能冒充精确等CPU截断实验或完整效能曲线。完整实际成本、35次GMRES停滞后LU回退和Mac迁移9项全部收费，主16项没有整项数值失败。
198: 
199: ## 11. 完整交付、成本与论文资格
200: 
201: 本轮科学队列全部结束。主证据池为113次二维机制反演、16项经典Windows比较、29次小型真矢量反演、原尺寸九臂及三个短物理校正；另保留Mac迁移前部分结果、所有中断与诊断。它们不是170个独立材料，不合并造成功率。所有方法原始预测、长轨迹、状态、失败与配置在同仓库Release无损保存，精简结果与图直接在Git浏览。
202: 
203: 远端真实使用i9 P核四lane和RTX4060 resident complex128 CUDA。经典主队列wall7952.01秒，原尺寸九臂whole-process wall3138.32秒，新增三臂wall1583.52秒；这些运行有重叠及不同范围，不能相加称项目历时。科学和分析已知CPU下界约Mac4.64小时、Windows10.38小时，归档另计；最终[总账](COST_LEDGER_EXTENSION.md)按平台、维度和费用范围去重，金额、能源和未测编排保持未知。
204: 
205: 独立审查核对了第一阶段封存与配对、原尺寸停止及current合同、经典误差分母与成本。最终修订明确：A+节省比例是整任务费用；80步比较引用含prior Φ；NIE/LSIE与匹配prior组分开；host CPU不是GPU速度。见[审查处理记录](audit/FINAL_ROOT_REVIEW.json)。
206: 
207: 现在适合形成一篇严谨的机制和失效论文初稿：噪声硬锁、辅助路径互补、resolvent传播的本构缺陷，以及相反的停止失败。尚不足以主张广泛优于成熟CIE的新算法。普通CIE解释了主要确认收益，A+有局部开发正例及相对原A的工程价值，三维/经典对象没有稳定新增恢复。下一步应将噪声和可信模型误差纳入可观测残差软约束，再验证是否真的改变共同失败区域；单纯更精确的Möbius内循环不够。

## research/research_max_struct_r9/CHAT_HANDOFF.md

1: # MAX-STRUCT-R9：最终研究交接
2: 
3: **本轮科学队列已全部结束。** R1–R8保持不变。本轮没有新增Fresnel实测或NN训练，不将继承证据记为新成果。
4: 
5: [研究网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r9/) · [中文报告](RESEARCH_REPORT_ZH.md) · [理论](THEORY_DERIVATIONS.md) · [英文初稿](PAPER_DRAFT.md) · [完整数组Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r9-20261106)
6: 
7: ## 最终判断
8: 
9: **A/A+保留真实的高对比开发正例，但尚未做到“普通问题普遍打平、困难问题稳定击败成熟CIE”。** A改变材料读出的空间残差度量；A+降低其部分费用。Möbius条件更新没有凸化整个Maxwell材料目标。比原A便宜与比普通CIE多恢复材料是不同主张。
10: 
11: 实际完成：第一阶段113次反演（47开发＋66确认）、经典Austria的16项Windows比较、小型全矢量三维29项、原3λ三维九臂，以及三起点各八步物理校正。Mac迁移前9个Austria终点及中断费用另存，不混入Windows性能。运行次数不是独立材料样本数。
12: 
13: ## 关键证据
14: 
15: |问题|实际结果|判断|
16: |---|---|---|
17: |历史b20d弯条|A+ NRMSE约.0406；Dual/Global约.390/.488|有用路径，属于公开开发对象|
18: |噪声硬DPCS|公开噪声球选秩后，六例single/Dual/Global/A+全部2/6→4/6|共享基线修复，不是A+独有|
19: |修复后|平均NRMSE .18666/.18882/.18766/.18552；CPU38.84/40.81/42.98/46.37秒|A+未多恢复对象、非最低成本|
20: |40ac开缝壳|完整current松弛.05807；弱阶段校正.05835，归属CPU19.63 vs28.10秒|CSI/WRI型路径互补；新六例仍2/6，组合仍4/6|
21: |Austria ε3.5|CIE .18513、A+ .18450、direct .60683、独立NIE .18644、LSIE .61526|CIE优势真实；A+基本打平|
22: |Austria ε4|CIE .47475、Global .57517、A+ .58300|A+未解决更强散射失败|
23: |A+工程收益|Austria两条件完整CPU比A少24.5%/35.5%|整任务含不同轨迹；非单独内循环测速|
24: |小三维29项|共同校正后A/A+/Global终点相近，direct更便宜|无稳定A+增量；真实矢量、小电尺寸|
25: |原3λ matched|CIE/Global/A+/LSIE χ相对L2为.42727/.42275/.42347/.92953|CIE明显优于LSIE；A+未进一步改善|
26: |共同八步校正|Global/A+ χ误差降至.24899/.24950；air为.97501|均有有效物理方向，无A+独有收益；不是等总成本比较|
27: 
28: 二维NRMSE与三维χ相对L2分母不同。经典NIE/LSIE为R=0参考，其余六臂同prior。全部候选按公开物理目标选择，GT仅封存后评分；公开论文几何不重新称盲测。
29: 
30: ## 最值得继续的公式和现象
31: 
32: 辅助电流J、e=L(m)J−D_mE、d=SJ−Y满足：
33: `F(m)−Y=d−S L(m)^{-1}e`。
34: 
35: 这是已知恒等式。新证据是原3λ matched辅助data RMS只有.017–.027，完整物理RMS却约.37；实际缺陷方向的传播增益约7–13（不同β范数不可直接排名，上端不是operator-norm上界）。Global/A+同范数结果近乎相同，小系数步没有控制resolvent传播的端口缺陷。
36: 
37: 共同八步校正确实改善两个端点，说明它们尚有有效原物理方向。但Austria ε3.5从预存80步继续到300步，材料NRMSE约.136恶化至.185–.190，single CIE的含prior Φ仅由.0043534降到.0043184。不能按GT改选80步，也不能将sqrt(2Φ)叫纯数据残差。
38: 
39: **下一轮优先：利用公开噪声和可信模型误差，区分需要释放的本构缺陷与不该继续拟合的端口残差。** 可研究可观测残差软约束和低罚权阶段的材料方向，不能默认更多A+步、更长corrector或更深网络能同时解决两种相反问题。严格递增loss变换保留局部极小，局部可逆换元保留驻点Hessian惯性。
40: 
41: ## 物理与统计边界
42: 
43: - ε9球Mie散射平面误差8/10/16网格为61.6%/55.5%/43.1%；固定16网格积分Green为43.4%。已知界面arithmetic/tensor为27.8%/25.4%，仍不够准确；GT界面未交给反演。
44: - 原3λ CM+RR全60源30/60网格差53.5%；source59积分后50.8%→36.0%，仍远大于10%噪声。matched是另一个有限模型数据，不是同数据算法提升或连续物理认证。
45: - 27000实未知、至多21600实观测。维数限制内部常秩点，不能单独证明box边界GT多解。论文相对ε指标受背景影响，故并列χ误差和完整色轴图。
46: - 经典corrector全部300 CAP；原尺寸九臂达到论文系数停止；短校正全部八步预算停止。均非全局收敛证书。
47: - 六个确认材料不支持广泛统计优胜；没有新实测、复tensor或NN收益结论。
48: 
49: ## 远端与费用
50: 
51: Windows二维四lane固定i9 P核逻辑0/2/4/6；三维4060使用resident complex128 CUDA/cuFFT，host支持P核12/14。队列全部结束，GPU已释放。
52: 
53: 同一普通CIE三维current＋material步2.975→1.021秒（2.915倍布局优化）；fine单源RHS CPU99.91→CUDA5.89秒（包含Krylov实现差异）；二维4096材料value＋gradient仅1.400→1.103秒，GPU冷setup34.64秒，所以四P核跑二维、GPU跑三维。不能外推统一全任务加速。
54: 
55: 科学/分析CPU已知下界约Mac4.64小时、Windows10.38小时；这是消耗，不是并行历时或GPU时间。最终归档另加入[总账](COST_LEDGER_EXTENSION.md)。LU/RHS/FFT按维度分账；失败、迁移、预计算和拒绝trial保留。未计量CPU、网络、LLM、能源和金额保持null。
56: 
57: ## 入口与论文资格
58: 
59: [经典二维](CLASSIC_BENCHMARKS.md) · [小三维](VECTOR3D_EXTENSION.md) · [原3λ九臂](delegated/original3d_feasibility/ORIGINAL_CM_NINE_REPORT_ZH.md) · [短物理校正](delegated/3d_extension/physical_refinement/REPORT_ZH.md) · [物理诊断](VECTOR_DISCRETIZATION_DIAGNOSTIC.md)
60: 
61: [合同](BENCHMARKS.md) · [先例](PRIOR_ART_AND_NOVELTY.md) · [复现](REPRODUCE.md) · [第一阶段审查](delegated/independent_phaseI_review/REVIEW_ZH.md) · [三维审查](delegated/independent_original3d_contract_review/REVIEW_ZH.md) · [经典审查](delegated/independent_classic_review/REVIEW_ZH.md)
62: 
63: 目前能写的是有原论文对象与真实矢量实验支撑的机制、基线修复和失效研究。广泛优于成熟CIE的算法论文仍缺：新独立对象上的稳定恢复增量、合格连续物理/实测支撑、等实际成本优势。本轮没有用期望发表的结论替代证据。

## research/research_max_struct_r10/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT R10：Green分母、能量归因与CIE补偿路径
2: 
3: 本轮在R9基础上实际完成更大范围的CIE/A+比较，并从公式提出、实现和检验低秩矩阵分母。R1–R9原始结果未改写。24材料×7方法主确认、3例结构修正、4个新材料的条件补偿、18条三维主反演和3条固定预算续跑均已完成封存；各队列独立评分，不能合并成一个确认成功率。
4: 
5: ## 最终判断的依据
6: 
7: **原A+仍有真实的局部材料恢复作用，但不是当前值得默认替换成熟CIE的方案。新Green矩阵分母提供了一条不同的有效路径，也出现明确反例。** 本轮的评价不要求所有对象获胜；真正限制主张的是：收益目前不能由公开指标提前稳定定位，部分高对比对象发生大幅退化，且平均完整任务CPU增加。
8: 
9: 最重要的新正例不在历史弯条上：一个新高对比环中，CIE/原A+的NRMSE约.423，新Green19降至.099，原始数据残差由13.7%降至.960%。同时，一个新高对比桥由CIE的.097恶化至Green19的.337；开缝壳中CIE约.111，而原A+和新分母都约.515。三个例子明确说明辅助材料规则改变了到达的区域，不能说明它普遍把原Maxwell目标凸化。
10: 
11: 多数普通对象的连续材料差异很小，完整物理corrector会把不同辅助路径带到近似终点。部分成功数差异只是NRMSE在.15附近跨线；不能把11/24对8/24描述成每一例都有显著退化。另一方面，上述环、桥、开缝壳的差异很大，确属需要解释的物理优化现象。
12: 
13: ## 1. 本轮究竟改了什么公式
14: 
15: 固定辅助current J，令Q=Einc+GJ、e=J−D_mQ。原A/A+使用
16: 
17: `a_i=1+(β−Gii)m_i`，`z_i=m_i/a_i`，`u=β e/a`，`R_A=(I+D_z Goff,block)u`。
18: 
19: Goff,block只包含互不重叠2×2块内非对角作用，不是完整G去掉对角线。A+的Möbius/MM主要降低这个既有材料读出的计算费用，并未解决完整非局部反馈。没有证据表明原物理自项Gii数值写错。
20: 
21: 本轮新实现为
22: 
23: `R_K=β[I+βD_m−D_mK]⁻¹e`，`K=PGP`。
24: 
25: P来自公开背景算子的强右子空间，不能读取目标GT；材料仍有全部576个独立实未知量。Green19只用G；Two-fold增加接收子空间。**新确认的实际秩是19和30，历史DEV才是19和28。** 早期文档把历史28沿用到新集是文字标签错误，实际basis、小LU、代码和运行一直使用新算子自动给出的30，阈值.1没有改变。
26: 
27: Woodbury把分母求解缩为19/30维小系统。材料导数完整保留分母随m变化的项；不是冻结分母后的近似梯度。current更新、FFT/APCS表示、初始化、prior、合法box和最终完整Maxwell验收均继承。小LU和额外材料内循环全部收费。
28: 
29: ### 已建立的数学结果
30: 
31: - **DERIVED：** 精确材料导数及Woodbury伴随，近似误差式 `R_K−R_G=−A_K⁻¹D_m(G−K)R_G`。静态小Green尾不保证小材料误差，材料乘法、遗漏耦合和resolvent共同决定作用。
32: - **PROVED_UNDER_ASSUMPTIONS：** 对本轮精确解析二维盘脉冲G，Im_H G是严格正的全角度辐射Gram。因此PGP分母在非负m、β上定性可逆。它不是冻结浮点矩阵的区间证书、统一条件数界、材料唯一性或全局收敛定理；不自动覆盖其他离散或3D。
33: - **DERIVED反例：** 即使正分母下的单变量m²/(1+m)²也在m>1/2时有负二阶导数。分式归一化不保证全局凸。对原损失施加单调变换也不能据此声称消掉错误吸引域；R9已有相关边界，本轮不重记原创。
34: - **EMPIRICAL：** 24个历史状态中，局部one-hop反馈谱半径<.168，而投影全域反馈后期18个状态>1，小分母minσ仍为.312–.688。局部耦合弱、Neumann级数不收敛、矩阵接近奇异是三件不同的事。
35: 
36: 详细推导、标准定理来源及独立代码审查见[THEORY_DERIVATIONS](THEORY_DERIVATIONS.md)、[prior-art](PRIOR_ART_AND_NOVELTY.md)和[审查记录](delegated/denominator_score/INDEPENDENT_REVIEW.md)。新颖性不建立在重新命名Woodbury、Möbius、Two-fold或耗散性上。
37: 
38: ## 2. 新24材料主确认：全部结果
39: 
40: 八家族各有low/high/mixed三个独立材料：环、开缝壳、多岛、曲条、梳状体、包覆体、穿孔星形、桥。low χ2.5，high χ12/18，mixed含两级；ROI host χ.1。全部576维未知，12源/48接收，1%期望复噪声，box[.02,20]。高contrast的12/18按既有家族规则分配；本设计不能单独识别contrast与家族的因果作用。
41: 
42: 同一有限离散生成/反演，不能称跨网格或实测。168条反演全部完成、无运行失败，先封存全部公开预测再统一读GT。实际每条3000辅助步；共同原物理corrector最多3000接受步/12000评估，无wall截断。140条由投影梯度准则停止、4条由相对函数变化停止、24条到迭代上限，不能将全部终态称为已收敛。
43: 
44: |方法|主材料成功/24|平均NRMSE|中位NRMSE|平均完整child CPU秒|角色|
45: |---|---:|---:|---:|---:|---|
46: |CIE-single|11|.16287|.13801|157.71|主要成熟基线|
47: |Dual-CIE|10|.16346|.13829|153.98|主要成熟基线|
48: |A+global|8|.18044|.14531|169.49|原one-hop材料改进|
49: |A+sourcewise|10|.18045|.14479|174.98|匹配Dual current的材料消融|
50: |Green19矩阵分母|8|.17669|.14638|187.66|本轮新候选|
51: |Two-fold30矩阵分母|8|.18099|.14544|192.53|接收子空间候选|
52: |Global-Dual|10|.16369|.13894|150.13|内部matched-current控制|
53: 
54: 主材料成功沿用已冻结NRMSE≤.15、共同χ>1的IoU≥.7及合法性，**不是研究停止门**。完整连续分数、旧Otsu口径、区域误差、停止原因和全对象配对见[原表](delegated/twod_score/scored/SCORES_FINAL.csv)。主support在GT读取前冻结，避免将混合材料低χ相从分母中删除。没有按GT选择候选或迭代。
55: 
56: |方法|low平均NRMSE / 成功8例|high平均NRMSE / 成功8例|mixed平均NRMSE / 成功8例|
57: |---|---|---|---|
58: |CIE|.16907 / 4|.21316 / 3|.10638 / 4|
59: |Dual-CIE|.16907 / 3|.21460 / 3|.10671 / 4|
60: |A+global|.16959 / 2|.26514 / 2|.10660 / 4|
61: |A+sourcewise|.16949 / 3|.26533 / 2|.10653 / 5|
62: |Green19|.16967 / 2|.25363 / 2|.10676 / 4|
63: |Two-fold30|.16984 / 2|.26646 / 2|.10669 / 4|
64: 
65: low和mixed的平均连续差很小，不能从成功数跳变推导大幅优劣。高对比有实质性的相反路径结果，值得逐例研究。对象配对bootstrap仅为这个设计的描述性区间，不是任意几何总体或多重比较校正后的显著性保证。
66: 
67: ### 三个决定性的新正反例
68: 
69: |对象|CIE|Dual|原A+global|Green19|Two-fold30|原物理现象|
70: |---|---:|---:|---:|---:|---:|---|
71: |0001 高χ12环|.42288|.43372|.42275|.09889|.10248|新分母把残差13.7%降至约.960%，进入不同可恢复区域|
72: |0004 高χ18开缝壳|.11109|.11170|.51514|.51514|.51514|CIE残差.968%，A类均58.98%；新分母没有修复原A的坏路径|
73: |0022 高χ18桥|.09741|.09625|.09608|.33700|.43720|CIE/A+残差约.974%，新分母仍16.85%/14.40%，出现反向损害|
74: 
75: 0004的A类在约482–580物理接受步后到投影梯度停止，不能把失败统一解释为预算用完。0022的Green19由函数变化停止，Two-fold30到3000上限；需要分别看停机类型。0001的矩阵分母到3000步时已达到材料/数据正例，若要继续追求其最终驻点，需要另计预算，不能把尚未收敛说成完全解决。
76: 
77: ![全部八个新高对比对象](FIGURES/twod_comparison/MATERIALS_high.png)
78: 
79: [全部low](FIGURES/twod_comparison/MATERIALS_low.png)、[全部mixed](FIGURES/twod_comparison/MATERIALS_mixed.png)、[成本与性能](FIGURES/twod_comparison/COST_PERFORMANCE.png)同时保留，没有只展示获胜形状。
80: 
81: ## 3. FFT-TSOM提供了什么，A+提供了什么
82: 
83: 本轮重用历史69次运行449状态，并对新24对象七臂的所有终点与同m末stage做独立能量审计。新诊断只为24份GT current补24次576阶LU/288源RHS；端点和联合状态全部复用已保存数组，没有为图反复解forward。
84: 
85: 这里“能量”指离散平方Euclidean/Frobenius范数，不是电磁储能。Green/S右空间的材料坐标投影不等于材料Jacobian可观测性；Two-fold顺序投影、数值交集和实际FFT/APCS affine空间分开记录。
86: 
87: 历史b20d同FFT空间下，Global→A+的材料误差能量19753.828→136.765，物理current误差1184.654→8.566。40ac则31785.863→32500.705、1062.402→1080.634。末mf12的辅助oracle不可表示能量分别仅.232和.205；不能用“给A+更多Fourier模式”解释b20d，也不能把小oracle尾当作40ac必然能恢复的保证。
88: 
89: 辅助current与材料导出的物理current必须分开。一个历史状态满足200.610=781.243+624.985−1205.618，其中最后一项是有符号cross。小辅助误差可能来自抵消，不能当作小物理current误差。物理corrector不受辅助affine空间约束，其结果低于辅助oracle floor不矛盾。
90: 
91: 新集复用20个(rank,mf)空间缓存，得到96个case-oracle，所有七臂同对象共享实际表示。全部168终点和168同m配对均有结果、没有缺失。完整能量百分比、绝对量和cross见[energy原表](delegated/twod_energy/scored/STATE_ENERGY.csv)及[全算子分带](delegated/twod_energy/scored/SINGLE_OPERATOR_ENERGY.csv)。不能因某带百分比升高就说该带绝对误差更大。
92: 
93: ## 4. 失败后做了什么修正
94: 
95: ### 4.1 先排除开发例的替代解释
96: 
97: 在公开b20d/40ac上，保留原Gii的Green19仍得到.0404/.0540；因此不能将改善归为删掉或修正了自项。只把旧one-hop的材料优化器换成相同8步L-BFGS-B，得到.3846/.6235，没有复制新矩阵分母的结果。错误传参导致的两次inner60试跑已中止、保留并计费，正式控制实际确为8步。
98: 
99: R9完整current低罚权方法也曾恢复40ac，因此矩阵分母不是项目首次恢复该对象，更不是已经胜过所有物理优化路线。普通磁场、Hodge或source encoding没有在本轮替代实现；没有将未运行路线记成贡献。
100: 
101: ### 4.2 输入/输出投影的单次结构修正
102: 
103: G的强右奇异空间是current输入空间。PGP额外删掉了场输出的(I−P)GP；GP则保留该强输入的完整输出。公开算子测得这部分Frobenius相对幅度16.83%，即约2.83%的强输入输出平方范数。三个末stage的残差变化为8.48%、3.76%、2.86%，两个版本的小系统均未在这些点接近奇异。
104: 
105: 这给出有根据的可检验修改，不能据小截断误差直接推导材料改善。本轮仅在已解封0001/0004/0022上把PGP改为GP，rank19、inner8、全部其他配置不变；结果属于**EXPOSED_POSTHOC_STRUCTURAL_DIAGNOSTIC**，不替换原24确认。PGP耗散定理不自动覆盖GP。
106: 
107: 三条实际修正已完成且封存后评分：GP的NRMSE分别为.10167、.51514、.43608，原PGP为.09889、.51514、.33700。环正例仍在，开缝壳没有改善，桥反而更差；三个完整child共701.64 CPU秒。桥的Phi从.01421降至.01040，但材料更差，GP到3000上限，PGP先由relative-F停止。不能把不同停止成本藏起来，也不能由GP的标准更佳Green逼近推导更佳反演。这个有根据的修正已被实际检验，没有继续调rank、β或阈值去救它。详见[三例完整表](delegated/right_projection_score/FINAL_COMPARISON.csv)。
108: 
109: ### 4.3 补偿CIE，而不是强行替换所有材料读出
110: 
111: CIE主确认中只有0001和0016未达到公开CN99噪声球；其余22例数据合格，但主材料成功只有11。这要求把数据未拟合的坏路径和数据合格的材料偏差分开。对后者，改优化器不能自动提供额外信息或保证选出GT。进一步利用已保存GT current核对：GT24/24均在各自实际观测的CN99球内；CIE数据合格但材料未达主描述门的11例，GT和估计都合法且同在该噪声球。这是具体的有限噪声双材料实例，不是精确多解、信息论下界或24例联合99%保证，soft prior仍可偏好其中一份。
112: 
113: 事后让CIE与Green19全部跑完、再按最低原Phi选择，平均NRMSE从.16287改善到.14960，但成功反降到10/24：增加0001，同时0000和0021因很小的材料变化失去原门槛成功。完整两臂CPU为8288.66秒、单CIE的2.19倍。没有免费丢弃未选臂，也没有用GT选best@2。这不是值得直接宣传的低成本改进。
114: 
115: 据此另冻结更保守且公开可执行的补偿：CIE已经位于噪声球内就原样返回；不合格才运行Green19并用同原Phi选候选。普通Dual-CIE同触发重试作为额外候选/成本对照。仅用四个全新high对象做最后独立确认，家族为环、开缝壳、桥、包覆体，seed2026110711，沿用家族contrast映射。没有因触发不足而追加困难对象。这个规则本身不是新的数学定理，价值取决于实际触发、额外恢复和费用。
116: 
117: **四对象确认现已完成，条件补偿得到一次新的独立恢复。** CN99白化范数门25.16542来自576个复观测、声明σ和99%卡方半径；不是用GT拟合阈值。两个对象触发，共四个fallback，全部选择封存后统一评分。
118: 
119: |全新对象|触发|CIE NRMSE|CIE→Dual|CIE→Green19|结论|
120: |---|---|---:|---:|---:|---|
121: |0000 χ12环|是|1.05000|1.05000|.10958|Green19恢复；IoU .04899→.71098，白化范数135.70→23.10|
122: |0001 χ18开缝壳|否|.12797|.12797|.12797|三部署完全同一CIE材料，fallback未运行|
123: |0002 χ18桥|否|.09071|.09071|.09071|三部署完全同一CIE材料，fallback未运行|
124: |0003 χ18包覆体|是|约.46018|约.46018|约.46018|全部失败，白化范数约702；没有按GT换方法或加预算|
125: 
126: 材料成功为CIE 2/4、CIE→Dual 2/4、CIE→Green19 **3/4**；平均NRMSE为.43222、.43222、.19711。三个部署视图的完整child CPU分别692.84、1015.55、1013.44秒；Green19补偿为单CIE的**1.46倍**，与同触发Dual重试成本接近。物理576阶LU为8054、11560、10225，新分母另有19阶小LU。三个视图共享CIE prefix，不能相加；实际整个实验八个科学inverse及一次零物理启动失败的unique child CPU为1336.14秒。
127: 
128: 这是“对未拟合数据的CIE增加一条有互补性的合法路径”的初步独立证据，比全体替换材料步更合适。它不会改善已经noise-compatible却材料不准的CIE终点；未触发时按构造原样返回CIE，而不是证明CIE已恢复GT。四对象只产生两个触发，不能据3/4对2/4作统计显著或总体成功率主张。环的效果在两批独立几何中出现，但尚不能预测所有环、壳或其他拓扑的收益。完整候选、未选择成本、失败和未运行臂见[独立审查](delegated/conditional_confirmation_score/REVIEW.md)。
129: 
130: ![四个最终确认对象的全部部署结果](delegated/conditional_confirmation_score/ALL4_DEPLOYMENT_MATERIALS.png)
131: 
132: ## 5. 三维和独立网格的边界
133: 
134: 本轮真正运行四个三维材料、六个观测条件共18主反演，使用12³材料、三分量dyadic Maxwell、60源及180接收分量；其中两对象为18³独立生成→12³反演。CIE/Global/A+三臂共用完整物理corrector。新低秩分母尚未移植3D，不能把原A+三维结果当作新分母的验证。
135: 
136: 源为外部真实电偶极，不是非法纵向平面波；非周期零填充FFT和CM+RR自项、同LU多RHS及伴随均继承实际全矢量实现。高ε9内部波长在粗网格约6.67格，细网格10格，gap/壳厚仅两粗格。18→12差异单列，不能靠同算子低残差宣称连续Maxwell准确。R9 Mie及原3λ电尺寸负证据仍然有效。
137: 
138: 18条主预测已全部封存后统一评分，无失败。下表是χ相对L2，不与二维幅度NRMSE混用；所有物理corrector都达到80接受步上限。
139: 
140: |材料/观测|CIE χ误差|Global控制|A+global χ误差|CIE/A+数据残差|
141: |---|---:|---:|---:|---|
142: |ε9立方体 matched12|.06082|.05486|.05758|2.979% / 2.978%|
143: |ε2壳/ε9芯 matched12|.49222|.48836|.49656|2.993% / 2.993%|
144: |ε9/ε2双体 matched12|.39200|.39929|.38306|3.082% / 3.069%|
145: |ε5空心体 matched12|.31816|.32140|.32124|3.004% / 3.005%|
146: |ε9立方体 18→12|.30701|.30372|.30005|3.173% / 3.173%|
147: |包覆体 18→12|.50825|.50199|.50806|2.998% / 2.999%|
148: 
149: 实心高介电体三臂都恢复较好，A+差别小且有正有负。包覆体与空心体虽接近噪声拟合，χ误差仍大；不能只看端口。跨网格的两份clean端口差分别24.80%和6.41%，超过3%测量噪声，且实心体材料误差由matched的约.06上升至约.30。这个离散效应远大于当前A+相对CIE的差别，不能宣传成连续3D优势。
150: 
151: 这里继承的选择目标Phi3=mean_s[||F_s−Y_s||²/||Y_s||²]、R0，无1/2；三臂一致，但不是生成iid σ的ML目标，也不同于二维全Y归一化Phi。真实σ白化残差另列，不能把源归一化冒称噪声白化。全部逐cell/区域指标见[三维原表](delegated/vector_score/scored/FINAL_SCORES.csv)。
152: 
153: 已保存场还给出更明确的解释：matched包覆体GT和CIE估计的白化残差平方分别10806.71和10754.59，均低于单case CN99界11043.23，虽然材料χ误差为49.2%；空心体也有相同的双合法材料现象。因此这些较大材料偏差不能全部归为“还没把数据拟合好”。反之，cross条件下，原coarse材料GT的残差平方分别749142和60516，远超同一噪声界。用这个coarse forward既要求返回原coarse材料、又要求拟合该fine数据到测量噪声内，两要求在本观测下直接冲突。不能把所有cross失败归罪于A+或CIE；也不能据此断言整个coarse材料类都不存在拟合解。核对只用保存场，新增forward为0，详见[噪声/模型分离](delegated/vector_score/scored/NOISE_MODEL_COMPATIBILITY.json)。
154: 
155: 双体80步后公开PG∞约1.3–1.9e−5，比其他组更明显未收敛。为区分预算与材料metric，本轮另从三个已封存m80出发，各追加最多80步同原物理corrector，三臂同样重置L-BFGS历史、不重跑aux；属于曝光后的预算诊断，不替换原18结果。其余已noise拟合的对象没有为追求分数盲目延长。
156: 
157: 追加三臂已全部完成：CIE的χ误差从.39200降至.35572，Global从.39929降至.36594，A+从.38306降至.35225；三者都从CN99噪声球外进入球内，真实相对数据残差约3.022%–3.029%。这直接表明原80步确有预算限制，但多算80步仍只带来部分材料改善，且A+相对CIE的差从.00894缩至.00347，不能归为A+特有收益。三者再次达到预算上限，PG∞仍5.4–9.1e−6；没有证明已经到同一个或不同的精确驻点。额外243次5184阶LU和29160复RHS、GPU整child496.69 CPU秒另计，原prefix不重计。见[追加评分](delegated/vector_continuation_score/SCORES.csv)及[材料前后图](delegated/vector_continuation_score/MATERIALS.png)。
158: 
159: ## 6. 成本、复现和论文资格
160: 
161: 二维主168任务实际完整Windows child CPU28475.421875秒，约7.91 CPU小时。四lane主对照和两lane候选各固定P核、BLAS1，并行wall不能相加称总历时。Green19平均完整CPU较CIE增加约19.0%，Two-fold约22.1%；这组数据不支持低成本的普遍成像改进。等CPU和等物理LU保存点另列，不能用相同3000辅助迭代替代公平费用比较。
162: 
163: GPU释放后已真正用于3D。对同一v08 CIE、同数据与预算，从FFT/GMRES改用库提供的retained dense LU，完整任务wall882.21→200.07秒，约4.41倍；物理corrector约5.05倍。第一次公开单点探针不足以外推整任务，这里另有完整任务实测。旧慢运行、中止、等待、setup、GT评分和标签修正均保留收费。此次优化不改精度/驱动/物理，不能计为A+算法收益。
164: 
165: 全部cost与unknown见[COST_LEDGER](COST_LEDGER.md)，复现入口见[REPRODUCE](REPRODUCE.md)。未提供硬件/模型工具账单费率，不能伪造货币费用。
166: 
167: **论文判断应落在受限的机制与算法证据上。** 当前可以写清楚：局部材料传播的能力与损害、两批独立高对比环的恢复、同表示下的正反路径、定性分母资格，以及真实成本。四对象补偿确认提供了可公开执行的简单算法正例，且与同成本量级的Dual重试有差别；但样本和触发数仍少，三维没有同类增量，closest-prior全文也有缺口。目前可形成诚实完整的机制论文初稿，尚不能将TAP/Inverse Problems级新算法优势视为已经充分建立。下一轮最值得推进的是同冻结补偿规则的多家族独立确认，以及公式中的`A_K⁻¹ D_m(G−K)R_G`如何预测有害/有利的材料方向；不是继续放大rank或β扫描。不要把这一轮称作“所有测试均需获胜才可发表”：真正缺的是可解释的适用条件、重复增量与费用论证。
168: 
169: 没有新增NIE、NN训练或Fresnel实测，相关已有结果继续引用R8/R9。Two-fold最接近文献的全文公式对照仍有访问缺口；已读摘要和明确文献来源保留。机制性论文草稿已提供，但期刊级算法优势不能由希望发表的结论倒推。
170: 
171: 最终根线程[科学审阅](FINAL_SCIENTIFIC_REVIEW.md)另列归因限制：最后四对象未运行Global-Dual fallback，不能把完整部署的增量唯一归为分母；主24则已有匹配global current控制。

## research/research_max_struct_r10/CHAT_HANDOFF.md

1: # MAX-STRUCT R10：最终研究交接
2: 
3: 2026-11-07。全部科学任务完成封存，旧R1–R9未覆盖。主要对比CIE、Dual-CIE与A+；Global-Dual只是matched-current内部消融，未新增NIE。
4: 
5: **结论：原A+本轮没有默认替换CIE的整体收益。新Green矩阵分母更适合作为数据未拟合时的补偿：最后四个全新对象中，CIE恢复2/4、同触发Dual重试2/4、Green19补偿3/4，完整CPU为CIE的1.46倍。** 这不要求所有例子都获胜；真正限制主张的是触发数少、反向失败明确，以及适用条件尚不能可靠预测。
6: 
7: ## 1. 全部24新二维材料 × 7方法
8: 
9: 八家族×low/high/mixed；576实未知、同数据、FFT/APCS、box、prior和共同原Maxwell corrector。168条全封存后读GT，无在线GT选根。
10: 
11: |方法|成功/24|平均NRMSE|平均完整child CPU秒|
12: |---|---:|---:|---:|
13: |CIE|11|.16287|157.71|
14: |Dual-CIE|10|.16346|153.98|
15: |A+global|8|.18044|169.49|
16: |A+sourcewise|10|.18045|174.98|
17: |Green19分母|8|.17669|187.66|
18: |Two-fold30分母|8|.18099|192.53|
19: 
20: low/mixed连续差很小，部分成功数只是跨过.15阈值。高对比有三个实质路径分叉：
21: 
22: - 环0001：CIE/A+约.423，新Green19 **.09889**；数据残差13.7%→.960%。
23: - 开缝壳0004：CIE **.11109**，A+/Green19均.51514。A类到低投影梯度坏拟合，非统一预算解释。
24: - 桥0022：CIE/A+约.097/.096，Green19 .337、Two-fold .437，新机制也会损害好路径。
25: 
26: 主材料描述门为预冻结NRMSE≤.15、固定χ>1 support IoU≥.7及合法性，不是研究停止门。140条corrector由梯度停止、4条函数变化停止、24条到上限。
27: 
28: ## 2. 分母、Two-fold和FFT能量归因
29: 
30: 原A+只把Gii放进逐cell分母，再做局部2×2块的一次传播。没有证据表明Gii写错。新实现：
31: 
32: `R_K = β [I+βD_m−D_mK]⁻¹ [J−D_m(Einc+GJ)]`，`K=PGP`。
33: 
34: 背景basis不读GT，576材料DOF不减少；Woodbury把分母缩为19/30维，完整材料导数保留。新集Two-fold实际30，历史DEV才28；修正了标签，没有改变basis或阈值。
35: 
36: 历史DEV保留Gii仍成功，仅换旧one-hop优化器未复制改善。PGP→GP的单次结构修正已跑三个曝光对象：环约.102、壳.515、桥更差至.436。“更佳Green逼近”不保证更佳材料。
37: 
38: 精确解析盘脉冲有限模型中，辐射Gram严格正定给出非负m/β下PGP分母定性可逆。不是浮点区间证书、条件数界、材料唯一性或全局凸化。Woodbury、Möbius、Two-fold及耗散理论不是本轮原创。
39: 
40: 重用69次历史结果449状态，并审计新24对象所有终点/同m末stage。七方法共享实际FFT/APCS容量；新末容量GT current不可表示能量占比仅.000225%–.014336%，却仍有大幅材料路径分叉。不能把A+收益归于更多Fourier模式，也不能用小oracle尾保证恢复。辅助/物理current、有符号cross分开；“能量”是离散平方范数，不是电磁储能；静态G/S模式不是材料Jacobian可观测性。
41: 
42: ## 3. 最后四个新材料：条件补偿独立确认
43: 
44: 生成前冻结：先完成CIE，只有真实σ白化范数>25.16542（单case CN99）才冷启动一个fallback，再按同原Phi选候选。Dual同触发重试控制额外候选/费用。
45: 
46: |对象|CIE|CIE→Dual|CIE→Green19|说明|
47: |---|---:|---:|---:|---|
48: |环|1.05000|1.05000|**.10958**|G19白化范数135.70→23.10，IoU .04899→.71098|
49: |开缝壳|.12797|.12797|.12797|未触发，保留CIE|
50: |桥|.09071|.09071|.09071|未触发，保留CIE|
51: |包覆体|约.46018|约.46018|约.46018|共同坏拟合，未加预算或换GT规则|
52: 
53: 成功2/4、2/4、3/4，完整部署child CPU692.84、1015.55、1013.44秒。部署共享prefix，不可加总；实际八次科学inverse及一次零物理启动失败全部收费。环改善在corrector前已出现，末stage NRMSE .19148，完整corrector后.10958；包含共同corrector贡献。包覆体三路均PGtol停止，非单纯迭代耗尽。
54: 
55: 未触发时按构造原样保留CIE，但不会修复数据已合格却材料不准的终点。主24里11例GT与CIE估计均在同一实际noise球、材料却不同；不是精确多解或恢复误差下界定理。只有四对象、两触发，不作统计显著或总体成功率主张，不与主24追溯合并。
56: 
57: ## 4. 真三维与4060
58: 
59: 四独立材料、六观测条件、18主反演：12³材料、三分量dyadic Green、60外部偶极源；两个对象另有18³→12³数据。三维比较原A+，新低秩分母未移植。
60: 
61: - ε9实心块：CIE/Global/A+ χ误差.0608/.0549/.0576。
62: - 包覆体：约.492/.488/.497。GT与CIE同在实际CN99球内，不能只归数据未拟合。
63: - 高低双体：.392/.399/.383；公平再加80物理步后.356/.366/.352，全部进入noise球，但仍到预算上限。预算有影响，A+与CIE差距反而缩小。
64: - 空心体：约.318/.321/.321。两跨网格clean端口差24.80%/6.41%大于3%噪声，不能称连续3D恢复优势。
65: 
66: GPU释放后真正使用4060 retained dense LU。同v08 CIE完整任务wall **882.21→200.07秒，4.41×**，端点残差接近。二维主要用i9 P核；旧慢运行及失败保留，不把后端提速归给A+。
67: 
68: ## 5. 论文与下一步
69: 
70: 可以形成“投影Green材料分母与CIE条件补偿”的机制论文初稿，具有正面恢复、定性资格、同表示正反路径和成本证据。尚不足以确立普遍更强的CIE替代算法或Maxwell凸化。
71: 
72: 下一轮优先固定现有补偿规则，做更多独立家族/物理条件确认；研究`A_K⁻¹D_m(G−K)R_G`如何预测方向收益/损害，而非扩大rank/β扫描。本轮未新增NN、NIE或Fresnel；R8/R9相关证据继续有效。Two-fold全文公式对照仍有访问缺口，不能宣称全球首次。
73: 
74: [完整中文报告](RESEARCH_REPORT_ZH.md) · [理论](THEORY_DERIVATIONS.md) · [四对象独立中文审查](delegated/conditional_confirmation_score/REVIEW_ZH.md) · [全部成本](COST_LEDGER.md) · [复现](REPRODUCE.md) · [英文论文草稿](PAPER_DRAFT.md)
75: 
76: 完整数组、模型状态、未选候选、长轨迹在同仓库[Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r10-20261107)分包，清单audit/PUBLICATION_MANIFEST.json。发布回执最终版以Git为准，无SHA验证。
77: 
78: 最终根线程[科学审阅](FINAL_SCIENTIFIC_REVIEW.md)另列归因限制：最后四对象未运行Global-Dual fallback，不能把完整部署的增量唯一归为分母；主24则已有匹配global current控制。

## research/research_max_struct_r11/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT-R11：Green分母的收益、破坏与修正
2: 
3: **最终状态：本轮科学运行与封存评分已完成，停止扩展。** 72个新二维材料对应432次比较，8个新三维设计对应50次比较；不同cohort和重复网格条件不合并制造成功率。代码、原始结果、失败与费用全部保留。按用户决定，交付分析建议后结束，不再设计或启动算法/实验。
4: 
5: ## 最终判断
6: 
7: **Green分母已有值得保留的条件成像优势：独立新环、以及跨网格的环与开缝壳，确实比CIE/Dual-CIE恢复得明显更好。** 这不是仅降低residual，也不要求所有对象全胜才有意义。与此同时，扩大48对象没有给出普通global控制之外的新成功，几个看似合理的分母修正也引入反例。因此应围绕“何时改变到有利路径、何时破坏路径”建立主张，不能把结果压成一句普遍有效或普遍无效。
8: 
9: 扩大到48个独立高对比/混合材料后，CIE/Dual/Global/G19/freezeA/contractedT分别恢复21/22/23/23/16/23个。Global、G19、contractedT成功的是完全相同的23个对象，新分母没有贡献Global控制之外的新增成功。只冻结非局部反馈的contractedT修复了新桥、梳状体上Green造成的大幅连续误差退化，比原Green少用约5.5%的有效反演CPU；但在共同实际成本前缀下，Dual仍为22/48、contractedT为20/48，平均NRMSE差的不确定区间跨零。完整表见[48对象报告](CONFIRMATION48_REPORT_ZH.md)。因此已经找到有用的局部修复，尚未找到值得默认替代Dual-CIE的分母设计。
10: 
11: 随后完成的8对象配对网格检验提供了更强的局部正证据：matched24与cross48→24条件下，原G19均6/8通过材料精度标准，CIE/Dual/Global均4/8，新增完整恢复为高对比环与开缝壳。形态图确实显示G19救回这两类，另外两个未达IoU门的多岛仍有大体形态，不能与灾难性错解混称“没成像”。与Dual的共同CPU前缀两条件也均4对6，但cross的新增对象换成环与混合多岛，不能称壳也同成本恢复。完整结果、费用与对象交换见[配对网格报告](CROSSGRID8_REPORT_ZH.md)。这支持有条件、可跨网格保留的收益，不支持用48组零专属增量否定它，也不支持用这8组覆盖48组反例。
12: 
13: Green材料分母有真实但有条件的恢复增量。它没有消除Maxwell非凸性，而是改变辅助current残差投向材料方向的方式；这一改变可以救回环，也可以把壳和桥导向错误区域。新16个独立高/低混合对象中，原Green材料分母恢复7/16，成熟CIE为4/16、Dual-CIE为5/16，平均幅度NRMSE分别为.14485、.23134、.21068；代价是Green实际CPU为CIE的1.269倍。完整逐对象结果及负图见[16对象报告](CONFIRMATION16_REPORT_ZH.md)。这些数字支持保留并继续检验Green，而不支持无条件替换CIE。
14: 
15: 必须把整条算法的差异与分母本身的增量分开：同global current的普通材料控制已为6/16。G19相对CIE增加环0000、壳0002/0003三个成功；相对Dual增加环0000与壳0002；相对Global仅增加环0000，没有反向成功回退。因此本组最清楚的分母专有效果是新高对比环，不能把三个增量全部归给Green。这个环上G19实际CPU54.48秒，低于CIE73.70及Dual64.79，但高于Global47.40；局部优势和全组平均成本应同时保留。
16: 
17: 本轮已实际尝试三种公式驱动修正：冻结完整分母、只冻结非局部feedback、将同一Green残差族接入current块。前两种保住开发环却没有修复开发壳/桥；第三种救回开发桥，但新16确认下降到5/16，并破坏一个原Green成功的新环。因此“两块一致”修正目前应被否定为默认替代，不用开发正例掩盖独立反例。
18: 
19: ## 1. 数学上究竟改变了什么
20: 
21: 固定辅助current，令
22: 
23: `e=J−D_mQ, A=I+βD_m−D_mK, w=A⁻¹e`。
24: 
25: 完整材料方向导数为
26: 
27: `dw[h]=−A⁻¹D_h[Q+(βI−K)w]`。
28: 
29: 分母同时改变空间残差metric和有效照明。有限rank K的metric改动既可增权也可减权，没有普遍的正定优越顺序。当前二维K为19维Galerkin截断，三维为公开source空间的rank120投影；方法名fullGreen的full仅指保留分母完整材料导数，绝不指已包含完整G。
30: 
31: 已保存10个开发状态的梯度拆分显示：抵消主要来自普通CIE已有的β收缩，而不是所有抵消都由新Green项造成。冻结分母的材料子问题即使为凸，也不是原移动目标的MM，更不等于Maxwell反演凸化。完整假设、推导与反例见[理论文件](THEORY_DERIVATIONS.md)。
32: 
33: current也改用Green后，均匀材料/full-G极限中的强传播通道从线性增长变为饱和，给出合理的修正动机；实际同状态检查却只见中等曲率改善及很大方向旋转，壳在条件数改善后仍失败。因此它解释“路径为什么能够改变”，尚不能从几何直接预测“哪个对象一定被救回”。
34: 
35: ## 2. 环、壳与桥的决定性干预
36: 
37: |已公开开发对象|原CIE|原Green材料分母|冻结A|冻结非局部T|两块Green|
38: |---|---:|---:|---:|---:|---:|
39: |环|.42288|.09889|.09985|.10339|.10517|
40: |开缝壳|.11109|约.51514|.51514|.52961|.51514|
41: |桥|.09741|.33700|.43375|.43710|.10088|
42: 
43: 表为幅度NRMSE，不能替代同时冻结的IoU及合法性判据。桥的修复已发生于辅助阶段，随后相同corrector的起点不同；不是同一起点使用更强corrector。补齐四格消融后，Green-current/普通材料单独用于环、壳、桥均失败，说明两块相互作用确实重要，但并不保证有利。
44: 
45: 新高对比环给出反例：原Green NRMSE=.1003、IoU=.7237；两块修正变为.5495、.2968。它足以否定把旧桥修复外推为通用方案。详细轨迹与全部曝光图见[轨迹分析](TRAJECTORY_ANALYSIS.md)。
46: 
47: 新增的完整物理线段诊断进一步定位了失败。只从已付费CIE/G19候选选取高Phi→低Phi方向，不用GT：错误环/壳分别有467/562个cell落在下界，投影梯度∞范数已约1e−8；但通向另一方法好端点的初始Phi导数分别为+.0151/+1.3385。预先固定的九点线段上，Phi从.00939/.17393先上升到采样峰.15283/1.59608，最后才降至约5e−5。桥也有相同方向的采样壁垒，但其梯度尚不能判收敛。峰值主要来自真实数据项，不能只责怪prior。这不是所有路径的能垒证明，却明确说明“同一坏终点再做局部下降”与“辅助阶段走另一条材料路径”是不同能力。
48: 
49: 同一原物理目标的二阶检查进一步区分了失败类型。环/CIE和壳/G19坏端点的自由face为62/14维，完整Hessian最小特征值为正，所有活动box乘子符号正确，投影梯度约1e−8；这与错误约束局部极小相符，但不是严格浮点证书。桥/G19仍有24个负face曲率方向及8个错号下界乘子，显然不能说它已经同样充分收敛。Green有时改变了进入哪个约束区域，有时也留下未完成的优化；不能把两者混成一种失败。
50: 
51: 基于真实Green的耗散符号，本轮又完成固定吸收型分母DEV：完整Green下current变换有统一增益上界，但可能压低近物理共振方向，因此不是凸化保证。固定β幅度、相位−i，比较有/无Green两臂的8次曝光反演后，Green版只救回壳（NRMSE .11063），同时丢掉环/桥/新环；K0版只恢复桥（.09721）。两臂各1/4，未发现默认替换价值。8次零执行失败、未触CPU上限，全部预测先封存再评分，花费759.28反演CPU秒；详细全成本另记。复CIE已有先例，不将复移位本身认作创新。详见[吸收型分母推导](ABSORPTIVE_DENOMINATOR_DERIVATION.md)及[全负结果表](delegated/diagnostics/absorptive_dev/FACTUAL_SUMMARY.md)。其后未继续扫描相位，按原冻结配置完成48与三维确认。
52: 
53: ## 3. 同成本、停止与材料误差的分开记账
54: 
55: 48主比较的有效科学child总CPU为49968.922秒。G19/Dual有效运行费用比约1.334，contractedT/Dual约1.261；有工程中断的case0022端到端费用未知，故上述比值不是包含所有中断的精确campaign倍数，另有47对象敏感性。共同实际成本B=min(两臂完整有效CPU)下，Dual/G19成功22/21，Dual/contractedT为22/20，双方仅按公共已付Phi选状态。48终点最优总误差不自动意味着同成本恢复优势。条件CIE→G19、CIE→Dual、CIE→Global也均23/48，G19没有超过付费控制。所有CPU、LU与原始连续指标单列。
56: 
57: 16对象的五臂全部同一Mac平台、同一允许CPU上限、相同材料空间/prior/观测和corrector。共同达到的全波LU预算下成功数仍为4/5/6/7/5；但辅助CPU不同，因此不能称同LU就是同总成本。240秒以上的允许预算对本组多数自然结束运行不再起作用，不用空闲填充制造“实际等成本”。
58: 
59: 额外后验预算审计为CIE和Dual-CIE中仍可继续且实际花费低于同对象Green的候选使用剩余CPU。公开停止状态及成本触发7次corrector重启；新增62.86 CPU秒，没有增加材料成功。由于L-BFGS历史重置、部分正常提前终止及原子操作尾部，此控制不是无缝延长，也仍须报告实际花费，而不能仅报统一上限。
60: 
61: 后验首次CN99命中检查也没有修复：65个命中状态的材料NRMSE全部比最终更差，成功标签0改善、19回退。已保存Phi-incumbents中没有一次先满足材料成功标准又在最终丢失。因此当前证据不支持“只需更早停止”解释这些失败；范围仅限已保存状态，不排除其他停止规则。
62: 
63: 65/80候选和GT都在同一实际观测的声明噪声球内，其中38个候选材料失败；另外15个候选在球外。前者不能据此称精确多解或Gaussian统计不可区分，后者也不能用信息不足代替算法解释。完整诊断严格区分物理Phi、数据噪声半径、prior与材料指标。
64: 
65: ## 4. 扩大确认与真正三维
66: 
67: 二维主确认48个独立材料、8种几何家族、high/mixed及χ=8/12/18的六臂288次反演已经完成。全部公共预测、144条件选择和384共同实际CPU比较先封存，再由root明确授权评分。逐对象连续误差、60/120/240/480 CPU前缀、已付LU轨迹与47对象敏感性均已交付。主算法没有因16对象结果改变。
68: 
69: 主队列完成132次后发生控制程序/Windows回执写入故障，打断同一case0022的四个在途attempt。已完成结果没有重跑，四个中断和三次零forward准备失败保留后另行恢复。原中断的完整CPU不能重建，按已知下界/未知null记账；成功重试的精确CPU不能冒充端到端支出。主48对象仍全部报告，另提供去掉该对象的47对象成本敏感性分析。[原中断清单](delegated/twod/independent48/CONTROL_ABORT_INVENTORY.json)与恢复回执保留，相关对象不作精确端到端同CPU主张。
70: 
71: 另在读取主48真值之前冻结并生成8个新独立材料，每个有matched24和cross48→24两条件，CIE/Dual/Global/G19四臂共64次。共同评分target为fine材料的2×2均值，两条件共享标准噪声draw、各自使用声明σ。64次在同四P核完成且零科学失败，全部公共选择封存后才授权评分；统计单位仍为8材料。G19环/壳完整恢复在两条件保留，细/粗清洁观测相对差为1.45%–5.30%，说明正例不只存在于同网格，但不是连续PDE误差证书。总科学CPU10921.40625秒，完整G19比Dual多20.8%/34.7%，共同前缀及全部失败另列。[完整配对合同](delegated/diagnostics/crossgrid8_protocol/PLAN_V2.md)。一次离线诊断检索误显示了私有生成metadata的单行inner字段，已记录流程偏差；线上worker、算法和公共选择没有使用该信息，也未改动冻结对象或重抽数据。
72: 
73: 三维新分母已经实际接入三分量dyadic Maxwell材料反演：12³个共享实材料未知量、每源5184复current自由度、60个外部偶极源、60个三分量接收位置。先做公开开发对象，再执行8个新设计×5臂，另有两个对象18→12跨网格数据的10次独立标记运行。采用CM+RR有限点模型与FP64原物理验收，不能称二维TM堆叠、完整连续Maxwell误差证明或2019作者设置的原样复现。新50次已经完整终态并在联合公共封存后授权解盲。8 matched对象的CIE/Dual/Global/fullGreen/contracted平均χ相对L2为.22844/.22608/.23127/.22780/.22774；平均IoU为.80869/.81248/.80382/.81335/.81259。fullGreen在5/8个对象的χ误差略低于Dual，但全组均值略高；主要形态均可辨认，没有二维环那样的大幅分母专属改进。两个cross对象的平均χ误差为.46483/.45439/.43356/.46714/.46319，Global最低。见[三维最终报告](VECTOR_CONFIRMATION_REPORT_ZH.md)。
74: 
75: 在新三维GT读取前作纯工程分流：最初Mac/4060各25任务；首臂冷启动计时促使尚未启动的c18十任务暂派4060，随后已付Mac暖态/历史current落盘计时表明可以更早完成，取消尚未执行的GPU追加队列，c18回到Mac，最终仍为25/25。每个对象全部五臂及matched/cross始终同host，原计划、变更、等待控制费用全部保留；原精度、算法和预算不变。跨平台秒数不混合成算法加速。Mac落盘只释放已归档历史场的RAM，未改变优化状态或数学计算。
76: 
77: GPU原队列15臂后，资源查询nvidia-smi超时导致c17的sourceDual在315接受步中断，额外320个5184阶LU保留计费。该臂同设置重跑后正常封存；查询重试避免再次丢失状态，但耗时明显。该臂完成后利用下一臂首次FW前既有目录保护自然退出，计数确认下一臂0新增FW，再由新进程使用30秒有效空闲查询缓存完成原九任务。这个变化只修资源协作，不是新算法或物理加速结论。受影响对象的端到端费用/前缀资格与7对象敏感性单列。全部50份预测及400个公共成本前缀经root逐项核对后才授权GT评分；没有新增科学调用。
78: 
79: 已完成的三维旧DEV12给出有限预算下的实际反例：v04 Global/Green的幅度NRMSE为.02120/.02038，但IoU为.8551/.8169；v09六臂均有高核心衰减，幅度NRMSE约.10、相对contrast L2约.49–.50。全臂均用完320 accepted corrector步，不能称已收敛或已证明三维优势。完整12臂、固定切片及成本见[真三维开发报告](VECTOR_DEV_REPORT_ZH.md)。新50确认已完成，算法未依据这两个DEV调参。40个matched终点全部CN99兼容，10个cross终点全部在球外；c13/c18干净12/18观测差达17.68%/3.53%，相对于3%名义噪声是实质模型风险。跨网格仍能看出空腔/双体等形态，但材料边界与数值偏差明显。它既不能全归咎于优化器，也不能证明粗模型不可拟合；所有50臂均到320 accepted-step预算终点，不能称充分收敛。
80: 
81: 二维IoU的support明确是材料超水平集`m>1`，不是相对真空背景的整个contrast support；全ROI未知且存在.1 host。三维另用冻结的有效材料超水平集`χ>.5`。跨维度成功数不混合。
82: 
83: 混合材料也不能只看以高相幅度归一化的总NRMSE：16组的低相真值为2.5，各方法low区域估计均值的对象中位数约2.426–2.467，但区域RMSE中位数仍为1.644–1.686。均值近似正确不等于低相空间细节恢复。保持原冻结成功标准，同时单独报告低/高相及背景误差，不临时调整门限制造更好结果。
84: 
85: 扩大48组也有同样限制：24个mixed对象中，Dual/G19/contractedT低相绝对χ-RMSE中位数为1.3210/1.3474/1.3331，未出现低相细节的清楚增量。更强的χ18纯高相组六方法全部0/8成功，Dual/G19/contractedT平均NRMSE为.23639/.29022/.21926；局部修复没有消除这一主要难点。χ18混合组三者均4/8。各contrast几何独立，不把该横截面比较冒充只改变contrast的因果实验。见[完整分相及六分组表](delegated/diagnostics/main48_final_review/mixed_phase_summary/SUMMARY.md)。
86: 
87: ## 5. 最终论文判断与交给用户的分析建议
88: 
89: 已成立的论文内容是：材料相关非局部metric的精确导数、作用边界、开发干预与独立反例、保留全部成本的确认比较。低秩惯性、resolvent预条件、冻结权重和变量变换本身均有成熟先例，不将其重新命名当原创理论。
90: 
91: 当前不宜提交“通用Green修正稳定超越CIE”这一主张，但绝不要求所有对象都胜出才有论文价值。48组限定了contractedT的局部修复与整体增量；16组独有环以及8组跨网格环/壳正例则支持有条件的真实优势。几何名称不能独自预测成败，平均分和二元门槛都不足以替代逐对象分析。已完成的真三维进一步限定范围：主要形态可恢复、个别连续指标有小收益，但没有观察到二维环/壳级别的分母专属突破；跨网格的共同模型差异突出。论文应围绕有证据的条件收益、失败机制与成本，而不是预设通用求解器胜出。本轮已经完成，后续仅由用户自行分析，不启动新的研究循环。
92: 
93: 建议保留原Full G19作为有条件有效的方法，以CIE/Dual-CIE为主要比较，Global-Dual只作隔离current更新影响的必要消融。不要默认推广freezeA、BothBlocks或吸收型分母；本轮反例已经说明它们会丢失原收益。contractedT值得作为解释桥/梳状体损伤的简化控制，尚不适合作为默认替代。
94: 
95: 接下来的分析重点应是已保存状态中的`Qeff=Q+(βI−K)w`、材料切向方向、活动box与进入corrector前的路径差异，而非仅看Green能量百分比或分母可逆性。二维同容量控制已经分离了FFT/APCS扩容；三维rank120远非完整可观测空间，不能仅靠提高rank来解释或保证材料恢复。跨网格c13的17.68%观测差提醒：应先明确模型与几何差异，再解释优化景观。以上是分析建议，不附带新的自动运行。
96: 
97: 方法边界见[METHODS_AND_ABLATIONS](METHODS_AND_ABLATIONS.md)，物理合同见[BENCHMARKS](BENCHMARKS.md)，实际成本见[COST_LEDGER](COST_LEDGER.md)，复现见[REPRODUCE](REPRODUCE.md)，英文工作稿见[PAPER_DRAFT](PAPER_DRAFT.md)。

## research/research_max_struct_r11/CHAT_HANDOFF.md

1: # MAX-STRUCT R11 · 最终交接
2: 
3: **本轮已完成并结束科学运行。** 72个新二维材料、432次比较，以及8个新三维设计、50次比较均已完整封存评分。不同cohort、同对象的跨网格条件不混成一组成功率。依用户要求，此后只交付分析与建议，不启动新算法、新对象、扫参或corrector扩展。R1–R10冻结证据不改写。
4: 
5: ## 首先纠正一句容易误解的话
6: 
7: **G19确实已经把所展示配对网格组的主要形态恢复出来了，环和开缝壳还有明显超越CIE/Dual的材料收益。** 未过严格IoU标准的多岛也不是“完全没成像”。论文不需要所有case全胜；需要把这些真实正例与48组、三维中的限制一起报告。
8: 
9: ## Full G19究竟是什么
10: 
11: 它是R10继承的实验材料读出，不是标准文献名称，也不是原A/A+的同义词。固定辅助current时，`e=J−D_m(Einc+GJ)`；普通CIE使用`(I+βD_m)⁻¹e`，G19使用`(I+βD_m−D_mK19)⁻¹e`。19指分母中公开固定的19维Green投影；Full指保留完整材料导数，**不指K19等于完整G**。
12: 
13: 原FFT-TSOM/current容量、576个独立实材料未知量、完整G正演及公共原物理corrector保持。三维对应fullGreen用rank120源空间，不能称G19。主要比较仍是CIE/Dual-CIE，Global-Dual只是隔离global current与材料分母影响的必要消融。[公式解释与方法表](METHODS_AND_ABLATIONS.md)。
14: 
15: ## 主要证据
16: 
17: |独立cohort|CIE|Dual-CIE|Global控制|原Full Green|结论|
18: |---|---:|---:|---:|---:|---|
19: |二维16材料|4/16|5/16|6/16|7/16|新高对比环是分母专属增量|
20: |二维48材料|21/48|22/48|23/48|23/48|后三者成功集合完全相同，未新增Global之外的成功|
21: |二维8材料matched|4/8|4/8|4/8|6/8|高对比环、开缝壳额外恢复|
22: |同8材料cross48→24|4/8|4/8|4/8|6/8|完整终点两项增量保留；不是新增8材料|
23: 
24: 表为原冻结材料精度标准，不能代替连续误差或图像阅读。配对8组的另两个多岛主要形态仍在。与Dual共同实际CPU界限下两grid仍4对6，但cross新增对象变为环与较早混合多岛；壳在Dual较短预算内没有完成同等形状恢复。完整G19 CPU比Dual多20.8%/34.7%。
25: 
26: 48组原G19比Dual多约33.4%有效CPU。contractedT修复新桥/梳状体的额外伤害，桥NRMSE由.36949降至.08307，形态明显改善，即使IoU .6762略低于原.7标准；但总体成功集合仍与Global相同，共同CPU前缀Dual/contractedT为22/20。冻结整个A明显更差；BothBlocks在开发桥上有效，却毁掉16组中的新环，不能升为默认方案。
27: 
28: ## 真三维已经完成
29: 
30: 12³、1728共享实材料未知量、每源5184复矢量current、60外部dipole照明；8设计×5法matched及2同对象×5法cross18→12，共50次。全部预测和400个公开成本前缀经root独立核对后才授权GT。
31: 
32: 8 matched平均χ相对L2：CIE .228436、Dual .226079、Global .231274、fullGreen-r120 .227803、contracted .227738。full在5/8对象的χ误差略低于Dual，但总体均值略高，主要形态与普通基线接近，没有复制二维环级别的跃升。2 cross平均分别.464833/.454391/.433562/.467138/.463195。
33: 
34: 40个matched终点全部CN99兼容，10个cross均在球外。两例粗/细干净观测差17.68%/3.53%，是实质物理/离散风险；cross图仍有形态，不能说完全没恢复，也不能把共同材料误差全部归给算法。50臂全到320 accepted-step预算，不称充分收敛；没有事后设置三维成功门。
35: 
36: [三维最终科学判断](VECTOR_CONFIRMATION_REPORT_ZH.md)。完整图集与复现位于delegated/vector；粗网格有限dyadic模型不冒充连续Maxwell精度证书或2019论文原样复现。
37: 
38: ## 公式给出的认识与没有给出的保证
39: 
40: `dw[h]=−A⁻¹D_h[Q+(βI−K)w]`说明分母同时改变残差metric和有效照明，可能旋转材料方向。10个真实开发状态中普通β收缩项大于新增Green项，不能把梯度抵消全部怪给Green。固定A材料块可以凸，却不构成移动目标的MM，更不等于逆问题凸化。
41: 
42: 原物理线段及Hessian把错误环/壳的近约束局部极小，与尚有负曲率/错号乘子的桥区分开来。Green能改变进入哪个区域，但目前没有一个仅凭几何名称、Green能量百分比或分母增益界就能可靠预测收益的规则。
43: 
44: 最值得用户继续分析的是已有`Qeff=Q+(βI−K)w`、材料方向/box活动与辅助stage2/3的差异。保留G19条件收益，以Dual作强基线；不要默认继续增rank、引入Two-fold或扩大反馈阶数。本轮没有启动下一轮研究。
45: 
46: ## 成本与流程边界
47: 
48: 全部失败和额外工作收费：主48有四个工程中断及三次零forward准备失败；原中断CPU部分未知，保留48主表与47成本敏感性。三维c17/Dual资源查询超时中断315步，额外320个5184阶LU计入，重试前缀只具attempt内资格，同时保留8主表与7对象敏感性。查询缓存是工程改动，不是科学算法速度。
49: 
50: 配对8组一次离线过宽检索误显示私有生成metadata的一行inner字段；未读GT数组或改变线上算法/选择，仍披露流程偏差，不声称全流程完美盲态。精确同算量、共同预算前缀、完整实际费用分别报告。不同平台/矩阵阶数不混算性能；货币费用没有费率，保持null。
51: 
52: ## 直接阅读
53: 
54: - [完整原始数组、轨迹与历史源码：40分包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r11-20261107)
55: - [完整中文报告与建议](RESEARCH_REPORT_ZH.md)
56: - [48对象完整比较](CONFIRMATION48_REPORT_ZH.md)
57: - [16对象正例与反例](CONFIRMATION16_REPORT_ZH.md)
58: - [配对网格报告及全图](CROSSGRID8_REPORT_ZH.md)
59: - [三维最终报告](VECTOR_CONFIRMATION_REPORT_ZH.md)
60: - [数学推导](THEORY_DERIVATIONS.md) · [文献与新颖性](PRIOR_ART_AND_NOVELTY.md)
61: - [英文论文草稿](PAPER_DRAFT.md) · [成本账本](COST_LEDGER.md) · [复现](REPRODUCE.md)
62: 
63: 论文最适合围绕“材料相关非局部残差metric的条件恢复收益与失败机制”，包含二维真实大幅正例及三维适用边界。证据不支持把它包装成通用凸化或所有维度更快的替代求解器；也不应因没有全胜而抹掉已经取得的成像增量。

## research/research_max_struct_r12/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT R12：双场反射、错误路径逃逸与材料选择
2: 
3: ## 2026-11-12 同起点配对补充
4: 
5: 新增44条受控运行，完整结果见[配对补充报告](supplement_paired/RESEARCH_REPORT_ZH.md)。它将原R11坏终点、R12同24对象的G19终点，以及四个冻结NN初值分别放进同对象/同起点/同whole-child CPU上限比较。
6: 
7: 最重要的新结果是：按原CN99球触发的4个G19球外终点，DR使3个显著改善材料且进球，G19合法重启/AP/direct未复制；完整24平均NRMSE从0.143891降至0.106261，新增DR费用约7.85%。0001 IoU0.659仍是有意义的主要结构恢复，原合取统计只新增2例，全部失败/20个旁路对象均保留。A组环上G19重启更好，壳/桥DR更好；NN组G19/DR质量接近，DR较省CPU，AP能解释部分NN收益但在0001严重退化。
8: 
9: **这是一组既有公开对象上的后续配对结果，不是新的盲确认。** G19是已有合法current投影重启、并非完整PR历史续跑；同预算上限不等于实际花费相同。原冻结结果与旧图在下文保留；新同起点对照应优先用于判断“DR是否比G19更善于救援”。
10: 
11: ## 科学判断
12: 
13: **本轮建立了两类不同的正面证据：DR在指定坏终点的再搜索和完整NN初值接续中有实际价值；G19在独立24对象的常规二维反演中有材料误差与结构恢复增量。** 前者主要是已暴露对象上的机制/接口研究，后者是冻结独立确认，不能混合成同一种统计。DR不能据此被宣布为全面替代CIE/G19，但也不能因为均匀初值整组排名没有胜出就抹去其条件优势。完整正反证据见[综合评价：DR、NN与G19](COMPREHENSIVE_EVALUATION_ZH.md)。
14: 
15: **这里的7/24只是旧定量合取条件，不代表只有7个对象有成像。** 按2026-11-12用户要求，对全24×8主图作统一的事后目视描述：DR有10个主要整体结构恢复、8个主成分恢复但弱结构不足，6个主要结构错误；Dual有19/24、G19与NN各20/24有S或P类可辨恢复。这是单人事后描述，不冒充预注册确认成功率。0009应称“主要高对比分支恢复，弱分支与幅值仍偏”；0007的G19、0013的多种改进法也有旧IoU门掩盖的结构增量。详见[结构恢复补充判读](STRUCTURAL_REVIEW_ZH.md)，原分数和原阈值统计保留。
16: 
17: 本轮更值得保留的结构认识是：**搜索能找到更好的材料，不等于原物理D+R目标会选择并保留它。** 两个高对比三维开发对象中，DR原始候选材料更准且处于声明噪声球内；共同corrector降低Phi，却增加数据误差、降低平滑prior并把候选推出噪声球。这不是“再多跑corrector一定更好”，也不是只看residual就能解释的现象。
18: 
19: 目前最诚实的论文定位是“反射可行性作为条件性的错误路径干预，以及搜索—物理可行性—材料选择的实证分离”。还不足以写成“普遍优于CIE的高对比Maxwell新算法”或“实现全局凸化”。不要求所有对象全胜，但独立样本上可重复的优势范围必须由证据界定，不能事后用GT挑出擅长对象。
20: 
21: ## 1. 本轮确实执行了什么
22: 
23: 用户原任务全文见[sources/user/R12_ORIGINAL_TASK.txt](sources/user/R12_ORIGINAL_TASK.txt)，原Pro包保留在[sources/pro/MAX_STRUCT_RESEARCH/](sources/pro/MAX_STRUCT_RESEARCH/)。旧R1–R11没有改写。两个Chat取得的是有界读取，不是无损完整导出；访问边界见[audit/SOURCE_ACCESS.md](audit/SOURCE_ACCESS.md)。
24: 
25: |研究组|实际工作|用途与边界|
26: |---|---|---|
27: |Pro源码重放|4对象direct/AP/DR与原后续controls|继承公开GT逐步诊断；不是新盲测，不把微小BLAS轨迹差异藏掉|
28: |成熟toy对照|4对象×5个CIE系方法×2种mode合同，共40运行|补足原Pro“只赢direct”的缺口；不同mode合同分别保留|
29: |R11错误终点|3对象×4臂，共12远端运行，另3本地开发运行|真实旧坏材料、同完整物理初态；预算前缀和超限完整端点分开|
30: |独立二维|6开发×8=48；24确认×8=192|全部确认预测封存后评分；24材料才是统计簇，不把方法/前缀当新材料|
31: |NN接口|4公开R8对象×4臂，共16运行|复用24,737参数冻结网络；不重新训练，不冒充独立确认|
32: |moderate矢量3D|2独立几何＋一个配对跨网格观测；18原运行＋18统一续跑|真三分量dyadic-DDA；较长共同corrector后差异消失|
33: |high矢量3D|2独立几何＋一个跨网格配对；4法×3观测=12运行|高χ8/混合2与8；原始/校正后分解，粗网格资格限制|
34: |实测与Austria|12原运行＋3完整Dual续跑＋3matched32诊断|真实Fresnel与Austria150→32；没有伪造实测pixelGT|
35: 
36: 上述跨网格观测、续跑和多方法不是额外独立材料。所有执行与失败成本见[COST_LEDGER.md](COST_LEDGER.md)，具体合同见[BENCHMARKS.md](BENCHMARKS.md)。
37: 
38: ## 2. 公式究竟利用了什么
39: 
40: 核心不是重新命名CIE残差，而是选择不同的搜索状态：完整总场Q与对比电流J。传播/接收集合A约束Q=Einc+GJ及观测球；材料集合C要求每cell的同一个实材料同时解释所有源、所有三维分量。固定尺度X=J/6，C投影化成共享斜率TLS小问题；A投影消元后复用H=I+36G*G和小端口系统。
41: 
42: AP和DR共享这两个投影。在同一个完整物理初态、ω=1时第一步严格相同；后续DR保留反射状态，允许经过不满足完整物理的中间对。这提供与原物理梯度下降不同的路线，却没有抹去原Phi的局部极小。若错误材料已经是合法且数据可行的物理解，纯反射同样会把它固定住。
43: 
44: [THEORY_DERIVATIONS.md](THEORY_DERIVATIONS.md)给出共享投影、条件shadow-gap→物理误差界及其假设。该界含材料resolvent，因此小shadow gap不自动说明端口精度，更不证明GT材料。iteration0的gap=0只是记录初态的占位值，不可用作集合一致证据。数值rank/边界容差明确公开。
45: 
46: 这些投影和DR数学工具已有先例。相近WRI、DCTMC、CSI、EM-ADMM和PADMM-CIE限制了新颖性主张；详见[PRIOR_ART_AND_NOVELTY.md](PRIOR_ART_AND_NOVELTY.md)。本轮贡献不是“首次使用非物理状态”或一般非凸全局收敛定理。
47: 
48: ## 3. 最清楚的正例：从已陷入的坏材料重新出发
49: 
50: R11环CIE、开缝壳G19、桥G19的保存材料原样作为共同起点。重新forward获得完整J/Q后，比较原物理续跑、固定分布平滑扰动多初值、AP和DR。
51: 
52: |幅度NRMSE，付费240 body-CPU内|续跑|平滑多初值|AP|DR|
53: |---|---:|---:|---:|---:|
54: |高对比环|.42288|.74888|.42288|.16585|
55: |开缝壳|.51514|.45633|.51514|.11377|
56: |桥|.43673|.43673|.35144|.10617|
57: 
58: 三者主要形态确实重现；DR的IoU约.647/.763/.695，仍有边界和宽度误差。按沿用旧描述门只有壳通过，不能因此把环和桥写成“完全没有恢复”。环与壳的旧点有近受约束局部极小证据，桥尚有负曲率/非收敛迹象，所以桥只支持路径干预，不单独证明严格局部极小逃逸。
59: 
60: 本地环审计显示Phi先从.009385升至第50步.014582，随后第100/200/400步降至.002002/.000167/.000065。AP保留旧坏incumbent。这个对照支持反射记忆对该路线有用，不能证明所有成功都需同样上升，也不能从三例推出普适盆地扩大。
61: 
62: ![真实旧坏终点全部对照](FIGURES/basin_remote_240.png)
63: 
64: 旧G19从均匀初值已经能恢复环，CIE能恢复壳/桥，本研究没有抹掉这些事实。完整结果及v1序列化造成预算超额的前缀更正见[BASIN_ESCAPE_REPORT_ZH.md](BASIN_ESCAPE_REPORT_ZH.md)。
65: 
66: ## 4. 独立24确认：收益没有成为通用优势
67: 
68: 24个独立对象覆盖levelsets、枝状、岛状、连续随机四族，高8/高18/高低混合各组合。固定24²、576个独立实材料未知量、12源48接收、1%复圆高斯噪声、公开box和原prior。八法有相同240 whole-child CPU上限、144proposal/96corrector窗口；记录实际支出和30/60/120/240前缀，不能把上限相同说成消耗必然相同。
69: 
70: 封存评分如下；[独立确认完整报告](delegated/confirmation/CONFIRMATION_RESULTS_ZH.md)提供逐对象、条件可用性、配对与全部失败：
71: 
72: |方法|平均幅度NRMSE|旧描述成功/24|
73: |---|---:|---:|
74: |成熟CIE|.202932|9|
75: |Dual-CIE|.185456|9|
76: |Global-Dual|.184007|9|
77: |A+|.171692|9|
78: |G19|.143891|9|
79: |DR|.184631|7|
80: |物理平滑多初值|.619361|4|
81: |冻结NN→direct|.149514|9|
82: 
83: 不能只说DR“没过gate”便停止分析，也不能说这组已证明所有结构无价值。它否定的是当前冻结DR作为普遍默认替代器的主张。反射从均匀场出发与从指定坏端点干预是不同实验；一些坏材料已经满足噪声球时，继续反射没有GT选择依据。
84: 
85: 新增次规则在24生成后、任何确认逆解/评分前登记：从已付费、合法、原始CN99球内候选中最小化同一R，Phi破平局。它没有增加forward，也没有改轨迹。八方法合计144/192个可用、48/192个不可用的次候选；单法最多24个，DR为18/24可用；不能删掉不可用对象后把均值降低包装成整体改善，也不能拿它悄悄替换主结果。最终需在相同对象交集比较主/次差异，并同时报告覆盖率。
86: 
87: ### 同成本前缀与配对不确定性
88: 
89: |whole-child CPU前缀|CIE|Dual|A+|G19|DR|NN→direct|
90: |---|---:|---:|---:|---:|---:|---:|
91: |30秒|.19411|.17677|.18169|.21718|.18551|.15044|
92: |60秒|.19876|.17942|.16787|.14988|.18396|.14881|
93: |120秒|.20149|.18404|.17014|.14411|.18408|.14951|
94: |240秒|.20293|.18546|.17169|.14389|.18463|.14951|
95: 
96: 均为全部24对象的平均幅度NRMSE，公开规则按Phi而非GT选择，所以材料误差可以随预算变差。G19早期setup较重，但60/120/240前缀已优于DR，不能把其最终优势全部解释为多用CPU。NN30秒已较强，然而它带有既有训练prior，不是同信息的纯优化创新。
97: 
98: DR−Dual平均差−.000825，按对象10000次bootstrap的95%区间[−.040982,+.053017]；不能据此宣称统计等价或优势。DR−G19为+.040740，区间[+.001934,+.095095]，本组数据支持G19平均材料误差较低；这是小样本、未校正多重比较的探索性区间，不是跨全部问题的定理。
99: 
100: DR17个旧描述未过中，11个在球内、6个在球外。11个球内对象均NRMSE≤.15，未过的是IoU≥.7这一项，relativeL2为.324–.600；应报告结构/边界不准，不能一概说完全没成像。原始DR proposal minPhi平均NRMSE.185814、最终.184631，12对象改善12退化；本组没有high3D那种原始进球→corrector出球的转移。因此不能把三维开发反例强套为独立二维的统一失败原因。
101: 
102: ![全部方法的CPU前缀](delegated/confirmation/FIGURES/offline_scoring_v1/confirm24_primary_CPU_PREFIXES.png)
103: 
104: ### 候选池是否隐藏了相同的好解
105: 
106: 另对三个R11端点的已记录付费候选池做OFFLINE_ORACLE诊断，未新求解：控制池best@K的NRMSE为.37454/.45202/.33477，仍高于DR在线返回的.16585/.11377/.10617。多初值确有被Phi丢弃的较好材料，但未追上DR。这加强了这些指定例子的搜索差异证据；只覆盖已物理审计的保存池，不声称所有未保存状态都差。
107: 
108: ## 5. 材料选择的问题：共同corrector为何会毁掉更好的结果
109: 
110: 两个high3D matched对象中，DR原始相对材料L2约.06163/.07831，corrector后变为.12978/.22724；CIE/Dual/AP最终也趋近相同材料。分解表明：
111: 
112: - 高χ对象D由1.0210e-4升至1.3577e-4，R由7.1732e-4降至5.9158e-4。
113: - 混合对象D由1.0034e-4升至1.3033e-4，R由3.5856e-4降至2.4505e-4。
114: - 原始DR均通过raw CN99，校正后两者均出球。
115: 
116: 因此Phi下降由prior改善支付，数据和材料同时变差。不能把同Phi校正器自动当作更可靠的最终裁判。这也是预注册次规则的依据。不过另外的AP例子中D与R都改善而材料仍变差，说明所有失败并不能统一归结为prior权重。
117: 
118: 三维高χ的6³网格内部cell相位大，8³→6³干净响应差约26.16%。这是真正矢量有限模型研究，尚非连续高对比Maxwell精度验证。moderate组在更长共同corrector后各法材料几乎一致，说明短预算优势可能是预算/收敛差异。详见[EXTERNAL_AND_VECTOR_REPORT_ZH.md](EXTERNAL_AND_VECTOR_REPORT_ZH.md)。
119: 
120: ## 6. NN、实测及经典对象没有被跳过
121: 
122: 冻结小NN的四对象接口中，uniform→direct、NN→direct、NN→DR、NN→CIE的平均幅度NRMSE约.4470/.1419/.0862/.3089。NN→DR的某些对象好于NN→direct，但两个对象轻微变差，完整CPU也增加；另一个对象CIE已解释相同收益。没有NN→AP单因素臂，不能把全部差异唯一归因于反射。独立24只确认NN→direct；不能宣称全部NN组合已经独立泛化验证。[NN报告](NEURAL_AND_SELECTION_REPORT_ZH.md)。
123: 
124: Fresnel用真实2001 single/twin原数据、空目标入射校准、源相关实际mask、711留出通道。多种算法与补足原生3000步的Dual最终接近相同corrector支配结果，留出残差约.0915/.0986。没有可信逐像素GT，因此不报告伪造NRMSE，也没有证明DR实测成像更好。
125: 
126: Austria原观测实际是150²→32²，旧64²元数据已保留并更正。粗GT残差.22526高于选用.12球，存在模型差异。随后用同粗GT、同预算/规则生成matched32，新GT残差.09951进球，DR仍差（最终NRMSE2.25391，残差.13490）；完整Dual也在corrector中从原始.60638退化至1.08143。因此物理模型误差确实存在，却不是全部失败原因。没有靠改阈值或反复换球挽救DR。
127: 
128: ## 7. 效率、失败和公开证据
129: 
130: 固定背景H只分解一次不等于免费：所有投影RHS、Green/接收乘法、物理验收、corrector、状态写出都计费。grouped v2在原toy固定投影workload的CPU为v1的约1/2.40，属于相同数学问题的工程优化，不是恢复优势或GPU速度结论。
131: 
132: 确认192子进程实际CPU为28147.796875秒，开发48有效子进程为7189.6875秒。初期13个路径求值启动失败合计.90625CPU秒，全在科学setup之前，源文件与失败回执保留。生成、运输、打包及后处理另记，内部/外部计时不重复累加。准备/阅读/部分旧timer缺失保持null。不能把并行CPU总和称作wall时间。
133: 
134: 远端固定独立P核、单BLAS线程；负载下遥测样本报告活跃核约1.5GHz，未中途改电源。五baseline显式BelowNormal，search继承controller；有一次有效优先级快照，但没有全程历史。故不宣称严格同调度下的wall吞吐优势。此次科学运行使用CPU，未宣称4060 CUDA已运行。
135: 
136: 完整源、模型、配置、原始状态、成功与失败由现有GitHub仓库及release分包交付；发布尚在准备时不会把未发布URL写成已可访问。复现入口见[REPRODUCE.md](REPRODUCE.md)。
137: 
138: ## 8. 能写什么论文，还缺什么
139: 
140: 可保留的论文主线：同一有限Maxwell模型下，反射在某些真实坏端点上产生AP/direct不具备的有效搜索路线；独立多样材料中这种能力没有成为通用优势；即使找到更准确的候选，同一物理D+R选择/校正也可能毁掉它。三种失败必须分别处理，不能用更低residual、更多迭代或更大的NN混淆。
141: 
142: 当前尚不支持：全局凸化、所有高对比对象优于CIE、连续高χ三维保证、实测材料优势、首次非物理状态/DR理论。论文资格仍缺最接近WRI/EM-ADMM/PADMM-CIE的合格数值比较、较高精度三维及机制适用区域的独立预测。未做完的这些对照应是明确限制，不以列举文献冒充已比较。
143: 
144: 下一轮最值得研究的公式是条件shadow-gap界中的resolvent放大，以及球内材料选择与soft D+R校正的差别。可先研究“何时保留观测可行候选、何时启动反射干预”的可观测判据；它必须在新对象上预注册，不能用本轮GT赢家事后拼成一个所谓全能混合算法。是否推进这个方向，应由逐对象确认诊断来决定，而不是本报告强行宣布新算法成功。
145: 
146: 英文初稿见[PAPER_DRAFT.md](PAPER_DRAFT.md)，科学审查见[audit/SCIENTIFIC_REVIEW.md](audit/SCIENTIFIC_REVIEW.md)。结果解释与发布状态分别记账；公开文件入口以CHAT_HANDOFF与ARTIFACTS为准。

## research/research_max_struct_r12/CHAT_HANDOFF.md

1: # MAX-STRUCT R12 — Chat 交接入口
2: 
3: 2026-11-12：原R12及完整归档已公开；新增44条同起点配对补充已封存评分，最新入口如下。
4: 
5: ## 最新配对补充：优先阅读
6: 
7: [中文报告](supplement_paired/RESEARCH_REPORT_ZH.md)、[GitHub直接图集](supplement_paired/GALLERY.md)、[费用](supplement_paired/COST_LEDGER.md)。A三坏终点×4，B原G19全24保留、4球外×4救援，C四冻结NN×4；共44条，先封存再评分。新图统一完整ID和whole-child CPU。
8: 
9: - B：G19→DR全24平均NRMSE0.106261，原G190.143891，增加约7.85%实际CPU；3/4触发对象材料显著改善且进球，2/4满足旧合取条件，0010仍失败。20旁路不运行，不是DR执行后无损。
10: - A：环G19重启0.114232优于DR0.140651；壳/桥DR0.112699/0.097269优于G190.515137/0.371690。当前比较是合法投影重启而非恢复原PR历史。
11: - C：NN→G19/DR均值0.086768/0.086094，实际CPU790.36/544.33；AP均值0.187543主要由0001失败造成，0022却能同样修复direct。secondary改变G19/DR均值排序，主/次选择规则均保留。
12: - corrector在新A与B救援中帮助精修，但C0004/0022仍损伤材料；旧3D反例依然成立，不一律删除或强制。
13: 
14: 本补充不新增材料/训练/算法，既有对象已暴露，不声称新独立确认。用户要求本组收尾，不自动开启新研究循环。
15: 
16: ## 综合评价入口
17: 
18: [DR、神经初值与G19综合评价](COMPREHENSIVE_EVALUATION_ZH.md)：DR坏终点干预是最强正面证据；NN→DR四开发对象平均NRMSE0.0862对NN-direct0.1419、NN-CIE0.3089，费用78.43对43.44/104.46 CPU秒；G19独立24平均0.1439对Dual0.1855且60CPU前缀已有正向。应分用途评价，不以全体全胜作为科研价值前提。
19: 
20: ## 最终科学判断
21: 
22: 双场反射（DR）是有用的指定坏终点干预，尚不是CIE/G19的通用替代器。它改变搜索路线，没有证明全局凸化。真正需要分开的是搜索失败、候选选择和模型误差。
23: 
24: - R11三个保存坏端点，同起点、同240 body-CPU上限：DR幅值NRMSE为环0.16585、壳0.11377、桥0.10617，优于对应direct/AP/平滑多起点返回值。已保存、已付费控制池的离线材料oracle仍不及DR在线结果；此结论不涵盖未保存状态或任意重启分布。
25: - 环/壳有近约束局部极小数值依据；桥尚未收敛且有负曲率，不能称三个严格局部极小逃逸。旧G19均匀初值能恢复环、CIE能恢复壳/桥，不能抹去这些历史结果。
26: - 24独立二维材料、8方法、192运行全部封存后离线评分。下表为全24主选择规则，非GT选根。
27: 
28: |方法|平均幅值NRMSE|旧描述成功/24|平均whole-child CPU秒|
29: |---|---:|---:|---:|
30: |cie_single|0.202932|9|136.609|
31: |cie_dual|0.185456|9|136.185|
32: |global_dual|0.184007|9|135.710|
33: |a_plus|0.171692|9|152.024|
34: |fullG19|0.143891|9|178.758|
35: |dr|0.184631|7|128.448|
36: |multistart|0.619361|4|219.280|
37: |nn_direct|0.149514|9|85.811|
38: 
39: DR−Dual配对均差−0.000825，95%未校正bootstrap区间[−0.040982,+0.053017]；DR−G19为+0.040740，[+0.001934,+0.095095]。不是等价证明或普遍优越结论。
40: 
41: DR的17个旧阈值未通过中，11个已进噪声球且NRMSE≤0.15，只是阈值IoU不足；另6个球外。24/24生成材料也在同球内，支持有限噪声下不同合法解释，不是exact非唯一性定理。次选择规则八法合计144/192可用，不可把条件均值当完整队列胜率。
42: 
43: ## 不把旧阈值未过等同无成像
44: 
45: 用户指出0009等例的IoU约0.5仍有明确结构。root逐图审视全24×8，新增[结构恢复事后判读](STRUCTURAL_REVIEW_ZH.md)：DR有10个整体主要结构恢复、8个主成分恢复/弱结构不足、6个主要结构错误；S+P为18/24，Dual为19/24、G19与NN各20/24。此为单人非盲事后解释，不改预注册统计，不得称新增确认成功率。0009主要高分支恢复、弱支和数值仍偏；0007 G19及0013多种对照的增量应得到承认。旧7/24不是DR全部成像能力的概括。
46: 
47: ## 三维、实测和NN
48: 
49: - 高对比真实矢量DDA两个matched开发对象：DR raw相对材料L2为0.06163/0.07831，corrector后0.12978/0.22724；D上升而R下降使Phi下降，候选从CN99球内到球外。独立二维未复制这一出球模式，不能泛化。
50: - moderate3D长corrector后各法接近；high跨网格干净响应差26.16%，尚不具连续高对比三维精度资格。
51: - Fresnel真实单/双介质，711留出通道，最终预测残差约0.0915/0.0986；没有像素GT材料准确率或DR实测优势。
52: - Austria实际150²→32²；matched32修正模型差异后DR仍失败，不能把全部失败归为inverse crime。
53: - 冻结24737参数NN四公开开发对象：NN→DR平均NRMSE0.0862、NN→direct0.1419，但增量集中个别对象且更贵。独立24只包含NN→direct，不宣称NN→DR独立泛化。
54: 
55: ## 入口与限制
56: 
57: - [中文完整报告](RESEARCH_REPORT_ZH.md)、[英文论文初稿](PAPER_DRAFT.md)、[理论](THEORY_DERIVATIONS.md)
58: - [独立24完整评分](delegated/confirmation/CONFIRMATION_RESULTS_ZH.md)、[全图集](GALLERY.html)、[研究首页](index.html)
59: - [科学批判审查](audit/SCIENTIFIC_REVIEW.md)、[成本](COST_LEDGER.md)、[失败](FAILURES.md)
60: - [归档与下载](ARTIFACTS.md)、[复现](REPRODUCE.md)
61: - [原任务及压缩后恢复](TASK_RESUME.md)、[原任务全文](sources/user/R12_ORIGINAL_TASK.txt)、[Pro原包](sources/pro/MAX_STRUCT_RESEARCH/README.md)
62: 
63: 审查者参与实现，不是外部盲评。论文可定位条件路径干预与搜索/选择分离；尚缺最近WRI/EM-ADMM/PADMM-CIE完整合格对照、更精细3D及预注册适用范围预测。不要把成熟DR工具或已知投影改称新数学。
64: 
65: ## 恢复与发布记录
66: 
67: 11月11日外盘I/O错误后掉线，部分晚期文稿/网页未持久化；12日恢复后由保存评分与原回执核对恢复。未重跑科学、未改变预测或评分规则。所有新增大文件只写外盘，内盘旧暂存不自动删除。当前发布未验收，不能将计划链接当已公开结果。
68: 
69: ## 2026-11-12最终发布
70: 
71: 完整原始证据已归档至[GitHub Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r12-20261111)，12个逻辑包、89分片，共11,416,183,356字节；逐个公开回读通过。新增综合评价与结构复核已纳入，不改变冻结原指标或在线规则。所有新增大文件位于外置盘。
72: 
73: 发布补充：新结果与完整原始证据Release已公开。验收时GitHub Actions/Pages官方重大故障导致新Pages仍排队；使用supplement_paired/GALLERY.md可直接在GitHub查看21图。Pages部署未冒称完成，状态回执与官方快照在补充audit。

## research/research_max_struct_r13/RESEARCH_REPORT_ZH.md

1: # MAX-STRUCT R13：有限材料事件、完整状态记忆与 Maxwell 分支可访问性
2: 
3: 本轮最终判断：**DR 是有实际价值的互补搜索器；G19 的恢复增量也真实存在。有限事件阈值解释了材料投影何时换分支，但还没有成为稳定优于普通 DR、简单重置或法向外推的新算法。** 最有价值的新认识是：同一个材料影子并不决定下一步恢复路径；完整辅助状态中的法向记忆和传播约束共同决定分支可访问性。它不等于把 Maxwell 反演凸化，也不保证跳向正确材料。
4: 
5: 研究已经完成真实代码、数值运行、封存、离线评分、二维/三维/实测图和成本归档。共 **962 份规范方法运行**，包括 917 份固定运行与 45 份条件救援；不是 962 个材料。新独立二维材料为 48 个，另有 3 个新三维几何的 4 个物理条件、2 个既有 Fresnel 实测对象和 Austria 的两种离散。固定 16 子集、不同方法、NN 初值与噪声条件不算独立材料。42 份早期物理导入失败和 4 份中断尝试另行收费、保留，不从记录中删掉。run-contract 的 1004 行包含 42 份恢复别名，即 962+42，别名不增加规范运行或独立材料数。
6: 
7: 入口：[全部图集](GALLERY.html)、[紧凑逐对象数据](analysis/review_tables/README.md)、[数学推导](THEORY_DERIVATIONS.md)、[方法合同](METHODS_AND_ABLATIONS.md)、[基线审查](BASELINE_AUDIT.md)、[英文论文初稿](PAPER_DRAFT.md)、[下载与复现](ARTIFACTS.md)。
8: 
9: ## 1. 本轮回答了什么
10: 
11: |问题|结论|证据等级|
12: |---|---|---|
13: |首次材料事件能否解析预测？|在有界、共享实材料的固定法向射线上可以；近似 AP 对的未来多步事件不能仅据该公式认证。|PROVED_UNDER_ASSUMPTIONS / DERIVED|
14: |跳过静止影子的整数步会增加搜索能力吗？|精确 AP 最近点对上只是同一 DR 轨迹的快进，不增加可到达分支。|PROVED_UNDER_ASSUMPTIONS|
15: |近似事件能否改变 basin？|能改变实际材料轨迹；开发对象出现明确正例，也存在极小状态变化引起相反结局的反例。|EMPIRICAL|
16: |阈值是否比一般扰动独有地改善恢复？|新 16 对象上在线事件、随机、倍增、无事件重置、记忆重置表现接近；未证明独有收益。|EMPIRICAL，主张未获支持|
17: |DR 是否能救援 G19？|可以，48 对象统一球外触发中有两个显著救回；事件也救回两个，但对象不同。|EMPIRICAL，条件优势|
18: |解除 Fourier/DPCS 是否自动更好？|不是；同 CPU 下完整 APCS/current 更差，且每步较贵。不能外推为更大空间本身有害。|EMPIRICAL|
19: |NN 是否有用？|冻结小网络明显改善 direct；DR/G19 上较小、方向有利但不确定。|EMPIRICAL|
20: |三维、实测支持新方法普遍领先吗？|不支持。已完成真实矢量反演和实测留出预测，但存在离散/照明误差，终点多接近。|EMPIRICAL / 物理资格有限|
21: 
22: ## 2. 来源、冻结与无 GT 在线合同
23: 
24: 继承 R4–R12、R12 paired supplement、Formula Atlas、FBG、P3L/Z-ALG 的真实实现。Pro ZIP 的报告、证明、代码、原始结果、协议、失败记录已抽取并审查；可见 Chat 缓存仅提供有限历史，没有获得两个会话的完整内部推理，未把预览当作授权或已验证事实。原任务及读取范围保存在 `sources/user/` 和来源审查中。
25: 
26: 算法冻结为 `3a91b35605634a18a731441693a46152eeb3363e`；发运 bundle 与后续三臂控制各有独立来源记录。全部新确认预测、公开候选选择及保留策略先封存，八组相关 seal 联合验收后才授权七组规范评分。GT 不进入 proposal、候选审计日程、在线选根或测试参数调整。[联合授权](audit/JOINT_OFFLINE_AUTHORIZATION.json) 与 [跨组输入核对](audit/CROSS_CAMPAIGN_PUBLIC_INTERFACES.json) 保留了这个顺序。
27: 
28: 原物理主规则始终是合法、预算内且已完成原 Maxwell 审计的最小 `Phi=D+R`。球内最小 R 与最小改变是另列的选择语义；不把它们偷换成同目标算法提升。primary、raw、corrected、secondary、付费前缀及预算外操作均有来源。每个对象的 contrast 幅度归一化连续误差、IoU、数据 residual、prior、球资格分别报告。**没有新增一个 IoU 门槛来定义所有成败。** 图像结构、孔洞、连接和幅度以事后描述给出，不包装成新预注册成功率。
29: 
30: ## 3. 公式带来的实际认识
31: 
32: 令 `X=J/s, s=6`，每单元跨照明（矢量时还跨三分量）共享实材料。材料投影可写为角度区间内的秩一实矩阵投影。对 `c=PC(a), d=a-c`，每单元存在 `c=u v^T, d=u_perp w^T`。设 `a0=||v||², b0=||w||², gamma=<v,w>`，沿 `c+t d` 的竞争角差满足
33: 
34: \[
35: \Delta(\phi,t)=(t^2b_0-a_0)\sin^2\phi+2t\gamma\sin\phi\cos\phi.
36: \]
37: 
38: 内点 `gamma=0` 时，首次竞争为 `t*=sqrt(a0/b0)`；边界由有符号二次方程决定。这个量回答“旧本构方向何时不再最优”，没有使用接收数据来承诺新分支是 GT。并列分支、零法向及材料区间条件均在数学报告中明确。
39: 
40: 精确 AP 对上，DR 在影子保持期间积累同一法向；整数快进与普通 DR 相同。近似对的传播缺陷 `e=PA(c)-a` 会被放大：一步差可界为 `(t+1)||e||`，但不是多步轨迹证书。状态扰动半径对应 Gram 扰动 `rho=2||Z||r+r²`，给出竞争 margin 的稳健性范围；若 r 只在跑完后测得，结论只是事后解释。
41: 
42: 本轮进一步给出带 resolvent 修正的 Maxwell 接收伴随联系。不能仅由“小 pair defect / gap”推断“材料压力近似正确”：实际环的 AP 对只有约 0.96% 相对 pair defect，方向修正却约为原方向的 41.7%；DR→A 的例子修正甚至约为 100–118%。因此可靠事件时间还缺少对 resolvent 放大和未来状态偏离的有效先验界。
43: 
44: 一个明确的抽象反例证明：两个状态可以有完全相同的材料投影影子，却在下一次反射后投向不同材料。这是投影动力学中的记忆机制，不是连续 Maxwell 多解定理。开发旧 B 的实际例子也显示：首次阈值仅差约 0.3% 的两个动作，初始材料影子最大差约 `9.64e-14`，完整状态相对差约 `0.0552%`，后续一条恢复、一条失败。它否定了“看材料跳变大小便足够预测收益”的解释。
45: 
46: ## 4. 新 48 对象：同实际 CPU 的主比较
47: 
48: 全部 576 实材料未知量，不用 shape 标签、GT 支撑或数字类别降维。36 matched，12 cross-grid；包括多岛、连接、孔洞、连续和高低混合 contrast。每方法 120 秒整进程 CPU 额度，proposal 72 秒、共同 corrector 48 秒；direct 使用完整额度。不可中断操作越限仍收费，但预算外候选不进入预算内结果。
49: 
50: |方法|平均 NRMSE|中位 NRMSE|平均 IoU|数据球内 /48|平均实际 CPU 秒|
51: |---|---:|---:|---:|---:|---:|
52: |CIE|0.19454|0.10279|0.69456|35|106.15|
53: |Dual-CIE|0.17985|0.09731|0.71592|37|106.26|
54: |G19|**0.16758**|0.10108|0.68577|39|119.08|
55: |G19 full APCS|0.29235|0.13189|0.52682|27|117.24|
56: |G19 full current|0.30985|0.27012|0.48188|24|115.86|
57: |DR|0.17672|0.09818|0.68726|39|**84.50**|
58: |整数事件|0.18149|0.09596|0.66314|38|86.33|
59: |同范数随机整数动作|0.19840|0.10040|0.65012|36|84.50|
60: |direct L-BFGS-B|0.49491|0.36858|0.28438|9|59.46|
61: 
62: 球内数是数据资格，**不是恢复成功数**。均值不替代配对图和逐例材料；所有 48 个都保留。[ECDF 与费用图](analysis/paper_figures/FIGURES/01_fresh48_ecdf_actual_cpu.png) 保留 direct 的完整长尾，另以插图放大密集区。
63: 
64: G19 比 CIE 的平均 NRMSE 低约 13.9%，比 Dual 低约 6.8%，但平均 IoU 并不更高；它有明确的材料收益，也有不同失败模式。DR 的平均 NRMSE 接近 Dual，实际 CPU 低约 20.5%，并能恢复一些 CIE/G19 失败的结构。这是值得保留的条件性竞争力，不要求它每个对象都胜出。
65: 
66: matched 36 对象的 CIE / Dual / G19 / DR 平均 NRMSE 为 0.21391 / 0.19437 / 0.17334 / 0.19979；cross-grid 12 对象为 0.13643 / 0.13628 / 0.15029 / 0.10751。分组同时改变物理/获取条件，不能把差异因果归于网格、噪声或孔径中某一个因素。36 个 matched 真材料全部在原 CN99 球内，12 个 cross-grid 的粗网格真材料全部在噪声球外：模型误差确实超过了噪声合同，未事后放宽 delta。
67: 
68: ### 正反例的结构解释
69: 
70: - `r13_fresh_0002`：DR/事件保留主要弯曲连通体，CIE/Dual/G19 输出碎片或错位环。G19 NRMSE 0.49750，DR 0.09363，事件 0.08562；属于实际恢复增量。
71: - `r13_fresh_0005`、`0031`：CIE/Dual/G19 保留孔洞或实心体，DR/事件退化为碎片。0031 的 CIE/Dual/G19 约 0.026–0.027，DR/事件约 0.386；反射没有自动优势。`0008` 则进一步区分 DR 与事件：普通 DR 保留有孔主体，整数事件塌缩为稀疏碎片，不能把两者合为一个失败结论。
72: - `0007`：direct 约 0.0772、IoU 0.852，而几类辅助方法输出错误分支。不能把 direct 的整体较差平均值说成它没有独特好例。
73: - `0029`：各主要方法的几何均错误；相似坏材料不能自动解释为数据不可辨识。其坏 Phi 仍远高于真材料参考，搜索和预算问题尚存在。
74: - `0042`：DR/full APCS 的嵌套环结构可较好，IoU 或幅度误差单指标不能抹掉可见结构；球外救援未解决其噪声资格。全 48 逐例说明见 [视觉审查表](analysis/visual_review/CASE_OBSERVATIONS.md)。
75: 
76: ### 基线充分性与完整 current
77: 
78: i9 上 CIE 46/48、Dual 48/48 完成原生 3000 步，G19 仅 1/48，其余被真实 CPU 截断；full 两臂各 2/48。故这是固定成本结论，不能称“充分运行的完整 G19 本质更差”。在新 GT 打开前另冻 M4 上相同 16 对象的 360 秒参考，CIE/Dual/G19 原生日程全部完成：平均 NRMSE 分别 0.14521 / 0.12464 / 0.13737，DR 0.11509。direct / physical GlobalDual / Full-J TRF 为 0.41081 / 0.36171 / 0.37670，均使用真实解析物理目标/导数，不能以一个弱 GN 冒充强基线。
79: 
80: M4 与 i9 的硬件和科学库版本不同，跨机差异不是纯预算消融。完整 current 从起步解除截断/硬锁同时改变 continuation、每步费用和优化几何；较小 profile 偏差并不保证更好的可访问轨迹。NFFT Eq3–6 的完整 current 适配已运行，但没有复现原作者全部优化细节，不能宣称击败成熟 NFFT-SOM。
81: 
82: ## 5. 事件、NN、A+：固定 16 对象机制审查
83: 
84: |方法|平均 NRMSE|范围说明|
85: |---|---:|---|
86: |CIE / Dual / G19|0.14920 / 0.12786 / 0.15008|同一预定 16 子集|
87: |原 A+|**0.10966**|实际平均 CPU 120.39 秒，越限收费；没有跑满 48|
88: |global-PR Dual-CIE|0.12620|不是 physical GlobalDual|
89: |DR / 整数事件 / 同范数随机|0.12313 / 0.12189 / 0.15973|整数匹配方向对照|
90: |AP / depth-one Anderson|0.36214 / 0.36155|普通 AP 局限仍在|
91: |AP128→DR|0.12413|准备成本已计入|
92: |在线分数事件 / 同范数随机|0.12352 / 0.12392|相同准备与候选审计|
93: |法向倍增 / 无事件重置 / 物理记忆重置|0.12290 / 0.11979 / 0.12146|阈值独有收益未分离|
94: |NN→DR / NN→G19 / NN→direct|0.11936 / 0.12656 / 0.27223|冻结 R8 小网络|
95: 
96: 在线事件相对 random/double/no-event/reset 的配对差均很小，描述性 bootstrap 区间含零。并非因为“未过显著性 gate”而否定方向；关键是最简单的对照已经解释了相近收益。旧 B 的强正例仍成立，但不足以推出新分布上的普遍阈值优势。
97: 
98: NN 只复用 24,737 参数的冻结权重，本轮没有训练。相对各自均匀初值，平均 NRMSE 差 DR 为 −0.00377，G19 为 −0.02352，direct 为 −0.13743；相应配对 95% bootstrap 区间为 [−0.01154,0.00197]、[−0.07238,0.00355]、[−0.24451,−0.04010]。NN 为 direct 提供实质帮助，DR/G19 上并非所有对象受益。训练先验是额外信息，旧训练费用未知，不能声称总体摊销净加速。
99: 
100: A+ 的这个 16 对象正面结果需要明确承认；它不是事件算法，也没有被 G19 或 DR 普遍取代。尚无完整 48、相同部署时间和外部物理分布的同等规模证据，故不能据该子集宣布 A+ 已成为最终最强方法。
101: 
102: ## 6. G19 坏终点救援：同起点，完整 48 分母
103: 
104: 按原 G19 primary 的公开 residual/delta>1 触发：9 个救援，39 个原样保留。每臂额外 120 秒 CPU，初始材料相同；合法重建内部状态及 fixed projector 均计费，未冒称恢复旧 PR 历史。
105: 
106: |部署策略|全 48 平均 NRMSE|全 48 累计 CPU 秒|相对原 G19|
107: |---|---:|---:|---|
108: |原 G19|0.16758|5716.00|基准|
109: |G19 合法重启|0.16762|6605.63|基本未改善|
110: |full-current G19 重启|0.16763|6607.55|基本未改善|
111: |direct 续跑|0.16769|6082.75|基本未改善|
112: |DR 救援|**0.15135**|6385.66|约 9.7% 误差下降，11.7% 额外 CPU|
113: |事件救援|**0.14947**|6349.20|约 10.8% 误差下降，11.1% 额外 CPU|
114: 
115: DR 显著救回 `0002`（0.49750→0.09113）和 `0024`（0.51817→0.14332）；事件救回 `0002` 和 `0017`，没有同样救回 `0024`。这是搜索路径互补。球内从 39 增至 41，仍不是 41 个材料成功。大部分改善来自两个对象，配对均值 bootstrap 区间接近/覆盖零；保留实际收益与对象差异，不制造强总体统计主张。图见 [完整分母与救援热图](analysis/paper_figures/FIGURES/02_conditional_all48_cost_trigger_changes.png)。
116: 
117: ### 一个必须披露的简单替代解释
118: 
119: 评分后增加了**不新增反演的事后组合诊断**：对 CIE、Dual、G19 各自球外对象，调用已经付费完成的“均匀起点 DR”结果，按原 Phi 选择，不看 GT。全部三种规则一起报告，不能只挑最好的。
120: 
121: CIE→uniform DR 全 48 平均 NRMSE 为 0.13377，CPU 加总 6197.83 秒；Dual→uniform DR 为 0.13317、6035.80 秒；G19→uniform DR 为 0.14803、6492.72 秒。它们没有新的条件队列实测历时，属于复用现成运行的 counterfactual，不是新的独立确认。
122: 
123: 这项反例说明：一部分部署收益可以由**已有求解器的互补失败模式**解释，不需要有限事件阈值。不能把 G19→事件包装成唯一最强的统一算法，也不能把事后更好的 portfolio 直接升格为已确认新方法。它为下一轮提供可冻结的简单对照，本轮不再调规则或重启算法循环。
124: 
125: ## 7. Corrector、驻点与候选选择
126: 
127: 原 corrector 后，新 48 的平均 raw→corrected NRMSE：CIE 0.19148→0.19454，Dual 0.18295→0.17985，G19 0.20908→0.16758，DR 0.20067→0.17672，事件 0.21822→0.18149。它平均上明显帮助 G19/DR/事件，也会损伤个别材料；Phi 下降不等于材料更准确。逐对象 raw→corrected 材料误差下降/上升数分别为 CIE 33/15、Dual 35/13、G19 39/9、DR 40/8、事件 41/7（各48，未改成功定义）。direct 的 raw 是初值，不能将其类似差值称作共同 corrector 损伤。
128: 
129: 把球内候选中“离 raw 最近”作为公开保留规则，平均反而较差：G19 0.19184、DR 0.20747、事件 0.22374。球内最小 R 也未稳定更好。因此本轮没有把 corrector 一律删除，亦没有把保留 raw 一律称为改进。保留候选集是必要证据，但选择它需要明确的材料偏好或更多信息。
130: 
131: 本轮另对 400 份固定终点作真实 box 梯度诊断，393 个不同材料，393 LU、7632 complex RHS；没有做新优化。多种错误材料具有很小投影梯度，saved Phi 与精确重算一致，支持“单纯续跑不够”的解释，但没有 Hessian/严格局部最小证书。0031 的 DR 坏终点投影梯度仍明显较大，不能把所有失败都归为同一种局部极小。
132: 
133: 付费候选池的 GT-best 只在封存后用于 OFFLINE_DIAGNOSTIC，绝不可部署。有限噪声球内存在材料相距较大的候选对，这说明当前观测加 prior 的选择存在歧义；在核心48中，43个对象有这样的已保存合格候选对；逐对象最大 RMS 材料距离的中位数为1.3416，最大3.5374，分别约为材料box宽度的6.71%和17.70%。这是有限噪声与软prior下的候选歧义，不是新的连续 PDE 非置换多解定理，不替代 Z-ALG 的既有结论。对于本身已满足两集合的错误材料，纯 AP/DR 固定点不会凭空得到 GT。
134: 
135: ## 8. 真三维矢量 Maxwell：已运行，结论有限
136: 
137: 复用 TensorDDA：6³=216 实材料，648 complex current 行/源，三个分量真正耦合、共享 scalar 材料；60 个外部偶极源，180 个接收通道。网格宽 0.24 lambda、h=0.04；使用 dyadic Green 与 CM+RR self，不是二维堆叠。3 个新几何含开缝框、非等高强耦合体、高低混合交叉支架，另做框的 8³→6³ cross-grid。
138: 
139: |几何条件|CIE NRMSE|Dual|DR|整数事件|定性说明|
140: |---|---:|---:|---:|---:|---|
141: |开缝框|0.28080|0.25368|0.28776|**0.22223**|事件误差较低，但空腔/开缝仍被填充；不是拓扑完全恢复|
142: |非等高耦合体|约0.09|约0.09|约0.09|约0.09|主要位置、组合形态保留；高 contrast 峰值明显低估|
143: |混合交叉支架|约0.07|约0.07|约0.07|约0.07|主要杆臂可见，细节和幅度不足|
144: |框 cross-grid|约0.227|约0.227|约0.227|约0.227|模型误差主导，方法终点接近|
145: 
146: 例如非等高体真 chi 峰值 8，重建约 3.1；IoU 约 0.778 的结果可称主要布局恢复，不能仅因未达 0.8 全部判失败，也不能忽略幅度偏差。32 页真实三维 voxel 图使用同相机、同 [0,8] 色标，全部方法展示，见 [三维图件](analysis/vector_volumes/FIGURE_MANIFEST.csv)。
147: 
148: 独立解析 Mie 的 εr=9 球，h=.02/.015 的远场差约 9.02%/7.06%，体平均场差约 19.76%/16.93%；几何细化有改进但不足以达到百分之一精度。跨网格框 clean mismatch 约 8.01%，超过 1% 噪声。全部 primary 噪声球比约 1.08–1.54，未通过纯噪声模型资格。正确定位是有限 DDA 下真实矢量反演与机制对照，不宣称连续高精度三维重建突破。
149: 
150: ## 9. Fresnel 实测与 Austria
151: 
152: Fresnel 原始 2001 实测 2 GHz 单/双介质柱，1053 训练通道（27×39），711 留出通道（新源 441、留出接收 270）。用空目标数据二阶多极照明标定，校准中位误差约 2.515%/2.408%；源相关 receiver mask、伴随均保留。15% 是模型 discrepancy 球，非已知 CN99 噪声球；接收 horn 与噪声协方差仍未完全合格。没有逐像素真值，因此没有伪造材料 NRMSE/IoU。
153: 
154: corrected primary 的训练误差约 8.76% / 9.14%，真正留出预测约 **9.15% / 9.86%**，各方法终点接近。DR raw 约 10.12% / 11.17% 留出误差；本场景 corrector 确有帮助。CIE 系 raw 的物理 incumbent 停在均匀初值，即使已运行辅助步并审计；这是这个 source-specific 适配的负结果，不能据此宣称击败原文成熟 Fresnel CIE。最终没有新的 DR/事件实测领先证据。[真实留出记录](analysis/fresnel_holdout/README.md) 与物理审查保留所有视图。
155: 
156: Austria 采用 400 MHz、2m ROI、20 平面波、40 接收、chi=3，分别做 150²→32² 与 matched32；不同于原 NIE 64²逆网格及其全部设置。matched 下 CIE/Dual 可以看到两个小体和环的主要布局，但幅度/环厚偏差大；G19 NRMSE 0.3983，CIE 0.4913，IoU 却分别约 0.55 与 0.693。DR/事件多碎片化。cross-grid 的粗真材料 residual 22.53%，超过 12% 诊断球，不能只怪求解器。未声称完整复现原论文恢复，更未声称 DR 在经典对象全面领先。
157: 
158: ## 10. 计算、失败与可复现边界
159: 
160: 主要计算使用 i9 P 核单 BLAS 线程，另用 Apple M4 参考；**本轮没有 CUDA 加速，也没有 GPU 训练**。i9 在线已知 CPU 为 80,241.94 秒（约 22.29 CPU 小时），Apple 在线 3,995.17 秒（约 1.11 CPU 小时）；它们不是并行历时，也不是同机速度比较。控制、生成、校验、评分、图、压缩、传输另列。不同维度 LU/RHS、固定背景投影/Gram/FFT 作用分开；不以材料 LU 少便宣称 DR 免费。
161: 
162: 详细阶段和未知费用见 [COST_LEDGER](COST_LEDGER.md)。早期失败、元数据/依赖恢复见 [FAILURES](FAILURES.md)。42 个物理任务最初在创建物理模型前因缺失继承 R1 namespace 失败，后来在独立 runtime 补齐原文件重跑，预算/算法/数据未改；原失败归档和恢复归档都提供。恢复别名不当成额外 42 个科学样本或重复计费。
163: 
164: 候选材料池、物理场与部分 DR/current 端点保存；没有每次 native PR 梯度/方向或每个候选完整 J/Q。因此不能宣称重建了所有优化历史。封存文件保留绝对路径/mtime 的原现场含义；在新路径复现应重新运行并封存，不能伪装原 seal 可直接验证另一份解压路径。
165: 
166: ## 11. 论文价值及下一轮最值得做的事
167: 
168: 可以写成一篇完整的**Maxwell 材料分支可访问性与反射状态记忆的机制研究初稿**，附条件救援、表示消融、NN 接口、corrector 和实际物理边界。已有价值不以全测试集全胜为条件：0002、0024、0017 等恢复增量、G19 的平均正面、A+ 固定 16 的较好表现，都值得保留。
169: 
170: 但尚不能把它包装成强 TAP 的广泛新算法突破：事件没有对简单控制显示独有增量，三维和实测没有稳定领先且物理误差较大。SIAM JIS 方向需要证明事件/完整状态机制相对已有非凸投影理论的不可替代贡献，或拿出可重复的阈值驱动效率收益；normal drift、RRR、局部双线性投影已有直接先例。Inverse Problems 方向尚缺新的连续 Maxwell 可辨识/稳定性定理，有限离散反例不能替代它。
171: 
172: 下一轮最值得冻结的最简单科学问题，是：**原物理残差触发的 CIE/Dual→DR 互补组合，能否在新对象上以较小额外费用稳定保留彼此的正例？** 这比继续增加事件阶数更有证据依据。理论上则应研究可观测且有 resolvent 控制的事件区间，而不是仅预测某个 cell 先跳。本轮未把这个后续建议当作已完成新算法，也不再启动新的调参循环。

## research/research_max_struct_r13/CHAT_HANDOFF.md

1: # MAX-STRUCT R13：最终 Chat 交接
2: 
3: 请先读 [中文完整报告](RESEARCH_REPORT_ZH.md)，再看 [综合图入口](GALLERY.html) 和 [全图集](GALLERY.html)。原始数组通过 [附件入口](ARTIFACTS.md) 提供。R4–R12 冻结结果没有改写。
4: 
5: ## 先给判断
6: 
7: **DR 与 G19 有互补的恢复能力；有限事件公式成立，但事件算法还没有显示出阈值独有的稳定增量。** A+ 固定16对象的正面也保留，不因本轮强调DR而隐藏。没有要求所有case全胜才承认价值，也没有把球内/低残差当GT恢复。
8: 
9: ## 完成规模与可信度
10: 
11: - 962份规范方法运行，48个新独立二维材料；重复方法/NN/子集不增加独立分母。
12: - 真三维：3个新几何、4个条件，216实材料、三分量dyadic DDA；不是2D堆叠。
13: - 实测：2个Fresnel对象，1053训练/711留出；Austria两种离散。
14: - 7组规范评分、8组相关seal先联合授权。预测/公开选择在新GT评分前冻结。
15: - 42次导入失败、4次partial中断保留。物理恢复用原R1依赖，未换算法/预算。
16: 
17: ## 必须继承的数字
18: 
19: 新48，i9，120秒整进程CPU额度：
20: 
21: |方法|均值NRMSE|实际均值CPU秒|
22: |---|---:|---:|
23: |CIE|.19454|106.15|
24: |Dual-CIE|.17985|106.26|
25: |G19|.16758|119.08|
26: |DR|.17672|84.50|
27: |整数事件|.18149|86.33|
28: |fullAPCS G19 / fullcurrent G19|.29235 / .30985|117.24 /115.86|
29: |direct|.49491|59.46|
30: 
31: G19大部分被CPU截断，不能称充分运行后的极限能力。另一个同16对象M4参考中CIE/Dual/G19全部完成原生日程，NRMSE为.14521/.12464/.13737，DR.11509。不要跨硬件拼速度。
32: 
33: 固定16：A+ .10966；DR .12313；在线事件 .12352、random .12392、double .12290、noevent .11979、memory reset .12146。没有阈值独有优势。NN→DR/G19/direct为.11936/.12656/.27223；主要改善direct，旧训练费用未知。
34: 
35: 注册G19条件救援：9触发39旁路，全48均值原.16758→DR.15135/事件.14947，全部费用+11.7%/+11.1%。DR救0002和0024，事件救0002和0017；收益集中，不能说每个坏终点都能救。
36: 
37: 事后诊断（不是新确认）：公开球外触发、复用已运行uniform DR，CIE/Dual/G19三种portfolio全48为.13377/.13317/.14803。较简单solver互补已经能解释部分收益，不可把G19-event叫唯一最强统一框架。
38: 
39: ## 理论到底新增什么
40: 
41: 1. 有界实材料normal-ray首次竞争角公式，区分内点/边界/并列。
42: 2. 精确AP整数快进同轨迹，无新搜索能力。
43: 3. 近似对需resolvent修正；小pair defect不保证小材料压力误差。
44: 4. 状态球竞争margin界，r未知时只有事后解释。
45: 5. same-shadow反例：材料相同不能决定DR下一步；老B的极近完整状态却不同终点是实际见证。
46: 
47: Normal drift、RRR、局部双线性投影已有先例，尤其Elser2021；不是本轮发明。有限离散结果不当成连续PDE唯一性/多解/凸化定理。
48: 
49: ## Corrector与物理边界
50: 
51: Corrector平均帮助G19/DR/事件；个别损伤。球内最小改变保raw并未改善均值，不建议一律停用corrector。保留raw/corrected/付费池很重要，但选择需要额外语义。
52: 
53: 3D开缝框事件NRMSE .222优于Dual .254，但空腔仍被填，不能称拓扑全恢复；另外两对象主要布局可见且峰值低估。Mie误差7–9%，crossgrid mismatch约8%，不声称1%精度连续三维资格。
54: 
55: Fresnel终点留出预测约9.15%/9.86%，方法接近，无新最终领先。Austria匹配下CIE/Dual保主要布局，G19材料误差较低，DR多碎片；不是原文完整复现。
56: 
57: ## 论文定位与下一步
58: 
59: 已经形成完整机制研究初稿。强TAP还缺更可靠3D/实测及新法恢复增量；SIAM JIS还缺相对成熟投影理论独有的预测/效率贡献；Inverse Problems尚无新连续Maxwell稳定性/唯一性定理。
60: 
61: 最值得下一轮冻结的是CIE/Dual与普通DR的残差触发互补，而不是继续扫事件过冲。理论应研究带resolvent误差控制的可观测事件区间。本轮未再开启这个新算法循环。
62: 
63: 代码、预算/名称区别、保存状态缺口、费用和复现均见 [METHODS_AND_ABLATIONS](METHODS_AND_ABLATIONS.md)、[BASELINE_AUDIT](BASELINE_AUDIT.md)、[COST_LEDGER](COST_LEDGER.md)、[REPRODUCE](REPRODUCE.md)。

## research/research_max_struct_r12/supplement_paired/RESEARCH_REPORT_ZH.md

1: # R12 配对补充：G19条件救援、NN接口与校正器
2: 
3: ## 最终判断
4: 
5: **DR最明确的新证据是作为G19球外终点的条件救援模块。** 24个既有对象全部纳入评价，按公开噪声球触发4个，DR在其中3个显著降低材料误差并进球；G19合法重启、AP、direct没有复制这些改善。另1个仍失败。完整24对象平均NRMSE从0.143891降至0.106261，约降低26.2%，新增DR救援CPU336.875秒，相对原G19整组费用约增加7.85%。这是既有公开对象上的后续配对结果，不能升级为新盲确认。
6: 
7: **DR并非所有坏终点上都胜过G19。** 环上G19合法重启的材料更准；壳和桥上DR明显更好。NN组中G19与DR最终质量接近，DR费用较低；AP能解释部分NN收益，却在一个对象上发生严重路径退化。因此不能把NN所有收益都归给反射，也不能因均值接近抹掉具体的路径差异。
8: 
9: 本次只执行用户批准的44条补充运行，无新材料、无新训练、无新分母/约束corrector。所有预测先封存，再统一评分。原R12数据、旧图、旧阈值和旧研究记录均保留。
10: 
11: ## 1. 对象、起点与预算已经配对
12: 
13: |组|实际对象|起点|方法与计费|
14: |---|---|---|---|
15: |A|r10_twod_0001/0004/0022：R11讨论过的真实旧坏终点|对应旧CIE/G19实际保存材料|G19合法重启、DR、AP、direct；各240 whole-child CPU上限|
16: |B|r12_confirm24_0000–0023，全24|原fullG19的primary材料|20球内原样保留；0001/0003/0010/0022四球外各运行四法，新增120 CPU上限|
17: |C|r8_test_0001/0004/0022/0027|同一个冻结NN输出|NN→direct/G19/AP/DR，各240 whole-child CPU上限|
18: 
19: 不同前缀的同尾号不属于同一个对象。A、C的proposal最多144秒、corrector最多96秒；B为72/48秒。direct可使用整个上限。所有方法保留同一物理Phi、box[.02,20]及原smooth/TV prior。
20: 
21: 新G19使用已有initialize_physical_current=True，把起点物理current最小二乘投影进原DPCS/Fourier affine空间，保留原global PR、双beta及full_material_readout。它重置PR历史，**不是原优化器无缝续跑**，也不是保留全部current的physical-complement变体。AP/DR初始材料、J、Q逐例完全一致，已由保存数组核对。
22: 
23: ## 2. A：同一个坏终点，谁能改善材料
24: 
25: 下表为原Phi在线选择的最终候选，NRMSE越低越好；没有GT选根。
26: 
27: |对象|G19重启|DR|AP|direct|
28: |---|---:|---:|---:|---:|
29: |环 r10_twod_0001|**0.114232**|0.140651|0.422880|0.422880|
30: |开缝壳 r10_twod_0004|0.515137|**0.112699**|0.515137|0.515137|
31: |桥 r10_twod_0022|0.371690|**0.097269**|0.369449|0.436727|
32: |均值|0.333686|**0.116873**|0.435822|0.458248|
33: |新增三条运行实际CPU秒|518.234|426.859|234.109|237.203|
34: 
35: 环中G19和DR都找回主要环结构，G19较准但费用也较高；不能把环再写成G19没有逃逸能力。开缝壳和桥的DR优势是真正的材料结构恢复，非单纯residual更低。桥原起点并非严格收敛局部极小，沿用既有资格限制。
36: 
37: 这里补齐了原实验缺失的G19同起点对照。它支持**条件性访问差异**，不证明所有G19热启动或无限预算失败。全三对象图保留GT、共同起点、四法raw与corrected。
38: 
39: ## 3. B：把救援连接到同一批24对象
40: 
41: 触发完全由原G19 primary的观测残差是否超过原CN99半径决定。未触发20对象未运行救援，保留它们是策略的旁路设计，不能当成DR运行后“无损”的实验结果。
42: 
43: |触发对象|原G19 NRMSE|G19重启|DR|AP/direct约|DR残差/δ|
44: |---|---:|---:|---:|---:|---:|
45: |0001|0.4318|0.4318|**0.1072**|0.4318|0.926|
46: |0003|0.4297|0.4297|**0.1191**|0.4297|0.922|
47: |0010|0.3658|0.3656|0.3656|0.3656|30.08|
48: |0022|0.2959|0.2959|**0.0283**|0.2959|0.916|
49: 
50: 0001恢复主要弯曲结构，IoU0.659仍有边界/宽度偏差；0003恢复主要高对比部分，低对比成分仍不完整；0022的连续材料恢复较好。三例进球及材料改善应正面报告。按原NRMSE≤.15且IoU≥.7合取条件，只有其中2例新增通过，原全24从9变11；它不等于只有两例有成像。0010仍是实质失败，保留完整图与费用。
51: 
52: 所有24对象主结果：原G19平均NRMSE0.143891；G19→G19约0.143883；G19→AP约0.143883；G19→direct约0.143883；**G19→DR为0.106261**。DR平均IoU从原0.601887增至0.676128。20个旁路对象和0010失败都在分母中。
53: 
54: 原G19整组实际external CPU4290.1875秒。新增四条救援运行CPU：G19重启355.000、DR336.875、AP239.328、direct62.500秒。因此DR策略总CPU4627.0625秒，原G19为4290.1875秒；其恢复增量不是免费。新G19重启使用了稍多于DR的实际CPU，仍未复制3例改善，这是较强的同成本竞争证据。
55: 
56: AP/direct实际花费较少，不能把共同上限说成相同实际支出。AP/DR均完成3000步且使用相同投影/审计合同，但投影乘子求解与后续corrector费用不同。未运行更长AP，不宣称任意等支出AP必定失败。
57: 
58: ## 4. C：NN、完整场接口与反射的贡献
59: 
60: |R8对象|NN→G19|NN→DR|NN→AP|NN→direct|
61: |---|---:|---:|---:|---:|
62: |0001|0.077571|0.077988|0.487912|0.076902|
63: |0004|0.082156|0.078824|0.080356|0.082944|
64: |0022|0.064368|0.067208|0.062564|0.289829|
65: |0027|0.122977|0.120354|0.119340|0.120818|
66: |均值|0.086768|0.086094|0.187543|0.142623|
67: |四条运行新增CPU秒|790.359|544.328|550.375|374.297|
68: 
69: - 0022：AP、DR、G19都能修复direct停留的坏材料。它证明扩大辅助状态搜索的价值，不能独归反射。
70: - 0001：AP进入明显错误结构；DR与G19/direct保留两岛结构。AP与DR费用接近（130.48对128.25秒），共同完整初场与投影也一致；因此这个差异不能简单解释为DR多花CPU。
71: - 0027：新G19合法热启动也能恢复结构。原NN→单beta CIE接口失败不能外推为所有CIE系初始化都破坏NN信息；本次G19同时改变了读出、global PR与current初始化，不能把修复唯一归因于投影。
72: - 0004：几种方法结果接近。没有理由把细小数值差写成重大恢复突破。
73: 
74: DR与G19主均值几乎相同，DR实际CPU低约31%；相对direct，DR材料均值低约39.6%、实际CPU高约45.4%。它是4个暴露开发对象的费用—质量观察，不是总体加速保证，也不是原R12局部机器费用的直接复刻。
75: 
76: secondary球内minR部署均值为G190.080104、DR0.082901；排序可随公开估计器变化，不能GT挑选有利规则。原NN共同预处理0.763496 CPU秒另计一次，旧训练费用未重测。本轮无训练和新推断，复用保存初值。
77: 
78: ## 5. Corrector：作用取决于到达的材料区域
79: 
80: 本次DR的raw→corrected平均NRMSE：A 0.200448→0.116873；B完整24 0.110176→0.106261；C 0.100916→0.086094。A三例与B三次有效救援中，反射先产生不同候选，corrector再完成有价值的局部精修；direct从原坏起点无法复制这一过程。
81: 
82: 但C0004从0.074037恶化至0.078824，C0022从0.055806恶化至0.067208，虽然Phi下降。C0022数据项和prior均下降，材料仍更差，不能把所有损伤都归结为prior权重。旧高对比三维DR raw0.06163/0.07832经原corrector变0.12978/0.22724、并被推出噪声球的反例也仍成立。
83: 
84: 因此本次不采用“一律删除corrector”或“一律强制校正”的结论。raw是已完成完整Maxwell评价的合法材料，corrector是局部D+R估计器。全部raw/primary/secondary及噪声球资格都保留，现有secondary缺失时明确回退primary。没有开发新约束corrector。
85: 
86: ## 6. 公式解释与不能越过的边界
87: 
88: 共享材料投影C与传播—观测投影A均未改变。AP/DR从同一物理初态属于C，omega=1时第一步相同；后续DR通过z←z+PA(2PCz−z)−PCz保留集合不相容状态。新配对显示这种记忆有时改变访问路径：B三例及NN0001不能由本次AP路径复制。
89: 
90: G19与DR的差别还包括辅助表示。G19的初始current投影相对损失A约0.818/0.956/0.845，B约0.806–0.980，C约0.437–0.716；原DPCS/Fourier族没有被偷偷扩展。不能将DR相对G19的全部收益唯一归为反射，也不能说G19数据/先验不同。
91: 
92: 新结果也说明“图像错误”不等于“该起点完全没有算法价值”。例如原独立组0022从均匀初值运行DR未恢复，而从G19保存材料重新建立完整物理J/Q后，上限120秒的救援流程得到NRMSE0.0283。二者的反射公式与物理观测相同，起点、阶段预算和实际费用不同，因此这只是互补线索，主要配对证据仍来自本补充；它支持把G19与DR视为互补路径，而不是将旧DR单独失败认定为其永远无法访问该材料。它不证明起点更接近GT，也不证明新策略已全局凸化。
93: 
94: 如果原物理材料已经在合法box和观测球内，从其完整物理状态启动纯AP/DR会成为可行固定点。因此本策略不会自动修复球内但材料不准确的解，20个未触发对象并非完美GT恢复。该限制是信息/选择与触发覆盖问题，不靠换评价阈值隐藏。
95: 
96: ## 7. 可信度、费用与论文含义
97: 
98: 44预测先封存，公开接线/预算验收无错误无警告，封存后经授权进行离线评分。源GT只用于离线评分，所有候选和前缀按已付费Phi或注册secondary选择。首次评分费用表误优先内部父CPU，v2改为外部完整父CPU；v1与修正记录均保留，预测未重跑。
99: 
100: 全部新增solver子进程CPU4669.46875秒，四核并行历时1240.1518秒；controller CPU5.765625秒另计。物理LU47716、complex RHS1042614；固定投影另有22个H Cholesky、756000迭代RHS，不能漏计或与物理LU混为相同动作。归档、失败启动/探测、评分和制图另见COST_LEDGER.md。
101: 
102: 这补上了“旧坏终点救援与常规G19比较脱节”的关键缺口，支持论文围绕**G19的观测触发条件救援，以及反射/表示/候选选择的分离**展开。它仍不是新独立确认、连续三维精度验证、对所有初始化的全局凸化或对所有CIE变体的优势。已完成补充后收尾，不再自动扩展新算法实验。
103: 
104: ## 入口
105: 
106: - [全部逐对象评分与成本](RESULTS/offline_v2/)
107: - [原始运行与候选](RESULTS/)
108: - [图集](GALLERY.html)
109: - [复现](REPRODUCE.md)
110: - [公开接线验收](audit/ACCEPTANCE_PUBLIC.json)
111: - [完整预测封存](ALL44_PREDICTIONS_SEALED.json)
112: 
113: ## 发布状态
114: 
115: 代码、逐对象数据、图和完整原始归档已提交既有GitHub，二进制资产已逐个匿名Range回读。可直接阅读[GitHub图集](GALLERY.md)，不依赖Pages。最终验收时GitHub官方报告Actions/Pages major_outage（2026-11-12新加坡04:48开始），Pages构建仍排队，因此没有把新Pages实际部署伪称验收完成。官方快照见audit/GITHUB_SERVICE_OUTAGE_SNAPSHOT.json；本地HTML与直接GitHub图像另行验收。

## research/research_max_struct_r8/NEURAL_WARM_STARTS.md

1: # R8：轻量NN初始化到底帮助了谁？
2: 
3: **已经实际训练并接入CIE、Dual-CIE、Global-Dual、sourcewise A+、global A+和直接原物理优化。最清楚的独立确认收益来自NN＋直接物理；简单接入A+经常损失这个收益。** 固定物理电流偏置产生部分修复，但没有成为稳定胜过直接物理的算法。所有结果来自先封存、后评分的876项预测，统计单位是42个材料对象。
4: 
5: ## 模型、数据和公开选择
6: 
7: 主模型为24,737参数两尺度CNN，16通道，十个物理/几何输入。1536训练、192验证对象独立生成，覆盖8类材料、contrast 2.5/6/12/18、复RMS噪声.001/.01/.03、6/12源及半/全孔径。输入含Born/DPCS复读出、场能量和几何，不含测试GT。CPU完成120epoch；独立验证选择权重，另训练同结构Born-only控制。代码和权重路径见[复现](REPRODUCE.md)，完整冻结条件见[协议](NEURAL_PROTOCOL.md)。
8: 
9: 42个新材料分为24个训练家族内、12个未训练家族、6个36²生成→24²反演对象。每个保留576实材料自由度；不以形状或GT support缩减。材料box与prior沿用公开合成合同。主630=42×15臂，额外120仅在预先指定8对象上做初始化消融，偏置扩展126=42×3。两套全部完成后才开测试GT；未因NN初始Phi较高而拒绝运行。
10: 
11: 所有最终端点取已支付原Maxwell Phi的incumbent。辅助方法最多3000步辅助＋3000步共同物理校正，直接方法最多6000步；另给同实际CPU/LU/RHS allowance下的已执行端点，不能把相同最大步数当相同成本。没有把GT用于初值、阶段、根或停止选择。
12: 
13: ## 独立确认结果
14: 
15: 下表是幅度归一NRMSE均值，括号是继承复合标签成功数；标签还要求IoU≥.7和合法性，不能只凭这个二值数判断连续材料质量。
16: 
17: |方法|家族内24|家族外12|跨网格6|
18: |---|---:|---:|---:|
19: |均匀＋直接物理|.24966 (2)|.40474 (1)|.19593 (2)|
20: |NN原始材料|.14690|.21898|.12052|
21: |NN＋直接物理|**.09647 (9)**|**.16245 (6)**|**.08431 (2)**|
22: |均匀＋Dual-CIE|.36843 (7)|.63825 (1)|.32821 (2)|
23: |NN＋Dual-CIE|.32822 (7)|.56779 (2)|.32835 (2)|
24: |NN＋Global-Dual|.28919 (8)|.60481 (1)|.32834 (2)|
25: |NN＋global A+|.29609 (7)|.59393 (2)|.32933 (2)|
26: |NN＋sourcewise A+|.32182 (7)|.59164 (2)|.32811 (2)|
27: |NN固定偏置＋Dual|.29007 (7)|.52299 (4)|.25988 (2)|
28: |NN固定偏置＋global A+|.31395 (5)|.51831 (4)|.22800 (2)|
29: |NN固定偏置＋sourcewise A+|.31293 (6)|.50160 (4)|.22833 (2)|
30: 
31: 完整15臂、固定8对象的m=.1/1、Born、DPCS、随机、两个NN和joint状态消融、按对象配对区间与停止原因在[原数字备忘录](delegated/warm_analysis/FINAL_NUMERIC_MEMO.md)。不能把固定8对象的Born-only NN局部优势拿来事后替换42对象的冻结网络。所有42对象及所有失败见[主图集](FIGURES/warm_material_all_42.pdf)、[原NN与后续物理图集](FIGURES/neural_initial_to_final_all42.pdf)、[偏置图集](FIGURES/warm_anchor_confirmation_all_materials.pdf)。
32: 
33: ## 数学解释与被否定的简单预测器
34: 
35: 令物理电流为`jp=jphys(mNN)`，硬数据投影为Pd，早期Fourier可调空间投影为Pu。投影初始化的差可拆成
36: 
37: `δjD = Jd − Pd jp`，`δjN = Pu jp − (I−Pd)jp`。
38: 
39: 两项在欧氏current范数下正交，但Maxwell本构缺陷是`(I−XG)(δjD+δjN)`，强耦合可以放大、混合这两个方向。因此NN材料很有用，并不意味着送入硬DPCS/低模态辅助族后仍保留其搜索轨迹。
40: 
41: 最小修正是固定`Ja=Jd+(I−Pd)jp`，只继续优化原Fourier系数。它不增加优化DOF，不额外求初始LU，但完整current存储、投影和所有后续求解仍计费；最后完整APCS空间与原法相同。8个验证对象的初始本构缺陷全部下降，材料却不都改善，已经否定该标量诊断作为独立成功预测器。确认集显示它可以部分修复家族外/跨网格退化，但仍远逊于NN直接物理。
42: 
43: 预定对象0027展示另一个反例：NN初Phi=.30965，global A+把校正前Phi降至.01565。然而NN直接物理在1814 accepted steps按投影梯度停机，得到NRMSE=.12082；global A+用满3000校正步仍为1.01370。这里直接支路实际少于3000步，所以不能仅用6000上限解释；仍不能据此证明错误吸引域不可逃逸。
44: 
45: ![预定家族外对象0027，全部使用相同合法色标](FIGURES/neural_initial_to_final_fixed_r8_test_0027.png)
46: 
47: 再保留两条完整支路并按最终Phi选择也不够。NN直接＋global A+的均值变成.07754/.34498/.22882：家族内改善，另外两总体恶化，CPU为单NN直接的2.38–3.17倍。Global-Dual组合出现相似现象，不是A+专有机制。此选择完全不读GT；它说明增加候选和降低数据目标本身不能保证更准确材料。
48: 
49: ## 成本与摊销
50: 
51: 相对均匀直接物理，NN直接的最终平均CPU增加1.06–2.74倍；均匀初值常较快停在坏材料。因而这一比较支持质量提升，**不支持最终时间加速**。同对象共同CPU/LU/RHS allowance下，五个硬件/总体层的平均误差仍较低，置信区间、未用完预算及稀疏快照间距完整报告。
52: 
53: 相对五种均匀和五种NN辅助方法，NN直接在全部五层的最终平均CPU和误差都更低。这十组细比较是在看主结果后加入的描述分析，没有重新调算法；小层区间常跨零，不宣称普遍统计胜出。同CPU也有反例：ID/i9的NN直接减NN Global-Dual平均误差差为+.00281，区间[−.06069,+.06401]。
54: 
55: 例如Mac ID部署CPU为NN直接32.03秒、均匀Dual62.12秒；家族外分别39.12/69.89秒。物理特征与推理都收费。主physics网络训练加训练/验证生成为402.558889 Mac累计CPU秒；Born-only和采集域修复训练另列全研究账本。相对均匀直接没有有限CPU摊销盈亏平衡点；相对较贵辅助法的条件性成本收益也不能免去训练、加载、候选验证或跨硬件价比。详见[成本对照](delegated/warm_direct_budget_analysis/REPORT.md)与[全账本](COST_LEDGER.md)。
56: 
57: 若明确改用**均匀Dual-CIE**作为成本参照，三组Mac平均每次节省约30.09/30.77/31.82 CPU秒。主网络及训练/验证生成402.56 CPU秒的条件性warm部署摊销为14/14/13次。它只按本次对象均值记CPU账，排除冷加载，不是经济盈亏或等质量保证；全研究中Born/domain训练仍全部收费。两机器之间没有指定CPU价格比，不给混合硬件经济摊销数。计算见[条件性摊销表](delegated/warm_direct_budget_analysis/CONDITIONAL_DEPLOYMENT_AMORTIZATION.json)。
58: 
59: ## 与A+互补还是替代？
60: 
61: 历史7b8d和b20d上NN＋两A+仍能恢复，NN＋Dual不能，属于已公开对象的条件互补。新42对象则主要是NN直接更好，不能继续把“NN＋A+”当默认组合。两种情况都保留；不存在一个从历史两例外推的通用关系。
62: 
63: Austria采集域另外训练了无Austria标签的轻量修复模型并实际下游求解，长校正仍没有稳定收益。Fresnel两圆柱也实际试过NN，没有解决测量模型问题；主网络没有3D输入/训练，未声称3D迁移。
64: 
65: 最值得继续的不是加深网络，而是处理**有噪声DPCS约束与学习物理状态之间的冲突**，并研究什么时候应保留直接物理轨迹。现有投影分解是可测诊断，尚不是能可靠预测材料恢复的规则。若新方法不能在冻结新对象、同实际费用下超过NN直接物理，就不能把一个更复杂接口称为新突破。

## research/research_max_struct_r8/BENCHMARKS.md

1: # R8 数据、物理与比较合同
2: 
3: 不同组不合并成一个“总成功率”。同一材料的多算法、初值、噪声或离散也不计为新的独立材料。完整来源分别见 PAPER_BENCHMARKS、PAPER_FOLLOWUPS、HISTORICAL_WARM_STARTS、FRESNEL_REPORT、VECTOR3D_REPORT 与 NEURAL_PROTOCOL。
4: 
5: | 组别 | 独立对象/证据性质 | 观测与材料 | 适用的比较 |
6: |---|---|---|---|
7: | 历史困难回归 | 四个早已公开对象；21个新端点、19个旧封存端点复用 | R6的24²实scalar-TM、固定旧Y、576实未知量 | 保留旧A收益、NN与sourcewise/global的关系；不能作新盲测 |
8: | 独立神经初始化确认 | 42个新材料：24同分布、12未训练家族、6跨网格 | 原24²反演；最后6个36²生成；6/12源、半/全孔径、三噪声等级 | 原750项固定比较，加126项明确记录的仿射电流机制扩展；全部预测封存后才开测试GT |
9: | 开发验证 | 固定前8个验证材料 | 与训练物理相同但独立材料种子 | 48原始轨迹、8最小偏置实验、24其余偏置对照；用于发展机制，不能加入确认分母 |
10: | NIE Austria | 一个经典几何、四个介电/频率条件；不是四种独立拓扑 | 2m ROI、20平面波、40接收、10%复RMS噪声；原150²生成/64²反演及32²适配 | 核对原论文物理合同、独立NIE实现、初值域偏移、长corrector；无额外prior的原合同单列 |
11: | Fresnel2001 | 实测单/双圆柱 | 真实空目标和散射场，36源、每源49有效接收；1053训练、711保留通道 | 校准与源相关接收算子、同数据五方法、固定NN迁移；无像素GT |
12: | Fresnel2005 | 实测FoamTwinDiel | 8GHz、18源×241接收；2496训练、1842保留；60²全未知材料 | 原NIE场景适配、五方法、独立NIE；另6+8GHz是新增信息实验 |
13: | 真正3D矢量 | 球、四立方体、包覆/空腔/有损体等小型对象；同对象多个算法不增样本数 | TensorDDA/dyadic Green，3分量场，独立8³/10³生成、6³/8³反演，60外部电偶极 | 实标量A+矢量推广、sourcewise/global、独立Mie与网格；不是2D堆叠 |
14: | 原3D电尺寸尝试 | 2019论文3λ域的实际缩网格尝试 | 16³生成、12³反演、60源；原论文60³/30³未完全复现 | 保存真实材料反演及截断成本；不能称收敛或原尺寸全复现 |
15: 
16: ## 合成主比较
17: 
18: 采用 exp(−iωt)、λ0=1、εb=μb=1，实contrast m=εr−1。背景外源/接收器均在未知ROI之外。ROI宽.96λ；继承历史disk-pulse自项和Green，圆盘半径.45h。这是冻结有限模型，不用新的等面积模型悄悄替换历史任务。接收48；外部线源与接收半径分别1.4λ、1.1λ。全域材料box[.02,20]，576独立实未知量；没有形状分类、字体或GT支撑缩减。
19: 
20: 主目标 `Phi=.5||F(m)−Y||²/||Y||²+R(m)`；同组所有方法使用相同smooth/TV prior。它是相对数据能量目标，不能叫按σ白化的似然。复噪声为 CN(0,σ²)，σ指复幅度RMS，即实/虚标准差σ/√2。声明的99%白化噪声球半径来自2p个实自由度的χ²分位数；跨网格另有模型误差，不能用纯噪声门认证其全误差。
21: 
22: 材料NRMSE为RMSE除生成前声明的contrast幅度；归一分母置于独立评分元数据，反演不读。阈值形态是 `m>.1+.5a`，预测由自身256bin Otsu生成；因为ROI内host=.1，形态mask不是相对真空的完整contrast support。连续光滑/嵌套材料的IoU只作形态诊断。沿用历史NRMSE≤.15且IoU≥.7的标签便于比较，同时逐例报告连续误差；标签不用于在线选择或停止实验。
23: 
24: 普通/NN主臂保持相同FFT-TSOM阶段mf3/6/10各500、mf12另1500，成熟辅助预算最多3000；共同original-Maxwell corrector最多3000 accepted。独立直接物理控制最多6000，明确多出上限，同时用已支付轨迹做相同CPU/LU/RHS预算比较。不能仅用最终端点把不同实际成本称相等。对NN本身不因初始Phi偏高而拒绝后续运行。
25: 
26: 原750项中630=42×15主臂，120是预定8对象的初始化/Born网络消融。用户要求的NN sourcewise/global完整因子比较在未读任何测试GT时补入，保留amendment时间。新增126=42×3仅改变固定初始APCS偏置；原权重、材料prior、算法参数与预算不变。此扩展由公开开发验证推导，不删除原750项，不混称最初一次性预注册。
27: 
28: NN由1536训练材料、192独立验证材料训练；原8家族和12个未训练家族对象的划分见生成器。模型和测试对象种子均保存。训练集中的不同源/噪声配置属于其所属材料对象，不跨split重复同一材料状态。额外paper-domain网络只用于明确的经典基准修复，不替换主确认的原权重。
29: 
30: ## 论文与实测的可比边界
31: 
32: Austria原合同采用空气初值、零额外prior、等面积圆盘和论文的source/rank/mode设置；这些与历史.96λ、非零host任务不同，不做混合平均。独立实现NIE意味着按原方程和停止规则重现，不等于获得作者软件或作者认证。32²早期80步corrector是探索，补充1000步及原64²结果单列，不能用短基线作最终胜负。
33: 
34: Fresnel首先校准空目标照明，逐源mask保持真实有效通道，使用相应adjoint；保留源与保留接收分别报告。少数名义几何/复gain诊断是OFFLINE_MODEL_DIAGNOSTIC，不送入主反演选根。原五方法继承[.02,20]下界，不能表示精确空气；独立NIE另用air0合同。发现此问题后完成同锚点单变量物理对照，再以[0,20]从共同初值重跑Dual-CIE与A+完整原预算；保留原五方法和全部25端点。修正了材料集合的两臂不能与旧五臂混算同合同算法收益，也没有像素GT指标。原始测量需从来源合法取得，发布派生结果及获取脚本；不把合成数据称实测。
35: 
36: ## 三维资格和不足
37: 
38: 3D状态及极化是真矢量；各向同性介电的A+共享一个体素材料参数，而非三个分量各选材料。自项为标量乘I3时实材料Möbius条件结构成立；复材料采用合法复参数优化，未套用实圆弧。原投影IDFT系数PR与曾用的SVD白化坐标分开报告，后者不能冒充原CIE轨迹。
39: 
40: Mie验证、独立生成网格和实际反演是三种证据。Mie测试的几何/电尺寸不同于所有成像例，不能把其误差直接贴到每个成像网格。小网格跨离散偏差仍显著，3λ原尺寸尝试被明确截断；未完成原论文60³/30³反演、Fresnel3D实测或3D NN训练，不伪装已完成。
41: 
42: ## 公平性与统计
43: 
44: 同一新材料的所有方法分配到同一硬件，绝对CPU按Mac/i9分层。线程wall之和不是项目历时。网络训练/数据生成、冷载入、特征、推理、两个完整候选支路及失败均收费。固定PR迭代、相同全波动作和相同实际CPU是不同对照；FFT/辅助材料工作不能藏在LU计数之外。
45: 
46: 主要分析按独立对象配对。小cohort给逐例和配对bootstrap区间，不把多算法端点数当样本量或夸大显著性。公开Phi组合先封存选择再评分，付费两支路均计入，best@K真值选择只允许明确离线诊断。不同组的NRMSE分母可能不同，三维原报告中的relative L2 NRMSE不得与主组幅度归一NRMSE直接合并。

## research/research_max_struct_r7/BENCHMARKS.md

1: # R7 数据、物理与评价合同
2: 
3: 本轮实际实验分为30个新二维材料主比较、6个重复对象初始化对照、公开历史开发诊断、Fresnel实测先导、真正矢量DDA反演及有限材料有理搜索。实际完成数以最终评分/费用回执为准；本合同不充当完成证据。
4: 
5: ## 新二维确认集
6: 
7: 独立生成30个材料，每个均保留576个独立实contrast未知量，不使用类别、轮廓、支撑或真材料初始化。六种解析几何各五个独立几何实现：曲条、偏心空腔、开缝壳、细弯连接结构、多岛、近距离胶囊体。几何参数只在生成/评分端；公开包没有family/mask/GT。随机种子2026110407，生成与反演是不同进程。算法开发仅用R4–R6已公开对象；新30对象在全部主预测封存后才评分。基线与原A先冻结，A+开发单独完成并在其新对象反演之前冻结；生成先于最终A+选择，不宣称全算法冻结后才生成数据。
8: 
9: 五个预定条件，每条件六个不同材料：
10: 
11: | 条件 | contrast幅度 | 复RMS相对噪声 | 源 | 源孔径 | 生成/反演网格 |
12: |---|---:|---:|---:|---|---|
13: |0|2.5|.001|12|全孔径|24/24|
14: |1|6|.01|12|全孔径|24/24|
15: |2|12|.001|6|全孔径|24/24|
16: |3|18|.001|12|半孔径|24/24|
17: |4|12|.005|12|全孔径|36/24|
18: 
19: 这些因素按条件捆绑；不能从每组6对象估计单因素的独立因果效应。每个对象的方法比较是配对的。六个额外初始值实验使用对象00/06/12/18/24/29的m0=1，对照原m0=.1，只是重复对象，不是新增独立材料。
20: 
21: 时间因子exp(-i omega t)，lambda0=1、k=2pi、mu=1，实正epsilon_r=1+m。物理对比电流密度j= m E，省略共同的-i omega epsilon_b因子。源为ROI之外的二维线源，接收在背景，完整近场散射E观测。ROI为.96lambda方形，源半径1.4lambda、接收1.1lambda；48个接收。24网格h=.04lambda，36网格h=.02667lambda；源极化为标量TM出平面Ez。无3D伪称。
22: 
23: 原disk-pulse有限VIE包含解析自项、Hankel1辐射，不是周期FFT边界。每cell圆盘半径.45h，填充率pi*.45^2；这不是连续满方形像素精确积分。36→24打破相同离散inverse crime，但不是独立连续PDE误差证书。几何按固定8×8子采样面积平均；实际细结构/高contrast网格误差由跨网格残差显露，不把代数线性残差当物理精度。
24: 
25: CN噪声E|e_rs|²=sigma²，实/虚各sigma²/2。sigma由生成端声明为relativeRMS*||Yclean||/sqrt(nr*ns)，公开sigma与99%噪声球；反演端不读clean。sigma不是离散模型失配的上界。
26: 
27: ## 方法和预算
28: 
29: 主比较：成熟sourcewise CIE、Dual-CIE、Global-Dual、原A-material、冻结A+、原物理L-BFGS-B。前五保留R6 beta/Fourier模式/PR/材料prior：mf3/6/10各500，mf12为1500，共3000辅助步；共同原物理corrector最多3000accepted。direct从同m0最多6000accepted，没有辅助算法开销。没有45秒wall截断。所有最终材料按同原Phi在线保留。
30: 
31: Phi=.5||F(m)-Y||F²/||Y||F²+R(m)。R为继承smooth=TV=1e-4，差分先除19.5，smoothedTV epsilon=.01。box=[.02,20]。这是历史归一化目标，不冒充sigma白化高斯似然；噪声合格性另列。各方法都用真全波LU物理评价；辅助传播使用真实有限Toeplitz零填充FFT，接收与物理corrector仍dense。
32: 
33: 相同辅助/校正iteration上限、相同CPU、相同LU或RHS是不同对照。每次已支付物理评价每10次保存incumbent材料，离线在共同CPU/LU/RHS范围取不超过该预算的最近保存检查点；因此是检查点分辨率的近似匹配，不能称逐时刻严格同步。float32预算材料仅用于离线评分；在线和最终材料/物理值float64。FFT、材料内循环和setup均收费。不同维度LU不合并作等成本证明。
34: 
35: ## 评分和统计
36: 
37: RMSE=||mhat-mref||/sqrt(576)，NRMSE=RMSE/预声明contrast幅度；support为生成几何面积分数>=.5。预测自身Otsu256阈值IoU，绝不GT调阈。连续误差为主；继承image成功为NRMSE<=.15且IoU>=.7且合法box。物理噪声99%合格与image成功分别报告；跨网格不能假设噪声球包含模型误差。所有失败、上限和未收敛保留分母。30个材料是统计簇；配对bootstrap按对象，不把源/初始值/迭代态当独立样本。每条件只有6个对象，区间可能宽，不能以未显著证明等价。
38: 
39: Fresnel的无逐像素GT合同、校准和留出通道见FRESNEL_AUDIT.md；三维材料/状态维度、自项、Mie几何误差、必要的体积分数修正见BENCHMARK_3D.md。二者各自先导结果不混入30个二维合成统计分母。