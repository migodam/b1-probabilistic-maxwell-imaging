# 原论文经典二维／三维案例扩展审计

状态：**READ_ONLY_AUDIT；R9扩展尚未运行（NOT_RUN）**。本次只核对原文、已有源码、冻结配置和成本回执；没有启动正演、反演、训练或新GT评分。原R9第一阶段113次反演及其 `COST_LEDGER.md`、`audit/ALL_COSTS.json` 保持冻结；后续经典案例必须另建phase2增量账本，历史复用不得重复计费。

## 1. 最小选择与结论边界

最小的**两个新运行对象**应是：2016论文的二维Austria **εr=3.5、400 MHz**，以及2019论文的三维**四立方体 εr=2/3/4/5、原3λ域**。它们都是原论文明确声称CIE/NIE优于LSIE的正例，完整几何可由正文确定；不能仅沿用形状后缩小电尺寸。二维εr=4保留为已知困难负例；若补做，也必须预先这样命名，不能把它当作原论文承诺成功的条件。750 MHz、εr=2 Austria是额外电尺寸控制，不是最小两例的必需追加。

R8已覆盖真实二维采集与150²合成数据，但主要九臂比较是32²反演；64²只完成εr=3.5的一个CIE臂和εr=4的五臂。R8三维主比较是**0.6λ、6³/8³**的明确适配；原3λ案例只有12³/16³ direct截断试算。因此，既不能说原论文3D高对比优势已复现，也不能说A/A+已在该原始合同下超越CIE。

## 2. 原始证据与论文本身的限定

- 2016：Zhong、Lambert、Lesselier、Chen，*A New Integral Equation Method to Solve Highly Nonlinear Inverse Scattering Problems*，IEEE TAP 64(5)，1788–1799，DOI [10.1109/TAP.2016.2535492](https://doi.org/10.1109/TAP.2016.2535492)。核对的是[作者上传的完整正文](https://www.researchgate.net/publication/296475773_A_New_Integral_Equation_Method_to_Solve_Highly_Nonlinear_Inverse_Scattering_Problems)，重点§IV及Fig.3–12；不是只读摘要。R8转录合同见下方索引。本文没有附作者可执行代码或同一噪声样本。
- 2019：Zhong、Xu，*Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*，J. Imaging 5(2),27，DOI [10.3390/jimaging5020027](https://doi.org/10.3390/jimaging5020027)。[开放完整正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)，§2.2–2.4、§3.1–3.3、Eq.(15)–(17)。本地可直接读 `research/research_fbg_r3_experiments/audit/literature_cie2019_body.txt` 和 `literature_cie2019.xml`（R9工作树内）。

2016的NIE就是2019所称CIE-I的直接前身；不能把两个名字当作两套互不相关的成熟算法。论文优势来自CIE方程与FFT-TSOM约束、分阶段β和模态数的联合条件，不是任意CIE权重、任意优化器或单纯几何轮廓的普遍优势。2019也明确显示过大的第二阶段current β会过早停。

## 3. 二维原合同

所有Austria条件共用：空气背景，二维标量电磁散射；2 m×2 m居中DoI；20个方位均匀覆盖整圆的平面波；40个均匀整圆接收点，接收圆半径3.75 m；正演MoM150×150，反演64×64；DPCS截断阶数L=15。两上方圆盘半径0.2 m，圆心(−0.3,0.6)、(0.3,0.6) m；下方圆环中心(0,−0.2) m，内外半径0.3/0.6 m。常规数据为10%加性高斯白噪声。初始化为空气材料、零APCS；后续阶段继承材料及APCS。sourcewise一步Polak–Ribière更新APCS，材料使用论文的修正对比度最小二乘更新。

| 条件 | εr／频率 | mf顺序 | β顺序 | 原文实际判断 |
|---|---|---|---|---|
| 最小二维正例 | 3.5／400 MHz | 5,7,10 | 6,2,0.5 | NIE找回Austria；相同mf的LSIE漏下方环；作者另试LSIE β型控制、CS-EB仍差 |
| 预定困难负例 | 4／400 MHz | 5,8,10 | 8,2,0.5 | 环内假目标；作者认为接近方法极限；无噪初跑也类似 |
| 可选电尺寸控制 | 2／750 MHz | 6,8,10 | 8,2,0.5 | 2 m域约5λ；NIE成功，LSIE只找上方两个盘 |
| 非高非线性速度控制 | 2／400 MHz | 5 | 6 | NIE与LSIE均可，NIE约10次而LSIE约50次；不宜代替高对比正例 |

二维APCS相对变化δ阈值为首轮0.01、后续0.001。εr=4的进一步噪声对照使用最后一轮更紧的1e−4，并比较10%、3%、无噪：噪声越低环内伪影越弱，不能概括成所有条件下完全失败，也不能称常规10%稳定成功。论文不提供R8统一的80步物理corrector合同。

原文还提到两个盘加有耗涂层矩形，但把细节转引2011 FFT-TSOM论文，本文没有完整尺寸／材料合同；**不选作最小对象**，不从示意图猜参数。Fresnel实测是另一数据合同，已有R8单独研究；不是本次经典合成对象的替代。

## 4. 三维原合同（不能缩放电尺寸）

共用合同：300 MHz，空气λ=1 m，3λ×3λ×3λ DoI。60个电偶极分布于半径3 m的三个圆，每圆20个等角位置：xy圆中心(0.2,0,−0.1)m、偶极z向；yz圆中心(0.1,0,−0.15)m、偶极x向；xz圆中心(−0.05,0.1,0)m、偶极y向。接收器与60源同位置，记录全部三分量，故每次完整数据为180接收分量×60源＝10800个复数。正演60³，反演30³，10%高斯白噪声。

初始化为背景材料与零APCS（§2.3）。L按噪声水平选择：保留高于噪声水平的主奇异值，**文章没有给三维固定整数L**；不能把R8的 `σ≥0.1σmax` 擅自称为作者唯一实现。APCS原始投影IFFT系数γ，每循环一次sourcewise PR更新；材料最小二乘。δ3D是每源系数相对变化平方均值的**平方根**，阈值1e−3，即 `sqrt(mean_l(||Δγ_l||²/||γ_old,l||²))<1e−3`。若内部用平方值，等价阈值应为1e−6；不能只给固定240循环后称论文收敛。

| 对象 | 几何／材料 | 原始stage合同 | 原文实际结论 |
|---|---|---|---|
| **最小三维高对比正例：四立方体** | 每个边长0.8λ；中心依次(−.7,.7,−.7)、(.7,−.7,−.7)、(−.7,−.7,.7)、(.7,.7,.7)λ；εr依次2、3、4、5 | mf=6,10；普通CIE current/material β均为6,1；有效relaxed第二阶段current β=3、material β=1 | LSIE只找底部两个低对比方块；CIE找齐四个。第二阶段current β=6提前停，不能把这个不良参数作为唯一Dual对照 |
| 复材料中等控制：涂层立方体 | 外边长2λ、εr=1.5+0.3i；内边长1λ、εr=2+0.8i | mf=6,10；β=6,1；有效relaxed第二阶段current β=6、material β=1 | LSIE/FFT-TSOM已能恢复；CIE质量相近而更快，非材料成功率优势 |
| 非最小：空心涂层立方体＋两外杆 | εr=2.5；正文未列Fig.11的所有尺寸 | mf6,10/β6,1留空腔伪影；mf4,8,10/β6,1,.5去伪影 | 需先从原图或前文献补全几何；不可直接把R8自定壳体尺寸称原合同 |

原文离散张量Green用引用[6]附录的有限差分积分方案；R8 Clausius–Mossotti＋radiation-reaction DDA不等于作者的同一格点离散。即使用30³/60³，也应称“原几何与采集合同比较，另用DDA离散”，直至实际核对其积分／自项定义。复材料A+尚无R8合法实现，因此四实材料立方体是当前A+最小公平对象；不能把实轴圆弧读出直接搬到复材料。

## 5. R8真实已做／未做

### 二维

R8合成数据确实用150²、原2 m域／原频率／20平面波／40接收点。脉冲使用等面积圆盘半径h/√π与解析圆盘自项，是方格积分近似；没有作者代码以证明完全一致。复高斯样本整体缩放至 `||noise||/||cleanY||=0.1`，不是原论文同一随机样本。四条件仍是同一公开Austria几何，不能算四种新盲测拓扑。

主九臂是32²、每stage最多300循环、随后最多80 accepted physical-corrector步。CIE/A/direct实材料逐像素自由χ∈[0,20]且有共同smooth+TV；独立paper_nie/lsie没有额外材料prior。后者是本地独立方程实现，不是作者成熟软件；两组结果也并非只差积分方程。各stage以原物理Phi选端点，另有common corrector，与论文只按APCS停止输出的路径不同。

| 原64²已完成内容 | 算法CPU秒 | wall秒 | 正演／反演及终止 |
|---|---:|---:|---|
| ε3.5 CIE single | 234.144440 | 251.853006 | 150²→64²；corrector到80步上限 |
| ε4 CIE single | 349.771393 | 445.180008 | 含复用已付aux CPU72.589199；重新启动L-BFGS history；80步上限 |
| ε4 Dual | 318.121192 | 441.251225 | 80步上限 |
| ε4 global A+ | 368.944895 | 432.489882 | 80步上限 |
| ε4 sourcewise A+ | 326.982497 | 377.664935 | 80步上限 |
| ε4 independent paper NIE | 316.105668 | 440.446420 | 无附加材料prior；80步上限 |

这些是对应RECEIPT的比较配方CPU，不应再与含同一已付aux的历史失败CPU重复相加。所有64²最终corrector都到上限；不称收敛。**ε3.5没有64²完整A/A+/Dual/LSIE强对照；ε4没有64²LSIE/direct完整配对。** 32² direct另延长1000步已做，不能把初期80步弱终点当最终强direct。原64²ε4的maxχ约19.4–20，GTχ=3；即使IoU不错也有严重幅度偏差，已有全尺度图。

R8已知真值场离线网格诊断：400 MHz ε3.5的clean 32→150误差约21.19%，64→150约9.96%，96→150约4.95%；150→200约3.71%。这支持区分离散误差与优化失败，不能把32²“没找回”直接解释为新机制失效。150²合成也不是连续真值。

### 三维

R8原系数PR240最终控制确实有三向场、60源、60三分量接收，LSIE/CIE/Dual/global/A/A+六臂。但域宽0.6λ，方块边0.16λ、中心按0.2缩放，源圆半径0.6λ、圆心也缩放；mf=1,2，初值χ=.25，材料χ∈[0,4]，固定ρ=10与quadratic空间prior，240循环＋最多300物理corrector，远非原mf6/10、空气初值、paper stopping。主6³inverse/8³synthetic；另8³/10³细化。先前exact-current及白化γ也是不同适配，最终原γ控制不能抹掉旧运行成本。

| R8三维已测项 | 实测算法CPU秒 | 适用范围 |
|---|---:|---|
| 0.6λ、6³原γ四方块 LSIE/CIE/Dual/global/A/A+ | 18.70 / 19.75 / 17.47 / 19.15 / 22.09 / 40.84 | 六臂relative-contrast L2均约0.46；A+没有显著优势且更贵；不外推3λ |
| 52个最终成功比较臂合计 | 1199.898052 | 仅在线成功臂；不含共享组装、合成、网格/Mie诊断、原3λ试算与失败 |
| 原3λ、12³inverse/16³synthetic direct截断试算 | 379.592774 CPU；395.867609 wall | 仅20初始步＋3个accepted corrector；不是完成原定30步，不是LSIE/CIE对照 |

原3λ试算峰值RSS356,040,704字节；实测12³/16³真值场偏差36.113%，终点数据残差27.632%，材料relative L2约0.244。算过而未合格，不是未测就声称资源阻塞。回执包含16³/12³共238009个FFT matvec、3160 primal RHS、1560 adjoint RHS、零failed RHS。它不能替代30³/60³原网格；也不能用它对原论文负面定论。

## 6. 后续最小公平实现范围（建议，NOT_RUN）

1. **二维150²→64² ε3.5**：优先补齐原论文意义上的 sourcewise LSIE/FFT-TSOM 与 NIE/CIE single；同mf、DPCS、初始化、噪声、原sourcewise权重和δ。再配有效Dual、ordinary global、A、global A+及sourcewise A+，以分开改权重、改current优化与材料读出的收益。不要用β=0重命名一个带额外prior的CIE黑箱就称成熟paper LSIE。
2. **三维原3λ四立方体、60³→30³**：先有全张量matrix-free LSIE/CIE源级PR和论文6/10 IFFT基、空气/零APCS初值，再配current β第二阶段3的有效Dual；A与A+共用相同物理、prior、stopping与current预算。ordinary global是global A/A+的必要优化控制。R8只实现小尺寸dense current路径，尚无可直接启动的原30³全套CIE runner。
3. 每条纸面复现基线先报告其**原stage输出／paper stopping**；若再统一加长物理corrector，单列before/after、实际成本与预算终点。论文材料域与A+实域若不能完全一致，应分别列“paper-like complex baseline”和“matched real-domain comparison”，不能把其差异隐藏成算法胜负。最小全实对象允许先固定所有臂同一实域，但必须将其标为论文条件适配。
4. 保留全部独立像素／体素材料自由度；不采用GT support、四块分区参数、已知几何拟合。公布统一全幅色轴、材料误差、残差、stage停因与费用；不能仅用轮廓IoU宣称高对比恢复。
5. 使用论文误差与R8误差的**独立标签**：2016 Eq.(16)与2019 Eq.(17)均是逐像素相对ε误差的RMS，后者同样有外层平方根。R8 relative-contrast L2、amplitude NRMSE、Otsu/fixed-threshold IoU另列，不混用阈值或分母。

原三维30³有81000个复vector-current分量；若直接dense复矩阵仅一个约105 GB，60³正演648000分量约6.7 TB。因此应复用已有matrix-free dyadic FFT，不能搬用R8小网格dense矩阵。这里只是存储算术，不是实际耗时预测。原论文“约25秒CPU/iteration、8 MATLAB workers、32GB”只是其机器与实现下的单迭代描述，不能折算为本机总CPU承诺；新阶段必须实测并保留失败成本。没有理由因预算自动改成0.6λ再叫原例。

## 7. 可直接核查的入口、配置与成本路径

以下R8路径以独立原工作树为准：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r8_worktree/research/research_max_struct_r8/`。R9同名拷贝可辅助读，但不替换原运行的来源声明。

| 文件（上述R8根目录下） | 用途／限制 |
|---|---|
| `delegated/paper_baselines/PAPER_CONTRACTS.json` | 两原论文逐case转录合同，含未完整指定项 |
| `PAPER_BENCHMARKS.md`、`PAPER_FOLLOWUPS.md` | 二维主表、64²/长corrector/网格偏差与限定 |
| `experiments/paper_benchmark.py` | build/solve/score主九臂；默认build非150，复用时须显式 `--forward 150 --inverse 64`；原64²九臂不曾因此自动完成 |
| `experiments/paper_original64.py` | 硬编码ε3.5、single一臂；不能当通用五臂runner |
| `experiments/paper_compare64.py` | 硬编码ε4、五臂 `--arm`；非任意case通用入口 |
| `experiments/paper_strong_corrector.py`、`paper_grid_diagnostic.py` | 既有长校正与离线网格诊断 |
| `src/paper_physics.py`、`src/paper_cie_fft.py`、`src/paper_nie.py` | 原采集/等面积积分、通用CIE、独立real NIE/LSIE；不能声称作者代码 |
| `delegated/paper_baselines/source_snapshot_final64_fast_LU/` | 最终64² GMRES≤1000→dense-LU工程修复快照；旧stagnation成本留存 |
| `RESULTS/paper_benchmarks/*/DATA_RECEIPT.json` | 数据生成原150²与各case成本；例如ε3.5总生成CPU2.172195秒、wall2.200880秒 |
| `RESULTS/paper_benchmarks/original64/*/*/RECEIPT.json` | 上表各64²实际CPU、wall、RHS、失败拒绝与stop |
| `RESULTS/paper_benchmarks/COST_LEDGER.json` | 已知CPU下界3566.388676秒；首个失败CPU未知；含长校正/失败的细项；不要加总比较配方费与复用费 |
| `VECTOR3D_REPORT.md`、`delegated/vector3d/FINAL_SUMMARY_ZH.md` | 全矢量适配与负证据总述 |
| `src/vector_cie3d.py`、`src/vector_pr3d.py` | 小dense-DDA全矢量、实域A+；prepare仅mf1/2、源位置硬编码缩放，不能只传width=3就复现原合同 |
| `experiments/vector_pr240_fourier.py` | 最终原γ六臂；**内部末段直接读已有truth评分**，不能原样用作新blind runner，应分离seal与score |
| `RESULTS/vector3d/FROZEN_CONFIG.json`、`pr240_fourier_coordinates/PR_FROZEN_CONFIG.json` | 主适配配置与PR240配置 |
| `RESULTS/vector3d/pr240_fourier_coordinates/PR_SCORES.json` | 原γ六臂实际CPU及材料误差；无需重算 |
| `RESULTS/vector3d/ONLINE_COST_LEDGER.json` | 52成功在线臂1199.898052秒及LU/RHS等，非完整grand total |
| `experiments/vector_original_size.py`、`src/vector_fft3d.py` | 原电尺寸FFT direct已有入口；非原30³CIE/LSIE入口 |
| `RESULTS/vector3d/original_size/FROZEN_CONFIG.json`、`COSTS.json`、`TRUNCATION.json` | 原3λ12³/16³实测、截断及物理不合格证据 |

R8 `paper_benchmark.py` 的OUT字符串写为小写 `results/paper_benchmarks`，实际历史文件常以 `RESULTS` 展示；当前Mac文件系统可解析，跨大小写敏感系统复现应固定目录映射，不能因路径不匹配默默重跑或覆盖。所有以上程序均是证据入口；本审计没有执行它们。

## 8. 原文公式抽取更正

初版审计误把2019 Eq.(16)、(17)均记为无外层根号，原因是本地正文纯文本抽取丢掉MathML `msqrt`。Root读原始PMC公式后指出；现已直接核查本地原始XML，两式均含 `msqrt`，以上§4、§6已修正。R8冻结材料不回写；R9数据包中的R8原始合同拷贝保留原样，并由 `data_public/classic_extension/PAPER_CONTRACT_CORRECTIONS.json` 显式覆盖错误的三维阈值元数据。此次错误没有进入R9三维反演；已启动的R9二维使用正确的二维δ平方阈值，未受影响。R8固定循环PR适配不应重新叙述为原论文δ停机。
