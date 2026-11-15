# Pro ZIP 完整清点与证据提取（机械审计）

2026-11-15。输入 `/Users/migodam/Downloads/MAX_STRUCT_research_bundle.zip`；仅解压一次至本目录 `extracted/MAX_STRUCT_research_bundle/`。共112个文件，62个文本文件均已读取/解析；全部2000条主分数行均参与字段与状态清点。PNG未独立目视；NPZ仅列容器成员名，未读数组。未运行包内脚本、未访问sealed_truth、未核SHA、未请求外部模型、未改历史文件。这里不判断证明有效性或新颖性。

## 1. 实际交付是什么

包内含8个public输入、4对象×A+/DR的8个diagnostic raw、0002 A+静态度量子集、2000行历史分数、19张原图、7篇研究文档、10个脚本、冻结physics.py，以及回放/24²几何/8²非线性/holonomy/五例离线真值正演的保存回执。不是一个新的48对象反演结果包。物理后端是scalar-TM disk-pulse VIE；没有新增3D或实测。

主分数筛选 `initialization=uniform,view=primary,solver!=direct` 得288行、每法48对象；逐行重算历史NRMSE≤.15且IoU≥.7，与METHOD_SUMMARY一致：CIE30、Dual32、G1930、DR29、MM36、Arc36。六法并集39；所有2000行run_status=COMPLETE。原图A+列仅MM，无Arc。

后验公共选择回放39/48，Arc→球外才G19→仍球外才DR，每次保留较小Phi。G19触发8例、DR7例，均已存独立uniform子运行。平均历史CPU相加145.55078125秒，Arc113.22265625秒；不是同120秒预算、热重启或控制器实测。无条件三法minPhi为38/48、316.80208333秒。

## 2. M0/M1 的确切代码定义

文档的M0/M1是拟议2×2干预的概念标签；包内实际只固定current为CIE metric，再切换material readout，尚未实现四组合实验。代码对应 `scripts/metric_geometry_pilot.py::geometry,onehop_matrix` 和 `scripts/nonlinear_map_validation.py::onehop_metric,exact_current,block_residual_jacobian`。

令 `L=I-diag(m)G`，本构缺陷 `r=J-diag(m)(E+GJ)`。`beta_c=1`、`beta_m=.2`。

- `p_c,i=beta_c/(1+beta_c*m_i)`；空间current metric `M0=diag(p_c²)`。逐源 `ws_s=1/sum_i|E_is|²`，`wd_s=1/sum_r|Y_rs|²`；current目标为半个 `sum_s wd_s||SJ_s-Y_s||² + ws_s r_s* M0 r_s`（乘整体常数无关）。24²的ws/wd来自原public观测，并未按新F(anchor)归一；8²来自各人工物理锚点。
- `d_i=1+(beta_m-G_ii)m_i`，`a_i=m_i/d_i`，`p_i=beta_m/d_i`。以不重叠2×2 cell构成四cell block，`Goff,b=G_bb-diag(diag(G_bb))`；`H_b=(I+diag(a_b)Goff,b)diag(p_b)`，`M1=H*H`。**不是全域Goff，也不是更换current物理后端。**
- 第三种self_diagonal为`diag(|p|²)`；无prior精确根处空间逐row正权重在共享实材料读出中抵消，故其局部映射与CIE一致到数值精度。
- M0、M1一般都依赖m。非线性A+材料最小二乘代码包含H(m)的完整局部导数；反向错配或M1/M1共同目标版本尚无实现/回执。文档要求共同目标随m变化时保留metric导数，不能把冻结metric迭代称为完整共同势算法。

理论提取：在实化空间，`K a=L U a`、`C h=diag(h)Q`、`D a=W S U a`，`Ac=K^T Mc K+D^T D`、`Bm=C^T Mm C`，`T=Bm^-1 C^T Mm K Ac^-1 K^T Mc C`。假设有限维可逆L、实共享材料、SPD度量、可逆Ac/Bm、精确无约束current及局部精确material、零prior、精确物理根、光滑内点或人为固定face。24²只在`.02+1e-8<m<20-1e-8`的自由坐标展开；不是完整box切锥，不是实际finite PR/MM/Arc程序的Jacobian。

Holonomy定义固定M（实测M=I）、固定rank(S)、`N=ker(S)`正交基、`J=J0+Na`、`K=LN`、`H=K*MK`，速度`a_h=H^-1K*MD_hQ`；曲率`H^-1[K_h*M(I-P)C_k-K_k*M(I-P)C_h]`。材料相关M需额外导数。文档明确精确current profiling清除闭环辅助记忆，未声称物理解空间增加。

## 3. 已存数值证据与成本

- 24²八个raw锚点×三readout共24行；另四个离线真值锚点×三共12行。均重新生成无噪声F(anchor)。CPU分别11.877513545、10.05558269秒，errors均空。A+ raw谱半径范围1.000000000007–1.000000145254；同Bc范数1.000000020416–1.000975401785。真值锚点A+谱1.000000019245–1.000000401644，同Bc范数1.000484881829–1.000849141377。0005 DR raw被说明为uniform .1初始候选，并非坏收敛终点。
- 单0002 A+代数回执：441自由材料；current normal residual1.9123e-15、直接块解/消元相差1.19–1.23e-14、native block动作误差1.4220e-16。这是保存回执，本次未重跑。
- 8²两个人工锚点×两readout×三epsilon=12行；高壳A+谱1.002166326544771，另±.001扰动4行。非线性+扩张CPU合计.201607568秒；12 LU、72 physical RHS、192 receiver-transpose RHS、48 correction RHS、24 analytic port Cholesky、1536 analytic port RHS、48 current port Cholesky/RHS、128四cell材料解。不含进程启动/导入，无网格收敛。
- Holonomy一个8²壳锚点、h=e27/k=e28、三epsilon正反矩形、每边RK4四步。曲率范数.002466510353；memory为2.47299e-7/2.22161e-8/2.46716e-9；receiver变化≤1.51e-23；同路径往返4.899e-19。CPU .170976205秒、418个48阶Cholesky、2496运输RHS，另18曲率RHS、128 span RHS、576×2016 SVD。数值rank576（相对阈值1e-6），stable rank8.02185。未做恢复。
- 五例离线真值正演保存表：0007/0016/0029/0037/0044，残差δ比分别.9763188993/.9341488167/.9682356738/.9850286640/.9237769000；5 LU/54 RHS、总wall .509907秒，所有matched-grid/in-bounds/inside=True。本次只读表，不读真值。其余共同失败0001/0004/0014/0034六法均球内；后两者跨网格，无纯噪声统计证书资格。

## 4. 未完成、失败与来源差异

无显式NOT_RUN字段；主分数2000行均COMPLETE，两套几何errors=[]。但以下仅建议，无实现或新结果：完整M0/M1四组合、实际PR/MM/Arc全状态Jacobian/干预、冻结规则的独立确认、新对象恢复、holonomy kick与同费用random nullspace等控制、网格收敛、3D/实测。holonomy建议的S=0平坦和exact-profile消记忆属于文中论证；脚本未执行这两个数值控制，也没有RK4步数收敛扫描（只三epsilon与往返）。不要把文档的全部建议计为已跑。

来源/自包含性差异：`RESEARCH_VERDICT_ZH.md`引用`literature.md`，FP-T文档引用`theory_sources/`，两者均不在包中；历史长链公式/原runner/完整轨迹没有完整收录，需父线程回原仓库核对。`docs/AUDIT_R13_ZH.md`“8个诊断对象”指原访问范围，包内仅4对象的raw子集；README已正确限定。旧`analysis/...`、`/workspace/...`是原会话路径，非包内可直接访问目录。原主图与dynamics8是不同运行，原图无Arc，不能串接为同一轨迹。

`CORRECTIONS_AND_SCOPE`说24²A+同Bc扩张约1.0002–1.001是粗略概括；实际raw最小值1.0000000204（0005 DR均匀锚点），应采用CSV准确范围。README/脚本大多记录线程控制，24²几何和代数实际threadpool_limits=2，8²为1；CPU来自Linux x86_64，不能直接当本机Apple Silicon性能。

可选truth脚本含SHA计算，当前任务明确禁止SHA，故如后续要用需由父线程遵守用户偏好先处理；本次未执行。打包PORTABILITY_CHECK自报PASS：18 CSV MATCH、4 JSON数值字段MATCH、默认4步骤returncode=0；这属于来源回执，非本次独立复现。

## 5. 源链接、文献与入口

冻结仓库commit `959c165e77935d49ad95152e0339fb21d9c416bc`；来源manifest提供固定commit URLs、release资产、原图main URLs。本次未联网验证可达，不能声称包本身已发布。所有提取URL见SOURCE_URLS.txt。

- 仓库：https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/959c165e77935d49ad95152e0339fb21d9c416bc
- 主分数固定源：https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/959c165e77935d49ad95152e0339fb21d9c416bc/research/research_max_struct_r13/supplement_aplus48/analysis/stack_comparison/SCOREBOARD_COMPACT.csv
- 本包没有新研究脚本/文档的发布URL，`sources/PROVENANCE_ZH.md`明确新脚本来自该轮会话，不能误标成历史公开算法。
- 引用包括van Leeuwen/Yang2025（10.1088/1361-6420/adab89）、Balduzzi2018、Lewis/Malick2008（10.1287/moor.1070.0291）、De Luca/Oriolo1997（10.1109/70.631239）、Shamir/Yomdin1988（10.1109/9.14412，仅摘要/元数据核验自报）、Sun等CC-CSI2017（10.1109/TAP.2017.2673758）。原PDF未打包，文献内容/新颖性须另核。

运行入口`run_reproductions.py`默认顺序selector→nonlinear→expansion→holonomy；`--with-geometry`加geometry→common metric→algebra。默认写reproduced，零外部下载；可选`--offline-archive`明确读取另供truth archive。依赖钉版本NumPy2.3.5/SciPy1.17.0/pandas2.2.3/threadpoolctl3.6.0；原环境Python3.12.14、Linux x86_64。FILE_MAP.csv含全部文件与读取范围，CLAIM_EVIDENCE.csv为声明—证据矩阵，MECHANICAL_CHECKS.json保存本次无实验计算的表格核对。
