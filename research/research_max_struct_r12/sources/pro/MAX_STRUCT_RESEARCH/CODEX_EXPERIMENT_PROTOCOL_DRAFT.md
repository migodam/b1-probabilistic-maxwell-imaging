# MAX-STRUCT 下一轮 Codex 实验合同草案

## 0. 任务与证据身份

本文件是下一轮可执行实验设计；**本次起草没有运行任何新增实验**。目标是检验新的 current–total-field 集合反射结构是否改变已有错误终点的可达性，同时明确它解决不了哪些材料选择问题。不能把这轮改成继续扫β、Green rank或高阶局部展开。

本轮交付应形成：可复现的实现、一组能否定具体机制的对照、24个独立二维对象上的冻结比较，以及少量资格明确的经典/实测/矢量3D锚点。不得以“某门槛没过”为理由只返回停止日志；失败后应解释哪个假设失效、完成必要对照，并形成有根据的保留、简化或转向判断。

本合同所指“严格”是指相同信息、明确变量与目标、完整费用及证据身份，不是要求所有对象都胜出。

### 已核对的历史来源

仓库：`migodam/b1-probabilistic-maxwell-imaging`。本报告审计快照为`24ab7e85d750a477ab337e03cdb1a77d5c9ff2c3`。实施前记录实际checkout commit；若更新，先列与下列来源相关的改变，不自动重跑全部历史。

- `research/research_max_struct_r4/ABC_VALIDATION_ZH.md`
- `research/research_max_struct_r4/THEORY_DERIVATIONS.md`
- `research/research_max_struct_r5/THEORY_DERIVATIONS.md`
- `research/research_max_struct_r6/RESEARCH_REPORT_ZH.md`
- `research/research_max_struct_r6/NIE_BASELINE_AUDIT.md`
- `research/research_max_struct_r7/THEORY_DERIVATIONS.md`
- `research/research_max_struct_r7/FAILURES.md`
- `research/research_max_struct_r8/CHAT_HANDOFF.md`
- `research/research_max_struct_r8/NEURAL_WARM_STARTS.md`
- `research/research_max_struct_r8/PAPER_BENCHMARKS.md`
- `research/research_max_struct_r8/FRESNEL_REPORT.md`
- `research/research_max_struct_r8/VECTOR3D_REPORT.md`
- `research/research_max_struct_r11/CHAT_HANDOFF.md`
- `research/research_max_struct_r11/RESEARCH_REPORT_ZH.md`
- `research/research_max_struct_r11/THEORY_DERIVATIONS.md`
- `research/research_max_struct_r11/FAILURES.md`
- `research/research_max_struct_r11/METHODS_AND_ABLATIONS.md`
- `research/research_max_struct_r11/CONFIRMATION48_REPORT_ZH.md`
- `research/research_max_struct_r11/CROSSGRID8_REPORT_ZH.md`
- `research/research_max_struct_r11/VECTOR_CONFIRMATION_REPORT_ZH.md`
- `research/maxwell_formula_atlas/06_DERIVED_OPTIMIZATIONS.md`
- `research/maxwell_formula_atlas/12_BILATERAL_MATERIAL_MODES_AND_FIXED_POINTS.md`
- `research/p3l_r12_bilateral_fixedpoint/TRANSVERSE_CONTINUATION.md`
- `research/research_fbg_r3_experiments/input/FBG_R3_Reproducibility/memos/agent_branches.md`

从这些确知报告/复现入口定位真实数组和checkpoint，不凭命名猜测R11端点文件路径。历史GT已经公开，只能作为暴露开发/机制证据；不能重新称盲测。

## 1. 不得重复作为本轮创新的内容

下列已有结论直接继承：完整resolvent与有限support Woodbury响应；R4纯端口full-current profile退化为原数据目标缩放；finite η 的WRI/profile metric；加权材料读出差式；material-compatible current tangent投影；source–receiver双边Jacobian/Hadamard Gram；R7 rank2/8/32有限动作负结果；P3L-R12固定/刷新弱坐标的横向corrector负结果。

不再安排“恒等式误差达到机器精度”的完整新实验。新代码只针对实际实现风险做少量单元资格，例如共轭/梯度符号、数据mask、共享实材料投影和矢量分量接线。已验证的旧核不反复验。

R11必须同时保留两类证据：16及配对8中的G19真实环/壳增量；48中Global/G19同一23成功集合与额外费用。三维matched多数方法形态可见但无稳定G19优势；cross-grid全部噪声球外且干净模型差异显著，不能全归给优化器。

## 2. 新候选的唯一主版本

令J为contrast current，Q为total field。取公开、固定的正尺度s，X=J/s，B=sG，T=sS。采用product norm ||ΔX||²+||ΔQ||²，并定义

\[
\mathcal A_\delta=\{(X,Q):Q=E^{inc}+BX,\ \|\mathcal W(TX-Y)\|\le\delta\},
\]
\[
\mathcal C_{box}=\{(X,Q):X_{i,:}=\alpha_iQ_{i,:},\ \alpha_i=m_i/s,\ m\in\mathcal B\}.
\]

\(\mathcal W\)必须是该数据合同真正使用的噪声白化及观测mask。若各源有不同缺失通道，按各源实际算子实现，不能沿用完整矩形数据的错误共用port Gram。相关跨源噪声需整体实化白化，或明确说明本轮只支持独立噪声。

主更新冻结为

\[
c=P_{\mathcal C}z,\quad a=P_{\mathcal A}(2c-z),\quad z^+=z+a-c.
\]

材料候选来自C shadow。A先C后是另一种变体，不能在运行中静默切换。普通AP采用相同两个投影、相同metric、同一初始状态，只去掉反射记忆。

每cell投影共享所有照明；三维各向同性时还必须共享三极化。固定实α后Q=(Q0+αX0)/(1+α²)，距离为(aa−2ccα+bbα²)/(1+α²)，比较合法二次驻点与box端点。稳定二次求根/2×2 eigensystem可作为工程选择，但不要重新引入复杂非线性内层。

**保持原估计问题的方式：** DR只是候选生成器。主比较最终仍按原完整Maxwell数据目标加原材料prior \(\Phi\) 选择，接同一个原物理corrector。C投影不含空间prior这一事实须公开；不能把其feasibility objective称为原MAP。要研究新的prior或直接noise-constrained估计器，放入第5节的独立信息/选择实验，不混入主算法净效应。

s及投影容差只允许在六个独立开发对象/公开均匀校准状态上冻结。默认用一个公开材料尺度，不能按测试对象真值幅度设置s。最多保留一个主版本；一个有限尺度敏感性实验只在历史机制对象上做，不展开成测试集扫参。

## 3. 预算、共同校正与选择规则

### 3.1 每种费用都要有准确名字

按矩阵阶数分别记录：

1. 原Maxwell factorization、forward RHS、adjoint RHS、tangent RHS；
2. 固定投影系统H的Cholesky/QR或preconditioner setup、H-RHS、Krylov iterations、G/G*及S/S* actions；
3. port eigendecomposition、μ求根、C投影、材料内层求值；
4. 物理审计、失败/拒绝候选、共同corrector、NN inference与训练/验证；
5. 硬件、精度、BLAS线程、进程CPU秒、wall span、峰值RAM/VRAM。

不能把DR中只有少量“原Maxwell LU”解释成全部线性求解很少。固定H反复求解虽然不需新的factorization，也不是免费。设置缓存可复用，但应同时给cold case、warm reuse及明确的摊销对象数。

共同主要资源坐标采用相同host/线程设置下的实际CPU秒；同时报告以上物理费用，不能把不同矩阵阶数的RHS简单相加当等flops。GPU实验各方法必须在同一目标上用同一host；CPU秒/GPU wall分别报告，不拼为单一加速比。

### 3.2 有限但可执行的默认预算

二维机制试验与24主测试默认每对象每臂给同一总预算B，建议在实施时将B冻结为同host上的240进程CPU秒；这只是计算额度，不是成功标准。若后端规模使该数不合适，可依据**公开校准任务**和已有R11成本统一换算一次，并记录理由。不能按测试表现逐对象改B。

默认将B的60%用于候选生成/原方法推进/必要审计，40%留给一个共同原物理corrector。终止早的普通方法如实保留实际费用，不靠重复相同已收敛步“耗尽预算”；另设真正使用剩余预算的direct多起点控制，避免把局部方法早停与全局方法新增算量混淆。

每臂在候选期只按已付费的原\(\Phi\)选一个合法incumbent（含初始点），之后用同一corrector目标、box、prior、梯度及停止准则。若增加多个候选corrector，必须给所有相关控制相同候选机会并从B扣除，不能只给新方法额外分支。

在公开成本前缀B/4、B/2、B报告incumbent，选择只能使用该前缀之前已计算的信息。最终全终点与共同前缀是两种结果，不能混用。只在所有对比臂都有对应完整物理审计的已付状态上配对；不靠插值虚构材料图。

physical audit频率依据固定投影步数/固定计算量冻结，全部收费。其目的是真物理候选选择，不是自动结束研究。算法迭代、投影误差、物理残差和科学成功是四个不同概念。

继承R11停止审计教训：同时记录未缩放原\(\Phi\)梯度、优化器实际缩放目标及其停止统计。box投影与目标尺度不交换，不能把有缩放目标的投影梯度容差简单除以常数，就宣称获得同一未缩放KKT证书。活动乘子符号、free-face曲率和预算终止分别报告。

## 4. M1：从R11同一个错误终点开始的决定性机制试验

### 对象与起点

只用三个公开历史终点：R11环/CIE坏终点、壳/G19坏终点、桥/G19坏终点。原报告指出前两者活动乘子与free-face Hessian符合近约束局部极小；桥还有负曲率/错号乘子，属于重要反控制。按报告定位和保存哈希，不再通过GT挑新“漂亮”锚点。

各臂从完全相同m_bad开始，统一完整物理解Jphys(m_bad),Qphys(m_bad)，不得给新方法更好的初值。R11其它方法已付较低\(\Phi\)终点可以作已知可达reference，但**只用于离线路径解释，不能输入DR提案**。这会直接区分“从均匀初值避开错误区域”与“从已进入的错误终点离开”。

### 四个必要臂

|臂|保留/改变|回答的问题|
|---|---|---|
|原物理继续优化|同m、同原目标、同corrector|是否其实只是原求解尚未完成？|
|AP|与DR完全相同A/C投影和metric，没有反射|新收益是否只来自换了joint投影/噪声球？|
|DR|唯一主反射版本|反射记忆是否改变错误端点后的可达区域？|
|同成本扰动＋原物理多起点|公开固定的低频/空间相关材料扰动，box合法；不读目标形状|收益是否只是额外搜索预算或任意restart？|

扰动控制的尺度与数量在公开开发对象冻结。若比较current-space扰动，应使用同一product norm和同阶费用；不能让随机控制只动几个无效坐标而DR用全空间。不要再把R7已否定的相同rank32射线包装成新的核心臂。

### 必须记录的机制量

保存完整DR状态z、两shadow及活动材料；AP同存两投影点。至少每次物理审计保存：原\(\Phi\)、数据项、prior、噪声球资格、box active set、相对初始材料距离、fixed-point gap、A可行误差、C可行误差，以及真实状态缺陷。

令J_aux来自当前shadow，m来自同一候选，已付物理求解给Jphys(m)。无需新resolvent solve即可计算

\[
d_{vis}=\mathcal W S[J_{aux}-J_{phys}(m)],
\quad r_{phys}=r_{aux}-d_{vis}.
\]

记录二者**绝对norm和相位交叉项**，不能只看能量百分比。对活动box上的材料步，报告真实物理梯度与可行方向关系；梯度接近零时不要用不稳定cosine充当方向证据。

### 结果怎样否定机制

- direct续跑就到同样好终点：该对象反驳“必须新basin机制”，更可能是停止/预算。
- AP与DR同样改善：joint projection/放弃硬DPCS已足够，不能声称反射记忆独占。
- 同成本扰动＋direct同样改善：属于通用restart收益，需承认Maxwell特定结构增量未成立。
- DR原\(\Phi\)明显下降但材料恶化：可达性可能改善，材料选择没有解决，转M2解释。
- DR gap不降、shadow漂移/周期样行为：准确固定点定理未到适用状态；保存状态分析原因，不把理论命题判错，也不称已发生混沌。
- DR只改善桥，不改善近KKT环/壳：证据支持加速未完成优化，尚不支持逃离错误局部极小。
- DR仅在过松投影容差下“成功”：先查真实物理误差；如果收益依赖数值不准确，应把它重新定义成近似proposal并给同近似控制，而不是宣布精确DR成功。

M1失败后不无限增加rank/步数。应根据对应反例选择一次有明确数学理由的修正（例如稳定投影、合理metric、约束处理），保留v1结果；修改版只能作为新的暴露开发版本，后面的未见数据继续独立。

## 5. M2：球内错误材料的no-go与材料选择/新增观测分流

从R11已公开报告中选择两个“原物理noise-ball合格但材料错误”的端点，另选一个同数据的不同合格材料作为候选对。用途是定义机制边界，不估计总体失败率。

若从精确物理状态启动且数据已在球内，z∈Aδ∩Cbox，**准确AP/DR都是固定点**。此结论直接由集合定义得到，不再做海量机器精度检查。它说明feasibility反射不可能自动替换一个已经可行的错误材料。

同时必须避免另一个错误推论：两个候选都在一次观测的CN99球内，不代表它们统计不可区分。计算候选干净响应的白化pair distance，并在清楚的实Gaussian白化约定下报告KL=.5||Δμ_white||²。若含模型/标定nuisance，应先说明同一nuisance已知、拟合或边缘化；不能混用。

只做两个分离实验：

**M2-P：改变材料选择信息。** 选一种在开发/训练数据上预先固定的材料prior或learned proposal；同一prior给DR、direct及必要CIE参考。对同一已生成候选池，比较旧\(\Phi\)选择、新prior选择和仅离线oracle容量。若oracle池内已有好材料而online选错，瓶颈是选择；若oracle也没有，瓶颈仍是proposal覆盖。不要用GT设计环洞/桥位置规则。

**M2-O：增加一组真实不同照明。** 基于当前候选对、现有物理和公开候选源集合，选择预测响应差最大的少量新照明，并与等测量成本的预定均匀/随机照明控制比较。所有算法取得相同新增数据；新源不能只是原照明乘非零常数。对合成任务由独立生成器给新观测；实际测量不可取得时只做设计预测，标NOT_ACQUIRED。

新增照明的收益归给新增信息，不能归成同数据优化器增益。若新观测仍区分不了候选，报告模型允许的材料/任务不确定性；有限候选集合不能自动称posterior，也不能宣称已证明全局多解。

## 6. M3：NN初值信息是否在接口中被保留

继承R8已训练模型/权重及训练来源，先不重新训练网络。选4个公开开发对象，其中包含R8已披露的0027类型反例：NN直接物理好，而经辅助管线后材料变差。确切数据文件必须由R8复现入口取得，不猜文件名。

共同NN材料mNN、共同观测与原目标，比较：

1. NN→direct物理；
2. NN→原Dual-CIE接口；
3. NN→已披露R8固定物理current偏置接口（作为旧修复控制，不称新发明）；
4. NN→DR，初始z用完整Jphys(mNN),Qphys(mNN)，避免先强行替换可见current分量。

只有在主假设需要时，增加一个soft/noise-consistent DPCS接口；该方法必须明确其改变的是current约束，不是暗中新增prior。无需把sourcewise/global A+及全部反馈变体都在4例上重复。

分开记录接口引入的hard-DPCS误差、Fourier/APCS截断误差、真实Maxwell defect \((I-D_mG)\delta J\)、receiver-visible defect以及第一次材料读出后的m变化。R8已有“初始defect降低仍未保证最终材料更好”的反例，所以“接口误差较小”只是机制变量，不能作为成功判据。

offline评分可分可见材料切向/弱切向，但必须报告误差绝对值，不用能量比例掩盖强方向变差。完整current tangent巨型QR不是本实验前置要求；优先用已有观测Jacobian及缓存action完成必要诊断。

若DR保留NN状态但与NN直接同结果/更慢，结论是接口正确性或稳健性，不是新成像能力。若反射反而破坏好初值，应记录该反例，使用公开物理incumbent保护；这不能靠GT回退。

## 7. 独立二维确认：6开发＋24锁定对象

### 7.1 数据结构

六个新开发对象与24测试对象完全分离。24测试采用4个预先冻结的非手工家族×6个独立随机几何，不再以少数手绘环/开缝壳代表泛化。建议家族为：

- 随机level-set/Boolean组合的多连通材料；
- 随机分支网络及窄连接结构；
- 随机多岛、多尺度混合相；
- 平滑分级/非二值随机场。

家族生成代码及分布事前冻结。每族6例覆盖公开高/中contrast及mixed材料，但几何独立时，不把不同contrast的横截面平均差解释成纯contrast因果。空目标/极小目标的处理规则事先规定，不根据方法难易重抽。

保留完整576独立材料自由度，不能把生成参数、target mask、家族标签或已知support交给求解器。求解器只得公开采集、原先验/box、Y、σ/协方差、mask和匿名ID。

24可拆两个12对象执行批次以便恢复/调度，**两个批次之间不看GT调参**。若资源只完成12例，应明确为pilot，不宣称论文级总体优势；完成哪一批由预定顺序决定，不看成功率。

### 7.2 冻结的主方法矩阵

主测试最多8臂，不带全部开发消融：

|主臂|必要性|
|---|---|
|CIE single|原文family最直接对照，保留真实single优势情形|
|Dual-CIE|成熟强基线，不能弱化停止或budget|
|Global-Dual|隔离global current搜索规则效应|
|A+|已有低成本nonlocal material readout|
|Full G19|已有条件成功的非局部matrix分母|
|冻结NN＋direct|R8已有最强独立初始化管线|
|同预算direct多起点|控制新增搜索算量，首入口包含均匀直接优化|
|DR主版本|本次唯一新主结构|

AP、freezeA、contractedT、BothBlocks、多个metric/scale及NN接口消融不进入24×全组合；它们只在机制对象或已存在历史结果中使用。

主物理对比明确标记初始化/信息差异。NN带训练prior，不能把其与无NN的方法比较说成纯优化器净效应；M1的同端点比较提供这一归因。NN模型不支持新采集mask时不得随意补零伪装同域推理。

若预算需要缩为6臂，事前将CIE single、A+仅保留在预指定外部/机制锚点，主比较保留Dual、Global、G19、NN-direct、direct多起点、DR；不能看24结果后删掉最强参考。

### 7.3 少量配对stress，不做全因子扫描

- 在24中事前按ID选4材料做fine生成→coarse反演，保留matched对。比较DR、Dual、G19三臂即可。
- 另事前选少量配对对象分别只改变一个因素：噪声、孔径或源数。每个因素2个对象足够作机制stress，不估计总体鲁棒性。
- initialization敏感性主要由M1/M3回答，不为每个24对象再乘十个random seed。

所有配对条件共享底层材料，不能计作新的独立材料。cross-grid目标用明确的粗网格平均/体积分数映射，干净预测间模型差异独立列出；名义noise-ball不包含该模型差异，不自动给“正确/错误”的同一阈值含义。

### 7.4 评分与统计

同时报告绝对χ-RMSE、relative material L2、继承的amplitude-normalized NRMSE、IoU/拓扑、分相/背景误差、真实data residual、noise-ball、prior及成本。低相误差不能被高相归一化掩盖。沿用历史阈值时注明定义；不要把toy的relative L2直接套成R11的NRMSE。

主要比较预先固定为DR vs Dual、DR vs同预算direct多起点；NN-direct/G19等用于强实践基线和适用边界。报告paired连续差、每例散点与全图、成功集合overlap，以及按材料簇bootstrap的不确定区间。24例只能支持有限强度统计，不把CI未跨零当万能gate，不把未显著解释成完全等效。

**任何低残差坏图都保留。** 若均值改善仅来自几例，则明确展示beneficial conditions；如果所有方法同一成功集合，不能用更小loss宣称新增恢复能力。

## 8. 经典原论文与Austria锚点

复用R8已有原文核对/实现，先检查实际输入合同，而不是从零写一个名叫NIE的弱替身。

R8已披露Austria原尺度：150²生成/64²反演、2m ROI、20平面波/40接收、400MHz、10%复噪声、εr=4。实施时与原论文实际图号、材料参数、noise定义和采集位置核对并记录适配差异；χ=εr−1，不能把εr=4写成χ=4。

最低必要比较为独立NIE实现、原CIE/Dual参考和DR；只在另一个原论文高contrast对象上加G19/A+作机制参考，无需全部8臂×所有paper cases。NIE已有实现不等于作者认证复现；若缺原代码或参数，仍可依据论文做披露重实现，不能冒称literal reproduction。

原论文3D设置按真实矢量模型处理，不能把二维TM堆叠或缩小若干倍的样本标为原论文已复现。scale-down可作为独立有限规模实验，标题和表格必须明示。

## 9. Fresnel真实测量：先量化校准，再比较求解器

继承R8真实2001/2005数据、空场多极照明、每源mask/adjoint、独立圆柱参考及FoamTwinDiel相关合同。原χ≥.02曾排除空气，R8已纠正为允许χ=0；不要回退到旧错误box。

只选两个在现有照明/校准范围内可解释的公开对象，主比较Dual、DR、一个direct/NN-direct参考。照明、校准nuisance、prior和缺失通道处理全部相同，不能让新方法额外拟合每源任意gain。

事前保留整组未用于拟合的照明/接收或频率作预测检验；使用6+8GHz共同反演是新增信息实验，与只用6GHz的same-data比较分开。若真实资料没有可信pixel GT，禁止输出假NRMSE；报告文档支持的尺寸/区域物性误差、heldout complex residual、校准残差、box触边与稳定性。

若所有方法同终点且heldout差，优先判断illumination/model/regularization共同风险，不能仅加迭代。若新方法训练/拟合校准后改善，需要同样的校准给基线；只在新方法中引入额外自由度的结果独立命名。

## 10. 三维矢量Maxwell：少量合格对照

沿用R11三分量dyadic后端、真实外部dipole、共享实各向同性材料，不重新把scalar kernel当vector。R11已到1728材料DOF，50臂均预算终止、cross clean差17.68%/3.53%；这些是当前边界，不用再跑同一50臂复述。

新投影实现的资格只针对新的代码连接：三极化/多照明共用m、复杂共轭、实际mask及scaled metric、原方程真实审计。完整dense H setup的时间/内存必须计入；若用FFT/Krylov近似投影，保留原dyadic自项、非周期边界及投影误差，不能继续套“准确投影固定点”结论而不注明误差。

独立物理控制复用已有Mie球/部分体素refinement，补一组对结论有意义的独立离散比较即可，不为了新表格重复成熟machine-precision恒等式。若该几何/contrast下模型误差仍大，继续完成有限DDA范围内的机制研究，同时明确连续成像主张尚不成立；不因为资格不足就把整个研究机械终止。

最小新反演是2个独立非手绘3D设计，每个matched与cross两条件，主比较Dual、DR；只在一个预指定对象加Global/fullGreen控制。全部各向同性材料参数和source/receiver配置公开，算法不接收GT support。该规模用于判断可实现性与failure mode，不统计三维普遍胜率。

如果3D只获得同图更低成本，应定位求解器加速；如果同成本主要材料更好，应再解释哪些可观测/极化方向不同。二维四cell同相权重的解释不能直接用于3D纵/横极化。

## 11. GT隔离、开发循环与非机械决策

生成器/评分器与在线solver分进程，测试包只含白名单公开输入。私有生成参数、GT、clean field、family label不进入solver日志或宽泛文件搜索。选根、incumbent、停止、下一候选与预算分配全部只用已公开输入和在线物理量。

每个arm预测/轨迹/配置/源commit先封存并哈希，再由独立scorer读取GT。不同方法的预测应在统一选择规则固定后再评分；不能看先跑方法的GT结果给后跑方法调设置。24分两批也保持这一规则。

历史公开对象允许用GT作开发分析，但必须注明“暴露开发”。若独立测试被读取后据其调法，这批从此转开发；保留原失败，不改名称，不用同一批再次主张新确认。新增小确认集是版本变化的验证，不是为了凑成功率无限抽样。

科学失败后的流程是：保存反例→检查物理/代码资格→判断否定哪条机制→一次有理由的简化或改动→用独立数据验证。数值bug可修复并按同协议重跑受影响任务，旧attempt费用保留；科学反例不能当bug删除。

没有“成功率必须提高X%否则停止所有工作”的自动gate。若主方法没有独立收益，应完成no-go/反例/适用条件报告，比较更简单AP/直接方法并转向材料选择或观测设计；这是完整研究结论，不是未完成。

## 12. 输出schema与最低交付

以下是**本轮建议新建的输出名**，不是声称R11已有同名文件：

- `PROTOCOL_FROZEN.json`：问题定义、cohort、方法、指标、预算、选择规则和hash。
- `RUN_MANIFEST.jsonl`：每次attempt来源、配置、硬件、状态、错误与费用。
- `TRACE.jsonl`及独立数组：必要checkpoint、两个shadow、原物理诊断和完整迭代状态。
- `PREDICTIONS.npz`：每臂所有可选合法候选、最后shadow、最终选中材料和索引，不含GT。
- `COST_LEDGER.csv`：每类线性求解、matrix action、setup、NN、失败和I/O的费用。
- `SCORES.csv`：解盲后材料/data/phase/topology指标。
- `MECHANISM_VERDICTS.md`：每个预测支持、否定、未到适用条件的理由。
- `RESEARCH_REPORT_ZH.md`、`REPRODUCE.md`、`CHAT_HANDOFF.md`：完整科学结论、复现、下一轮入口。

每条run记录至少包含：

```json
{
  "experiment_id": "M1_or_M2_or_M3_or_LOCKED2D_or_EXTERNAL",
  "case_id": "opaque_id",
  "material_cluster_id": "shared_across_grid_noise_aperture_pairs",
  "cohort_role": "exposed_dev_or_new_dev_or_locked_test_or_real",
  "source_commit": "exact_sha",
  "protocol_hash": "exact_hash",
  "method_id": "frozen_method_name",
  "init_hash": "same_anchor_or_declared_NN",
  "physics": {
    "dimension": 2,
    "state_dofs": 576,
    "material_dofs": 576,
    "time_convention": "exp(-iwt)",
    "observation_contract_id": "exact_public_contract",
    "noise_and_mask_id": "exact_public_contract",
    "matched_or_cross_grid": "declared"
  },
  "selection": {
    "criterion": "original_full_Maxwell_Phi_with_original_prior",
    "chosen_candidate_id": "public_online_id",
    "truth_used_for_selection": false
  },
  "cost": {
    "cpu_seconds": 0.0,
    "wall_seconds": 0.0,
    "maxwell_factorizations_by_order": {},
    "maxwell_rhs_by_order_and_type": {},
    "fixed_projection_factorizations_by_order": {},
    "fixed_projection_rhs_by_order": {},
    "operator_actions_by_type": {},
    "nn_training_and_inference_seconds": {},
    "unknown_cost_fields": []
  },
  "termination": "actual_status_not_imputed_convergence",
  "prediction_hash": "hash_before_scoring",
  "gt_access": "scorer_only_after_prediction_seal",
  "failure_or_limitation": []
}
```

schema示例中的零仅是字段类型示意；未测费用必须null＋原因，不能复制零。材料指标在独立scorer输出，不能预先写进在线run作为选择输入。source/receiver/current/material维度分开，实际3D字段必须填真实值，不沿用示例576。

最低科学交付必须包含全部材料图、正负例、同预算连续曲线、成功集合overlap、one-to-one配对费用、初始至最终NN信息变化、模型误差与算法错误分层，以及“新结构到底比AP/通用restart多了什么”的明确答复。没有优势时应给出最强反例与下一项可执行理论修正，不只输出一句FAIL。

