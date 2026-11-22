## Maxwell 实验手册

查已冻结的数据、观测身份、物理算子、设备选择、基线与错题，按同一口径接续和归档实验。[项目手册](experiment_handbook/README.md) · [数据与算子目录](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/experiment_handbook/catalog.html)

<!-- Z_ALG_OPTICAL_INERTIA_PUBLICATION -->
## Z-ALG · 光学定理与材料对比度下界

在 96² 有限模型中，用实部和光学定理约束检验六件新强散射材料。2% 噪声时六件均有严格排除，10% 时有三件，两件显示严格多源优势。报告也给出 4×4 区域检验和 Fresnel 实测的适用边界。

[中文报告与图表](z_alg_optical_inertia/index.html) · [数学论证](z_alg_optical_inertia/theory.html) · [Fresnel 实测](z_alg_optical_inertia/fresnel.html) · [ChatGPT 研究交接](z_alg_optical_inertia/handoff.html)
<!-- /Z_ALG_OPTICAL_INERTIA_PUBLICATION -->

<!-- Z_ALG_OPTICAL_INERTIA_SUBMISSION -->
## Z-ALG · 有限孔径散射数据的介电峰值认证

英文投稿稿、原生 96² 分块对照、配对噪声曲线与严格证书，保留定位和实测的失败范围。

[中文定稿报告](z_alg_optical_inertia/submission/index.html) · [英文论文 PDF](z_alg_optical_inertia/submission/manuscript_ip.pdf) · [复现与下载](z_alg_optical_inertia/submission/reproduce.html) · [研究交接](z_alg_optical_inertia/submission/handoff.html)
<!-- /Z_ALG_OPTICAL_INERTIA_SUBMISSION -->

<!-- Z_ALG_OPTICAL_INERTIA_VALIDATION -->
## Z-ALG 光学惯性 · 独立验证修订版

既有 SRAP、MAX-STRUCT 材料迁移，Manitoba 实测误差、校准转移及各自适用范围。

[中文验证报告](z_alg_optical_inertia/validation/index.html) · [修订英文稿 PDF](z_alg_optical_inertia/validation/manuscript_ip.pdf) · [复现与下载](z_alg_optical_inertia/validation/reproduce.html) · [研究交接](z_alg_optical_inertia/validation/handoff.html)
<!-- /Z_ALG_OPTICAL_INERTIA_VALIDATION -->

<!-- Z_ALG_SECANT_TASKS_PUBLICATION -->
## 固定 Green 下的材料峰值排除与镜像停测

在冻结的二维离散物理模型下，旧环形观测排除了 576 个像素的介电对比度全都不超过 5 的情形。该材料范围内的最低预报误差下界约为 0.014695，超过原噪声半径 0.010517。无需重建整图即可断定至少一个像素超过 5，位置仍未知。另用镜像排除定理回放历史合成读数，152 份符合预设误差上限的 SRAP 视图中，逐组 428 次、全批 1,517 次和事后同布局最短固定前缀 926 次，均排除同一批 150 份视图的镜像歧义。

[多源排除证书](z_alg_secant_tasks/multiple-illuminations.html) · [中文报告与图表](z_alg_secant_tasks/index.html) · [ChatGPT 研究交接](z_alg_secant_tasks/handoff.html)
<!-- /Z_ALG_SECANT_TASKS_PUBLICATION -->

<!-- Z_ALG_HOST_FOLD_PUBLICATION -->
## 已知非均匀宿主中的材料折叠

一个固定二维标量 TM 模型，30 个已知宿主单元与 6 个未知材料单元。沿两组指定三通道复读数之间的直线，整个材料盒的合法解数依次为两份、一份、零份。报告展示折叠、全阶正演与带噪材料范围的证据。

[中文报告与图表](z_alg_host_fold/index.html) · [ChatGPT 研究交接](z_alg_host_fold/handoff.html)
<!-- /Z_ALG_HOST_FOLD_PUBLICATION -->

<!-- Z_ALG_R16_PUBLICATION -->
## Z-ALG R16 · 在电流域选方向，能否改善 576 格材料重建？

用 FFT-TSOM 更新多源感应电流，每轮再独立更新 576 个材料格。报告给出电流空间可表示的材料方向，并用历史二维 Maxwell 观测检验这些方向是否改善材料恢复。

[中文报告与图表](research_z_alg_r16/index.html) · [ChatGPT 研究交接](research_z_alg_r16/handoff.html)
<!-- /Z_ALG_R16_PUBLICATION -->

<!-- Z_ALG_R15_PUBLICATION -->
## Z-ALG R15 · 多测几个频率，还是多换几个角度？

2–30 GHz、576 个未知材料像素。相同的 64 个读数，增加发射角能降低整图误差，换频有时更能分开相邻缺陷。报告给出两种宽带 MIMO 优劣反转的材料、独立选频实验和 Fresnel 实测对照。

[中文报告与图表](research_z_alg_r15/index.html) · [ChatGPT 研究交接](research_z_alg_r15/handoff.html)
<!-- /Z_ALG_R15_PUBLICATION -->

<!-- Z_ALG_R12_PUBLICATION -->
## Z-ALG R14 · 换频后还有多少份材料解释同一组读数？

用 62 份已有带噪观测估计材料范围。在给定噪声合同下，18 个自然案例加测双频 ±2% 后，每案仍有两份相差超过 1 个 RMS 单位的材料能够解释同一组旧新读数。

[中文报告与图表](research_z_alg_r14/index.html) · [ChatGPT 研究交接](research_z_alg_r14/handoff.html)

## Z-ALG R12 · 换频率、换照射位置，能解开哪些材料歧义？

给定已知局部分支，相对噪声和相对带宽可以预估错误极小留下的概率。两个带噪仿真验证了风险随带宽下降。另在 576 维仿真中，近场看到了原距离漏掉的材料差，低到高累计拟合同份宽频读数也改善了三例旧坏起点的材料恢复。

[中文报告与图表](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r12/) · [ChatGPT 研究交接](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r12/handoff.html) · [交接原文](research_z_alg_r12/NEXT_PRO_HANDOFF.md)
<!-- /Z_ALG_R12_PUBLICATION -->

<!-- Z_ALG_R10_PUBLICATION -->
## Z-ALG R10 · 再照一次，排除一类材料解释

二维全波仿真中，同一组 22 束旧读数容许 24 个独立未知材料的两份完整解释。新增一束后，第 23 个材料的所有可行值都落在 [9.28, 9.66]。同功率随机束和只按几何选的束也排除了远支。在四片六材料区域上，降维能以直接六维法约 1/5 到 1/48 的 CPU 时间给出相同的严格排除。另有一份错误材料能通过全部 190 个两项检验，只有共同使用整批读数才被识破。

[中文报告与图表](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r10/) · [下一轮研究交接](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r10/handoff.html) · [交接原文](research_z_alg_r10/NEXT_PRO_HANDOFF.md)
<!-- /Z_ALG_R10_PUBLICATION -->

<!-- Z_ALG_R7_PUBLICATION -->
## Z-ALG R7 · 三次照射估计指定区域的材料值

在两套二维电磁仿真里，十个发射器配合照三次，就能估计指定方形区域的材料值。背景材料整体比模型高 1% 时，32 个样本中只有 14 个的估计范围还包含真值。两套布局各单测一次背景，再调整照射方式，32 个样本的范围都重新包含真值。

[中文报告与图表](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r7/) · [下一轮理论交接](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r7/handoff.html) · [交接原文](research_z_alg_r7/NEXT_PRO_HANDOFF.md)
<!-- /Z_ALG_R7_PUBLICATION -->

<!-- Z_ALG_R5_PUBLICATION -->
## Z-ALG R5：认证歧义与物理主动测量

CONDITIONAL GO：固定精确观测的全box材料覆盖与一源分支分离已认证；全噪声可行集覆盖、总体恢复优势仍未成立。原失败和后续Codex修订分开报告。

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r5_measurement/) · [数学证书](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r5_measurement/theory.html) · [完整复现包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r5-measurement-20261118)
<!-- /Z_ALG_R5_PUBLICATION -->

<!-- MAX_STRUCT_R3_PUBLICATION -->
## MAX-STRUCT-R3：公式驱动的机制研究循环

噪声硬锁代价、被动单元圆弧与更新秩层级；真实边界停滞和3D适用性控制。11个同一公开对象的开发控制未建立新增成像优势，全部失败、成本和后续方向保留。

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r3/) · [给Chat的入口](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r3/CHAT_HANDOFF.md) · [数学推导](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r3/THEORY_DERIVATIONS.md)
<!-- /MAX_STRUCT_R3_PUBLICATION -->

<!-- MAX_STRUCT_R2_PUBLICATION -->
## MAX-STRUCT-R2：结构审计、新盲测与新增照明

三新对象四方法：成熟CIE/global均2/3，精确profile1/3，未建立超越成熟CIE的新成像机制。满容量事后控制仍失败；三维新增照明只支持公开一参数族中的相容估计。代码、原始结果、独立审查与全部费用保留。

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r2/) · [给Chat的入口](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r2/CHAT_HANDOFF.md) · [完整代码/证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_max_struct_r2)
<!-- /MAX_STRUCT_R2_PUBLICATION -->

<!-- MAX_STRUCT_R1_PUBLICATION -->
## MAX-STRUCT-R1：张量同场族与机制边界

合法张量材料在有限照明下保持同场；新增接收与旧源混合不增加信息，独立新照明仅部分例子超过噪声尺度。M表示未完成资格、T非空证书未建立，未声称新成像优势。代码、原始数组、失败与成本全部保留。

[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r1/) · [Chat原始Markdown](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r1/RESEARCH_REPORT_ZH.md) · [代码/证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_max_struct_r1) · [完整复现包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r1-20261102)
<!-- /MAX_STRUCT_R1_PUBLICATION -->

<!-- Z_ALG_R3_DIGITS -->
## Z-ALG-R3-DIGITS · 576 材料九数字盲成像与基线审计

[中文报告](research_z_alg_r3_digits/RESEARCH_REPORT_ZH.md) · [九数字真实图](research_z_alg_r3_digits/figures/digits_1to9_reconstruction.png) · [代码与复现](research_z_alg_r3_digits/REPRODUCE.md)

九数字全部运行，清楚形状恢复有限；Native实现、成熟TRF及跨session物理区间已审计。保留失败、GT后诊断、噪声14/27与完整费用；不声明普遍Full-J优势或九数字成功。

<!-- Z_ALG_R3_CURRENT -->
## Research Z-ALG-R3：代数候选与真正的盲成像

12个材料全部未知的非线性VIE盲成像。新四例中，完整物理残差在线选出的兼容性候选组合通过3/4个数据门，强倒数Full-J对照2/4；总wall没有下降。局部Q6/Q20分支搜索未证明成像增益，失败保留。一般优势与独立论文新颖性仍HOLD。

[中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r3/) · [latent与特征子空间推导](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r3/theory/LATENT_AND_CHARACTERISTIC_SUBSPACES.html) · [原始Markdown](research_z_alg_r3/RESEARCH_REPORT_ZH.md) · [代码、图像与完整复现包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r3-20261030)
<!-- /Z_ALG_R3_CURRENT -->

<!-- P3LR10_CURRENT -->
## P3L-R10-BI：双边传递、条件物理复用与材料覆盖边界

[阅读完整中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r10_bilateral_structure/) · [140条公式集](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/maxwell_formula_atlas/) · [无损证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr10-bilateral-20261030)

固定形状的320次新幅度查询全部通过误差门；实际32节点冷形状费用节省51.82%。新材料对象的覆盖优势未成立，Full66 posterior仍未认证，强TAP贡献HOLD。完整源码、失败、独立QA与论文草稿已保存。
<!-- /P3LR10_CURRENT -->

<!-- Z_ALG_R2_CURRENT -->
## Research Z-ALG-R2：行列式商代数、材料解枚举与盲成像

四材料六维商与六材料20维商求解器已实现。六材料新盲测返回2/1/1个合法根，固定64起点对照为1/1/1；双方包含GT，额外材料相距1.15879 RMS。三目标单批次0.1044秒对1.5180秒，仅支持本协议。48单元成像及SVD噪声失败独立保留。53,951完整复RHS；全球新颖性与严格全根认证仍HOLD。

[中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r2/) · [15组推导公式](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r2/FORMULA_COMPENDIUM.html) · [原始Markdown](research_z_alg_r2/RESEARCH_REPORT_ZH.md) · [阅读包与完整代码数据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r2-20261030)
<!-- /Z_ALG_R2_CURRENT -->

<!-- P3LR9_CURRENT -->
## P3L-R9-S：材料可能性体积、稳定子空间与长期公式集

32个全新parent先封存预测再评分。固定方法354/512个方向预算产生不同且数据相容的材料，但精细GT覆盖仅6/32（anchor3/32），高概率覆盖未建立。相对动态方法整链wall节省2.75%，未达到20%门。固定host精确ROI消元的匹配物理核节时58.65%；不冒称全系统收益。正确低维积分与全部失败独立保留。

[中文完整报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r9_stabilized_subspace/) · [长期维护的122条推导公式](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/maxwell_formula_atlas/) · [匿名Markdown](p3l_r9_stabilized_subspace/RESEARCH_REPORT_ZH.md) · [证据与复现](p3l_r9_stabilized_subspace/EVIDENCE_PACKAGE.md)
<!-- /P3LR9_CURRENT -->

<!-- P3LR4_CURRENT -->
## P3L-R4 已完成：风险成像与完整条件概率分别裁决

[中文报告](p3l_r4_conditional_transport/RESEARCH_REPORT_ZH.md) · [报告与实际图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r4_conditional_transport/) · [GPT Pro交接](p3l_r4_conditional_transport/GPT_PRO_HANDOFF/START_HERE.md)

512+256个全新对象：主要材料风险动作相对R2 CNN改善13.10%/11.36%，相对更强CNN改善5.36%/4.90%；普通Gaussian融合相当。Full-66 proper score有相对进步，但完整后验仍未校准。固定核的数值Maxwell决定重放wall下降66.17%；动态非完整路径未建立独立优势。原始负结果与完整费用保留。
<!-- /P3LR4_CURRENT -->

- [Research T：Maxwell基础、拓扑与补偿机制](research_t_20261015/START_HERE.md) · [实际成像图库](research_t_20261015/)

# 研究总览

[网页总览](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/) · [GitHub主页](https://github.com/migodam)

2026-10-15两条独立路线：[Research A 多材料世界发现](research_a_20261015/START_HERE.md) · [Research B 主动测量与候选修正](research_b/START_HERE.md)。各自的原始结果、费用和负结论独立保留。

2026-10-15 Research G：[统一材料下一解搜索](research_g_20261015/START_HERE.md) · [实际图库](research_g_20261015/MULTISOLUTION_GALLERY/index.html)。新增一个近门槛合法边界结构，Best-of-K无改善，保留负结果；与T的物理工作独立计费。

## 主项目

|项目|内容|入口|
|---|---|---|
|B1|Maxwell物理参考、连续材料posterior与概率正确性|[最新E0](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/b1_next/e0/) / [仓库](https://github.com/migodam/b1-probabilistic-maxwell-imaging)|
|C1|统一材料生成与物理proposal研究|[仓库](https://github.com/migodam/c1-physics-state-generative-imaging)|
|A17|独立Maxwell电流–材料反馈理论与论文|[仓库](https://github.com/migodam/a17-maxwell-certified-exchange)|
|Full-wave self-calibration|相位、校准与早期研究资料；仍有未结事项|[仓库](https://github.com/migodam/fullwave-self-calibration)|

## 历史实验档案

归档是研究阶段管理，不改变科学结论。所有原URL、commit、Release和license保留。

|仓库|研究主题|状态|
|---|---|---|
|[a13-shared-physical-flow](https://github.com/migodam/a13-shared-physical-flow)|共享物理latent可行性|历史快照，原证据与负结果保留|
|[a14-task-relevant-posterior-coordinates](https://github.com/migodam/a14-task-relevant-posterior-coordinates)|任务相关posterior坐标|历史快照，原证据与负结果保留|
|[a15-fullwave-posterior-ambiguity](https://github.com/migodam/a15-fullwave-posterior-ambiguity)|全波posterior歧义|历史快照，原证据与负结果保留|
|[a14-pairtsom-maxwell](https://github.com/migodam/a14-pairtsom-maxwell)|PairTSOM与材料更新|历史快照，原证据与负结果保留|
|[a16-maxwell-current-selection](https://github.com/migodam/a16-maxwell-current-selection)|电流空间选择|历史快照，原证据与负结果保留|
|[a18-nn-native-maxwell-stage0](https://github.com/migodam/a18-nn-native-maxwell-stage0)|NN-native初始研究|历史快照，原证据与负结果保留|
|[a18-b-maxwell-feedback-closure](https://github.com/migodam/a18-b-maxwell-feedback-closure)|反馈闭包|历史快照，原证据与负结果保留|
|[a18-physics-native-current-proposal](https://github.com/migodam/a18-physics-native-current-proposal)|current proposal集成|历史快照，原证据与负结果保留|
|[a19-task-adaptive-feedback](https://github.com/migodam/a19-task-adaptive-feedback)|任务自适应反馈|历史快照，原证据与负结果保留|
|[a20-opm-imaging](https://github.com/migodam/a20-opm-imaging)|OPM成像|历史快照，原证据与负结果保留|
|[a23-fast-material-imaging](https://github.com/migodam/a23-fast-material-imaging)|快速材料成像|历史快照，原证据与负结果保留|
|[a24-source-anchored-om-feedback](https://github.com/migodam/a24-source-anchored-om-feedback)|source-anchored反馈|历史快照，原证据与负结果保留|
|[a25-adaptive-om-feedback](https://github.com/migodam/a25-adaptive-om-feedback)|自适应OM反馈|历史快照，原证据与负结果保留|
|[a26-om-feedback-prior-completion](https://github.com/migodam/a26-om-feedback-prior-completion)|prior completion|历史快照，原证据与负结果保留|
|[a27-feedback-failure-attribution](https://github.com/migodam/a27-feedback-failure-attribution)|反馈失败归因|历史快照，原证据与负结果保留|
|[a28-maxwell-material-information](https://github.com/migodam/a28-maxwell-material-information)|材料信息诊断|历史快照，原证据与负结果保留|
|[a29-material-hypothesis-coverage](https://github.com/migodam/a29-material-hypothesis-coverage)|材料假设coverage|历史快照，原证据与负结果保留|
|[a30-maxwell-observable-material-manifold](https://github.com/migodam/a30-maxwell-observable-material-manifold)|可观测材料manifold|历史快照，原证据与负结果保留|
|[a32-material-risk-directed-feedback](https://github.com/migodam/a32-material-risk-directed-feedback)|连续材料风险反馈|历史快照，原证据与负结果保留|

后续B1、C1分别复用本项目仓库，阶段放在目录/Release，不再为每轮新建repo。[发布规则](../RESEARCH_PUBLICATION_POLICY.md)。


## Research W P0 — 2026-10-22

[执行入口](research_w_p0_20261022/START_HERE.md) · [中文报告](research_w_p0_20261022/RESEARCH_REPORT_ZH.md) · [图文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_p0_20261022/)

一维解析回归PASS；二维前向网格资格未过，P0-B为UNRESOLVED。没有合格墙体歧义见证、搜索优势或P1结论。

## Research W continuation — 2026-10-23

[图文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_continuation/) · [完整中文报告](research_w_continuation/RESEARCH_W_CONTINUATION_REPORT_ZH.md) · [复现与原始数组](research_w_continuation/REPRODUCE.md)

P0-Q通过，限定为QUALIFIED_DISCRETE_FULLWAVE；原两墙体8项完整观测细化差0.0065–0.0332。P0-B完成3/4组四方法搜索，最小cal距离19.91，未达到1，未发现合格歧义对；搜索另有0.1最小特征未严格执行的实施缺口。joint与P1未运行；不把搜索失败说成唯一性。旧P0及其结论保持原样。


## Research W P0-R → P1 — 2026-10-24

[中文报告与实际图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_p0r_to_p1/) · [研究报告](research_w_p0r_to_p1/RESEARCH_REPORT_ZH.md) · [概率材料交接包](research_w_p0r_to_p1/MaterialHypothesisPackage/EXAMPLE.json)

从同合同真值回放与FWI诊断推进到条件固定支撑的joint补偿见证（d=0.798716，目标位移0.1λ₀），随后独立新观测的proper两墙概率参考和四格消融。原minimum_feature仍未定义，在线anchor未达标；多假设有受控平方风险收益，progressive与压缩没有独立端到端优势。历史P0/continuation及其负结果保留。

- [Research W-SLAM](research_w_slam/)：已完成的物理、概率与成本结果；限制及负结果保留。
