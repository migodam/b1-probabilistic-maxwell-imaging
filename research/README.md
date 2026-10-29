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
