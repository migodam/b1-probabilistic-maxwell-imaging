# Research A 上游原始证据核对

日期2026-10-15。只读核对；新增物理0、执行测试0、SHA检查0；没有修改历史结果/账本或B工作区。所有结论来自下述当前本地实文件，未把旧计划当已运行结果。详细机器表见 EVIDENCE.json，A23全部gate行见 A23_ALL_GATE_ROWS.csv。

## 科研影响与可复用资产

|资产|实际版本与证据|可复用内容|必须保留的限制|
|---|---|---|---|
|A17 closure_final|无独立Git；快照在Gaussian/A17/CLOSURE_R1/remote_snapshots/closure_final_closed_v1/source；Exchange baseline记录56097653|verified_exchange.py的纯FP64 DirectionAction、online_tolerance.py、closure_state_v2.py、材料注入/Schur/finite-gain数学|这是冻结GN步收益，不能充当候选full nonlinear数据相容性验收。历史wrapper会自动hash；本任务没有调用。|
|A23 actualpilot|a47cda424dfb55792e69d952a25f54ce9a2a5613；A23_RESEARCH_REPORT_ZH.md和A23_GATE_DECISION.json|src/a23/physics.py/encoder.py与src/a20/backend.py；全3456实材料、相对反馈/极化率归因、真实失败|H1/H2/H3及vanillaIBS2control全部FAIL；STOP_NO_EXPANSION。旧EXECUTIVE_DECISION是实验前建议，不能替代实际负结果。|
|A24 actualimaging|ca80ca2a29867e1b03589061dd5c17688c9acbe0；A24_EXECUTION_REPORT.md和results/b0/GATE_DECISION.json|src/a24/physics.py Physics2D(24x24,8源32接收)，vendor/a2/physics.py，relative_feedback.py、material_updates.py|G1PASS但G2FAIL/G3HOLD/G5质量未达；CURRENT_COMPRESSED_FEEDBACK_NO_GO。理论目录的GO仅是运行小实验。|
|C1-P1|C1当前83056881e89329212b16f5be47a4d0dbc4465943中的p1/，历史P1原样|2D bounded normalized density、likelihood接口、成熟distillation和失败轨迹|P0是13dipole intensity问题，非本次Mie/A30。State/memory没有稳定增量，未达到完整posterior。|
|C1-P2|同8305688的p2/；15真实.pt本地存在|4D bounded MDN、FP64 numpylogq/sample、full/diag/Gaussian对照；不必重训|仅冻结径向proxy目标，event好不等于条件材料CDF好。不能把rawmodel当fullMaxwell合法候选。|
|B1×C1P2|B1 7ed4faef9bc3dd214d4de689c2203cc383f7495f，research/joint_p2/|相同prior/decoder/噪声、已支付responsecache、acceptedcontinuousreference、Miefullsolver、candidateforward数组|径向校准台；本轮不得借用剩余541次。普通ring编辑仍不可补偿，r5取消已被否定。|

## A17真实正面和反面

- 十个完整额外对象的冻结GN gap中位下降50.7%；96/99额外计划cell有审计，95改善1相同，1失败+2未跑保留。72noisecell为71改善1no-op。
- 三个完成非线性pair材料复误差改善10.84%、0.22%、7.16%，child-wall比2.68、1.74、1.46，全部更慢。第四对象仅17-update安全前缀，无独立receiverbaseline完成pair。
- paired repair仍输同秩Krylov；宽域from-scratch选择中位receiver风险比1.0947。总成本19,962.060秒；不能把局部保真/费用换成全局多解或更快成像。
- 用户要求的closure_final materialized目录存在：Gaussian/A17/CLOSURE_R1/publication/package_final_v1/public_reader_materialized_v1/closure_final/。主要可读报告Gaussian/A17/CLOSURE_R1/A17_CLOSURE_REPORT_ZH.md；closure源快照不是独立git，不能伪造新commit。

## A23全部失败Gate

|Gate|实际证据|结论|
|---|---|---|
|H1|cached与original transfer相同但完整cold仅降0.35–0.86%；sourceanchored约9.7–10.4%transfererror；directadjoint精确但更慢。阈值为同≤5%质量下降≥20%。|FAIL|
|H2|四暴露父对象，nominal两噪声全部不qualify；2001/2003/shell多数材料非法，asymmetric合法却无增益。|FAIL|
|vanillaIBS2control|2003原始NRMSE改善约26.6/23.3%，仍物理拒绝，没有第二家族，不能以raw低误差宣布成功。|FAIL|
|H3|唯一运行full合法prediction为2014，材料NRMSE.81237>共同.60目标；无合格成本点。|FAIL|
最终204预测79OK/125物理拒绝；fullvalidation1RUN/19REJECTED/184NOT_RUN。复噪声Re/Im各σ²/2，不能套P2的σ约定。四历史场景是统计单位，不是独立新holdout。全tangent数据行rank1536，不能把稳定截断或32chart的complement当full材料新信息。

## A24实际成像失效

224条件136OK、4FAILED、84共同初始化失败。失败是TV-prox KKT cap，不能抹掉或说Maxwellidentity失败。R16有六个已观测transfer>.05的反例(最高.0535177)；task20/32，不完整图像也不能撤销真实反例。G3四parent完整八条件均缺失，官方均值null；弱背景完整子组F改善3.97936%、R3.74175%，未达10%。没有方法达到共同质量，因此无matchedcost收益。3D与NN均NOT_RUN。详见EVIDENCE.json保存全部反例。

## C1生成器的正负结果和可用性

P1 P0 full-distill三seed事件MAE .01031/.01487/.01205，均值.012409、KL.13208；diagjoint .050736/.996777。U0 fullcov不稳定，diagonal-distill .042439优于full-distill .042692；P0state/memory一seed改善两恶化，memory事件误差三seed均坏。U0float64logistic还有端点舍入负例。教师监督是额外信息，不能称同监督架构的纯优势。

P2采用4D材料latent、不输入topology/state/memory。Full/diagdistill accepted事件MAE .013215/.013028，paired差+.000187，95%区间[-.000580,.000880]；无full事件优势。原13阶latentKL .216/.275有full改进，但非accepted/非materialKL。accepted重要事件内gridCDF误差.14176/.16177，都远未到.01。15模型本地.pt真实存在；模型只负责proposal，候选仍需新fullMaxwell验证。P1 p0/u0目录各36/24个.pt（含观察选择副本），P2为15个，共详CHECKPOINT_INVENTORY.json。没有打开旧lockedtruth。

## P2可继承和不能外推的机制

固定实心ε2.5+.1i与中心r=.25R候选，noise-distance raw1.52965→epsilon1.29522→epsilon+boundedgain.90017；合法参数ε2.466617+.130009i,g1.053221。两频共享gain后残差2.46794，材料距离仍.12846。它是best-found可行witness不是全局min。
ringh.08单频补偿残差4.45054、两频8.02208，仍不相容；central+ring也失败。小radius raw/projected slope2.9778/2.9747，不支持普遍r5取消。实际P23459/4000已封存，本ResearchA新账本独立。

## 资产缺口与立即适用的边界

没有发现ResearchA目标必需的C1P2checkpoint缺失；不应为补模型重跑P1/P2。A17closure没有独立gitrevision，只有源快照和记录baseline；报告应如实标注。A23部分早期源码snapshot缺失属于历史缺口，不应重跑掩盖。SOM/FFTTSOM、CIE-I/WTTSOM不是A24已实现公平baseline：A24明确记录前者prior/optimizer不同、后者NOT_IMPLEMENTED。若本轮不能真正适配，可报告NOT_RUN，不应命名普通多起点为TSOM复现。A30由父线程/另一原生worker负责，不在本只读子任务内重复读取或运行。

主Codex决定ResearchA理论、实验方向与最终有效性。此表只冻结可核对的先验事实；0新增物理。
