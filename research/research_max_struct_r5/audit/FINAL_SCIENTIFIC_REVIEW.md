# R5 最终独立科学审查

2026-11-03。Sol high；只读主报告、推导、计划、实现及已经落盘的诊断/离线评分/谱/成本结果。没有读取 GT 材料数组，没有执行 forward、Jacobian、HVP、优化、恒等式测试或新科学实验；没有改主报告或诊断程序。本审查是给父线程的独立批评，最终科学判断与发表决定仍由父线程负责。

## 审查结论

**没有发现应阻止当前“公开开发机制研究／有限干预及负结果”交付的科学错误。** 位置消融、单一自然阶段切换、块谱和 vector-pair 结果足以支持报告已经限定的作用位置与适用边界；它们不支持普适恢复、稳定成本优势、渐近盆地定理、未知材料三维成像或首次发现标准公式。报告大体保持了这些边界。下列两处措辞建议在最终版本澄清，无需新物理。

## 具体问题与严重性

### P2：把负梯度方向说清，避免将梯度符号说反

位置：`RESEARCH_REPORT_ZH.md:60`；相关推导 `THEORY_DERIVATIONS.md:62–65`。

报告称“普通CIE的这一个分量恒指向提高contrast”。公式中的 `a_r` 本身是非正**梯度**分量，真正指向 contrast 增加的是 `-a_r` 这一负梯度／下降方向。读者若按梯度向量理解现有句子，会得到反号。建议改为“该分量的负梯度指向提高contrast”。这不改变公式、结果或机制边界；报告已正确说明这一分量也出现在成功轨迹，不能判失败。已通知父线程。

### P3：首封文档演化披露正确，但建议直接链接审计记录

位置：`RESEARCH_REPORT_ZH.md:89`；`audit/POST_SEAL_DOCUMENT_EVOLUTION.json:4–13`；`audit/THEORY_DERIVATIONS_AT_FIRST_SEAL.md`。

现有披露已经足以否定“发布时全部首封路径未变”的错误读法：主报告说明前六节首封版本保存在 audit，后续附节不是预注册内容；演化 JSON 明确原始 seal 包含 THEORY、首封记录6480 bytes、7–9节在评分后增加、首封前缀是在该记录时间重建保留而非声称首封当时已复制，并明确最终 THEORY 故意比首封新。审查只作直接 byte-prefix 比较，确认保存副本为6480 bytes，且是当前 THEORY 未变的完整前缀；未作 hash 检查。

建议主报告89行直接链接 `audit/POST_SEAL_DOCUMENT_EVOLUTION.json` 和保存副本，而非仅泛指“audit”及复现／账本。不要在后续网页、README或发布说明中重新写“all sealed paths unchanged”。本次读取范围内没有发现主报告作此错误声明。此项是可追溯性改善，不需要改历史 seal，也不使既有预测失去封存资格。

## 位置消融与阶段切换

- 旧弯条四臂表与 `RESULTS/OFFLINE_FACTORIAL_SCORES.json` 对应：普通／普通 `.389182533`，one-hop／one-hop `.041672481`，current-only `.558513628`，material-only `.041480425`；IoU与物理残差亦一致。后两臂是新R5，前两臂为已付费历史参考，报告22行有明确说明。
- 该干预支持“这一实例的恢复不需要 current one-hop，完整材料 readout 位置足以复制该正例”。它没有隔离材料 readout 内所有机制，也没有证明非局部项在其他设置中必要或充分。报告24、73、79、95行没有把它夸大成普适算法优势。
- `cie_stage_handoff.py:172–226` 根据原生 Fourier frequency keys 搬运系数并在自然新阶段清空 PR 历史，符合原策略；两份 `HANDOFF_WIRING.json` 的重建 J 差均为0。不是只搬运 m 后丢弃辅助状态的伪阶段消融。
- 两 handoff final NRMSE 分别 `.388613370`、`.558737600`，均 FAIL；corrector均到1000步上限、没有收敛资格。主报告28行和 THEORY107行把结论限定于这个切换点，并明确不能推出每一步都需要传播或任意切换时间失败，边界恰当。
- 三个 extension 与离线评分一致：42cb `.766322569`、63ab `.404823109`、环 `.072525…`；两个失败仍为失败。报告42行明确跨平台历史参考不等于同执行对照，63ab较低材料误差伴更差数据残差，未用于GT在线择优。无新独立 blind cases。

## Walsh 谱与异质材料

`THEORY_DERIVATIONS.md:84–95` 与 `experiments/analyze_block_modes.py:16–59` 一致。四个单位 Walsh 模态在给定正方块次序下对应 `2ge+gd,-gd,-gd,-2ge+gd`。对 uniform m，`a`是标量，固定正交 Walsh 使算子正规，`|1+a lambda|²` 确为作用于 **u** 的残差权重；它不应混作整个从原始 r 到 R 的权重（后者还有 self scaling），主报告明确以残差滤波解释并没有拿它替代物理 Jacobian。

对异质 m，真实 modal operator 是 `H diag(a) Goff Hᵀ=(H diag(a) Hᵀ)diag(lambda)`；实现47行直接保留实际 block，不把 homogeneous 公式逐cell应用。modal cross 使用真实保存的 u、Lu，且另列 diagonal/offdiagonal 项，因此没有因忽略模态混合而伪造总体 cross 符号。

从保存 JSON 只读汇总，one-hop四对象×三阶段最后检查点共12态的 constant-energy fraction 为 `.709493301–.834385757`，全部 constant cross正、其他三模态总cross负，与主报告71行一致。它们是同一批状态的评分后解释；`MODAL_ANALYSIS_PLAN.md:3`、modal JSON classification和主报告89行均明确 post-score development。不是新独立预测，也不等于存在成功分类器。

## 三维 vector-pair 反例

`THEORY_DERIVATIONS.md:99–101`、主报告81行与实际配置/脚本一致。两体素沿x，入射传播方向均为+y，极化分别x/z，因而两波都严格横向；“longitudinal”是相对体素连线的x方向。没有非法的入射纵向电磁平面波。

脚本 `vector_pair_mechanism.py:55–75` 使用真实 `model.G` 和 scalar self density，正确保留体积缩放；没有把未乘体积的 dipole kernel 当密度 Green。对相同各向同性 susceptibility m=2，正反空间模态乘纵向／两个横向给出六模态，`.2`下同相纵向权重 `1.342836868`、同相横向 `.856192834` 与 JSON一致。静电 `diag(2,-1,-1)` 仅解释符号，没有替换实际retarded场或权重。

这一例子确实反驳“二维 scalar 同相增强可无条件搬到所有矢量同相极化”的外推；不是关于三维反演恢复的反例。完整全波解只有已知材料两cell的1LU/2RHS；显示的权重是 residual operator 的谱，并非完整物理场／接收数据的信息谱。在物理根处本构残差仍为零，权重表本身不证明某种实际残差轨迹已出现。主报告和 SPECTRUM scope 没有跨越该边界。

## 原物理曲率、GN、prior与约束

先前审查提出的 prior 口径已经修正。保存端点 JSON / `MECHANISM_SUMMARY.json` 现在分列 `gn_data_directional_curvature`、`prior_directional_curvature`、`exact_data_directional_curvature`、`residual_feedback_directional`、`exact_Phi_directional_curvature`，未继续把 full-Phi-minus-data-GN 叫纯 residual feedback。ordinary/self 解析零 transport 也已置零并将夹角标null，保留 raw roundoff norm。

在stage2最后one-hop检查点，两个失败例 full-Phi方向曲率分别 `-5.8452491e-9`、`-3.7922429e-9`；旧弯条／环分别 `8.4085610e-10`、`1.9994025e-8`，所述符号一致，且精确 prior曲率为非零正量、量级另列。非零 projected gradient 仍存在，所以不能称极小点、吸引盆地或负曲率逃逸已经成立；主报告46行明确排除这些说法。

这些几何在 `m_before`（i500最终readout前）沿该实际readout步计算，既不是corrected final prediction的曲率，也不是整个box切锥上的最小曲率。MECHANISM_SUMMARY warning已经声明。主报告46行“最后辅助阶段的检查点”准确；95行的“受约束物理曲率的区分”宜按此有限含义读，不能改写成已求得受约束 Hessian 谱或切锥证书。下一轮105行明确提出尚待检验自由面／切锥，未冒称已经完成。

## 成本、样本、物理资格与新颖性

账本分组之和与主报告一致：`7717 LU / 183686 complex RHS`，其中一个LU为6维vector系统，其余为576维2D系统；task wall总和 `290.670004584s`，vector约`.273934s`，二维约`290.396071s`。报告明确task wall不是并行历时或LLM／文献／绘图全成本；货币／tokens null、失败1LU/24RHS收费、历史引用不重计。没有把LU数量直接当跨维等成本，也没有等wall胜出主张。

174 states、120同态readout counterfactuals中含6个历史固定态的新增诊断；MECHANISM_SUMMARY的168 actual steps与96 counterfactual只覆盖冻结replay，两个分母不可混用。主报告7行说明168不是独立样本，62行把34／19次材料与Phi相反变化当相关描述，无显著性夸大。

新颖性账本保留了既有 variable/fixed metric、CIE、CSI、Schur/profile及饱和先例，明确关键文献访问缺口；“未核实精确先例”没有被写成“首次”。标准乘积微分、Walsh、dyadic及双边导数均不主张原创，最终创新候选聚焦可复现的具体干预与反例。这个保守定位与已读证据相容；本次审查不是一次新的穷尽文献检索，也不提供IEEE TAP接受资格或标准公式首创裁决。

## 交付范围

可保留当前机制研究／开发负结果定位。最终发布前建议完成上述负梯度措辞和首封演化直链两项编辑，并保持相关文档在父线程完成后链接可达。审查时 `TRAJECTORY_ANALYSIS.md`、`RESEARCH_LOOP.md`、`FAILURES.md`尚未落盘；已知由父线程另写，故不把暂缺当科学错误，也未等待或替其编写。主报告未由本审查修改。
