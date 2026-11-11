# R12 科学批判审查

审查对象为本轮已完成分支、根报告和 `PAPER_DRAFT.md`。没有新增反演/forward/实验，没有自行对确认24做另一套scoring；已读baseline worker完成的官方paired/bootstrap及追加保存Yclean诊断。审查提出的basin pool范围问题，后由root明确授权追加0FW封存后离线oracle诊断，见第3节。该审查与执行/评分分开，但审查者参与过core、NN与high3D实现，**不是外部盲评**；最终数学、机制及novelty裁决仍归root。

## 1. 当前总判断

1. **窄化后的endpoint escape主张有实际证据。** 环/壳从相同保存坏材料完整物理初始化，direct继续停留，AP未给同样改善，DR进入较低原Phi和更好材料区域；原端点的近局部极小数值证据使这一干预有意义。桥存在负曲率/未收敛，只能作为endpoint/path干预，不是第三个已认证局部极小逃逸。
2. **不能升级为默认求解器普适优越。** 成熟toy shell控制、NN逐对象反例、较长corrector后3D趋同、Austria matched32中DR仍失败，以及确认初步primary均值都反对该升级。不同初始化下老CIE/G19已有好解不能从历史中删去。
3. **最清晰的新实证是搜索与选择分离。** high3D中原Phi所选DR raw既准确又物理CN99可行，corrector下降Phi却牺牲data、材料并出球。另有D/R同时改善而材料变差的真实反例。这是对当前estimator的有限经验边界，不是梯度错误或任意prior都失败的证明。
4. **文稿大部分限制已有披露，但仍有具体应改的全称/归因。** 以下先列明确错误，再列结论必须保留的限定。修正这些句子不要求新实验。

## 2. 需修正文稿的具体项

|优先级|审阅版本中的句子/位置|直接证据与问题|建议修正|
|---|---|---|---|
|P1|[BENCHMARKS.md](../BENCHMARKS.md) 首段：“所有材料误差都在相应预测封存后离线计算”|原Pro重放脚本在 `run_case` 中已计算 direct `stats`，每50步记录 `norm(cm-truth)`，最后才写npz；[重放receipt](../experiments/pro_reproduction/REPLAY_RECEIPT.json) 明确 `original exposed toy recording, not online selection`。它没有GT选候选，但也不是先seal后首次评分。|将全称限定到新solver cohorts，另写Pro重放保留原暴露diagnostic truth logging；与paper introduction现有正确限定一致。|
|P1|[PAPER_DRAFT.md](../PAPER_DRAFT.md) §3.1：“All comparisons within a cohort use … identical … initial material”|NN四对象cohort含uniform与NN初态；独立24 `nn_direct`也不从uniform.1开始。相同的是逐数据集算子/数据/box/原物理R，不是全部arm初态或信息先验。|写“for each dataset operators/data/box/physical R are shared; initialization and learned information are declared per arm; AP/DR endpoint controls share the same complete physical initialization”。|
|P1→受限范围已补证|[BASIN_ESCAPE_REPORT_ZH.md](../BASIN_ESCAPE_REPORT_ZH.md) 首段：“控制没有找到的更好材料区域”|原报告仅有Phi-selected估计，不能推出整个visited region。审查后root授权的saved-pool oracle表明：实际保存且bodyCPU≤240的九个control池的材料oracle仍均差于对应DR online；新证据见第3节。|可加强到“这些实际保存的cap内有限pool未含同样好材料”；仍不能扩为全部未保存inner或其他restart分布从未到过。oracle不得变成在线预测。|
|P2|[BENCHMARKS.md](../BENCHMARKS.md) 确认行、[COST_LEDGER.md](../COST_LEDGER.md) 确认行、[BASIN_ESCAPE_REPORT_ZH.md](../BASIN_ESCAPE_REPORT_ZH.md) 末句仍有“尚未完成/封存、仍在运行”|[CURRENT_EXECUTION_HANDOFF](../delegated/confirmation/CURRENT_EXECUTION_HANDOFF.md) 和 [CONFIRMATION_EXECUTION_COSTS](../delegated/confirmation/CONFIRMATION_EXECUTION_COSTS.json) 已记录192确认+48dev terminal sealed、0科学执行失败及28147.796875确认childCPU。只有正式评分/整合还在进行。|区分“execution+seal完成”与“官方统计/最终科学判断待完成”。全局CPU仍可因其他缺失保持null，但不能说确认算法CPU总表不存在。|
|P2|[PAPER_DRAFT.md](../PAPER_DRAFT.md) §6 “queue experienced low active-core frequency”|已有的是离散时刻Windows counter读数，近start标签还曾跨stage race后修正；并非逐task频率积分或原因识别。|更精确写“processor telemetry reported low frequency at a recorded under-load snapshot”；保留raw counter/source/time，不能归因某电源/温度原因或据此校正算法时间。|
|P2|[PAPER_DRAFT.md](../PAPER_DRAFT.md) §5.1 shell成熟控制更准确|事实数值成立，但Pro只有200 corrector/direct600 allowance，而新增成熟控制使用最多3000 corrector；toy current族与两套Fourier适配也不同。|保留为“成熟方法可以恢复该对象，反驳全面失败/普适优越”的控制；不要写成隔离算法机制或等compute优胜。|

这些是本次所读版本的审阅记录；若root已修相应句子，应标resolved而不是回改旧产物或重跑科学。

## 3. basinescape与controls究竟比较了什么

### 支持的部分

[basin scores](../RESULTS/basin_remote/OFFLINE_SCORES.json) 的240-body-CPU图与数值一致：DR NRMSE环/壳/桥为0.165853/0.113775/0.106175，direct为0.422880/0.515137/0.436729。环/壳不是由换初态、GT初始化或另一个算法好端点制造；已保存原坏材料重新做完整物理lift。固定AP/DR投影和omega1的第一更新恒等式也排除了首步初始化偷换。

局部环路径的Phi先升后降，以及原stationary endpoint的direct续跑停留，支持有限的非单调路径干预。没有必要把它包装成连续路径barrier高度、所有局部极小被消灭或全局收敛证明。不同BLAS的有限非凸轨迹不同并不自动推翻这条有限证据。

### 必须保留的限制

- multistart只覆盖事前固定的smooth Gaussian RMS1/4/10、24CPU/restart分布。不能将它写成所有restart、所有初始化或全部成熟CIE的下界。
- 四法得到同一预算**上限**，没有都花满同样实际CPU：环direct完整只0.21875CPU，壳0.0625CPU；它们自然停止。比较的是给定solver在该上限内的已报告估计，不是相同objective evaluation数/相同工作量。
- v1的projection archive在两阶段之间写出。远端DR环/壳full为263.09375/261.78125bodyCPU，不能称240CPU端点。已公布240前缀是合法已付候选；其NRMSE分别0.165853/0.113775，full分别0.144608/0.112089。两组不能混用残差/图像形成一个所谓同预算结果。
- [预算更正](../audit/BASIN_BUDGET_ACCOUNTING.md) 已诚实保留body vs whole-child区别。独立24采用whole-child，不追溯赋予basin旧时钟同等语义。
- AP/DR同A/C能隔离反射更新；DR与CIE还改变current族、辅助目标、材料动作及费用。不能将所有差异都归因反射。
- 环/壳“近受约束局部极小数值依据”不是严格连续Maxwell局部最优证明；桥尤其不能并列称同样的local-minimum escape。

### 新增的已授权0FW pool核查

root接受上述量词审查后，明确授权对已封存 `prediction.npz/candidates` 与history逐行配对，离线读原R10 truth。未重放算法，未读取projection全场；12/12池行号完全对应，online最小Phi材料与原sealed240前缀逐元素相同。只纳入bodyCPU≤240，完整超时候选仍排除。

|对象|control池中最小NRMSE|DR online / DR池oracle|控制selection确实遗漏的较好点|
|---|---:|---:|---|
|环|.374538（multistart）|.165853 / .165853|multistart online .748880，而池oracle .374538，仍差于DR|
|壳|.452016（multistart）|.113775 / .113775|multistart online .456326，池oracle .452016|
|桥|.334774（AP）|.106175 / .106175|direct/multistart online约.43673，池oracle .336028；AP oracle .334774|

因此已记录有限池支持DR在这些指定起点的访问优势，不只是control最终Phi选择偏差；同时control pool确有更好材料被Phi遗漏，不能完全消去selection效应。三例DR online恰等于各自记录池材料oracle是该次有限事实，不能当新的GT-aware方法。

证据：[oracle报告](../RESULTS/basin_remote/SAVED_POOL_ORACLE_DIAGNOSTIC_ZH.md)、[逐候选范围/online对照](../RESULTS/basin_remote/SAVED_POOL_ORACLE_DIAGNOSTIC.json)、[CPU前缀](../RESULTS/basin_remote/SAVED_POOL_CPU_PREFIXES.json)、[best@K](../RESULTS/basin_remote/SAVED_POOL_BEST_AT_K.json)、[0FW receipt](../RESULTS/basin_remote/SAVED_POOL_DIAGNOSTIC_RECEIPT.json)。K包含重复已付callback，不是独立draw数；未保存的inner或未审计shadow没有重建。

## 4. 原Phi、corrector与secondary：证据比单一叙事更复杂

high3D的 [保存状态分解](../delegated/vector_high/RAW_CORRECTOR_DECOMPOSITION.json) 是强且可复核的有限反例：

|matched DR|raw→corrected D|raw→corrected R|raw→corrected材料relativeL2|raw→corrected物理残差/δ|
|---|---|---|---|---|
|chi8|1.02102e-4→1.35767e-4|7.17317e-4→5.91582e-4|.061631→.129783|.990166→1.143895|
|mixed2/8|1.00340e-4→1.30326e-4|3.58556e-4→2.45045e-4|.078315→.227239|.983948→1.125702|

这里是R下降压过D上升的soft-objective tradeoff，不是“拟合更好所以材料变差”。但mixed AP又有D/R**同时下降**、材料.216416→.227260的反例；Austria fullDual则D大降、R升、材料变差。不能把全部失败归到同一种prior权重方向或反射不会找路。

Primary守原Phi使历史比较可解释，却不把Phi变成GT证书。新noise-feasible minR只能在已有pool中选，不是可行集上全局求解。其注册记录为24数据生成后、确认inverse0及GT评分0；不能说所有设计都在生成前冻结。

不同方法pool覆盖不一致：mature auxiliary常仅physical stage endpoints，DR每50步物理审计，corrector通常记录每个已付value-gradient。即使八方法使用同一secondary规则，收益仍是“search+pool coverage+selector”的组合。须报告availability及共同可用配对子集；把144个available当144个success，或将available-only均值与primary全部24比较，都会造成实质统计错误。

raw优势也不能按事后时间免费截断：high3D DR raw选中2000/2400/50步，对应14.27/24.91/.57CPU，但实际完整proposal花34.73/40.74/35.97CPU。raw-only CPU ceilings只是已可取得候选的视图，不等于已冻结early-stop算法。

## 5. 噪声球、模型误差与实测资格

- synthetic matched24/NN/high3D的raw iid CN99有明确sigma合同。Austria10%**固定范数**噪声的.12球、Fresnel.15模型差异球均不是同一种CN99统计证书。根外部报告已经正确区分，paper/abstract不可重新混成“noise-certified recovery”。
- Austria原150²→32²，coarse GT残差.2252565在.12之外；不同拟合材料能进入.12，因此既不能把可行性解释成正确材料，也不能说整个coarse可行集为空。matched32 GT残差.0995131入球但DR仍失败，排除了“全部只是网格差异”的单一解释。
- high3D同几何8³/6³ clean差异.261567远超1%噪声。它说明这对离散模型差异很大，不是从一个配对证明整个方法不可能或连续PDE模型错误。
- Fresnel留出预测近.0915/.0986是实际预测证据；无pixelGT、绝对幅值/喇叭响应/协方差未分离，峰值chi9–13不应写成材料介电常数正确。已有空场和圆柱参考误差诊断不替代实测材料资格。
- Austria64²旧metadata已被修正为150²，原合同作为历史保留。公开机器可读入口须同时指向amendment，避免读者仅取旧CONTRACT又得到64²。
- A上数据球可达不等于A∩C物理可行。某一次小gap也不等于physical fit；core迭代0的gap占位已在本次**仅docs**修正中明确排除。完整条件与rank容差见 [最终理论notes](FINAL_THEORY_REVIEW_NOTES.md)。

## 6. NN prior与DPCS归因

冻结NN输出与历史public initializer完全一致、无新训练是有用的provenance证据。四对象NN-DR相对NN-direct的均值改善主要由0022驱动；0001/0004轻微变差，0027几乎相同。没有独立泛化结论，更没有“任意错误prior可修复”证据。

尤其不能从0027原CIE坏材料直接推出DPCS单一责任：NN-CIE保留 `initialize_physical_current=False`，起始current依原DPCS/零Fourier系数；NN-DR则从完整物理NN current/field出发，后续材料动作及current族不同。本组没有NN→AP，也没有仅切换DPCS而保持所有其他动作/完整物理warm current的单因素比较。安全结论是**这些指定接口的观测差异**，不是已经定位DPCS不可逆信息丢失机制。

NN含监督训练信息；不能用相同physical R一词掩盖它与无NN控制的信息先验差异。独立24测试的是NN→direct，未测试NN→DR在新24分布的效果；开发的NN组合不能被确认NN-direct成绩追溯验证。

## 7. true3D的分辨率与代表性

完整3D向量coupling、cell-major三分量和真实216材料优化是实际实现，不是spectral proxy。但是样本与物理资格有限：

- moderate为两个指定DEV几何，high为两个seeded几何；各自另有一个同几何cross-grid配对，不能每组三观测都计三独立对象。
- high8的内部cell相位约1.885/1.414，mixed coarse实际峰值7.80691。6³/8³不是收敛网格序列，不足以证明continuum精度或高index工业规模资格。
- moderate短corrector的AP/DR材料优势在同一较长共同预算后消失；high matched最终四法也近乎一致，DR总费用明显高。不能只展示raw或短corrector图支持“最终3D优势”。
- inherited projected-IFFT range SVD白化、PR/current adaptations与本文box wrapper必须披露；不是literal2019 reproduction。

## 8. 成熟baseline与原创性

成熟toy控制保留原PR/读出但有明确10² mode schedule、box、dense-G action及no-prior适配，两套均保留，这是合理的受限control。不能称原论文作者实现的完整复现，或以固定3000 auxiliary步骤宣称与3000projection等费用。

完整Fresnel/Austria Dual追加控制纠正了初轮300总PR的欠预算；这些追加是在旧GT结果暴露后请求的诊断，保留原负结果且不当新确认是正确做法。literal DPCS15和R9 noise-adaptive policy属于不同冻结合同，不能看GT后替换以救结果。

原理层的DR/TLS/非物理中间态/Schur消元不是新数学。DCTMC主要primary文本仅检查abstract/introduction；PADMM-CIE/EM-ADMM部分全文访问未完成，也没有WRI/PADMM-CIE数值qualified baseline。现有证据支持一个application-level controlled study及机制/selection负结果，不足以宣布最先提出非物理搜索或优于这些closest methods。

如果独立24只是与Dual近似而落后实际G19，论文的贡献应聚焦“指定坏终点干预与selection边界”，不能通过反复扩大辅助目标/挑best-arm组合把其重命名为普适新机制。cascade/portfolio是新的estimator，需另行冻结与评价；当前paper已正确把它留为未证实方向。

## 9. 计时、优先级与统计解读

外部GetProcessTimes优先于重叠inner timers，body/import/archives和不同维度LU/RHS分别记录，是当前正确做法。总CPU不是wall，给定CPU ceiling不是相同objective calls或flops；early numerical stop不应被当作算法偷偷少拿预算。

已保存 [EXECUTION_PRIORITY_ENVIRONMENT.json](../delegated/confirmation/EXECUTION_PRIORITY_ENVIRONMENT.json) 的一次实际快照：controller和DR/multistart/NN为Normal，示例mature baseline children为BelowNormal。代码与此一致。这是**实际记录到的调度差异**，不应写成全部priority相同；但单snapshot不代表每个task全过程，也不能定量归因速度或材料差异。单P-core affinity/CPU-time caps仍有真实计账意义，不因此失效。

频率snapshot仅反映其时刻的reported fields；结束空闲读数不是8core负载实验。既不能用它证明固定同频，也不能凭低reported MHz推出特定硬件故障。paper当前不作strict identical-scheduling wall-throughput主张是合适的。

全局总账尚有prepare/审阅/旧probe等缺失可继续null；但确认执行总额已可填28147.796875CPU，startup失败.90625与dev/transport另计，不能继续将确认执行本身写成未知或算法运行中。

正式确认统计应至少保持：逐material配对；dev6不混入confirm24；primary全24；secondary availability、共同可用配对子集与缺失分开；bootstrap按material对象而非source/pixel/重复prefix抽样；固定success只描述、不作停止；per-family与多比较/小样本限制显式披露。DR与Dual均值相近不是统计等价证明，DR/G19均值差也要与paired不确定性并列。

## 10. 确认24官方结果补审（同一正式评分，不重算GT）

已读 [官方中文报告](../delegated/confirmation/CONFIRMATION_RESULTS_ZH.md)、[aggregate](../delegated/confirmation/RESULTS/offline_scoring_v1/aggregate.json)、[pairs](../delegated/confirmation/RESULTS/offline_scoring_v1/pairs.json)、[评分源冻结副本](../delegated/confirmation/RESULTS/offline_scoring_v1/frozen_inputs/score_confirmation_source.py) 和 [scoring receipt](../delegated/confirmation/RESULTS/offline_scoring_v1/SCORING_RECEIPT.json)。评分分开dev6/confirm24、primary/secondary，按材料对象配对bootstrap10000/seed20261112；NA不补0，采用全材料roster重抽、双侧可用连续差异。该实现与披露一致；percentile95%未多重比较校正，不能包装成普遍族级显著性。

primary全24 mean amplitudeNRMSE为DR.184631、Dual.185456、实际G19.143891；旧描述成功分别7/24、9/24、9/24。官方paired结果（负号有利DR）：

|primary差异|均值ΔNRMSE|95% material-bootstrap CI|可支持解释|
|---|---:|---|---|
|DR−Dual|−.000825|[−.040982,+.053017]|未建立均值优势，也未做等价性证明|
|DR−G19|+.040740|[+.001934,+.095095]|该cohort和该未校正paired区间支持DR更差；不外推所有物理模型|

新增secondary DR可用18/24，available-only mean .098258、旧描述成功仍7/24；Dual19/24、G1920/24。不能拿.098258与primary全24的.184631声称“规则将误差减半”，因为六个没有noise-feasible池候选的对象不在前者均值内。未取得候选本身是输出，不能事后填primary或删掉。

官方 [GT_saved_noise_qualification](../delegated/confirmation/RESULTS/offline_scoring_v1/GT_saved_noise_qualification.json) 表明24个已保存clean truth均在本次CN99球内。DR有11个球内旧描述未过、6个球外；但这11个**全是NRMSE≤.15而IoU<.7**，平均NRMSE .117711、relativeL2 .324088–.600206。应称“球内但旧超水平集/细节阈值未过”，不能写成11个完全没有成像。

这些11例同时有合法GT与不同合法DR材料都满足同一有限噪声球，支持该次noise-radius相容候选的诊断，不是exact nonuniqueness或连续Maxwell不可辨识定理。具体字段见 [DR_failure_continuous_diagnosis](../delegated/confirmation/RESULTS/offline_scoring_v1/DR_failure_continuous_diagnosis.json)。

确认raw proposal minPhi均值.185814，primary .184631；逐对象NRMSE12降12升。noise状态为18 inside→inside、6 outside→outside，**没有复制high3D matched的进球→corrector出球模式**。因此本cohort并非“DR普遍先找对、再被同corrector毁掉”；原始搜索/有限pool问题仍存在。查看 [raw/final并排表](../delegated/confirmation/RESULTS/offline_scoring_v1/DR_raw_primary_secondary_comparison.json) 与 [DR all views对Dual/G19](../delegated/confirmation/RESULTS/offline_scoring_v1/DR_all_views_vs_Dual_G19.json)，不能只借high3D图为全24解释。

paper不将DR替代G19为默认solver的当前判断符合这些数据。摘要仍应把mean/CI、有限matched distribution、secondary availability与失败分开；最强贡献保留为“特定端点干预与估计器边界”，而不是确认了普适新盆地机制。

## 11. 本次实际改动与未做事项

仅在 [THEORY_DERIVATIONS.md](../THEORY_DERIVATIONS.md) 第一更新段后增加iteration0 gap/A-state占位声明；没有改core、算法、预算或旧预测。重大数学/原创性争议仍提交root裁决。

root另外明确授权的两个后处理均已完成：

- [RAW_VS_CORRECTOR.png](../delegated/vector_high/RAW_VS_CORRECTOR.png)：固定upper-central z=3、统一chi0..8、GT/原minPhi DRraw/DRcorrected/CIEcorrected，注D/R/物理CN99；不是按GT选raw。已目视修正标题重叠，两次图版成本留于 [revision costs](../delegated/vector_high/RAW_VS_CORRECTOR_REVISION_COSTS.json)。0新FW/优化。
- 第3节basin saved-pool oracle：仅已封存有限pool、原评分指标、0新FW/优化；独立oracle产物标签，未改online材料或原seal。

## 12. 最终两根报告复核与resolved列表

复核 [RESEARCH_REPORT_ZH.md](../RESEARCH_REPORT_ZH.md) 与更新后的 [PAPER_DRAFT.md](../PAPER_DRAFT.md)，并逐项对照同一官方aggregate/pairs/raw诊断。**没有发现需要推翻主科学判断的重大数值冲突**：当前“特定坏终点可干预，但不能普遍取代成熟G19；搜索与选择需要区分”的有限定位与完成证据一致。

|复核项|最终文稿与真实证据|状态|
|---|---|---|
|确认primary均值/描述成功|DR .184631/7，Dual .185456/9，G19 .143891/9，均为全24|核对一致|
|paired不确定性|DR−Dual CI跨0，不能宣称等价；DR−G19未校正percentile CI为正，明确小样本/多比较限制|核对一致|
|预算前缀|60CPU G19 .1498817、DR .1839617；30/60/120/240表与aggregate一致|核对一致；是paid-prefix，不是相同flops/严格同调度|
|17个旧描述未过|11球内、6球外；11球内均NRMSE≤.15、仅IoU<.7，relativeL2 .324–.600|两稿已保留，未写17个完全没成像|
|有限噪声相容性|官方保存Yclean核查24/24 GT在球内；11例另有不同合法DR材料同球内|只支持该次噪声半径诊断；两稿没有宣称exact nonuniqueness，最终必须继续保持|
|raw→primary|18 inside→inside、6 outside→outside，材料12改善12退化，mean .185814→.184631|中文已明确不是high3D毁raw模式的普遍2D复制；英文没有将其泛化|
|3D corrector反例|matched DR D升/R降/Phi降、材料误差升并出CN99；另有AP D/R都降仍材料变差|两稿一致，无梯度错误或单一prior归因|
|NN信息与初始化|paper已撤销“所有arm初态相同”全称，区分learned initializer与physical R|主P1已修；下面Only用词仍可更精确|
|oracle与online|新basin finite-pool oracle仍劣于DR online，未当预测或省掉未审计范围|中文已整合，范围正确|

### 已解决的审查项

- BENCHMARKS已明确Pro保留逐步GT diagnostic logging、不是新seal后首次评分；不再对所有代码作全称。
- paper §3.1已按arm区分NN初态与原physical R，撤销相同初态全称。
- basin主张已从“控制从未找到”改为公开Phi规则返回估计；新增root授权0FW pool核查后，仅加强到实际保存的cap内pool。
- root成本表已取得确认28147.796875、dev7189.6875、失败.90625，未把它们与phase/secondary重复加；全局其他未知仍null合理。
- THEORY迭代0占位声明已仅改文档；paper明确不把占位gap用于界。
- 英文计时段已改成reported under-load telemetry，并保留一次priority snapshot及非严格同调度限定。

### 本次最终读到的剩余文字修正（不要求新科学）

1. **明确计数单位：** 中文报告“各方法有144个可用/48个不可用的次候选”应为“八方法合计144/192可用、48/192不可用”；每个方法只有24，DR为18/24。这是明确的量词/分母笔误，不改变官方评分。
2. **初态限定词：** paper §3.1的“Only the saved-endpoint … share the same complete physical initialization”宜改“In particular …”：high3D AP/DR也共享完整物理初态。已修的核心意思是NN与uniform不可全称相同，不必反向否定其他matched AP/DR。
3. **中文遥测措辞：** “实测活跃核约1.5GHz”宜与英文一致写为“一个under-load遥测样本报告约1.5GHz”，避免将离散counter样本升级为全过程silicon frequency测量。
4. **滚动状态：** BENCHMARKS等局部稿仍有“当前评分中”的旧句。正式评分已完成，应在发布前统一为“192 sealed+官方评分完成，科学审查/归档状态另列”。不要继续把已取得的确认执行费用写unknown；也不要因评分完成就说近邻benchmark或publication验收已完成。

以上文稿改动若由root完成，可直接将本表标resolved；不修改已封存预测、selector规则或官方统计。当前没有依据宣称全局凸化、exact材料不可辨识、普适DR优势、连续高index3D精度，或新增secondary已解决确认失败。论文接受/最终novelty仍是root保留判断。

## Root 2026-11-12恢复后的验收

已回读根文稿：144/192计数、In particular初态限定、遥测样本、确认评分完成状态均修正。root独立从192逐对象记录重新求均值/成功数/成本，与官方aggregate一致（ROOT_CONFIRMATION_READBACK_20261112.json）。高3D图和oracle产物复核完整。原审查记录不删除；本段记录处置，非新增科学。G19/A+为项目已有方法，不统称成熟文献算法。

## 用户要求的结构判读补充复核

2026-11-12：root亲看四家族全部主图，新增S/P/F事后描述，保留旧统计。复核者检查192标签计数与0009双尺度、0007、0013图，无重大越界。root采纳限定：0009是强对比主干恢复、弱分支部分恢复；0007G19骨架增量仍有尺寸/幅度偏差；0013相对CIE的五岛恢复并非DR独有。18/24始终是单人事后S或P对象数，不称新的确认成功率。原始定量结果、在线选择和算法预算均未改变。

## 综合评价复核

COMPREHENSIVE_EVALUATION_ZH.md的坏端点、NN三臂均值与bodyCPU、G19配对区间及30/60CPU前缀均对照已保存正式结果。未见实质数值/归因错误。应保留“已记录bodyCPU”，不把历史训练/共享预处理/未计导入自动算入。NN0022 CIE也恢复、0027仅本接口严重退化、缺NN→AP/current初态消融、G19区间未多比较校正、S/P事后非盲等限定均保留。root接受此复核；不新增实验。
