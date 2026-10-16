# LXZ proxy sampler：只读概率正确性审阅

2026-10-16。读取`src/lxz/proxy.py`、冻结LXZ protocol以及P2 material/response/likelihood/C1 proposal实现。没有执行sampler、controls、geometry检查或物理任务，没有改源代码。结论针对静态数学形式和声明的proxy numerical target，不是有限链收敛证明。

## 总体结论

**在修复Jacobian paired-index后，没有发现当前local/global/hybrid/tempered接受率公式改变声明target的错误。** 原geometry结果不能支持Jacobian/Fisher机制解释；父会话已保留旧数组、单独计费重跑geometry repair，不以静默覆盖掩盖。有限链、mode discovery、CDF、成本标签仍需以下限制。

审阅快照行号：prior mapping31–43；geometry61–80；normalized global density82–105；MCMC132–163；SMC165–185；保存/计费217–260。

## 1. 原prior的Gaussian坐标与frozen pCN：形式正确

原prior在bounded `x` 上uniform。令`x_j=a_j+w_j Φ(u_j)`，Jacobian为`∏w_j φ(u_j)`，故`p0(x(u))|dx/du|=φ4(u)`。后验在u坐标为`π(u)∝exp(ell(x(u)))φ4(u)`；不是说材料后验Gaussian，也不是使用C1的logit mixture坐标作standard-Gaussian prior。

普通local使用`u′=sqrt(1−β²)u+βη`。Frozen geometry用`R=V diag(sqrt(1−b_j²)) Vᵀ`、`B=V diag(b_j)Vᵀ`；R对称且`RRᵀ+BBᵀ=I`。于是standard Gaussian下联合`(u,u′)`的block covariance在交换两者后不变，proposal对原prior可逆。该kernel接受率仅需`min(1,exp(ell′−ell))`；target prior/Jacobian确实被prior-reversible proposal抵消，不是漏项。

Geometry只在链开始、使用本次观测及64-point初始化中选择的状态构造，之后冻结。依赖初始化/观测使初始law不是stationary，却不破坏每条固定kernel的target不变性。25% burnin不足以构成stationarity证明。若今后逐状态更新R/B而不加proposal density ratio，这个结论将失效。

## 2. 已发现并由父会话修复的geometry索引错误

原rows顺序是`u−he0,u+he0,u−he1,u+he1,…`，原列差误用了`f[2*j+1]−f[2*j−1]`：j0减去最后一行，其余混相邻plus方向。当前行71已改为`f[2*j+1]−f[2*j]`，与中心差分顺序一致。

这会改变A/H/eigenvectors/rates和混合效率，原结果不能叫真实Jacobian/Fisher优化。即便使用错误H，R/B仍是frozen symmetric prior-preserving矩阵，因此这主要是**物理几何解释/算法实现错误，不是已证明的posterior invariant-target偏差**。修复后表现仍需实际repair数据，不能由静态检查预设更好。

Round2的gain Schur项使用local untruncated Gaussian-gain近似，anchor response由对称扰动均值近似。它只是proposal metric，不是proper bounded-gain marginal target的精确Hessian/Fisher。代码最终target仍用带logdet与截断归一化的likelihood，因此该近似本身不改posterior。

## 3. C1全局proposal与MH ratio：形式正确

`q_x=.9 q_MDN(x|y)+.1 p0(x)`。MDN使用完整四维Cholesky Gaussian log density、exact bounded-logit Jacobian，并通过现有helper校正global categorical export roundoff；constant中的`−2log(2π)`正确对应D=4。Sampler的prior分支与MDN分支混合生成law和logpdf一致。

在当前uniform x-prior下，independence MH接受率为`min(1,exp(ell′−ell+logq(x)−logq(x′)))`，当前offset方向正确。独立prior proposal则q=p0，接受率只用likelihood ratio。这里logq是**dx density**；如果未来改成radius/volume prior或换支持，不能继续省略prior ratio。

Hybrid的global概率固定.25，两个component kernels各自π-invariant，固定凸组合仍π-invariant；不需要额外selector ratio。若selector变成状态/acceptance-history依赖，必须重新证明或补齐选择概率比；当前finite-state controls里的负例不能自动覆盖未来自适应版本。

所有方法甚至local/prior都预生成一批q draws，只计wall不调physics。丢弃未用draws不改变有效global proposals的law，但构成强制额外setup，应披露其对CPU公平比较的影响，不能说local纯算法本来需要C1 acquisition。

## 4. Tempered SMC：权重与mutation形式正确，finite-N限制保留

目标序列`πβ∝p0 exp(β ell)`。当前权重增量`(βnew−βold)ell`、归一化、ESS触发stratified resampling，以及prior-reversible local proposal的`βnew(ell′−ell)`接受率，均与该序列一致。未resample时保留旧weights再做target-invariant mutation是合法SMC步骤，不应因mutation改变状态就重新乘似然。

64个prior初始点加63轮64粒子mutation恰为4096 target evaluations；β在第6轮到1，checkpoint1024及4096都是β=1。这个schedule/粒子数可能漏mode；finite-N normalized event估计不无偏、不保证全mode coverage，历史importance ESS不等于MCMC ESS或结构发现证书。

## 5. Budget与成本：数量基本一致，异常和shared setup需单列

每arm自己的64个prior likelihood重复计费；geometry再计8个derivative response/likelihood，剩余数补足4096。Checkpoint按包含setup的target.calls触发，因此不同arm保留的chain states略不同，属于同调用预算而非同样本数预算。Reject states全部保留，没有accept-only selection偏差。

纯Response使用已付费bank，无forward solver。`MaxwellTarget.__call__`先计attempt，异常加failures；geometry_setup先计8但没有相同try/failure计数。且COST_LEDGER只在整arm完成后追加，异常中断arm若没有外层failure receipt，其已用proxy attempts可能不在completed-row总和中。当前完成结果不能抹去旧坏geometry成本；父会话说明repair两轮另计131072 proxy attempts、full physics0。

Per-arm `inclusive_seconds`不含每observation共享的C1 cold load，后者在PROPOSAL metadata和整轮COMPLETE wall中另列。主比较须披露cold/amortized model/teacher/setup成本；不能仅以arm内部时间宣称端到端快。当前代码只做per-arm硬上限，全stage wall/proxy保护须由父会话统一账本/进程控制负责。

## 6. 指标标签：不改target，但可能误导科学解读

- `significant_event_recall`用reference≥.05、model≥.005；它等于当前false-exclusion补集式coverage，**不是P2 model≥.01的mass-aware discovery指标**。跨报告比较必须显式区分阈值。
- `false_confidence`目前是`max(model_mass)−max(reference_mass)`的signed peak excess，不是P2的“max model≥.90且reference至少两重要事件”boolean flag。应在报告叫peak-probability excess，或另算原flag，不能直接解释成错误自信发生率。
- `unique_material_geometry`实际对`[r,h,ε′,ε″]`四坐标取唯一数，包含ε变化；它不是geometry-only数量、更不是分离mode数量。
- Event/continuous ESS遇到constant chain会返回None，不能把None删掉后仅平均其余大ESS，暗示所有事件均混合。IPS FFT估计也是finite-chain诊断，没有R-hat/多链全局mode证明。
- 四观测均已暴露，proposal未用当前reference拟合，但C1 teacher/setup是历史额外监督。它不是zero-supervision传统sampling的相同训练信息比较。
- CDF指标是保存grid上的最大差，对PILOT C连续积分作比较；不是continuous supremum。参照mean/SD是physical coordinates，略取整prior widths只影响报告scale，不改变target，但应避免声称精确相对误差认证。
- 当前MCMC从64-point best likelihood开始，SMC保留64点集合；同初始化信息/调用预算不等于同初始distribution，early checkpoint差异包含warm-start影响。四seed和finite误差必须分列。

## 最终审阅意见

可继续以修复后的kernel研究**同一已声明proxy target**的有限预算误差与mixing，不能由这次静态审阅、toy controls或高ESS升级成完整posterior/continuum/无漏mode保证。E3若要求两个分离材料modes，现有bins/事件质量/多图不足以满足；必须另有mode/路径barrier证据。未发现新的必然MH/SMC目标偏差；仍需保留上述geometry失败、有限链偏差风险、reference误差与成本/指标语义限制。

