# Research B：补偿感知主动 Maxwell 测量的实际实验

2026-10-15 · 独立Research B分支 · CPU · 原始A/P2资产只读

[在线入口](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_b/) · [完整证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-b-20261015) · [数据与费用索引](START_HERE.md)


本轮结论：**继续研究候选补全与新观测后的corrector；停止在当前径向菜单上扩大复杂选测策略。**

1. **仍难分的世界：**小中心孔与实心可沿合法介电常数和gain补偿保持接近；固定625池本身也遗漏大量接近truth的材料。三步后多数有限兼容集合仍跨结构family。
2. **有用的新测量：**在本台架中，高频尤其3.2通常比repeat1.6、2.4或单独校准更有材料选择价值；这不是任意目标/硬件的频率推荐。校准能通过限制gain间接提供材料信息。
3. **新策略的优势：**32-parent比较没有证实robust、材料风险或EIG胜过强固定3.2。普通Fisher首步与固定3.2重合。四个真实A交接案例中robust比固定3.2好，但差异主要来自一个parent，不能升级为总体优势。
4. **MAE是否下降：**固定池三次3.2使Selected-MAE从0.162279降至0.090885；真实A四例加入连续先验多起点corrector后从0.125480降至0.036170。这两组不是同一人群；后者混合了新信息与额外重建计算，不能单独归因于选测。
5. **真正共享补偿贡献：**full-Mie端点支持部分跨频gain冲突，但很多对照冲突很小。不能把全部频谱收益都归给gain共享；原raw-conflict指标倾向校准且材料结果较差。
6. **遗漏candidate的失败：**删去中心孔family后29/30相关noise episodes未报警。低残差或唯一剩余候选都不能保证没有漏世界。
7. **是否继续测量：**不为消耗余额继续campaign；对仍难分的小孔保留歧义。当前优先改善候选支持，并保留简单固定策略作强对照。
8. **潜在TAP贡献：**可复算的跨频补偿见证、候选缺漏失败机制和实际A→B→corrector闭环是研究素材。复杂策略收益、一般非同心3D、真实仪器优势与强原创性尚未建立。


## 1. 从物理问题到可检验任务

一个散射观测可以由几种不同的材料世界解释。新增频率的价值不只在散射变化大，而在候选还能否使用**自身同一组合法材料/校准参数**解释全部历史。不同世界允许不同gain；同一世界不允许每频任意重置gain。我们先用已保存响应排序，再计算新频率响应，在独立模拟观测后更新候选，最后按真实材料而非熵评分。

范围是三维矢量Mie的**同心、径向、各向同性球体台架**，不是任意3D成像。R=1；四维连续系数prior为两个几何系数[-3,1]、Re ε∈[2.2,2.8]、Im ε∈[.04,.16]上的均匀测度。材料decoder自然含实心、中心孔、环形孔及组合；inactive几何对应真实材料原子，积分保留其prior质量。topology标签仅在事后分析，不输入候选生成/优化。

观测为16复通道，packing为Re16后Im16。每frequency采用固定solid锚点的12% RMS噪声，不能由truth/candidate大小调整。每个episode的一个实gain来自截断N(1,.15²)、支持[.7,1.3]，跨历史持续不变。不同episode独立gain。宽带ε常数是本台架假设，未经真实色散验证。

材料主指标是R=1球内体积加权的Re/Im相对介电常数**分量L1和**，s_ref=1；在两材料壳层边界并集上精确积分。它对极小孔体积不敏感，因此另报事件/结构选择，不能用MAE小代替正确结构。

## 2. 两种信息语义与算法

Mode 1输入full-validated有限材料集合，不赋posterior权重。每个材料在完整历史上profile一个有界共享gain，用最小joint χ²选材料；raw/robust使用full缓存响应。robust在当前best+9兼容材料中最多48个筛选点、材料L1≥.02的所有跨/同事件对上，比较合法gain区间内最难区分的未来白化距离。它是有限集合启发式，**不是严格连续minimax regret**。整个625池始终保存。

Mode 2另用原continuous prior、原likelihood和响应proxy的有限条件近似（2048样本、128预测draw、96初始decision材料）计算EIG/材料风险。正确更新为joint-history likelihood减old-history likelihood，gain只积分一次；不得独立乘两个已边缘化gain的likelihood。Fisher先堆叠旧新Jacobian，再消去同一nuisance。采用参考挑action再对固定池profile的流水线仍是**reference-assisted**，不能称纯candidate-only。

没有训练NN，没有ROM。A17的实际墙钟与A23/A24失败不支持默认加速；本轮直接小型full solver即可。候选数量、profile residual、design prior和概率后验分别标注，绝不转换成residual softmax posterior。

## 3. E0：真正缓存排序，新增Maxwell=0

11项源包理论/概率检查通过后，对七个已暴露P2案例运行了实际排名与新独立噪声模拟。原1.6已在posterior中，不能作为新证据；repeat使用新noise ID。原2.4未用于当前C参考，所以可加入。2.4极化/接收子集、analytic校准与no-measurement也列入inventory；缺缓存动作明确不可用，未偷偷补求解。

| 旧暴露案例 | repeat1.6 预测风险下降 | 2.4 预测风险下降 | 校准 预测风险下降 |
| --- | --- | --- | --- |
| solid | 0.016637 | 0.032561 | 0.019930 |
| small_cavity | 0.019053 | 0.029510 | 0.022591 |
| large_cavity | 0.014245 | 0.030506 | 0.013998 |
| annular_cavity | 0.045243 | 0.065187 | 0.019764 |
| combined | 0.040926 | 0.063276 | 0.015613 |
| prior_boundary | 0.017355 | 0.027318 | 0.013688 |
| gain_negative_control | 0.033116 | 0.054580 | 0.013884 |

这是有限连续proxy参考的预测风险下降，不是实测MAE。七例中raw、Fisher、joint Fisher、EIG、材料风险和robust均选2.4，conflict均选校准。Mode1七例的单次realized平均MAE却由0.124276变为0.132117；暴露小面板不证明主动测量获益，也不据此删难例。E0随机基线共享同一菜单draw、uniform一步等于预设2.4，均保留限制。

## 4. E1/E2：受控新物理与最多三步闭环

先冻结六频菜单{1.2,1.6 repeat,2.0,2.4,2.8,3.2}与方法家族，再生成32个新parent，每parent两份噪声。新增四频full bank、独立审计、truth与noise锚共3016次attempt。初始13/17参考事件差最大.000277；这只支持初始观测。

E1是一新增block；E2复用同32对象，是**暴露后续分析**。E2最多三block，uniform固定1.2→2.4→3.2；random在六频均匀抽取；固定六频全部报告。外部void-volume几何传感器SD=.01、校准各计一个新block，不免费，也不假装来自原散射。测量block只是相同dwell/并行通道的归一化成本，没有硬件秒/能量报价。

| 策略 | E1 固定625池选择 MAE | E2 固定625池选择 MAE | E2 有限参考辅助选择 MAE† | E2 新增blocks均值 |
| --- | --- | --- | --- | --- |
| none | 0.162279 | 0.162279 | 0.104046 | 0.000000 |
| fixed_1.2 | 0.149763 | 0.140054 | 0.085696 | 3.000000 |
| fixed_1.6 | 0.143048 | 0.133089 | 0.066354 | 3.000000 |
| fixed_2 | 0.150269 | 0.114810 | 0.062181 | 3.000000 |
| fixed_2.4 | 0.122534 | 0.106367 | 0.059581 | 3.000000 |
| fixed_2.8 | 0.123681 | 0.107033 | 0.052725 | 3.000000 |
| fixed_3.2 | 0.101982 | 0.090885 | 0.039812 | 3.000000 |
| random | 0.133483 | 0.105370 | 0.055309 | 3.000000 |
| uniform | 0.149763 | 0.097503 | 0.053684 | 3.000000 |
| raw | 0.102027 | 0.093068 | 0.043414 | 3.000000 |
| classical_Fisher | 0.101982 | 0.090885 | 0.039812 | 3.000000 |
| joint_projected_Fisher | 0.101982 | 0.089865 | 0.040140 | 3.000000 |
| event_EIG | 0.104847 | 0.090107 | 0.042696 | 3.000000 |
| expected_material_risk | 0.101982 | 0.093775 | 0.040103 | 3.000000 |
| robust | 0.103714 | 0.093842 | 0.044747 | 3.000000 |
| conflict_diagnostic | 0.141574 | 0.113446 | 0.049028 | 3.000000 |
| calibration_only | 0.141574 | 0.134254 | 0.084007 | 3.000000 |
| calibrate_then_risk | 0.141574 | 0.094999 | 0.038684 | 3.000000 |
| geometry_only | 0.130434 | 0.125114 | 0.070882 | 3.000000 |
| risk_cost_stop | 未运行 | 0.093687 | 0.040103 | 2.984375 |
| robust_abstain | 未运行 | 0.094046 | 0.045372 | 2.984375 |
| two_step_bundle | 未运行 | 0.095791 | 0.047853 | 2.000000 |

†这些列是已执行算法的实际truth评分，但参考分布/效用数值并未得到完整收敛保证。`mode1_candidate_only`是输出选择池的字段名，EIG/材料风险行的action仍需要Mode2参考。

E2 risk-cost-stop平均2.984375 blocks，几乎未节省测量；二步calibration bundle有明确成本，但不能以两步结果和三步结果直接宣称等成本获胜。固定3.2是E1观察到的强描述性基线，E2前保留，不能冒称事先未知的独立最优频率。

按parent先平均两噪声，再2000次固定seed26101569 bootstrap，差为策略−baseline，负数较好：

| 选择池/信息来源 | 策略 | MAE差：策略−固定3.2 | parent配对95%区间 |
| --- | --- | --- | --- |
| mode1_candidate_only | event_EIG | -0.000778 | [-0.005217, +0.003528] |
| mode1_candidate_only | expected_material_risk | +0.002890 | [-0.000581, +0.007754] |
| mode1_candidate_only | robust | +0.002957 | [-0.002879, +0.010341] |
| mode2_reference_selection | event_EIG | +0.002884 | [+0.000838, +0.005081] |
| mode2_reference_selection | expected_material_risk | +0.000291 | [-0.001564, +0.002216] |
| mode2_reference_selection | robust | +0.004935 | [+0.000635, +0.010134] |

有限32parent和多重描述性比较不支持总体支配。共享prediction setup与portfolio计时不能换算成各策略独立运行速度。所有逐步/逐parent的结构、regret、coverage与成本在CSV/Release，未只保留最优案例。

## 5. 参考数值问题：不能称完整full posterior

固定四parent×两噪声、E1/E2和两政策标签对应16条唯一history。13/17比较中12条未稳，事件差最大.090781；2048样本对17阶事件差最大.186853、已选candidate期望L1风险差.016276。17阶本身也未全部收敛。

一次独立physical-coordinate自适应积分保留同prior Jacobian、inactive原子与shared gain：16条全部完成，只有2条达到原relative2e-6规则，14条在512 cells后仍UNRESOLVED。adaptive对17阶最大事件差.075959、对2048样本.186789。重复频率inverse-variance压缩与原目标身份误差约2e-14。没有放宽精度、21阶/全案例追加或重选获胜policy。

因此，本轮Mode2是独立有限proxy比较，**不是可信完整Maxwell posterior，也不是严格上界**。实际full-Mie truth和固定池选择的MAE仍可复算，负的复杂策略比较不依赖把这些参考概率当真。没有用未收敛的概率数值宣布概率校准成功。

## 6. 补偿机制与应该保留的歧义

12个固定small/large/annular结构、四新频的bounded ε/g joint优化端点经full-Mie验证。小孔r=.08的joint白化距离约.039–.142；large r=.25在2.8/3.2约3.225/3.176；ring旧频已很容易分，不能当作破解原歧义的成功。每个优化值只是best-found上界。大孔3.2同ε端点gain冲突平方约.589；许多其他组接近零。校准本身direct material derivative为零，却通过限制gain给出joint Fisher信息，说明补偿分析有物理意义，但单一冲突分数并非材料任务收益。

又用6次full调用核对**同一个**小孔ε/g世界在全菜单的近似相容性，避免拿不同频率的不同近对证明序贯不可能。对固定两个Gaussian世界及最多三次新增block，链式KL/Pinsker的等两世界设计先验分类错误率下界约.439742。它不是posterior质量或全材料空间定理；只是固定数值物理/噪声/菜单中的可检查见证。完整假设与字段见[SHARED_COMPENSATION_TEST.md](SHARED_COMPENSATION_TEST.md)。

## 7. 候选不完整与实际保留的世界

625池Best-of-K平均.062224，Coverage@K(MAE≤.05)=.40625；只重加权或选择，oracle下限不变。删central/annular/combined支持后重新排名，robust动作分别改变33/13/9个noise episodes。删annular后robust MAE约.28291。central漏解29/30、combined漏解15/16相关episodes未触发残差报警，sentinel不能作coverage证明。

三步后固定3.2/robust/risk的启发式兼容K范围为2–27/1–35/3–29；192条episode-policy中160条仍跨family。这不是192个独立对象。完整625材料、共享gain、profile χ²、history与mask已导出到`REMAINING_WORLDS.npz`，概率为none，未知支持质量为null；不只交一个点估计。

## 8. E3：真实Research A包与corrector，已执行

按lexicographic预先固定confirm_00..03，27个真实full-validated候选。原native schema、execution profile与B schema用version adapter衔接，保留材料数组/单位/solver/observation identity；gain-applied旧prediction先恢复未乘gain的F，防止重复乘。B自己扩展已授权菜单，不改A目录。新噪声与truth只供模拟/评分。

135候选未来响应、24truth、54保护几何ε-corrector端点、24三个连续prior起点的独立4D corrector端点、3参考端点核对，共240次attempt。新材料必须同时通过原1.6与所选新频full validation和原A兼容阈值才安装；原候选始终保留。12个prior起点均产生可用端点，其中2次optimizer未达其收敛标志，保留失败状态；通过full/data检查不等于找到优化全局最小。

| 真实A交接方法 | Selected-MAE | Best-of-K MAE | Coverage@K ≤.05 | 平均K |
| --- | --- | --- | --- | --- |
| A_base_only_no_measurement | 0.125480 | 0.125480 | 0% | 1.00 |
| A_fixed_candidate_profile_no_measurement | 0.125480 | 0.110881 | 0% | 6.75 |
| A_robust_action_then_profile | 0.119000 | 0.110881 | 0% | 6.75 |
| A_fixed3.2_then_profile | 0.146238 | 0.110881 | 0% | 6.75 |
| A_reference_action_then_profile | 0.146238 | 0.110881 | 0% | 6.75 |
| A_robust_then_protected_eps_corrector | 0.064836 | 0.057371 | 50% | 11.75 |
| A_robust_then_independent_prior_multistart_corrector | 0.036170 | 0.031589 | 100% | 9.75 |
| declared_continuous_P2_reference_comparator | 0.027106 | 0.019557 | 100% | 96.00 |
| A_base_only_same_action_eps_corrector | 0.089090 | 0.078574 | 0% | 2.00 |

参考comparator†未作完整双频端点验证、不安装A、参考UNRESOLVED。`A_base_only_same_action_eps_corrector`复用base对应的已付费端点，是同action的候选容量对照；选action仍看完整A集合，不能当作纯base主动策略。

robust选3.2/2.8/3.2/2.8；固定3.2的四例均值较差主要由confirm_03贡献。只有四例，不作统计胜利。protected安装20个、prior-multistart安装12个新fullvalidated端点；二者分别与原池合并，不把不同更新结果混成同一方法。新增ε/geometry候选可以改善Best-of-K，固定池选择不可以。prior-multistart的Selected结构只对2/4；小材料MAE不代表全部结构正确。

尚缺“同容量corrector但不新增测量”的四例对照，故改善不能完全拆成信息/优化/初始化收益。本轮不为补故事继续花费剩余额度，也不把B独立corrector冒称Research A已接收并重跑。

## 9. 费用、失败与科研去向

新full-Mie尝试3284/4096，失败物理调用0，GPU0；剩余812不消费。E1/E2受控wall 370.413/7200秒；全部已入账数值任务CPU 373.482秒。最终值见持续追加账本与[RESOURCE_SUMMARY.md](RESOURCE_SUMMARY.md)，文档/上传时间不冒充物理执行时间。一次adaptive报告字段KeyError计为失败job，45秒已算入，导出修复不重算积分。

32-parent策略结论支持“新增测量有用”，不支持“复杂选测更优”；四个真实A案例支持候选补全值得做更公平的后续试验，尚未证明普遍优越。下一步优先有界候选发现/修正与保持歧义；当前径向任务不继续扩展复杂policy或ROM。

## 10. 文献与新颖性边界

源包已存在主动电磁/Bayesian design近邻，本次只定向核对五篇，没有重做宽泛综述。Fouda/Teixeira的[作者全文](https://arxiv.org/html/1401.1092)已有adaptive sensing与uniform比较；Sandhu等的[正式论文](https://doi.org/10.1109/ACCESS.2023.3288997)讨论EM imaging EIG设计；Bartuska等的[作者摘要](https://arxiv.org/abs/2112.06794)已有nuisance-aware design；Naghshvar/Javidi的[作者摘要](https://arxiv.org/abs/1203.4626)已有active sequential hypothesis testing与成本。摘要阅读不冒充全文复核。

标准Schur、Bayes risk、复合假设、KL链式法则和Pinsker不是本轮新定理。最终新颖性结论保持未建立；可再研究的是Maxwell-specific补偿和漏候选机制是否在更广问题中产生真实、付费值得的材料收益。历史A17/A23/A24/A30的负结果没有被本实验覆盖或推翻。
