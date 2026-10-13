# B1 × C1：保留、解释并传递 Maxwell 材料歧义

实验阶段编号20261013；公开交付日期2026-10-14。研究支持范围为同心分层球的三维矢量Maxwell校准台，材料由四维连续prior生成。代码、原始概率、成像package、负结果与费用一并交付。

**核心结果是一个真正能运行的多解释概率成像接口。** 同一份散射测量可以同时保留实心与中心空腔材料族；新增独立几何传感器后，可以合法改变它们的概率而不重写Maxwell数据。完整生成器posterior精度尚未完成，完整协方差也没有在所有重要事件指标上胜出。当前建议继续受限歧义benchmark与下游handoff，优先研究主动测量，而不是扩大网络或预设唯一重建。

## 1. 逆散射为什么不应只交一张图

散射测量看到的是材料经过Maxwell方程后的响应。不同结构可能通过介电常数、几何或仪器calibration的变化产生相近响应。最小residual算法会选一个endpoint，但选中一个并不意味着其他材料世界的概率为零；连续图像的一点density也不是“这张图有多少概率”。

贝叶斯逆散射已经有多种相容形状家族、未知材料联合推断与evidence积分的先例。Carpio–Iakunin–Stadler的二维Helmholtz/holography研究说明MAP附近Gaussian近似可以遗漏观测相容形状，因此“多张散射重建”本身不是原创发现。[原文§6–8](https://arxiv.org/pdf/2003.09318)

条件density的高斯混合密度网络（MDN）与完整协方差也有既有方法基础；温度序列蒙特卡洛（SMC）、重要性修正和独立context重新加权同样不能单独当作新的Maxwell理论。本研究的可验证问题更窄：**在同一物理、prior、noise与nuisance条件下，分开几何/材料/gain补偿，算清材料事件质量，并让下游保留这些解释及其不确定性。** 文献检索429与全文缺口保留，未完成全球原创性认证。[主张账本](LITERATURE_AND_CLAIM_LEDGER.md)

## 2. 共同材料世界与观测target

生成器输入没有孔洞数量或拓扑标签。它生成连续`x=(zc,zs,ε′,ε″)`：zc、zs独立均匀于[−3,1]，ε′∈[2.2,2.8]，ε″∈[.04,.16]。两个互不重叠径向bump形成`u(s)`；u<0取**严格真空**，其他区间共享复ε。解析roots形成合法的分层球，复用已存在B1solver。

这个表示自然产生实心、中心空腔、环形空腔与组合，但支持仍是径向同心材料；它不是开放的一般三维形状先验。zc或zs∈[−1,1]时对应空腔几何处于未激活状态，各有prior质量.5，因此实心材料的先验质量.25。原系数均匀先验被正确推到半径上，含零半径原子；不能悄悄换成radius-uniform。

冻结R=1，八个角度20°到160°、两种偏振，每频16个复数测量通道，实部在前、虚部在后的实化排列。kR=1.6为主频，2.4为预先指定附加频率。每实通道噪声标准差为固定solidε=2.5+.1i响应实排列后的均方根的.12倍；没有为了“成功”随候选材料增大noise。real/imag各自方差σ²，complex总variance则为2σ²。

未知gain是一个**实数幅度**，propertruncatedN(1,.15²)，范围[.7,1.3]，两频共享。边缘化保留log determinant与截断归一化常数，不用最优 gain 的残差替代evidence。复介电常数在测量频带内不变跨有限频带是校准台假设；没有未知相位、现实色散或复杂仪器模型认证。

A/B/C/D四个信息条件分别是已知ε/gain、未知ε且gain 已知、未知ε+gain，以及C加第二频率。方法比较只在同target内；D增加通道，原始 log evidence的数值也不能直接与C当作无单位methodscore。

后验为`p(x|y)∝L(y|x)p0(x)`，材料事件质量为相应区域积分。十二个event是四中心空腔半径区间与三个环形空腔半宽区间的组合，互斥且覆盖声明support。材料事件划分用于输出与评价，不作为generator的强制专家标签。

## 3. Full Maxwell补偿witness与负结果

固定solidε=2.5+.1i、gain1作为参考，与r=.25R中心空腔比较。在单频1.6下原始白化距离 **1.529652**，仅ε补偿 **1.295221**，ε+gain补偿 **.900170**。补偿endpointε≈2.466617+.130009i、gain≈1.053221；prior内部，gain标准化位置约.355。材料体积距离仍.12846，结构差异没有因测量接近而消失。

两频共享gain时该补偿后的白化距离升至 **2.467944**。小空腔r=.08R则单频.05684→.03559，双频仍.12598→.10750；额外frequency不能消除一切弱可辨识方向。

环形h=.08R、中心加环形、先验边界案例是清晰负结果。单频环形raw13.60，即使ε+gain补偿仍4.45；双频仍8.02，多endpoint需要碰ε边界。低residual的best-found解不一定有高prior邻域质量，也不是evidence。全部参数、边界位置与full-Mie 实际动作终点保留，没有调大noise或改prior制造多解释。[物理补偿报告](AMBIGUITY_AND_COMPENSATION.md)

局部白化Jacobian投影减少部分response幅度，但r=.02–.08R拟合raw/projected斜率 **2.978/2.975**，未支持普遍r³主导项完全抵消后留下r⁵机制。明确停止该通用幂次主张；它与大空腔合法补偿witness是不同问题。[LOCAL_MECHANISM_STATUS.json](results/LOCAL_MECHANISM_STATUS.json)

## 4. 概率reference先算对，再谈生成器

共同两频响应数据组每频625物理节点；连续 Chebyshev 响应近似与full-Mie 数值目标各有identity，节点不是prior离散化。112独立先验抽样点×2freq的full核验给最大response相对误差2.986e−4、whitened差.014871。这是有限点经验性证据，没有全支持一致误差界。

七个开发集 C 条件案例的9/13/17连续积分统一evidence、means、SD、CDF与quantiles：CDF采用截断区间直接积分，quantiles直接反解，不再用prior-bin展宽。13/17最大CDF差.00011715，quantile/先验区间宽度差.00010308，mean差3.6703e−5、标准差相对差.00023042。

独立系数坐标自适应积分初次未达到内部2e−6精度，保留失败；改实际半径/半宽坐标并保留同priorJacobian后，七例5/9rules全部收敛，事件质量与17阶张量积分最大差4.80e−6。训练 teacher9阶均值误差除以先验宽度 .001433未达.001，保留旧数据再细化13阶，改善至7.816e−5。没有通过放宽标准或重跑truth遮盖数值缺口。

large-cavity和solid两case另用可计算提议密度 q、full-minus-proxy 残差控制变量，各384/192full-Mie 调用。large-cavity最大event质量改变量.00015128、最大残差蒙特卡洛标准误4.583e−5；solid相应.00010187、6.114e−5。窄bootstrap只代表已访问residual的MC误差，quadrature与模型误差另列；不能称全域或现实概率认证。[reference报告](POSTERIOR_REFERENCE.md)

128个新观测另做固定65点的直接连续积分 CDF。最初节点重采样CDF从连续精度主张中撤下；原截断积分有.00265042非单调缺陷，原数组与失败保留。最终改为对**同一连续密度**的相邻互斥区间分别积分正质量，再累加CDF，单调性缺陷为0。这不是单调回归，不假设区间内概率均匀，也不从网格推算分位数。

不同积分划分与原13阶规则的最大log evidence差.00998733、事件质量差.00550557，新旧CDF最大差.01476339，均完整披露；它们是同一exacttarget的不同数值近似，不误称evidence归一化错误。前八例13/17partitionCDF最大差.00097375、条件CDF差.00414378；直接反演分位数仍使用独立原计算。全体观测和所有阈值没有一致精度认证，生成器比较保持经验性固定网格范围。[最终CDF审计](results/FRESH_CONTINUOUS_PARTITION_AUDIT.json)

仅核对前八例不足以代表全部观测。最终已对128个新观测全部复核event、均值和SD：原13阶的观测12、39、48失败，升级17阶并对照21阶后均通过，无未解决项；三例CDF也重新按17阶positivepartition积分。主参考改为 [FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz](results/FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz)、[ACCEPTED_REFERENCE_AUDIT.json](results/ACCEPTED_REFERENCE_AUDIT.json)，原13阶和中间partition结果保留。该修正新增physics0，未修改模型、测量、prior或noise；CDF细化与直接quantiles仍限声明的子集，不据此认证所有阈值或现实Maxwell。

## 5. 同一份测量中的实际多材料解释

预注册large-cavity/Cmeasurement的full-Mie 残差修正后验给：

| 不重叠材料事件 | 后验概率质量 |
|---|---:|
| 较大中心空腔r∈(.20,rmax]、无ring | 67.7226% |
| 实心 | 16.2540% |
| 中等中心空腔r∈(.10,.20]、无ring | 7.1828% |
| 无中心空腔、小ringh∈(0,.06] | 4.1259% |
| 其余已建模事件 | 4.7148% |

前四个事件累计后验质量 **95.2852%**，没有因为展示四个结果重归一化到100%。另一个truth为solid的观测仍给solid40.7667%、较大中心空腔41.0889%、中等中心空腔11.6300%。这正是measurement不能确定唯一材料时应保留的ambiguity。

这些事件不自动等于分离modes；可能包含continuousridge及相邻材料族。每个package提供事件内材料距离的代表样本、ε实部/损耗截面、条件means/SD/covariance、合法完整正向预测、gain 条件后验、材料距离与噪声尺度下的观测距离、先验与后验 odds、全十二事件质量和事件内部带权样本。full-Mie 修正后的质量与proxy 条件带权样本是不同数值近似，明确标注，未宣称exact完整后验样本。未知模型外mass为`null/unquantified`。[large-cavityJSON](results/packages/large_cavity.json)、[solidJSON](results/packages/solid.json)、[图库](index.html)

后验预测检查按完整权重抽取1024个材料/gain样本，使用固定噪声下残差平方和；solid与large-cavity的tail统计约.7425/.7481，gain负控制约.09256。原uniform-list抽样会过度使用补充rare-eventparticles，已修正并保留原诊断。tail统计用于检查观测相容性，不是“图像正确的概率”，也不认证真实仪器模型；计算复用proxy及已支付的fullfield核验。

## 6. 哪些概率主要由prior或nuisance决定

不同prior不混入主方法排名。另设独立敏感性 target，在保持相同响应与噪声时改非零几何部分为radius-uniform、保留零几何原子质量，并独立实际几何坐标自适应积分。solid观测的大central-eventproxy质量从系数先验 **.410802**变为半径先验 **.147149**；large-cavity相应 **.677377→.262458**。这说明空腔半径事件质量相当先验敏感，不能只展示某组prior的候选当现实确定性。

正规 gain 先验的原始标准差从.15分别改.10/.20，仍保持[.7,1.3]归一化常数与相同 σ，solid大central-event质量约.398059/.416090；large-cavity约.689427/.672597。这里gainprior变化比geometryprior变化小，并不意味着所有cases或现实gain模型都不重要。[独立sensitivity数据](results/PRIOR_NUISANCE_SENSITIVITY_REFINED.json)

单频C加第二频率的D消融，large-cavity的solidproxy质量.162591→.000431、大central质量.677377→.966140；solidcase的solid质量.407769→.833822。这是附加measurement有价值的数值证据，不是新算法胜出，也没有对D做全支持完整后验certification。

gain负控制的真实gain=1.2，而A/B仍固定gain=1，专门模拟错误校准。它不是“已知真gain”的公平信息对照；C/D的正规有界先验包含1.2，说明为何不能把gain误设吸收到材料而装作唯一结构已经恢复。

## 7. C1共同target上的正面结果与不足

C1-P1的二维vector-DDA 强度测量结论仅作方法起点，不直接移用其checkpoint或成绩。P2为共同B1四维复场测量目标重新适配归一化有界 MDN：八个混合成分、完整/对角协方差、两个96单元隐层、三个固定随机种子、共同3000次参数更新。teacher与bank是额外监督，setup费用披露；没有拓扑标签、state/memory或大型GPU训练。

64个新材料对象×两份独立噪声在全部15个模型冻结后生成。full/diagonaldistillation平均事件质量 MAE **.013215/.013028**；按对象配对 bootstrap差.0001871，95%区间[−.00057984,.00088033]跨零。因此不宣称完整协方差事件质量全面更好。

旧13阶node-expectation的latentKL **.215723/.274824**、paired差−.059101及区间[−.067108,−.050940]保留为历史经验近似，没有把三例reference修正改称已重算完整KL。decoder的inactive等价纤维使latentKL不等于材料场KL；不能以该改善数值声明材料fieldKL相同改善，两者也远高于.02远期目标。[KL链式分解](THEORY_HANDOFF.md)

两蒸馏模型的十二event重要质量覆盖率1、错误排除/过度自信0，平均Top-K 的参考后验覆盖率约.97764/.97481。单高斯事件质量 MAE.027550；truth-pair 联合full/diagonal约.048227/.038368，表现退化，但监督信息少于distillation，不能把差异全归网络架构。

最终正区间CDF比较给full/diagonal平均边际网格误差 **.098906/.101614**，重要事件内条件CDF误差 **.141761/.161768**，明显未达.01研究目标。事件概率表现良好并不代表事件内部连续材料不确定性算对。Full的条件CDF有所改善，配对差−.0200067的95%区间[−.0255638,−.0148288]；边际差区间跨零。参考仍只有已披露的积分差异与子集细化证据，不能升级为现实概率校准。[方法报告](METHOD_COMPARISON.md)

三个蒸馏提议在同large-cavityC、seed26101321各64full-Mie 调用做合法`.9MDN+.1prior`SNIS：修正后的事件质量 MAEfull .012752、diag .025983、Gaussian .018627，ESS36.6/29.5/23.9。仍有零访问小事件和宽蒙特卡洛区间；单case不证明普遍优越性。通用64粒子 SMC两个seeds均到beta1、295/301calls，但重要事件质量变化大，高ESS不保证无漏模式。神经teacher额外setup与传统在线预算不同，不能直接宣布NN更准确更快。

## 8. 独立下游信息的具体demonstration

下游使用B1参考产生的带权样本和材料数据包；不是把C1原始生成分布直接交接后取得这些regret。它验证可信后验信息的传递价值，不意味着C1已达到同样的连续材料或下游精度。C1-only情景更新对照尚未运行。

在既有development顺序第一个`solid`case，独立几何传感器读数c=−.00628881、SD=.01。Maxwell 后验的solid **40.7667%**、largecentral **41.0889%**，合法contextupdate后变为solid **72.9739%**、largecentral **9.0147%**，mediumcentral仍 **14.0581%**。sensor带噪声，因此不该把其余候选删成零。

旧MAPbestendpoint在较大中心空腔事件9，只交一张图使下游无法恢复已经丢掉的solid支持。Top-K或完整权重都能改选solid。它证明的是handoff必须保留候选与uncertainty，不是证明重新算精确BayesianMAP会错误。

全128个新观测以十二eventzero-oneloss评价，平均后验决策 regret：single保存样本中密度最高的候选 **.247820**，Top-K代表 **.063901**，完整带权后验 **.00007225**。paired完整带权后验减单一候选差 **−.247748**，按对象计算的95%区间 **[−.302452,−.193408]**。实际0–1loss分别.609375/.437500/.367188；参考后验最优策略有限truthloss.375，不把sample方法的随机实际损失顺序解释成超过贝叶斯最优。

Top-K+条件带权样本+remainder与全部带权样本相同，是保留同支持/weights的表示identity。simpleTop-K代表忽略event内情景似然平均，因而保留.0639regret。没有新信息的常数似然不改变原weights。七例与fresh完整传感器噪声、reference、decisions和bootstrap均保存，未挑案例noise让truth概率每次上升。[下游报告](DOWNSTREAM_DECISION_TEST.md)

严格的交接关系为`p_j(c)∝p_j E[p(c|χ)|A_j,y]`，分母包含未展示事件。代表图不能一般替代这个条件期望；[THEORY_HANDOFF.md](THEORY_HANDOFF.md)给出条件期望、TV/regret误差界和latent/materialKL链式分解的证明。这些是已知概率关系的推导，不作为新原创定理；实验未取得对所有情景的uniform误差上界。

## 9. 费用与可复现交付

新增材料–频率attempt **3459/4000**，剩余541；mechanism/reference/fresh/reserve分别86/1847/758/768。失败1为pre-kernel最大半径setup，已append-only补收；successfulkernels3458、RHS96488、kernelwall12.544s。完整controlledwall包括teacher、training、CDF、失败与重启，远低于7200s保护，精确最终receipt见resourceledger。旧B1384/384与C1-P1原样保留。

全controlledjob计费与uniqueBEGIN不同：并行零物理job的旧start/endglobalcountdelta错归了其他owner的reservations，原row保留，通过独立attributioncorrection修复；总3459从未减少或reset。GPUtraining0、旧lockedtruth0、DeepSeek0。必要数学/物理检查已执行，无SHA或额外smokecampaign。[资源](RESOURCE_SUMMARY.md)、[失败](FAILURES.md)

交付包含可浏览Top-K图库、概率JSON/NPZ、完整带权后验、条件更新接口、code/config、三seedmodelweights、continuousreference、全部rawfullfields、job/physicsledgers和科学裁决。[MATERIAL_FIELDS.npz](results/MATERIAL_FIELDS.npz)还给出径向每个位置的真空概率、材料实部/损耗均值和SD；范围限于声明的球对称支持与有限带权样本近似。现有B1/C1仓库延续使用，历史结果保留，不为每组experiment新建repo。

## 10. 十个科研问题的裁决

1. **是否找到不同材料、近似等价散射？** 找到一个受限r=.25R的full-Mieε+gain补偿witness，且独立posterior确实保留solid/centralcavity显著质量；没有证明严格等价或一般三维nonuniqueness解决。
2. **等价性来自何处？** ε和realgain分别吸收部分geometrysecant；附加frequency限制补偿，smallcavity仍弱可辨识。ring/boundarycases是补偿失败负结果。
3. **算法能否保留候选？** 可信事件参考与未经物理修正的蒸馏密度能在预注册材料事件分辨率保留重要候选；有限SMC/SNIS与连续CDF仍失败，不能保证所有材料后验模态。
4. **C1完整协方差仍有价值？** 对共同target改善latent 密度和部分条件 CDF；事件质量 MAE与diagonal相当。保留它作为proposal，不升为全面赢家。
5. **比sampling更准更快？** 目前不能作一般结论：teacher/setup额外监督、预算不同，full-Mie 修正成本明确但MC不小；可信低维积分仍最有用。
6. **何时Top-K足够？** 当声明event充分细且其内部context相关变化可忽略时可作摘要；本几何context明确需要条件样本与remainder，完整带权样本更稳妥。
7. **独立context改善判断吗？** 在声明模拟传感器与固定案例集合下有明显交接决策 regret收益；不保证所有现实context独立或有效。
8. **概率保证是什么？** 明确归一化且明确的概率目标、Jacobian、合法marginal/weights/MH与有限数值误差；没有完整后验的一致误差界/现实概率校准、模型外mass量化或all-mode保证。correctMH仅指数值目标。
9. **是否够IEEE TAP 强原创贡献？** 尚未建立。Bayesian/MDN/SMC等有先例；当前窄贡献候选是受限完整 Maxwell补偿witness、event参考与下游交接的共同实证。通用r⁵机制被当前证据否定，不冒称新定理。
10. **下一步选什么？** 继续受限材料歧义基准，优先面向材料事件的主动测量与reference质量；保留full/diagonal蒸馏提议。停止盲目扩大NN、通用主导项抵消与神经网络全面优越性，非同心大campaign必须新scope/授权，当前NOT_RUN。

独立报告的结论是：Maxwell 概率成像可以交出“哪些不同材料世界仍可能、为何难区分、各有多少事件质量”的可运行package，并允许真实独立信息继续更新；它不能因生成更多样本或发布页面而宣称完成完整材料后验或逆问题唯一性。最终机器可读裁决见 [GATE_DECISION.json](GATE_DECISION.json)。
