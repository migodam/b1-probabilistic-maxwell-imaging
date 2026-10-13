# 共同 target 的推断方法比较

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

P2认真比较的是同一prior、decoder、noise、gain likelihood和measurement下的概率近似。已知ε/gain与双频条件属于信息消融，不能与未知ε/gain条件直接作算法优越性排名。完整fresh方法数据、seed和paired intervals在 [C1实际方法摘要](https://github.com/migodam/c1-physics-state-generative-imaging/tree/main/p2/results/summary)。

## 同teacher的raw条件密度

五种CPU density：八组件full-covariance MDN、八组件diagonal MDN、单Gaussian，以及full/diagonal的truth-pair joint训练对照。后两个不使用posterior teacher，监督信息不同。网络两个96-unit hidden layers、三个固定seeds、共同3000updates、batch128、Adam1e−3；所有checkpoint在fresh truths生成前冻结。训练bodywall38.814s、CPU37.969s，没有Maxwell或GPU调用。

下表是64fresh parents×两noise、三个seed的均值；两noise先parent合并再汇总。主eventMAE来自最终acceptedreference（全128event/mean/SD核对，三例升级17并对照21）；[REFINED_FRESH_EVENT_COMPARISON.json](results/REFINED_FRESH_EVENT_COMPARISON.json)保留每个seed、事件和bootstrap。其余KL、Top-K、Brier等列为旧13阶评价的历史经验近似，未宣称同时全部重新积分。共同exacttarget不变，raw模型频率仍是proposal分布。

| 方法 | event MAE | latent KL(P∥Q) | Top-K reference mass | 重要质量recall | false exclusion率 | false confidence率 | truth event Brier |
|---|---:|---:|---:|---:|---:|---:|---:|
| Full distillation | .013215 | .215723 | .977643 | 1 | 0 | 0 | .576698 |
| Diagonal distillation | .013028 | .274824 | .974813 | 1 | 0 | 0 | .576670 |
| Single Gaussian distillation | .027550 | .671771 | .974902 | 1 | 0 | 0 | .599932 |
| Full joint | .048227 | 20.029026 | .972161 | .999852 | 0 | 0 | .669180 |
| Diagonal joint | .038368 | 5.730530 | .947100 | .991701 | .033854 | .020833 | .650995 |

Full covariance并未在重要事件质量上明确超过diagonal。最终acceptedreference的paired parent bootstrap（2000次、seed26101329）给full−diagonal eventMAE差 **.000187101**，95%区间 **[−.000579842,.000880333]**。旧13阶node-expectation的latentKL差−.0591013及区间[−.0671084,−.0509402]、Top-K参考覆盖差.0028301及区间[.0011308,.0045764]仍保留为历史经验近似，未重算为完整材料KL。它们不支持“full covariance在所有材料概率指标上更好”。

KL是在四维latent `dx`中的数值densityKL。inactive coefficient fibers使decoder非单射，latentKL不能冒称材料场KL。所有KL仍显著高于用户远期.02目标。90%各坐标credible coverage的均值不是90%联合region coverage；不同生成器Brier/logscore值与finite-MC错误都保留。teacher resampling造成的节点CDF量化误差另列，原node-CDF表不得用作continuousCDF主张。

实际continuous comparison采用同冻结模型、同8192次draws、原prior的正区间partition积分CDF。Full/diagonal/Gaussian的三seed平均marginalmaterialgridCDF最大误差 **.098906/.101614/.173646**，重要eventconditionalgridCDF **.141761/.161768/.265015**，都明显未达.01目标；同prior下mean/SD偏差也远超前八例参考积分误差。[C1continuous partition comparison](https://github.com/migodam/c1-physics-state-generative-imaging/tree/main/p2/results/continuous_accepted)保存seed、parentbootstrap与缺失条件事件诊断。

Full的conditionalCDF有所改善，acceptedcomparison的full−diag差−.0200067，95%parent区间[−.0255638,−.0148288]；marginal差−.0027088的区间[−.0072868,.0019521]跨零。但两者未完整恢复continuousposterior。原截断CDF的.00265042非单调失败保留；最终同密度positivepartition单调缺陷0。全部128event/mean/SD已复核，三例细化后无未解决项；CDF和directquantile的误差检查仍限声明子集/固定grid，没有uniformfull-Mie认证。旧partition与原13阶的log evidence差.00998733、event差.00550557及旧新CDF差.01476339作为历史诊断保留，不以其代替新acceptedreference的全域误差界。

Posterior distillation较joint truth-pair训练有明显数值改善，但teacher包含更多监督。full/diagonal distillation用**相同teacher**，才是covariance结构的公平比较。不能把teacher成本消失后声称纯网络架构突破。

## 真full-likelihood修正：匹配64调用

冻结seed26101321、同一large-cavity/C observation，三种distilled proposals各64次full-Mie调用。proposal `.9 MDN + .1 original prior`，SNIS正确使用完整density、propergain marginal与full-Mie likelihood，所有失败、RHS与wall计费。

| 方法 | corrected eventMAE | ESS /64 | log evidence | 最大event绝对差 |
|---|---:|---:|---:|---:|
| Full distill + correction | .012752 | 36.628 | 49.086955 | .062072 |
| Diagonal distill + correction | .025983 | 29.463 | 49.094008 | .155896 |
| Gaussian distill + correction | .018627 | 23.873 | 49.484416 | .090630 |

相同observation的full residual reference logevidence为49.216529。三种都发现了reference的重要事件0、6、9，但还有未访问的小事件和宽MonteCarlo区间；64call单case单seed只能给预算明确的诊断，不能独立证明full优于diagonal、全模态覆盖、校准或普遍速度优势。低residualranking并未获得evidence身份。原始arrays、proposal branchcounts、intervals和成本见 [CORRECTED_PROPOSALS.json](results/CORRECTED_PROPOSALS.json)。

## Generic tempered SMC

64particles、两个seeds的genericSMC从原prior出发，自适应temperature并用full-likelihood delayed acceptance纠正。两run到达beta1，调用295/301次，full target身份成立；末端ESS64/44.64并不能认证事件质量。event0概率.296875/.086544、event6 .046875/.182959，与full reference .162540/.071828存在明显seed差。两run的normalizing constant误差、稀有事件漏访和全部粒子见 [GENERIC_SMC.json](results/GENERIC_SMC.json)。

SMC在线物理预算与64-call neural correction不同，且neural拥有额外paidbank/teacher/training。该实验不支持“神经网络比传统Bayesian sampling更准确更便宜”的一般结论；它说明低粒子预算下合法采样仍不等于可靠mass估计。共同小维度reference仍是概率准确性的基准。

## 冷/暖时间必须披露范围

三种density按轮换顺序各五次重放，使用64新proposaldraws及准确logq，未重复physics。full/diagonal/Gaussian cold medians约 **1.446/1.685/1.357ms**，warm约 **.216/.218/.182ms**。cold包含模型reload和conditionexport，但未flushOScache；warm复用已conditioneddensity并重新采样/求logq。[PROPOSAL_TIMING.json](results/PROPOSAL_TIMING.json) 保存每次顺序与数值。

这只是proposalacquisition，未含teacher、bank、training、posteriorreference或full-likelihoodcorrection。三种64-callcorrected proposals合计192full calls、bodywall18.865s；该wall不能和毫秒proposal时间混成完整成像成本。所有setup/amortized费用见 [RESOURCE_SUMMARY.md](RESOURCE_SUMMARY.md)。没有重复fullSMC冷/暖sweep，记NOT_RUN；不能用一次warmcache结果冒称完整成本优势。

当前方法选择是保留full与diagonaldistillation为便携proposal，两者都需要独立reference及在重要场景下fullcorrection。单Gaussian是认真简单baseline；停止把physicsstate/memory或更大NN当作默认下一步，也停止基于本受限实验宣称NNposterior全面完成。
