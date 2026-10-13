# B1 × C1 Joint P2：假设、结果与路线记录

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

日期2026-10-13。独立阶段编号`B1-C1-JOINT-P2-20261013`。B1 baseline1cd474c、C1 baseline9f2c518仅为readonly起点；旧384/384与P1账本封存。独立worktrees、唯一B1账本、4000新材料–频率attempt与7200受控wall上限由人类批准，没有追加另一campaign。

## 共同概率定义

**假设**：用同一continuousmaterialprior和propernuisancelikelihood，能够公平比较physicalreference、Bayesiansampling与C1density，并把多解释传递给下游。

**改变与基线**：decoder采用两个互不重叠radialbumps和复ε；拓扑是事后event，generator不接收H。原coefficientprior保留inactivefibers/geometryatoms，对照错误radiusuniformmeasure、shared/independentgain与highESSmissedregion。

**结果与成本**：七组零物理E0及六个normalizeddensity/Jacobianchecks通过；必要反例被检测。高ESS漏约半质量仍成立。无Maxwell新调用的概率资格并不等于fullposteriorcalibration。

**决定**：冻结properprior、noise、gain与eventpartition，再做fullphysicswitness。停止把profileminimum或networkconfidence当posterior。

## 物理补偿与附加频率

**假设**：结构secant有一部分被ε/gain吸收，另一部分可由额外frequency恢复。

**单一主要改变**：固定参考solid及noise，只释放ε，再释放同一个boundedgain；加2.4频率时仍共享材料/gain。匹配基线是同geometry的rawsecant，不改prior或noise。

**结果**：r=.25中心cavity的singlefreqdistance1.52965→1.29522→.90017，合法gain1.0532；twofreqε+gain仍2.46794。ring/combination/boundarystress不能充分补偿；小cavity仍弱可辨识。

**成本与故障**：全部qualification、selectedfullendpoints、重复与pre-kernelfailure预约计费；错误最大半径floatroundoff保留并追加补收。mechanism实际用量由唯一ledger给出，不把paidcachehit重复当newphysics。

**决定**：有明确受限ambiguitywitness，继续eventposterior；activefrequency有价值，不把所有结构称等价。

## 连续reference修正

**假设**：低维posteriorprobability应先由独立连续reference算对，surrogatesmallfielderror仍须在重要事件验证。

**改变**：同target下从coefficientadaptive换physicalcoordinateadaptive，保留priorJacobian；teacher9阶mean不满足标准则细化13阶；full-minus-proxyresidual以knownq独立修正。

**结果**：七个developmentCcontinuousCDF/means/SD/quantilescoherent；physicaladaptive七例收敛，eventcrossdiff<4.80e−6。large-cavity的fullcorrectionevent0=.162540、event6=.071828、event9=.677226，Top4=.952852。solidtruth测量仍给solid.407667与largecentral.410889。

**成本**：共同两频bank1250node，独立112材料×2freq224calls；两个residualreference384+192calls。原adaptive失败、teacher9和partialCDFwall保留；修改积分不新增truth或physics。完整stage用量见ledger。

**决定**：继续受限probabilitybenchmark与handoff，不称uniformfull-CDF或continuum认证。rawnodeCDF退到diagnostic，不用它给连续生成精度排名。

## 共同比较与fullcorrection

**假设**：C1-P1的fullcovariance+distillation能迁移到共同四维B1target，并改善density与材料事件质量。

**改变与基线**：同teacher、同updates、三seeds的full/diagonalMDN与singleGaussian；jointtruthpair保留为不同监督对照。另用generictemperedSMC与matched64call合法importancecorrection，所有fullcost披露。

**結果**：freshfull/diageventMAE=.013215/.013028，paired差区间跨零；latentKL=.215723/.274824，full改善约.05910。重要事件recall1，falseexclusion/confidence0只针对冻结12events。jointtraining退化；64particleSMCseed变化；64callcorrection的单case结果未认证全mode。

**成本**：15models×3000updates，CPUtrainingbodywall38.814s；fresh64parents×2noise模型冻结后生成，truth两freq128calls；SMC295+301calls；matchedcorrection64×3=192calls。teacher监督与banksetup不是免费。

**决定**：保留两种distilledproposal；停止fullcovariance事件质量全面优胜、NN已替代可信reference或sampling的主张。不扩展GRU/Transformer。

## 下游独立信息

**假设**：完整weightedposterior与event内uncertainty能利用真实独立context，代表点Top-K或单MAP可能丢掉必要材料族。

**设计**：冻结模拟geometrysensor`c=normalized void volume+N(0,.01²)`；预注册stream，不挑noise。wholeweights按independentlikelihood重新加权；zero-oneeventdecisionloss与referenceBayesregret，freshparent2000bootstrap。richpackage保留wholeconditionalparticles与remainder，等于wholeupdate是表示identity，不是新算法优势。

**实际结果**：contextposterior更窄，初始2048celladaptive不能收敛；保留故障和sensorvalues，再用context-aware8192cellpartition同target同容差完成七例及128freshobs。首个solid的p0 .407667→.729739、p9 .410889→.090147，sensornoise没有重抽。freshwholeweighted−singleparticle regret−.247748，parent95%区间[−.302452,−.193408]；Top-K代表仍有.06390regret，而wholeweighted约.00007225。

**成本与决定**：context新增physics0；fullcorrecteddemos复用paidresidualarrays，其他freshcontextreference为continuousproxy。丰富package等于wholeparticle是保留相同条件样本/remainder的identity，不能当独立算法优越性。继续保留多材料世界与下游hand-off，singleMAP仅为highest-densityparticle/多起点best-found，不称globalMAP。

## 连续CDF与局部机制的最终反证

**结果**：原nodeCDF撤下continuous主张；原截断CDF的非单调缺陷.00265042公开保留。改为同一density的相邻正区间积分，最终marginal/conditionalCDF单调缺陷0，不做binuniform展宽或binquantile。full/diagconditionalgridCDF .141761/.161768仍明显未达.01。partition与原13阶log evidence差.00998733、eventmass差.00550557、新旧CDF差.01476339，前八例13/17验证保持empirical，未宣称全panel认证。局部raw/projectedsecant斜率2.978/2.975，不支持普遍r⁵leadingcancellation。另设radiusuniformprior、独立5/9smoothcoordinate积分，largecentral-event概率由.677377变.262458，说明真prior敏感。

**决定**：不追加noise、不靠更多NN掩盖continuousuncertainty误差；停止通用leadingcancellation与NN全面posterior优胜主张。保留boundedbenchmark，优先event-awareactivesensing。剩余授权不要求耗尽，非同心大campaign仍NOT_RUN。

最终的数值表、具体sensor示例和研究裁决以 [DOWNSTREAM_DECISION_TEST.md](DOWNSTREAM_DECISION_TEST.md)、[GATE_DECISION.md](GATE_DECISION.md) 与最终ledger为准；此日志保留路线变化及失败原因，不把developmentreuse改称independentvalidation。

## 全fresh参考的最终细化

全128event/mean/SD13/17检查发现原13阶观测12/39/48未通过mean/SD准则；升级17并交叉21后无未解决项，三例positiveCDF也重算17阶。最终主参考为acceptedNPZ/audit，模型、truth、prior与noise均未变。新三seedfull/diageventMAE .01321508/.01302798，paired差.00018710、95%区间[−.00057984,.00088033]，结论仍无明确事件质量优势。AcceptedconditionalCDF .141761/.161768仍未达.01；旧13node-expectationlatentKL保留历史近似，不改称完整materialKL。

Unmatched-parenthesis、stringparentIDs引发NaN摘要、C1metadataNameError都是零物理执行/报告故障，原记录与wall保留；已完成CDF数组由resume复用，不重新积分掩盖错误。阶段physicalattempt仍3459；公开交付2026-10-14，最终wall由完整jobledger刷新。
