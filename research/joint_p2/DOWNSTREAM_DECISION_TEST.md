# 独立情景信息如何重新判断材料世界

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

此实验交接的是posteriorpackage，而不是让下游凭图片猜真值。原Maxwell测量不修改；下游得到一个独立、带噪声geometrysensor：`c=normalized vacuum volume+η`，η∼N(0,.01²)。sensor可以为负，因为读数加Gaussiannoise；没有把它截成无噪声truth。给定材料，sensor与scatteringnoise条件独立，更新合法likelihood为`p(c|material)`。

精确更新需要每个event内部的条件似然期望，而不是只用代表点的likelihood；分母还包括未显示event。[THEORY_HANDOFF.md](THEORY_HANDOFF.md)给出该关系、TV与decision-regret界及其前提。该数学界未获全contextuniform误差常数，因此实验使用实际reference和有限bootstrap，不把理论形式冒称已认证保证。

truth只用于生成模拟sensor与评价，不输入posterior inference。development七例的实际sensor值及原stream冻结于 [CONTEXT_SENSOR_MANIFEST.json](results/CONTEXT_SENSOR_MANIFEST.json)，seed26101328；freshpanel使用另一个固定substream261013216。第一次contextreference精度不足时，保留原sensor值与失败，不挑noise重新生成。两套结果均复用已支付responsebank及fullcorrectionpoints，context阶段新增physics0。

**下游实验使用B1 reference-derived的wholeweightedparticles/package，而非未经修正的C1生成样本。** 因此下述低regret验证了概率信息交接机制，不能解释成C1已实现相同下游精度；C1-onlycontextupdate对照尚NOT_RUN，已披露的C1条件CDF误差仍然存在。

## 第一预注册案例：实心与中心空腔都不能先删

按development既有顺序，第一个案例是`solid`。单频complexscattering、未知ε和properboundedgain的full-Mieresidualreference，在增加context前同时给：

| 材料事件 | Maxwell posterior | 加独立sensor后 |
|---|---:|---:|
| 实心event0 | .407667 | .729739 |
| 中心半径(.10,.20]、无环形event6 | .116300 | .140581 |
| 中心半径(.20,rmax]、无环形event9 | .410889 | .090147 |

该sensor读数c=−.006288811645、SD=.01，真实归一化voidvolume0只作模拟标签。sensor更支持低voidvolume，但仍不给实心100%：约9.01%大中心空腔与14.06%中等空腔继续保留。wholeweightedpackage更新与独立reference最大event差约8.46e−5；这不是精确fullposterior采样，错误与target分别声明。

只交一张旧MAP时，已找到bestendpoint落在event9；单点handoff没有另一种材料支持，无法因sensor把event0重新引入。Top-K、wholeweights和richpackage都选event0，zero-oneposteriorregret0；旧MAP的regret .639592、该模拟truth的0–1loss1。该结果说明**丢掉候选支持**会妨碍下游更新，不是证明重新运行完整BayesianMAP也会失败。

第二个`large_cavity`读数c=.004698922145、真实voidvolume.015625，sensornoise并不有利于真值。event9从.677226下降到.550074，event0从.162540上升到.290861，event6从.071828升到.140660；全部策略仍选event9。此例保留真实noisycontext能加强错误方向的信息，而不挑选每次都让truth概率上升的结果。

## 全freshpanel的决策比较

64freshparents×两noise共128observations。损失为冻结十二事件的zero-one分类；如果contextposterior为p，选择eventa的Bayesregret是`max(p)−p[a]`。它评价当前事件decision，不等于PSNR或一般应用utility。parent配对bootstrap2000、seed26101329，两个noise不伪增独立parents。

| 下游传递方式 | 平均posteriorregret | 95%parent区间 | 实际0–1loss |
|---|---:|---|---:|
| 单一highest-densityparticle | .247820 | [.197895,.303787] | .609375 |
| Top-K代表材料+eventmass | .063901 | [.036934,.094624] | .437500 |
| 完整weightedposterior | .00007225 | [0,.00020143] | .367188 |
| Top-K+条件样本+remainder | .00007225 | [0,.00020143] | .367188 |
| 独立contextreference的最优event | 0 | [0,0] | .375000 |
| 不使用新增信息的最大posterior事件 | .088193 | [.050239,.131336] | .421875 |

pairedwholeweighted−singleparticle regret差为 **−.247748**，95%区间 **[−.302452,−.193408]**。fresh“MAP”具体是保存的posterior节点中latentdensity最高的particle，不是全域连续MAP；development则为八个起点的best-foundendpoint。必须保留这个限制，不能把整个比较升级为精确MAP的数学下界。

richpackage保留了全conditionalparticles与remainder，因此与wholeparticleupdate一致是**表示identity**，不是另一算法独立胜出。完整weighted方法在有限truthpanel的实际loss稍低于referenceoptimal，不表示超过Bayesoptimal：它们occasionally选不同事件，有限truth抽样可使realizedloss顺序改变。regret仍以独立reference为准。

Top-K点代表的sensorlikelihood不等于event内部条件likelihood平均；由此产生.0639regret，而不是probabilitytarget改变。如果只展示Top-K且不传未显示质量，就不能对新context合法重归一化成完整posterior。用户需要可靠下游推断时，应交wholeweightedparticles，或Top-K+conditional samples+remainder。

## Reference、费用和界限

初始2048celladaptivecontextreference失败，之后采用context-awarephysicalpartition、8192cellcap与原2e−6内部准则；prior、likelihood与sensor未变。fullcorrecteddemos重新加权已支付full-minus-proxyresidual，其bootstrap仅包含该residualMC；freshreference属于连续proxytarget。全部原始decisions、sensorvalues、updatedeventmass、ESS和bootstrap在 [DOWNSTREAM.json](results/DOWNSTREAM.json)。

常数contextlikelihood原样保留weights；没有新信息不能消除ambiguity。若现实context与y并非条件独立，需要`p(c|material,y)`；现在的模拟independencesensor不能证明任何任意现实上下文都适合独立重新加权。模型外materialworlds仍null/unquantified，更多weightedparticles也不会自动扩展support。
