# 共同 posterior reference 与数值误差

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

共同 target 是原 coefficient-uniform prior、固定八角度两偏振 vector-Mie response、固定 real-channel noise，以及 proper bounded shared gain。所有 probability 都针对该声明支持；unknown/outside-model material worlds 没有被人为赋成零质量。

## 物理与 surrogate 分开

Full numerical target 使用既有 B1 matching solver，nmax=14、kR=1.6/2.4。七个 stress materials 在两频分别做 n10/n14/n18 与独立 admittance recursion：n14/n18 response 差在保存精度为零，recursion 最大相对差 4.44e−16，passive absorption 与 boundary residual 保留在 [PHYSICS_QUALIFICATION.json](results/PHYSICS_QUALIFICATION.json)。这是限定参数样本的资格核对，不是 continuum Maxwell 认证。

共同 response bank 以物理坐标 `(r,h,ε′,ε″)` 每轴五点、每频 625 点构造，共 1250 个材料–频率节点。节点用于连续 Chebyshev response approximation，**不是离散材料 prior**。已支付重复几何通过严格 cache 复用并记 provenance。两种 target 分别标识为 `vector-Mie-matching-n14-radial-v1` 与 `tensor-Chebyshev-five-per-physical-axis-radial-v1`。

112 个独立 prior probes、两频共224次 full solves 给最大 response 相对误差 2.986e−4、最大 whitened field discrepancy .014871。该有限点结果在 [INDEPENDENT_FORWARD_AUDIT.json](results/INDEPENDENT_FORWARD_AUDIT.json) 与原始 [arrays](results/INDEPENDENT_FORWARD_ARRAYS.npz)；它不能推出全支持 uniform bound，尤其不能自动认证重要事件 odds 或 CDF。

## 连续积分和同测度摘要

先按十二个 disjoint material events 划分原 prior，分别 tensor integration，再合并 evidence、moments 和条件 moments。对 C 的七个 development observations 使用9/13/17阶。CDF 通过截断区间直接积分，quantiles 通过求解该 CDF 反演，几何 zero atoms 保留。任何新摘要不做 prior-bin 展宽。

13/17 的七例最大差为：event mass 1.3501e−5，CDF .00011715，mean/prior width 3.6703e−5，quantile/prior width .00010308，SD relative .00023042。原注册比较标准分别是 .01、.01、.001、.001、.01，未放宽标准。[PILOT_REFERENCE.json](results/PILOT_REFERENCE.json) 同时保留三阶结果、条件事件 CDF/SD/quantiles、variance–quantile coherence 检查。

第一次 independent adaptive integration 在 coefficient 坐标使用512 cells仍未满足内部2e−6精度；记录未删除。改为实际几何坐标，保留**完全相同 prior 的 Jacobian 与零几何原子**，5/9 paired rules 在152–408cells收敛。七例与tensor17最大 event difference 4.80e−6、mean/priorwidth 2.14e−6、SD relative8.49e−6。[ADAPTIVE_PHYSICAL_REFERENCE.json](results/ADAPTIVE_PHYSICAL_REFERENCE.json) 保存原失败、leaf partition、cross-comparison与非认证状态。此改动改变积分坐标，不改变概率 target。

## Teacher 与 fresh parents

Train128、validation32、fresh64 parents，各两份独立 noise。fresh truths在全部15个checkpoint冻结后生成；历史 development panels仍是暴露资料。teacher 是 continuous proxy posterior 的分层 resampling，inactive coefficient 从正确条件 prior抽样，不使用 topology label。

最初 train order9 的前八例 mean error/priorwidth最大 .001433 超出 .001，因此保留原数据后用order13细化。与17阶比较，前八例最大event error2.196e−5、mean fraction7.816e−5；validation相应最大event4.584e−6、mean fraction1.768e−6。teacher增加了监督信息，其bank、truth和积分成本均披露，不与单一truth pair训练混成“相同监督”。

Fresh128 observations最初采用order13、仅前八例17阶复核；该检查不能外推为全部观测通过。最终对全部128例执行event/mean/SD13/17比较，观测12、39、48在原13阶mean/SD标准下失败，升级17阶并以21阶核对后均通过，无未解决项。模型、truth、prior与likelihood不变，新增physics0；原13阶结果与失败保留于 [ALL_FRESH_REFERENCE_REFINEMENT.json](results/ALL_FRESH_REFERENCE_REFINEMENT.json)。node-resampled empirical CDF仍只能作离散诊断，不能替代continuous参考。

专门的direct continuous fresh audit已完成128观测、四个实际材料坐标、固定65点grid与重要event条件CDF。前八例13/17directCDF最大差.00097279、quantile/priorwidth最大差.00025562、coherence均通过。然而原截断CDF有最大 **.00265042** 非单调缺陷，高于roundoff；[FRESH_CONTINUOUS_AUDIT.json](results/FRESH_CONTINUOUS_AUDIT.json)及原数组保持为失败证据。

最终CDF改用**同一个连续密度**的相邻、互不重叠区间积分，累加正区间质量。这不是单调回归，不假设区间内均匀density，也不从grid插值binquantiles。128观测的marginal/conditionalCDF单调性缺陷都为0，normalization端点差4.33e−15。各坐标partition的全支持evidence与原13阶估计最大log差 **.00998733**，event质量差 **.00550557**；新旧CDF最大差 **.01476339**，差异完整披露，不把它误称evidence归一化错误。前八例13/17partitionCDF最大差.00097375，conditional差.00414378。direct反演quantiles仍沿用前八例原独立计算，不从partitionCDF推算。

主比较使用 [FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz](results/FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz) 与 [ACCEPTED_REFERENCE_AUDIT.json](results/ACCEPTED_REFERENCE_AUDIT.json)。上述三例的positivepartitionCDF也升级17阶；全128marginal/conditionalgrid单调缺陷0。event/mean/SD完成全部观测复核，CDF13/17比较仍为固定前八例加三例细化；directinversequantiles仍限前八例。原partition的log evidence/event差和旧CDF差是历史数值诊断，不冒称新reference的uniformbound。不同rules始终近似同一exacttarget，原node/truncated/partition结果保留。最终C1主continuous结果在`p2/results/continuous_accepted`，新增physics0。

## Posterior predictive consistency

七个package按actualwholeposteriorweights抽取1024个material/gainsamples，统计fixed-noise白化残差平方和的posteriorpredictivetail。solid约.7425、large-cavity .7481、annular .3615、gain-negative-control .09256；结果及原uniform-list错误诊断保留在package的`posterior_predictive_consistency`。这里forward采用proxy并经独立fullfield有限点核验；不是fullposterior/真实仪器modeladequacy证书，也不将tailprobability称为“图像正确的概率”。新增physics0。

## Full-Mie residual reference

large-cavity 与 solid 两个预注册 observation 使用同一proper likelihood；每个event分别32和16次独立full solves，总384+192调用。事件proposal为`.95 adaptive-cell + .05 conditional prior`，原latent `dx`密度明确可算。

估计量是已收敛proxy event integral，加上full-minus-proxy importance residual，而非profile minimum。large-cavity最大event质量变化 .00015128、最大residualMC SE4.583e−5；solid相应 .00010187和6.114e−5。2000次bootstrap区间只包含residualMC，不包含cubature、全域solver、missed mode或现实模型误差。[FULL_MIE_CORRECTED_REFERENCE.json](results/FULL_MIE_CORRECTED_REFERENCE.json) 与两份原始NPZ保留全部full/proxyfields、proposal density、residual及weights。

Full corrected conditional moments与量化事件质量来自该control-variate估计。package内高数量conditional samples仍来自proxy条件分布再乘corrected事件质量；另一份付费full conditional importance particles数量更少。它们是同target的不同数值近似，不能称为exact full-Mie draws。proxy continuous quantiles与corrected moments必须分别标target，不组合成已认证的full-CDF摘要。

## 能保证与不能保证

积分、proposal与likelihood都明确定义；有限物理误差、积分细化与MC误差分层记录。经验可靠性足以支持受限事件质量 demonstration，但没有 uniform posterior bound、真实仪器calibration、无限维support coverage、一般三维材料正确性或必然无漏模式的保证。条件误差区间很窄不等于模型外不确定性很小。
