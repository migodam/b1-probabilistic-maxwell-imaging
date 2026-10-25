# E5：有限先验条件下的主动测量对照

E5保存了24个条件、七种策略共**168/168条结果**，即12个父世界各两个初始噪声repeat。TaskEVSI在本轮32次预测样本的平均预期目标风险上较低；实际未来观测后的风险和真实误差配对区间仍跨零，不能据此宣布普遍优越、统一最优采集准则或新的任务导向设计理论。结果仅限`FINITE_PRIOR_ENUM4_V1`，最终科学判断仍由父线程完成。

## 选择发生在未来观测之前

当前信念来自t1：TX0/RX8单个复数观测。四个动作是在source0/1/2/3之一取得新的full16RX观测，噪声独立于t1。即使再次选择source0，其RX8也是**新的独立测量**，使用不同因子ID；原t1因子保留，未重复计入同一个观测。

Random采用seed2026102506的条件子流；Uniform按预定round-robin。EIG和TaskEVSI对每个条件使用32次当前后验预测抽样，候选动作之间共享同一组假想世界和标准化噪声以减小比较方差。EIG最大化预期KL信息量，TaskEVSI最大化目标位置Euclidean平方损失Bayes风险的预期下降。两者目标不同，预期风险也不等于某次实际定位误差。

动作选择先保存到`ACTION_SELECTION.json`。之后离线生成器才读取父世界索引，并用seed2026102507的独立噪声生成实际未来观测；选择函数不接收该索引或未来观测。每种策略的后验仅使用唯一t1因子和自己选中的唯一新测量因子。真实误差保存在独立offline文件，未用于选动作。

三个局部对照在**当前早期后验MAP**处求原生Maxwell Jacobian，不在真值处求。A为两个墙材料列，B为两个目标位置列，噪声whitening为每实/虚分量σ=.01。SingularValue使用公开物理尺度[4.5,.25,2,1]后的完整四列J最小非零奇异值平方；秩阈值和原物理信息矩阵均保存。没有伪逆或jitter来补齐不可辨方向。

## 三种风险量必须分开

下表均先平均父世界内两个repeat，再平均12父世界；单位为λ₀²。

| 策略 | 32样本预测的动作后期望风险 | 实际未来观测后的后验Bayes风险 | 离线实际GT平方误差 |
|---|---:|---:|---:|
| Random | .26194 | .24540 | .41206 |
| Uniform | .26011 | .22874 | .41393 |
| EIG | .21266 | .25835 | .34792 |
| SingularValue | .24975 | .27753 | .37372 |
| PrincipalAngleSum | .25849 | .27691 | .35634 |
| GeneralizedRetention | .25849 | .27691 | .35634 |
| TaskEVSI | .19876 | .21706 | .32674 |

第一列是选动作时的有限MC期望估计；第二列是独立未来数据实现后的条件后验风险；第三列是该实现的目标后验均值相对离线真位置的平方误差。一个实现里三者不应被强行视为相同数值，也不能用第一列取代第三列证明实际收益。

所有策略都有24/24保存行，没有删掉不利条件。粒子/MC抽样和两个repeat不是新父世界，实际独立分母仍是12。

## 配对区间保留了未解决的不确定性

现成汇总用父世界配对bootstrap2000次、seed2026102510；差值为TaskEVSI减基线，负值表示本指标较低。下列区间都是描述性的，不是新的策略优越性通过门。

| 比较 | 指标 | 配对均值差 | 现成bootstrap区间 |
|---|---|---:|---|
| TaskEVSI−EIG | 预测期望风险 | −.01390 | [−.02361,−.00605] |
| TaskEVSI−EIG | 实现后后验风险 | −.04130 | [−.09486,.01047] |
| TaskEVSI−EIG | 实际GT平方误差 | −.02118 | [−.08024,.02223] |
| TaskEVSI−Uniform | 预测期望风险 | −.06135 | [−.08276,−.04018] |
| TaskEVSI−Uniform | 实现后后验风险 | −.01168 | [−.09362,.07377] |
| TaskEVSI−Uniform | 实际GT平方误差 | −.08719 | [−.31286,.15429] |

TaskEVSI用来选动作的也是这32个预测样本；第一列的负差带有选择内的乐观性，不是独立验证。保存的48条条件MC不确定性记录（24条件×两个基线）中，配对MC标准误范围为0至.05085，均值.02299。该条件MCSE不是经过选择校正的置信区间，更不能消除跨父世界和实际未来噪声的不确定性。实际风险/GT误差的上述区间均跨零，保留这一事实。

## 裸角度与无先验retention是同一几何量

本轮PrincipalAngleSum使用墙响应子空间与目标响应子空间夹角的`Σsin²θ`；GeneralizedRetention使用同一无先验投影后的retention特征值之和。在这里的响应子空间定义下，它们代数等价。96个候选动作审核中，两个score最大绝对差仅4.44e−16，24条件动作完全相同，结果行也相同。这不是两个独立成功的科学检验，也不是retention或任务设计的一项新普适性发现。先验辅助、不同参数度量或其他非线性准则的推广不由本轮证明。

实际动作计数（source0/1/2/3）为：Random 6/4/7/7；Uniform 6/6/6/6；EIG 4/18/1/1；SingularValue 3/12/4/5；裸角度与retention均5/8/6/5；TaskEVSI 11/11/1/1。这些是当前条件模型中的已实现选择频率，不是普适最优孔径分布。

## 费用：完整bank与局部J均已付费

E5实际新增**480 RHS、24 FE**，每个条件一次新的full4TX原生J（4 forward+16 transpose RHS），三个局部策略共享这一实际求导；没有借用假定可得的旧LU。E5 wrapper真实CPU为**174.407221秒**、wall176.376806秒、子进程最大RSS1110573056 bytes。

局部策略独立运行的选择费用包含该条件完整J和审核，平均约5.626 CPU秒；Random/Uniform选择本身没有被虚构收费为J。EIG/TaskEVSI共享有限表的32样本更新，其standalone选择CPU均值约.03673秒，但这不包含把Maxwell模型制成全表的设置费，也不能当成端到端速度优势。

继承bank为**11800 RHS、2950 FE、9288.023040 CPU秒**（含失败001），wall9496.344884秒。E2E3实际账本只计一次，E5/各standalone策略比较仍显式归属它。即使Random/Uniform不靠全表选动作，当前统一精确有限后验的推断/评价仍依赖这个共同bank；依赖早期精确MAP的局部策略同样继承设置工作。不能把共享设置或原生J藏掉后声称某策略总体更省。各standalone费用列是替代方案归属，不可相加成真实多次收费。

## 证据路径与报告边界

- [168条在线结果](results/e5/finite_sensing/ONLINE_RESULTS.json)、[完整终态与范围](results/e5/finite_sensing/MANIFEST.json)、[保存原始策略行](results/pilot_summary_tables/E5_ALL168_POLICY_ROWS.csv)。
- [父世界指标](results/pilot_summary_tables/E5_PARENT_METRICS.csv)、[Task/EIG/Uniform配对区间](results/pilot_summary_tables/E5_TASK_EIG_UNIFORM_BOOTSTRAP.csv)、[32样本条件MC不确定性](results/pilot_summary_tables/E5_SAVED32_MC_PAIRED_UNCERTAINTY.csv)。
- [独立离线GT汇总](results/e5/finite_sensing/offline/PARENT_GROUPED_GT_SCORES.json)、[E5 wrapper receipt](results/jobs/active-sensing-001/RECEIPT.json)、[最终bank费用](results/e2e3/finite_bank/online/JOB_COST_FINAL.json)。每个`parent_XX_repeat_R/`还保留动作选择、原生J与后验数组；实际未来噪声/观测保留在独立offline文件。

这次文档工作只读取保存数据和已完成的汇总，没有新求解、重新选动作、重新采样或创造新的实验支持。条件有限先验、已知目标本征参数和校准仪器的范围保持原样；连续概率、普适采集优势与最终科学接受均未由此被动计算建立。
