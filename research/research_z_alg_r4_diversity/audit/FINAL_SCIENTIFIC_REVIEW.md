# R4-DIVERSITY 最终独立科学审查

状态：`REVIEW_COMPLETE / GLOBALDUAL_SUPERIORITY_NOT_SUPPORTED`。日期：2026-10-30。主代理保留最终新颖性与科学裁决；此前 `INDEPENDENT_REVIEW.md` 未修改。

## 独立结论

本次完整冻结实验**不支持 GlobalDual 在低或中材料增量下具有实用材料精度优势，也不支持综合成本优势**。CPU72个配对观测中，GlobalDual 对 native TRF、reciprocal TRF、matched data-LBFGS 的 practical NRMSE 胜数均为0；CUDA18个同数据子集上，对 matched data-LBFGS 全部为 practical ties，对两种TRF也无 practical 胜。少数联合门增加和更低个别数据残差不能改写此结论。

这是本次固定网格、固定生成器、预算和算法配置下的否定结果，不是 GlobalDual 在所有问题上无用的定理。标准 resolvent、源空间加权和root-local统计性质仍可用于解释，但不构成新 Maxwell 理论、根/后验覆盖或概率解集结论。

## 审查范围与封存状态

本次只读取实际 `RESULTS_SUMMARY.json`、`audit/CLOSURE.json`、`THEORETICAL_EXPECTATION.md`，以及完成的 `analysis/ANALYSIS.json`、`analysis/COST_LEDGER.csv` 与 `analysis/COARSE_PATH_ATTRIBUTION.json`。从已有score字段定向摘取分层指标、outlier和门槛反例；没有重跑物理、比较实验、公式验证、hash、文献搜索或访问私有随机性。没有读取private truth文件。

CLOSURE记录360 planned predictions全部complete，failed=0、NOT_RUN=0，关闭时GT评分尚未开始，inverse source/budgets未变、future reruns=NONE。评分记录360个有限材料endpoint，评分新增物理LU=0。**complete是预测/传输关闭状态，不表示优化器全部收敛**；大量有效结果为BUDGET_STOP，具体停止与已付费工作必须与分数同行保留。仅有有限材料数组也不证明KKT、物理根或GT恢复。

## 连续精度与联合门须分开

以下CPU clean strata每项12个形状，NRMSE=RMSE/真材料增量：

| 材料增量 | native TRF mean NRMSE | reciprocal TRF | matched data-LBFGS | GlobalDual | joint PASS（native/reciprocal/data/GD） |
|---|---:|---:|---:|---:|---|
| .4 | .164351 | .161534 | .174929 | .172135 | 3 / 3 / 3 / 4 |
| 2.5 | .126124 | .123805 | .669848 | .676464 | 9 / 10 / 6 / 7 |
| 12 | .413788 | .409363 | .419220 | .435658 | 0 / 0 / 0 / 0 |

低clean下GD均值比matched data略好(.172135 vs .174929)，但比两种TRF更差，且所有配对差异按预声明practical规则没有GD胜。中clean下GD median=.139530与native=.140871接近，**均值却为.676464**，不能用median略低掩盖灾难失败。noisy的中增量也同样：GD mean=.679420，native=.125732，reciprocal=.125196。

低clean多出的一个joint PASS来自 `curved_bar`：GD NRMSE=.147557，native=.155997，reciprocal=.152105，matched data=.151402；各自Otsu IoU均>.7且数据残差都过门。GD确实跨过冻结NRMSE=.15门，这个PASS应该报告；但相对native的连续改善只有约.00844、小于practical tie绝对下限.02，不能称显著或实用精度胜。与data比较亦如此。中clean的GD7/12还低于native9/12、reciprocal10/12，不能挑data6→GD7来概括低中增量优越。

全部72个CPU变体的joint计数为native23、reciprocal25、data18、GD21；CUDA18变体为8、8、6、6。高增量所有method均失败：CPU每method 0/24、CUDA每method 0/6。这些都是固定门的描述性计数，不是独立Bernoulli试验成功概率估计。

## 预声明 practical 配对结果

| 后端 | 比较（GD相对对手） | GD胜 / tie / GD负 |
|---|---|---|
| CPU72 | native TRF | 0 / 59 / 13 |
| CPU72 | reciprocal TRF | 0 / 54 / 18 |
| CPU72 | matched data-LBFGS | 0 / 69 / 3 |
| CUDA18 | native TRF | 0 / 16 / 2 |
| CUDA18 | reciprocal TRF | 0 / 16 / 2 |
| CUDA18 | matched data-LBFGS | 0 / 18 / 0 |

这些数来自既有scorer，不是新增比较运行。practical tie为`max(.02,.10*min(NRMSEs))`，不等于统计等价或统计显著性判定。CPU12个形状均有六个相关contrast/noise变体；CUDA只有预选3形状且是CPU子集，不增加独立样本量。现有shape-cluster摘要亦无GD胜：CPU对两TRF均0胜/8tie/4负，对matched data为12tie；CUDA对两TRF均0胜/2tie/1负，对matched data为3tie。后端不合并，不给显著性结论。

## Outlier 是科学结果，不应删去

中增量2.5的 `ring` 和 `curved_bar` 在CPU的clean/noisy四条记录是主要严重失败。以clean为例：

| 形状 | native NRMSE / 残差 | matched data NRMSE / 残差 | GD NRMSE / 残差 |
|---|---|---|---|
| curved_bar | .08742 / .0001253 | 3.34856 / .006093 | 3.40425 / .005753 |
| ring | .07759 / .00007667 | 3.31519 / .006939 | 3.39719 / .006731 |

GD较小的数据残差伴随更差材料RMSE；它没有修复匹配first-order基线的成像失败。两种TRF在这些目标通过joint门，而GD/data未通过。noisy同类失败也存在。相同ring的CUDA结果仍失败（GD clean NRMSE3.31982，native .07759），不能把该重复子集记录视为额外独立反例。

这些现象支持“有限预算优化路径与残差度量可给不同材料误差”，但没有严格局部极小、不同根、局部满秩、真实物理等价解或全局非唯一性的证明。约.005–.007的残差不是exact alternative solution，且本轮数据门未过。

**不能把这两个outlier归因成Gram proposal造成的独立失败。** 已有phase-best记录表明Data/GD在共同原数据stage8就已偏离：ring NRMSE约2.619、残差.1624，curved_bar NRMSE约2.796、残差.1264；同一阶段两种TRF约为.188/.199、残差.0732/.0772，随后恢复。差距因此先于GD专属的final24第一phase出现。准确归因是：本次Data-LBFGS coarse路径遗留的坏材料未被后续GD proposal/data correction修复；既不能说Gram单独制造了全部失败，也不能说它已解决共享coarse路径问题。该分析来自既有候选，未新增优化运行。

小的对象收益必须保留：除低clean curved_bar的PASS外，中增量digit8 clean的GD NRMSE=.143778、IoU=.79195，joint PASS；matched data=.153173、IoU=.74172，joint FAIL；native/reciprocal在该clean对象已PASS。digit8 noisy还有更强的门槛例外：GD=.142593、IoU=.78082为PASS，native=.152122、IoU=.74667以及data=.151988、IoU=.73203为FAIL，reciprocal为PASS。GD相对native/data的这些连续NRMSE差仍小于.02 practical tie下限。**“中digit8 native已PASS”须限定clean，不能覆盖noisy。** 因此结论是“不支持整体实用材料精度优势”，并非“每个对象、每种次指标都毫无改善”。

## 成本裁决

CPU低clean native平均5.659s、1855 full complex RHS，GD37.670s、16678 RHS；中clean native12.305s、3502 RHS，GD45.024s、20021 RHS。GD用更多完整物理迭代抵消了每次source-adjoint相对Full-J的便宜，因此**不能由单次导数RHS更少推出端到端更快**。

全CPU每method平均wall为native19.616s、reciprocal44.860s、data35.934s、GD40.036s；对应平均full complex RHS为4239.7、9425.8、16200、17617。GD确实比reciprocal在总体wall较少，但reciprocal的材料精度比较没有GD practical胜；此事实不构成统一效率优势。GD低clean也比matched data略省wall(37.670 vs41.067s)，属于应保留的局部正成本差异，而非跨strata收益。

CUDA每method平均wall为native13.131s、reciprocal43.287s、data19.455s、GD25.633s；GD对matched data全practical tie却更费平均wall/RHS。它是CPU-host SciPy加complex128 CUDA物理核、双worker共享同GPU的实际策略结果。CPU与CUDA样本组成、并行竞争和运行时刻不同，不据此宣称硬件speedup或精度等同性。

ledger明确保留全部360 inference记录：CPU288共119660 LU、3418740 full complex RHS；CUDA72共48321 LU、1316052 RHS。12×12 Gram小系统RHS分别740304、373152，不能再当作full-system RHS重复加总。生成36 LU/432 RHS、CPU/GPU旧公开数据DEV、导入失败和launch/setup未知开销另列，不能混入method均值，也不能把UNKNOWN当零。求解器wall、worker总时间、operation timing互有包含，不能相加；并行任务wall求和不是整批实际经过时间。

## 理论新增内容的审查边界

`THEORETICAL_EXPECTATION.md` 新增实际固定非零背景b的展开是必要纠正：本批控制m=b+A a，A→0时Q趋于Q_b、Gram趋于H_b，而非必然趋于真空E0/H0。该展开只作离线受控家族解释，未把隐藏b交给solver。inverse-Gram扰动界明确依赖最小特征值及可逆性；低名义A不自动保证弱传播或稳定源Gram。

实际已付费generation的Born/full描述量进一步限制标签：A=.4的相对差范围约.209–.291（median .248），A=2.5约.872–1.222，A=12约3.744–5.468。这不是新的精度测试，也不说明Born反演曾被执行；它说明本次“低”是材料增量分层，不能改称已认证weak-scattering/Born regime。弱材料极限只能作为有条件理论预期，而非本批低增量结果的已验证渐近假设。

可逆source mixing下raw GD目标不变的恒等式成立，属于源基底坐标不变性，没有增加信息。文中正确保留reference normalization/正则相对权重的例外，不能把raw恒等式扩大为任意实现不变性或统计最优性。

小噪声root-local协方差论证在A满列秩、内点、无正则、局部线性化、同方差各向同性噪声且估计器线性无偏的条件下成立：`L_W-L_0=D,DA=0`使交叉项消失，协方差差为`σ² D D^T`。这只是标准加权最小二乘边界；没有证明本批576维GT Jacobian良态满秩，也不能解释为本次有限预算/active-bound/正则算法的性能定理。现有实际结果与“GD无普遍优势保证”相容，但不会反向证明上述未验证假设成立。

## 允许的最终研究表述

本轮完整有限研究表明：在所冻结的低/中/高材料增量、12形状及配对噪声设计中，GlobalDual没有建立相对matched data-LBFGS或成熟TRF的实用材料精度优势；部分门槛通过与微小局部差异存在，但伴随额外成本与严重中增量失败。native TRF在这批数据中给出较有利的精度/成本结果，reciprocal参数化带来不同精度和显著时间代价；均不能据此推成普遍推荐定理。

因此关闭“GlobalDual非高对比度优越”的本轮经验主张为`NOT_SUPPORTED`，保留有限模型、有限预算、paired clusters和后端范围。公式新颖性、完整文献优先权及论文最终定位由主代理决定；本独立审查不授予新理论、全局覆盖、概率解集或posterior结论。
