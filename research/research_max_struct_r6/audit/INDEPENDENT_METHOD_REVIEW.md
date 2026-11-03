# R6 独立方法与公平性审查

2026-11-03；Sol high。只读方法合同、公开清单、源码、localized diff及已存在的接线回执；没有读取 sealed_truth、GT材料/掩膜、clean Y或评分，没有运行新的physics、优化、action测试或恒等式测试，没有改正在执行的主队列。

## 当前结论

**未发现要求中断84任务主队列的算法、单位、归一化或公平性错误。** 目前可以称“成熟sourcewise CIE/Dual-CIE、Global-Dual控制与A-material在共同扩大日程/共同材料空间/共同corrector额度下的比较”。还不能称任何臂已收敛、已穷尽算法能力、A获普适优势，或同日程就是等实际算量。收益、泛化和FFT整体吞吐等待完整结果、封存及独立评分后判定。

## 公开输入、物理与单位

`BENCHMARK_PROTOCOL.md:7–9,31`、`configs/primary_manifest.json` 与 `data_public/cases.json` 的计数一致：21个case，每个4臂，共84任务；15个 `DEV_ALL_DISCLOSED`，6个 `NEW_CONFIRMATION_UNSCORED`。21对象不是84独立材料，更不能把源列、重复噪声或阶段检查点计作新材料。新6对象必须和15公开开发对象分列，不能把已披露形状混成fresh blind率。

manifest四臂对同case使用同public_files。求解入口 `run_case.py:39–56,105–114` 不加载case truth或labels，拒绝私有payload键与冲突的公开值，只用公开算子/观测；固定检查576 cells、12 sources、48 receivers、24×24。新6的公开记录使用已有Z-ALG-R3的算子，历史对象使用各原公开算子；共同之处是同模型/维度/算法合同，不是所有case必须共享相同源布局或同一Y。

`MaxwellModel` 直接包装冻结disk-pulse VIE：`j=diag(m)(E0+Gj), F=Sj`；m为实contrast，而不是epsilon或电流J混用。新观测生成回执明确相同离散、6次forward、6LU/72complex RHS；没有跨网格/连续PDE或实测资格。新噪声公开合同为 `E|noise|²=sigma²`、sigma=`.001||Yclean||/sqrt(48*12)`，与CN .1% RMS一致。生成回执不提供真值给solver。

公开case合同的host=.1、box[.02,20]、amplitude12、smooth/TV各1e-4、epsilon.01、difference width19.5，与硬编码冻结Regularizer及benchmark common相同。求解不使用amplitude12作目标缩放；它是后续材料RMSE的公开分母。跨case noise/discrepancy必须按各公开sigma与原合同处理，不能把不同历史噪声合同平均成同一绝对sigma。

## source Omega、prior与原Phi

`benchmark_solver.py:58–76` 一次性冻结 `wd_s=1/sum_r|Y_rs|²`、`ws_s=1/sum_i|Einc_is|²`，四臂共享同一源坐标中的实正权重。data/source penalty normalizations与最终原物理目标是不同对象；最终统一

```text
Phi(m)=0.5 ||F(m)-Y||_F²/||Y||_F² + R(m).
```

`run_case.py:144–146` 把physical_sigma置为 `||Y||`，和 `MaxwellModel.physical_score/data_value_gradient` 相同；没有将噪声sigma或辅助sourcewise wd混进最终Phi。`R(m)`对所有臂都是同一个真实材料prior，不是对变换材料另加ridge。mature `_readout` 直接用J/Q/ws；global标量readout与A seed先用Ls source whitening再用unit weights，是相同共同实材料目标的不同算术接线，不是逐源独立材料。

所有臂仍有576独立实材料参数；扩大mf只扩大辅助current的Fourier坐标，不新增材料自由度、GT信息或照明。各beta与材料box共同保持ordinary CIE合法正分母分支。

## 四臂与系数运输

`arm_config:66–80` 的实际配置符合协议：

| 臂 | PR current | material readout | current beta | material beta |
|---|---|---|---|---|
| cie_single | mature sourcewise | mature scalar | 1,.2,.2,.2 | 1,.2,.2,.2 |
| cie_dual | mature sourcewise | mature scalar | 1,1,1,1 | 1,.2,.2,.2 |
| global_dual | global PR | mature shared-real scalar | 1,1,1,1 | 1,.2,.2,.2 |
| a_material | ordinary global PR | R4 exact one-hop material | 1,1,1,1 | 1,.2,.2,.2 |

CIE“single beta”是每阶段current/material相同，不是全程固定某一个beta，也不是另行调参的最优CIE。Dual sourcewise与GlobalDual的PR比率/line search不同，这是被明确命名的方法差异；A应主要与同global-current控制比较作用位置，同时接受成熟sourcewise基线的整体比较，不能只挑较弱臂呈现。

四臂在同一次调用里执行 mf3/6/10各500 + mf12 1500。继承的 `fourier_apcs` 保留native Fourier coefficient metric和DPCS投影；各阶段按frequency keys零填充原A，未重设材料/电流初值。PR仅在自然阶段开始重置，符合各成熟实现。`coefficient_tol=0` 不会因正的系数变化值提前停止，当前没有wallcap；异常或浮点失败仍必须按实际记录。

mf12包含全部576 Fourier坐标，但APCS仍是DPCS正交补；系数表示冗余，DPCS null actions不应说成新增可观测材料信息。相同终点空间也不意味着固定3000步达到同一个最优解或相同收敛程度。

## 成熟基线、corrector与观察器

`audit/mature_cie_observation.diff` 和 `mature_cie_observed.py` 的数值更新、_readout、basis transport、目标和停止公式与FBG成熟原型一致。新增pr初值、restarted、stage copies、callback与计时均为记录；没有利用观察结果反馈调参。`OBSERVATION_EQUIVALENCE.json` 只是已执行的小尺寸观察等价证据，不是576网格的科学qualification或完整基线能力证明；本审查没有重跑该测试。

`refine_observed:164–198` 保留R4的scale1e6、L-BFGS-B、box、ftol1e-12、scaled gtol.01、maxls30、maxfun=4*maximum，统一maximum3000。snapshot callback只读取已经付费的最后state/value/gradient，并用相同x核对，没有新forward/adjoint。原始gtol.01是在1e6缩放目标上，所以不能把它误报成未缩放Phi梯度阈值.01。

`PaidStateObserver.state/physical_score` 的额外score调用在匹配m的缓存上返回，不改变物理state或LU；prior/数组复制/序列化有CPU和内存开销，纳入observed run成本。validation_records与原candidate数组/目标被入口核对。corrector从初始+全部真实阶段端点的最低原Phi出发；observer保留所有已付费corrector objective evaluations的最低原Phi，与R4 best-incumbent逻辑一致，允许保存line-search已求值的更优试点，不是GT择根。最终导出其已付费prediction，没有为finalization重跑forward。

accepted_history是优化器接受迭代，evaluation_history含line-search评估，两者不是相同费用单位；3000上限是接受迭代上限，不是3000次LU、3000次梯度或保证正好执行满额度。必须保留实际stop/iterations/nfev与projected-stationarity的可用信息，再讨论是否“足够strong”。

## FFT独立工程路径

静态检查 `fft_tsom_ops.py` 与 `cie_fft.py`，没有发现action或adjoint接线错误：

- `NativeFourierAPCS.apply` 是orthonormal IFFT后 `(I-VdVdᴴ)` 投影，adjoint先同投影再orthonormal FFT并取同mask。frequency顺序与成熟fourier_apcs一致；自然运输按key零填充，未换系数内积。
- `DiskPulseToeplitzGreen` 以实际affine Cartesian坐标、原disk-pulse Hankel核及原self项构造所有有限lags。每轴padding≥2n−1，负lag按padding取模，输入零填充并截取原box，保留有限非周期卷积；不是周期Green。adjoint对完整embedding谱取共轭，输入/截取互为伴随，次序正确。
- current固定材料时 `p=beta/(1+beta*m)`，`H=P(I-XG)U`；`H_adjoint` 实现 `Uᴴ(I-GᴴX)Pᴴ`，未交换非对易的m与G。`T=P(JD-X(E+GJD))` 与原标量CIE等价。
- `cie_fft:225–249` 保留sourcewise时的per-source PR+比率、非下降重置与per-source line alpha；A仍global。sourcewise material readout直接传ws（268–270），A走原one-hop readout，不会因FFT把sourcewise baseline偷换成global更新。
- U/GU/SU/H均改为actions与adjoints；S、DPCS low-rank投影、A的局部block数据及最终dense physical LU仍保留。可以称“辅助current matrix-free FFT接线”，不能称整个成像pipeline matrix-free或O(N log N)。

FFT只省略历史APCS的roundoff SVD重构；这是已披露的浮点扰动，不保证跨几千步的非线性路径逐位一致。现有 `FFT_WIRING_RECEIPT` 仅检查已付费单步/action/adjoint，不能拿其小步wall宣布整体speedup；两固定case的完整同日程dense/FFT运行及终点差异仍需实际记录。

## 非阻塞记录限制与后续报告要求

1. **每步辅助算术观察的粒度不完全一致。** mature callback提供operation_counts；R2 global与R5 A callback没有该字段，`run_case.py:131–134` 会在相应HISTORY的observed_cost记录null。物理LU/RHS、CPU/wall和最终各臂auxiliary cost仍存在，不影响算法或当前公平性；不要声称已经有四臂同粒度的逐步辅助算术ledger，也不需要为此中断/重跑科学队列。
2. **FFT sourcewise的`material_readout_source_whitening_products`字段是继承计数标签。** `cie_fft.py:268–275` 的sourcewise分支直接用ws，未真的做J@Ls/Q@Ls，却仍将该字段加2。这个辅助label不能当两次实际whitening矩阵乘的测量值；可在成本解释中注明直接加权分支，物理LU/RHS和外层wall不受影响。无需为修标签重跑physics。
3. 最终账本要包含生成6LU/72RHS、已经存在的观察等价与FFT wiring等研究准备/失败成本，并与84任务主科学预算分开；嵌套aux/corrector/observer timers不能重复相加。queue elapsed与sum_task_wall不同；四并发单线程observed wall包含竞争，CPU是有用补充。
4. 当前`load_public`的路径黑名单不含字面`sealed_truth`，但已冻结manifest没有任何该路径，payload key限制也存在；本审查未发现GT暴露。权限隔离不能仅依赖此黑名单，保持冻结public manifest与独立评分端即可。此为未来入口防护的限定，不构成本次GT泄漏证据或新gate。
5. 主任务尚在执行，本审查不读取partial scoring来宣判A成功/失败，也不因某条baseline路径较慢就提前改budget。完成/失败/停止及剩余资格由实际回执决定。

## 可采用的准确表述

“本轮在21对象（15已公开开发、6独立新合成对象）的同材料空间和统一扩展日程下，比较成熟sourcewise CIE、sourcewise Dual-CIE、Global-Dual和A-material；各自接同一较长原物理corrector，并报告实际费用与停止原因。该协议增强了基线预算，不预设任何方法已收敛，也不以固定迭代预算代替等实际成本。FFT是另列的两case工程实现比较，科学恢复结论仍依完整原物理结果与独立封存后评分。”

## 追加审查：公开material_bounds冲突修复

2026-11-03，读取 `audit/PUBLIC_LOADER_REPAIR.json`，比较初始 `audit/source_freeze/experiments/`、修复后 `audit/source_freeze_after_loader_repair/experiments/` 与当前两个loader。没有读取GT或新增physics。

**无新的算法／公平性阻塞。** 原新对象第一任务在合并public文件时，发现legacy operator metadata的material_bounds `[.5,20]` 与公开observation／冻结numerical box `[.02,20]` 不同，于forward之前失败。旧失败 `RESULTS/primary/obs_9775e164fc894239__cie_single/RECEIPT.json` 保留，physical_model_cost为null，报错发生于load_public而非MaxwellModel构造；修复回执明确0LU/0RHS与`.497679s`失败process wall。不要把null成本误说成已完成一个零误差物理任务。

实际两份loader diff只有material_bounds键的later-public-file precedence例外：`run_case.py:52–56`、`run_fft_case.py:56–60`。G/E/S/Y等其他冲突仍抛错，私有payload检查仍在。冻结manifest文件顺序先operator后observation，因而此次确实采用较后的公开observation metadata，而非任意改变numerical box。MaxwellModel不使用该metadata来构造材料约束；arm_config与corrector仍从未改的冻结配置使用 `[.02,20]`，prior、source normalization、beta、日程、初始化及选择均不变。修复是公开输入适配，不能描述成根据成败调整算法或物理数据。

只读比较确认修复后两个loader与 `source_freeze_after_loader_repair` 副本一致；未作hash检查。`primary_final_manifest.json` 保留原84个task_id与全部case/arm/cohort/public_files，仅首个失败任务的output_dir改到 `primary_retry`，其余83条仍指原路径；24条新对象任务与completion manifest完全一致，60个已完成历史任务也继续指向原路径。旧失败目录没有被覆盖，新旧执行来源可明确区分。失败task自身receipt计时`.382559s`，修复审计记录的`.497679s`是外层process wall，两者不能重复相加。

最终计数须以final manifest的84个科学任务为分母，另收费旧失败尝试；不能把“原60+completion24+旧失败1”误报成85个独立方法任务。原seal/freeze说明也应连同本修复记录发布，不能笼统称所有执行文件从首冻起未变。此项不要求重跑前60个已完成任务或暂停剩余队列；后续待完整评分与报告再作结果层面的独立审查。
