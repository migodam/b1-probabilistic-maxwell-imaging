# A/A+ 与成熟CIE：高对比度基线只读审计

审计结论：**A在部分历史高对比度对象上确有受控制的材料路径收益；A+主要降低同一材料读出的计算量。现有高对比度结果不支持普遍材料优势，也没有发现足以推翻R8主表的核心CIE求解器代码错误。** 更具体的问题是含噪DPCS硬约束、早期current表示及校正前候选选择；这些有实际量化证据，但不能合并成“已解释全部失败”。科学结论由父线程裁决。

本审计只读取R4–R8源码、既有评分和已付费回执。R8原test现为 **PUBLIC DEVELOPMENT**，以下任何案例都不能重新称盲测。读取630行R8主评分及630份完整回执，缺失0；另读取R7的180行、历史40行及R3控制。仅对4个公开receiver矩阵做小型SVD及噪声代数统计；**新forward/adjoint/LU/反演/训练/GT文件读取均为0**，不发布、不改求解器。原始完整数组来自只读的相邻 `r8_worktree`；本轮只写本目录。

## 1. A与A+究竟改变什么

记辅助current为J，总场Q=E+GJ，本构缺陷r=J−diag(m)Q。成熟CIE使用逐cell标量变换 `βr/(1+βm)`；A将其替换为固定2×2邻域的一次传播：

```text
b_i = β − G_ii
z_i = m_i/(1+b_i m_i)
u_i = β r_i/(1+b_i m_i)
R_A = (I + diag(z) G_off) u
loss_A = .5 ρ tr(R_A* R_A Ω_s) + R(m)
```

这里完整Maxwell forward没有被一次传播截断；一次传播只用于辅助材料损失。A的材料梯度包含变换对m的导数。A+保留同一个loss_A、同prior、box、DPCS/Fourier和最终原物理corrector，以Möbius条件结构替换通用L-BFGS材料内层：每次最多2个four-color sweep，候选为状态项二次驻点、box端点和旧值，用完整prior-MM评价，再做完整固定current损失检查与回退。它不是正则化条件全局最小，不保证原Phi或GT误差下降，更不是增加完整物理阶数。

配对必须正确：

| 比较 | 可主要隔离的差别 |
|---|---|
| A vs A+，均global PR | 同一个onehop材料损失的内层求解方式/工作量 |
| A/A+ vs Global-Dual | 同global current下的材料读出损失/更新 |
| sourcewise A+ vs Dual-CIE | 同逐源PR下的材料读出损失/更新 |
| global A+ vs Dual-CIE | 同时改变材料读出与global/sourcewise current，不能单独归因A |

代码入口：R8 `src/cie_fft.py`、`src/benchmark_solver.py`，R7 `src/aplus_readout.py`，R4 `src/cie_onehop.py`。R8 common box[.02,20]、host.1、rtol.001；mf3/6/10各500＋mf12×1500；material beta1/.2/.2/.2，Dual/Global/A/A+ current beta固定1，single随material beta。辅助停止阈值0，完整3000步；原物理corrector最多3000，direct最多6000。最终GT指标定义/预算不能跨历史合同混合。

## 2. 正控制与必须保留的负控制

R4旧弯条b20d的同数据开发控制如下。不同原版本的有限数值差别保留在原表，不把所有历史行当同一预算或独立新对象。

| 控制 | 幅度NRMSE | IoU | 说明 |
|---|---:|---:|---|
| global普通＋共同corrector | .38918 | .42268 | 强机制对照 |
| 仅diagonal self | .38942 | .42268 | 自项本身没有复制收益 |
| 一次邻域交换 | .04167 | .93976 | 旧正例 |
| 等幅π反相邻域 | .55870 | .25781 | 相同拓扑/幅度，反相失败 |
| sourcewise普通扩到8000辅助步 | .55767 | .25581 | 增加普通模态/步数未恢复此例 |
| 32² onehop / global | .05307 / .05310 | 均.93548 | 网格变更后比较优势消失 |

R5位置消融更直接：普通global current＋onehop材料恢复(.04148/.93976)，onehop current＋普通材料不恢复(.55851/.25781)。这定位了旧正例中的材料更新作用；不能归为新观测、更多材料容量或current传播必需。自然阶段handoff、平移分块平均均没有普遍保住收益。32网格同时改离散、块物理宽度和噪声实例，不是严格PDE收敛结论。

R7原30对象A/A+均21/30材料通过，A+−A平均NRMSE约+.000604、配对区间跨0；A+降低已观测费用，但没有预注册非劣margin。后选家族20例A+/Dual/Global/single为12/9/9/8，曲条9/8/8/5，壳体3/1/1/3；壳体成功集与single相同。该条件正证据不能删除，也不能并入原30或推广到全部高contrast。

按相同幅度≥12口径单列R7的18对象，A/A+各11/18、Dual/Global各9/18、single10/18；A+平均NRMSE.28700、Global.35865。这保留了有限高contrast正证据，但其采集/噪声/家族合同不同于R8高26，不能混分母或称相互矛盾。见 `R7_HIGH_CONTRAST_SUMMARY.csv`。

R8旧四对象回归仍保留两条弯条正例，而42cb和40ac共同失败。NN接global/sourcewise A+保住旧bar，但NN+Dual没有自动复制；这是互补的路径证据，不是新42对象确认。完整历史逐行值见 `HISTORICAL_CONTROL_SCORES.csv`；R4 controls见 `R4_CONTROLS.csv`。

## 3. R8统一初值：高contrast分层的实际结果

这里只看uniform，避免把NN先验混入A读出比较。高组按既有离散幅度12/18分组，共26对象；低组幅度2.5/6共16。该分层是描述性的：噪声、孔径、家族、生成网格也随条件变化，不能把差异解释为contrast的独立因果效应。成功沿用原幅度NRMSE≤.15、Otsu IoU≥.7及合法box；CN99另列。

| 方法 | 全42成功 | 全42平均NRMSE | 高26成功 | 高26平均NRMSE | 高26 corrector CAP/失败 |
|---|---:|---:|---:|---:|---:|
| single CIE |13|.41783|6|.45692|13|
| Dual-CIE |10|.43978|3|.49293|12|
| Global-Dual |12|.42958|5|.47994|12|
| A |12|.42231|5|.46258|13|
| global A+ |11|.42368|4|.46489|11|
| sourcewise A+ |12|.42356|5|.46507|10|
| direct uniform（不同预算） |5|.28629|0|.34062|1|

关键配对：高组global A+−A平均NRMSE **+.002307**；A+没有独有成功，A有1例。A+−Global平均−.015050，但中位只有−.0000713，A+无独有成功、Global有1例。因此均值改善不能写成稳定高contrast成功优势。sourcewise A+−Dual平均−.027857，中位+.0000501，有2个标签增益，其中0006为大幅连续差异、0029主要是IoU边界；仍不是逐例占优。

两个具体反例防止误读：

- **0006，contrast18、noise.001：** Dual NRMSE .42249失败；Global .05012、A .04989、global A+ .05041、sourcewise A+ .04957均恢复。global A+对Dual的改善不能全部算材料读出；普通Global已能恢复，但sourcewise A+也展示了相对逐源Dual的路径变化。
- **0029，contrast18：** A/Global/global A+的NRMSE分别.11728/.11817/.11861，都较低；IoU .70588/.70588/.69903，A+仅因原IoU门失败。不能把这一标签差当灾难性连续材料损失，也不能为保住成功数改阈值。

A+内层工作减少是真实且可分硬件读出的：Mac22对象A/A+平均aux CPU63.18/26.52s，材料inner evaluations总2008761/329257；XINAN20对象127.56/68.92s、1272701/299747。完整管线均值Mac105.27/70.78s、XINAN234.38/173.27s。它们是共享主机既有费用，不是隔离硬件基准，也不证明材料非劣。见 `R8_A_READOUT_COSTS.csv`、`R8_MATCHED_PAIRS.csv`。

## 4. DPCS硬锁与噪声：具体严重到什么程度

R8的四个公开算子只改变source数/Tx孔径，receiver矩阵S相同。rtol=.001均保留 **rank17**；最弱保留奇异值/最大值=.002538169，首个排除=.000905282，保留条件数 **393.9847**。这不是根据对象GT回选的秩。

辅助current固定 `V_r*J=Σ_r⁻¹U_r*Y`。声明独立圆复Gaussian噪声每复通道方差σ²时，固定S下噪声电流平方范数期望为 `Ns σ² Σ_i≤r 1/s_i²`。本次只用公开S/Y/σ计算其平方根相对于已观测JD范数的比例：

| 相对复RMS噪声 | 对象数 | min / median / max |
|---|---:|---:|
|.001|22|.0431 / .0659 / .0690|
|.01|8|.4369 / .5024 / .6002|
|.03|12|.8260 / .9109 / 1.0822|

**这些是固定S的预期RMS量，不是当前实际噪声实现的误差，不是JD/材料GT误差，也不能直接决定新秩。** 分母JD本身含噪，超过1不矛盾。它说明把rtol=.001称“噪声很小所以无碍”缺乏依据；不能反过来断言rank17必错或已经知道最佳截断。

R3在公开b20d真材料处已有更强的固定锚点证据：保留rank17，DPCS噪声current mismatch范数.481806，精确最小source-weighted状态代价.0289064；声明噪声期望.0226618、99%上界.0424996。约束最小能量形式为 `b*(D H⁻¹ H⁻* D*)⁻¹b`。它不是最终物理corrector材料误差下界，也不能拿不同m处的原阶段state loss作“解释百分比”。clean-native仍NRMSE.38967、residual.10782，strong-release仍.49048/.10066；**噪声硬锁不是全部失败的充分解释**。证据见 `R3_HARDLOCK_EVIDENCE.json`、`R3_DEVELOPMENT_CONTROL_SCORES.csv`。

## 5. 辅助端点、选择与corrector分别做了什么

R8所有uniform辅助路径完整运行3000步，当前全局coefficient_tol=0，并非旧45秒wall截断基线。每个initial/stage端点均用原Maxwell Phi验证；选最小Phi再接corrector，末端导出所有已付费物理调用的Phi incumbent。未按GT取best@K。

高26中，Dual/Global各有19例选择最后stage3、5例stage2、2例stage1；A为20/3/3；global A+为19/4/3。uniform均未选stage−1，不能声称它们完全没用辅助。NN高组Dual/Global/global A+/sourcewise A+分别有4/4/3/2例选stage−1，这些例子的最终改善主要来自保留NN材料再物理校正，不能算辅助current恢复贡献。

每支高组26例corrector都进一步降低Phi；global A+的选中pre-Phi/final-Phi中位比约30.48。但uniform global A+有15例报告数值成功，其中11例材料标签仍失败；Dual14例数值成功中12例材料失败。NN direct高组26例均数值成功，仍13例材料失败。PGTOL/relative-F/CAP不等于材料证书。

**0027是现成的选择反例：** 同一个NN初值Phi=.309653；global A+选择校正前Phi约.01565的辅助候选，最终3000步达到上限，Phi=.000463504、NRMSE1.013697。NN direct从同一初值按PGTOL结束，Phi=.000428659、NRMSE.120818、IoU.84967。较低pre-Phi没保住更好corrector终点。这不证明某个全局盆地、障碍或新选择器；两种支路预算/成本不同，不能用GT在线改选。

原paid wide表在高26的NN global A+中有13例满足“选中aux pre-Phi低于NNdirect初值、但最终Phi和材料都比NNdirect差”。这是已发生的选择/轨迹限制，不能推出固定偏置或保留current一定修复。NN高组global/sourcewise A+分别16/26最终NRMSE劣于其NN初值，而NN direct为4/26。普通Dual/Global也出现同类问题。具体逐行候选值、最后辅助值、选择stage和corrector停止见 `R8_AUXILIARY_CANDIDATE_SCORES.csv`、`R8_SELECTION_AND_CORRECTOR_630.csv`、`SELECTION_HIGH_CONTRAST_SUMMARY.json`。

## 6. 已证实的问题、已排除的bug与未验证假设

| 分类 | 已有证据 | 不允许推出 |
|---|---|---|
| 算法可行族限制 | 含噪DPCS固定；低阶Fourier和零APCS会改变初值current；定量投影诊断已保存 | 核心梯度有bug，或降低rank一定成功 |
| 选择/停止限制 | pre-Phi排序不保证corrector后排序；CAP和错误材料的数值停止均真实存在 | 已认证局部极小、充分预算上限，或单纯再加步数必能救回 |
| 对照混杂 | global A+ vs逐源Dual同时改变current；single有独立正例 | 挑Dual最弱结果代表全部成熟CIE |
| 历史报告错误（已记录） | 旧FBG E3 scorer通用NRMSE文字“8或18”与实际amplitude12不符；汇总cap/stop遗漏原receipt预算stop | R8评分公式或当前solver同样错误 |
| 历史基线时长限制 | FBG E3 single9/12、Dual11/12 current被wall截断；不是成熟性能天花板 | 用旧短端点制造A优势；也不能据此否定R8完整3000步 |
| 被排除的代码怀疑 | run_case初始state_dof看似固定，但finally用实际n×Ns覆盖；630最终receipt均一致 | 存在六源数值按12源运行的bug |
| 本审计未证实 | 没有新的梯度、adjoint或核心求解器bug证据；未做数值rerun | 把“没找到”称形式化无bug证明 |
| 待验证假设 | 噪声兼容DPCS、状态传递或更好支路选择可能有帮助 | 已证明任一修复能恢复高contrast、或可按GT调rank/β |

R4旧正例的离散/分块敏感性、R3去噪失败和R8 Global正控制同时约束解释。不以一个理由替换全部失败。Fresnel物理背景修正属于另一个测量合同，不移植来指控本次合成数据中匹配的[.02,20]背景。

## 7. 建议冻结的六个最小开发案例（本审计未运行）

| 固定案例 | 角色 | 必须保留的比较/问题 |
|---|---|---|
|r8_test_0000|正常control，contrast2.5/noise.001|single、Dual、Global、A、A+均应保持；防止修复破坏正常例 |
|obs_b20d5a55677d4eca|旧A-positive|保留global、self、phaseflip、material-only的历史控制，不重复包装成新独立对象 |
|obs_42cb233e4f224bb3|共享失败|不能只保留旧正例；长corrector与GT可行Phi证据仍有效 |
|r8_test_0003|3% noisy shared failure|按公开noise=.03中的最小ID固定；量化噪声约束问题，不按新GT挑rank |
|r8_test_0006|current模式混杂control|Global普通已恢复；同时保留sourcewise A+，区分current与material作用 |
|r8_test_0027|NN状态/选择counterexample|固定原NN与direct分支，检验是否真的保住原物理支路，而非只降低pre-Phi |

均为已暴露DEVELOPMENT；后四个R8案例的public文件和冻结对照列表已写入 `FIXED_MINIMAL_DEV_CASES.json`。它们是交给父线程的固定建议，不是本轮执行授权；不建议以现有GT挑一个新在线秩、阈值、β或停止点。若另立后续实验，应先固定一个明确机制及其可失败预测，而非把六案例变成搜索网格。

## 证据入口与执行边界

- `EVIDENCE_LEDGER.json`：逐项证据等级和原路径。
- `R8_PAIR_SUMMARY.csv` / `R8_UNIFORM_STRATIFIED.csv`：高/低contrast、幅度、噪声、cohort及硬件分层；没有跨硬件加速推断。
- `DPCS_PUBLIC_OPERATOR_SPECTRA.json` / `DPCS_CASE_NOISE_SEVERITY.csv`：全部4算子和42对象的谱/噪声量。
- `R8_SELECTION_AND_CORRECTOR_630.csv` / `R8_AUXILIARY_CANDIDATE_SCORES.csv`：已付费选择/校正轨迹入口，不是新GT评分。
- `R4_CONTROLS.csv` / `R7_ENDPOINTS_180.csv` / `HISTORICAL_CONTROL_SCORES.csv`：历史正负控制，不混分母。
- `ANALYSIS_RECEIPT.json`：只读分析费用及零物理调用声明。小型读取/整理另有未独立计时开销，不声称完全闭合所有人工审阅成本。

本目录的两个脚本只读CSV/JSON/公开S/Y/σ，没有导入任何Maxwell求解器或GT评分器；没有新网络、solver、GT在线选择、SHA核验或publication。审计中一次路径glob展开和一次临时摘要脚本语法错误被修正，均没有物理调用，也未改变原结果。

## 父线程后续四例的历史基线可复用性核对

父线程随后选定b20d、40ac、0000、0027开展独立DEV噪声/秩比较；这是父线程的新执行，本审计未运行。对b20d/40ac，旧uniform普通single/Dual/Global均来自R6显式Fourier与APCS roundoff-SVD路径；它们与R8同m0=.1、box/prior、rtol=.001、3000aux及3000corrector cap，但不与R8 nativeFFT逐步实现完全同版。R6 Dual的actual corrector为b20d1598、40ac2599步，Global为2003/2607步，均PGTOL；不能把nominalcap当实际相同工作。R8同版uniform sourcewise A+两对象都存在，R7同FFT global A+也存在，但它们不是ordinary fixed_001对照；NN普通Dual也不能替代uniform。

**因此，为严格只改noise/rank的R9普通CIE比较，b20d/40ac应新运行同版uniform fixed_001对照。** 旧端点保留作历史正负reference，不混作严格配对或新计费。确切绝对receipt/prediction路径、配置、数值schedule、source entrypoint与actual stop均在 `HISTORICAL_BASELINE_REUSE.json`，已及时通知父线程。
