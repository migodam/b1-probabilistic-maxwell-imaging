# R6 独立结果与报告审查

2026-11-03；Sol high。只读 `RESEARCH_REPORT_ZH.md`、主/FFT封存评分JSON、GT_FEASIBLE_WITNESS、matched工程成本、aggregate成本/曲线记账、推导与基线说明；不读取GT数组、不新增forward/adjoint/优化/action测试、不修改主报告。最终原创性及论文价值裁决由父线程负责。

## 总体判断

主结果、分层成功率、连续指标、A成本负担和FFT有限工程收益与保存评分一致。当前证据支持“两个已公开历史弯条的额外恢复、六个新对象未增加确认优势、FFT组合可运行且两例observed CPU降低”的有限结论。不能扩大为稳定胜过成熟Dual-CIE、相同完整轨迹免费加速、六个新对象上的FFT优势或任何唯一性/连续PDE/三维成像结论；主报告大体保持了这些边界。

发现一处需修正的事实措辞，以及一处建议补足的停止状态披露。均不需要新科学实验。

## P2：新对象生成与最终算法冻结的时间顺序写反

位置：审查版 `RESEARCH_REPORT_ZH.md:5`，称“六个新对象在算法冻结后生成”。

现有回执给出的先后是：`audit/GENERATION_RECEIPT.json:3` 的生成时间为 `2026-11-03T13:12:07.903326+00:00`；最终 `configs/primary.json:3` 及首冻副本的freeze时间为 `13:17:24.787862+00:00`。因此这些证据不支持“生成在算法冻结之后”。算法在首次求解/新GT访问前冻结、预测在评分前封存是另一件事，不能混写。

建议改为“六个新对象由独立生成端保管GT，算法在求解/新GT访问前冻结，所有预测封存后才评分”，或者提供另一个确实早于13:12的算法冻结证据再保留原句。已即时通知父线程。这个文档时间顺序错误本身不构成solver读GT或GT选择证据，也不要求重跑已封存结果。

## P2建议：明确历史42cb两条基线也耗尽3000步

位置：`RESEARCH_REPORT_ZH.md:41–48,73`。

报告准确写出薄环四臂达到3000上限及全84任务的76 PGTOL／6 iterlimit／2 ftol，但没有点名余下两个iterlimit来自历史42cb的Dual-CIE、Global-Dual。两者最终projected gradient inf分别 `3.8393516988e-6`、`1.6564335503e-7`，corrector均3000步；不能把它们当已到相同stationarity程度的基线。

建议在42cb段或73行补一句：“历史42cb的Dual/Global也达到3000步上限，其预算影响同样未排除。” 此处不会改变“42cb四臂均失败”或A没有确认优势的判断，但能防止读者把历史困难基线的有限预算失败误读为算法能力上限。主报告没有直接声称这些任务已收敛；这是建议补足，不是科学结果失效。

## 分层结果与失败类型：已核对一致

`RESULTS/scoring_primary/SCORES.json` 有84个固定预测、21个独立材料，四臂各21；FFT的4条是两个历史观察复跑，不增加独立材料。主报告表7–11的成功数一致：

| stratum | CIE | Dual-CIE | Global-Dual | A |
|---|---:|---:|---:|---:|
| DEV9_DIGITS | 8/9 | 9/9 | 8/9 | 8/9 |
| DEV6_DIFFICULT | 1/6 | 1/6 | 1/6 | 3/6 |
| NEW6_CONFIRMATION | 4/6 | 4/6 | 4/6 | 4/6 |

这里A的8/9是九数字成功数，**不是digit8或digit9失败**。实际未过联合标志的是digit5：A NRMSE约 `.1141759`、IoU `.666666667`，material/data通过但IoU未过；Global同例IoU `.693693694`，也是image gate未过；Dual `.709090909`刚过该既有阈值。CIE的digit7则NRMSE约 `.4514690`、relative residual约 `.1811428`，属于明显材料和数据失败。主报告13、117行区分正确，没有把IoU边界失败渲染成完全没有成像出该数字。

两历史额外恢复是FBG curved_bar `7b8d` 和R2旧弯条 `b20d`：A NRMSE分别 `.0468592`、`.0412751`，data residual `.0001535`、`.0009662`，联合通过；Dual与Global均失败。Global控制说明两例不能只归因为sourcewise改global PR。它们都是已经公开开发对象，不是两次新盲测突破；报告46行有明确限定。

R4已公开42cb、63ab仍失败。63ab的A NRMSE `.4047504` 低于Dual `.5818590`，但data residual `.1113841` 高于Dual `.0983922`；报告正确保留GT误差与物理失配的冲突。四臂FBG ring均停在约 `.4248` NRMSE，同样没有以单一形状名称虚构普遍收益。

新六对象各臂都是同四对象通过、同薄环/开缝壳失败；A对Dual NRMSE是3低3高，对Global是1低5高。全21对象对三个baseline分别都是13低8高，不能把这些小差/相同对象的成对多臂比较当独立大样本显著性；主报告63行已避免此错误。

## 可行真值目标见证与停止原因

GT_FEASIBLE_WITNESS来自独立评分端已有sealed m/Yclean/noise，不重新生成或forward；extra_LU/RHS均0。它只作为post-seal可行目标见证，不能回送在线候选选择。保存见证的物理模型/算子来源与主求解一致。

- 偏心薄环：GT Phi=`8.701694621e-6`，预测Phi=`.00484444–.00494972`，data residual为GT噪声残差的 `102.2328–103.3320` 倍。主报告约数及“未达到已存在低目标可行区”的解释正确。
- 开缝壳：GT Phi=`8.553296885e-6`，预测Phi=`.00381288–.00406722`，residual ratio=`91.1052–94.1193`，同样远离该可行目标。
- 所有四个成功新对象的四臂Phi均略低于GT Phi；公开prior和噪声下GT不必是regularized optimum。报告75行没有把低Phi当更高材料正确性的证书。

薄环四臂PG inf=`2.7672e-6–5.9015e-6`且3000步耗尽，不能给渐近收敛或差局部极小证明。开缝壳四臂按PGTOL停止，保存incumbent PG inf=`8.8162e-9–9.8355e-9`，因此“更清楚的差局部驻点证据”是合理的有限数值说法；仍不是严格局部极小证书。可行GT目标更低排除了把高失配端点直接解释成已恢复等数据非唯一解，**不证明数据唯一或排除同时存在非唯一性**。主报告71–75行正确保留该区别。

全84停止：76 projected-gradient条件、6 iterlimit、2 relative-ftol，和主报告73一致。所有辅助schedule都执行3000步；corrector的实际接受步/评估数不同。共同3000上限是统一额度，不是等实际算量、都执行满额度或“所有baseline已经足够收敛”。

## FFT费用、终点与信息边界

`scoring_fft/SCORES.json` 及 `analysis/FFT_ENGINEERING_MATCHED.json` 与主报告93–102表一致。四条都是同预指定两观察、同stage/box/prior/corrector配置；已有方法审查确认真实padded free-space action/adjoint与sourcewise PR，非只用IFFT生成dense列。

| 对象/臂 | CPU dense→FFT/s | NRMSE dense→FFT | LU dense→FFT |
|---|---:|---:|---:|
| digit1 Dual | 36.5124→34.1387 | .0965611→.0959243 |1320→1493|
| digit1 A |59.0188→49.7582|.0959817→.0964369|1298→1730|
| b20d Dual |39.4786→38.0157|.5107190→.3895708|1619→1974|
| b20d A |70.5818→54.4442|.0412751→.0404943|1342→1410|

A total CPU下降15.69%、22.86%，而aux CPU分别 `41.7637→35.5088`、`53.5261→43.9684`。两条A的corrector CPU也不同，且物理LU增加，所以整体数字不能解释成保持完全相同路径的纯辅助FFT speedup。报告已同时给aux和total费用并披露终点不同。FFT-A比同FFT-Dual的total CPU分别多45.76%、43.21%，约46%／43%正确。

b20d Dual FFT的Phi从`.00509610`升到`.00583240`，NRMSE反而从`.510719`降到`.389571`，两者仍联合失败。该反例确实说明objective与材料误差不同，也不能把新终点改善归因为FFT增加物理信息。S、DPCS和最终物理LU仍dense；报告没有整个pipeline O(NlogN)或大网格吞吐结论。

## 费用字段、维度与曲线

主评分直接求和：`157418 LU`、`3773976 actual_complex_RHS`，其中source_RHS=`1889016`、adjoint_RHS=`1884960`；FFT为 `6607 LU / 158376 complex RHS`。这些是576维complex线性系统与多源RHS；不是把2D Fourier transforms、局部4cell block RHS或系数个数混进物理RHS。全部84个stage3 receiver_rank为17，mf12 APCS rank559、材料576实未知量，与THEORY67–69一致。

主CPU四臂和分别为897.221917、932.812807、956.568711、1504.322188s；outer auxiliary CPU和分别394.984693、401.973627、466.369227、1049.638833s，corrector和492.120711、520.698389、480.513116、444.907471s。主报告106–113的四舍五入正确；A相对Dual总CPU约1.613倍。内部execution_receipt与外层aux timer有很小范围差别，不应因两层timer末位不同而误判算术错误，亦不应将其重复计费。

主队列两段elapsed和 `1140.370765s`，中间修复间隔另列；FFT queue elapsed `55.359759s`，都不是sum_task_wall。旧pre-forward失败、生成和接线检查的成本应随完整账本发布。审查时COST_LEDGER.md尚未落盘；其数据源AGGREGATE_COSTS已存在，不据临时缺文档否定科学结果。

`CURVE_ACCOUNTING.json` 明确原Phi只来自已付费host/阶段物理验收与corrector评估，没有每个aux update的完整物理Phi。CPU曲线的stage与corrector记录有不同起点，短offset未记录，标null而非编造对齐；主报告79行的“只连接已付费点”正确。解读曲线时应保留这个时钟范围限制，不把每段连接线当中间真实forward。

## 理论、原创新颖性与附文状态

推导的P(m)r乘积微分、原Maxwell全LU、固定source metric、native Fourier/DPCS affine约束和padded Toeplitz adjoint都与实现一致。A相对Global控制只有材料块公式改变，相对mature Dual还包含PR聚合差异；THEORY77行明确分开，因果归因没有把两者混为一次传播净效应。局部传播可能无法改变跨腔全局模式是报告125行标明的下一轮假说，当前数据尚未直接测量它，不能升级成失败定理。

已有FFT-TSOM/NIE/CIE身份被保留；缺少独立NIE重现不应说成其无效或不存在。当前报告把有限工程组合与历史弯条经验收益分开，没有给标准公式首创资格；本审查没有重新做文献穷尽检索，最终新颖性裁决仍由父线程完成。

非阻塞附文时效：`FFT_TSOM_INTEGRATION.md:3,78–90` 仍是运行器准备/等待父线程完整结果的历史接线说明，而主报告已完成4完整任务。建议追加完成状态及指向最终matched结果，保留早期wiring回执范围，不让读者误认为完整FFT任务仍NOT_RUN。它不影响已保存的四个完成回执与评分。

## 交付建议

修正生成/冻结时间顺序一句，补充42cb两条3000-step基线的预算边界，再保持当前“历史补偿候选、有额外成本、六新对象未增加优势、FFT两例工程收益”的结论即可；没有理由为本审查新增physics、重评分或扩大预算。

## 最终文档复核：上述三项均已解决

2026-11-03，父线程完成编辑后只读复核当前版本，保留以上原发现记录，不新增实验：

1. **生成/冻结时间顺序：RESOLVED。** `RESEARCH_REPORT_ZH.md:5` 现在明确独立生成、GT隔离下冻结求解，且直接披露生成回执早于最终config冻结约5分钟，不再声称反向时序。`BENCHMARKS.md:3` 同步说明；与已核对回执一致。
2. **42cb基线迭代上限：RESOLVED。** 主报告73行明确另两个cap是历史42cb的Dual/Global，并给出其约 `3.84e-6/1.66e-7` 的最终投影梯度数值；薄环四条cap和开缝壳PGTOL的区分仍保留。未把这些有限预算失败称为全体baseline已收敛。
3. **FFT附文完成状态：RESOLVED。** `FFT_TSOM_INTEGRATION.md:3,90` 已写完整两例四任务已执行、封存评分，指向最终matched结果，给出A两例 `59.02→49.76`、`70.58→54.44` 秒，并保留早期接线PASS不等于长轨迹完全一致的限制。

最终 `COST_LEDGER.md` 已落盘，分别保留主84、4FFT、生成、wiring、旧失败、评分与非物理工作。新增总调用算术为 `157418+6607+6+8+12=164051 LU`、`3773976+158376+72+96+36=3932556 complex RHS`；其中12LU/36RHS来自16维观察器小模型，其余物理为576维，不能视作等成本。旧失败的task与process计时明确重叠，queue elapsed、修复间隔与task wall和也被分开。此补充总账与主报告的主科学批次／FFT批次分别计数不冲突。

当前审查无未解决的事实性修正要求。最终有限科学结论保持：A是有成本的历史困难对象补偿候选；六个新对象未确认普遍优势；两例FFT实现有有限observed CPU收益并保留不同终点。最终新颖性与论文裁决仍由父线程负责。
