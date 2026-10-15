# E2：已保存 full-Maxwell 数组的失败分析

本文件只分析已付费数组，不调用 solver、不新增 RHS，不更改原 E2、A30 或失败结果。统计单位为4个 parent，各 weak/ring 两条件。全部是历史已暴露 development；不是独立验证。

## 主要结果

**没有找到比普通 multi-start 更多的结构解族；有限编辑补偿没有把任何原先不相容的编辑修正成新相容轨迹。Best-of-4 在所有具有合法候选的三个条件上与普通 multi-start、weak/null、无补偿均相同；其余5/8条件没有合法解，不能删掉后报告总体成功。**

有限编辑 archive 计数为`[0, 0, 2, 0, 1, 0, 1, 0]`；multi-start 和 weak/null 都为`[0, 0, 1, 0, 1, 0, 1, 0]`。有限编辑多出的一个是**材料场 RMS 去重**意义的额外成员；同样成员在无补偿的原始编辑中已经存在。它不满足预注册的至少6个结构 mask 改变单元，不能称新结构族，更不能归功于 nonlinear compensation。

|条件|base白化残差|合法去重材料数 finite / multistart / weak / no-comp|finite结构族数|三方法及no-comp相同的Best4复材料MAE|
|---|---:|---|---:|---:|
|p250065_weak_a1_snr20|24.1987|0 / 0 / 0 / 0|0|NOT_DEFINED|
|p250065_ring_a1_snr20|27.0654|0 / 0 / 0 / 0|0|NOT_DEFINED|
|p250069_weak_a1_snr20|20.9292|2 / 1 / 1 / 2|1|0.04351438|
|p250069_ring_a1_snr20|26.3769|0 / 0 / 0 / 0|0|NOT_DEFINED|
|p250073_weak_a1_snr20|20.8919|1 / 1 / 1 / 1|1|0.03477679|
|p250073_ring_a1_snr20|25.1159|0 / 0 / 0 / 0|0|NOT_DEFINED|
|p250077_weak_a1_snr20|22.9950|1 / 1 / 1 / 1|1|0.04845289|
|p250077_ring_a1_snr20|26.3348|0 / 0 / 0 / 0|0|NOT_DEFINED|

相容阈值固定为sqrt(χ²_512,.95)=23.785449004；原 sqrt(2)/σ、8源32接收器、固定gain1、材料box均未放宽。base只有3/8满足相容性，所有ring条件都失败；finite-history只在1/8具有第二个去重材料，结构族仍均不超过1。四个parent没有一个在两个背景条件都具有合法候选，所以全parent双条件 Best4 均为NOT_DEFINED。不得将仅weak成功的三个parent写成完整paired endpoint。

## 额外材料到底是什么

唯一额外有限编辑 archive 成员见`RESULT_SUMMARY.json/new_archive_members`。它来自 p250069_weak 的第三条轨迹，初始 full residual 已经低于相容阈值，无需补偿。该材料与base的 field-RMS距离达到.03去重标准，但阈值结构 mask 改变小于6格。其复材料MAE比base更差；TopK离线择优仍选择base。于是“多一个材料成员”不等于“多发现一个重要结构”，也不等于“成像改善”。no-comp保留同一成员，是直接消融反证。

三个有合法候选条件的 Best4 描述性均值为0.0422480213，各方法相同，改善0。这个数仅描述三个已成功条件；5个无合法解的MAE没有用0、无穷或base偷偷填补。Coverage@4按原离线field-distance≤.1定义，见完整表及JSON；该coverage门槛与材料去重.03、结构6格是不同定义。

必须保留一个去重策略诊断：如果离线查看**全部合法trial、尚未按.03去重**，multi-start在p250069_weak的一个trial使MAE从.04351438到.04308716，改善.9818%；finite-history在p250077_weak的一个原始编辑从.04845289到.04803630，改善.8598%。它们被在线first-member材料去重并入base邻域，因此没有改变已冻结archive的Best4。这些微小的within-family改善不构成新结构或finite-compensation优势；前者反而提醒不能把“archive无改善”误说成“所有物理trial的材料完全没有改善”。未经预注册的truth择优只是事后诊断，不能回写在线筛选规则。

## 修正器确实改变了材料，但未完成目标

共200个小材料QP，191个报告收敛，9个达到迭代上限。full-forward修正产生真实非零材料变化，并非接口空转。有限编辑32条轨迹的初始/最终 residual 中位数分别为58.8559/36.6712，配对残差下降中位比例11.03%，但 invalid→valid 轨迹数为0/32。protected edit cells 最大缺陷0，优化器没有偷偷把洞填回去。

普通 multi-start 和 weak/null 的全局起点更远，经过两个修正步骤后大幅降低 residual；部分进入相容阈值，但都落在base的.03材料邻域内。它们没有发现额外结构。这支持“当前局部修正/表示预算不足以产生可验证新分支”的解释；不能从中推导这些方法的全局最佳能力。weak/null使用的是N3处材料tangent，且有限步还经过材料box投影，投影与有限非线性都会破坏线性弱方向性质；它不是任意receiver-dark current，更不是已经验证的有限材料null族。

Round1为16实参数低频DCT，Round2为32实参数multiscale Haar；只运行每轨迹两次 nonlinear correction。材料是576格的完整χ，源共享；chart不是一般高维材料空间的完备表示。A30旧sparse64的35.69%仍是truth-fit容量，不能套用到本轮16/32在线优化。五个无解条件的base本身不相容，说明本轮连共同粗重建前置条件也未充分满足。缺少可行base、有限参数表示、有限迭代和物理不可补偿这几种原因，在此预算内没有被完全分离。

## parent汇总与成本

`RESULT_SUMMARY.json/parents`先平均每parent的weak/ring。finite合法材料数的parent均值依次0/1/.5/.5；对照0/.5/.5/.5；finite唯一+.5parent差对应上文无结构收益的同一个原始编辑。结构族parent均值各方法一致0/.5/.5/.5。全8条件算法新增结构族数0；与multi-start相比新增0；best4收益0或未定义。没有用32条轨迹当32个独立parent，也没有为4个暴露parent编造确认性置信区间。

累计328 /384新full-state attempts、9214 /9216complex LU RHS，包含首次失败audit的1state/91RHS；受控执行wall 18.12985s。剩56state但只有2RHS，无法支持下一个8源full state。无需继续付费探针。GPU0，旧P2额度0，旧A30 ledger不改。

成功audit的实际finite secant、Schur、dual、材料mass-adjoint和c=χe误差均<1e−8，见ACTUAL_BACKEND_AUDIT.json。首次dual小二阶量相对归一化误差1.91668e−8失败保留；parent改为可解析1%扰动的dual试验，容差不变，不把失败抹掉。为保留91RHS重复费用，省略两个原拟定no-comp端点（最初两个条件各少一次）；完整运行no-comp30states而不是32。原PROTOCOL仍保留原计划，实际偏差以AUDIT_REVISION和ledger为准。五方法/全预算speed优越性没有建立。

本分析脚本另记`E2_POSTPROCESS_JOB.json`，0physics/0RHS；其耗时不混入已冻结的18.13s求解stage。

## 科学判断边界

本轮否定的是**已实施低维受限corrector和编辑组合在该8条件下的新增结构/Best4优势**。它没有证明完整合法材料空间不存在其他相容结构，没有证明A30所有局部法一定失败，也没有证明多源TM逆问题全局唯一。没有概率认证；candidate-only字段应为null。A3 learned policy未训练，有限history parent重用未证明收益。没有比较TSOM/ROM、反馈特化或learned correction，不能为它们归因。

E3未授权；这轮也没有提供由二维成功机制进入三维的科学理由。优先停止当前补偿优越性主张，保留失败作为反证。若未来另获预算，应先解决强粗重建与真正可行base，再用同预算独立起点和扩展表示分离优化瓶颈；不能靠扩noise/freegain或加大网络制造多解。最终机制裁决由parent Codex负责。

## 发布前修正：独立 offline process

实际完成的原 E2 首次执行是在全部搜索、proposal、corrector 和数据筛选冻结后，才在同一个进程读取 `OFFLINE_TARGETS` 做离线评分。没有发现真值进入在线生成或排序，但它没有满足“独立离线进程”的实现要求。这是执行偏差，不能把原实验描述成当时已经进程隔离。

现已将评分提取为 `research_a.e2_offline`。新在线 runner 只写入冻结搜索数组和history，再启动新的 Python 子进程；在线进程不打开真值文件。对已经保存的8条件运行了一次独立 scorer，新增物理0、RHS0、solver构造0，未加载a30物理adapter模块。32行METRICS JSON的数值及结构逐项完全相同，全部图像/prediction/truth数组值和键完全相同，详见 `OFFLINE_PROCESS_EQUIVALENCE.json`。

原同进程评分的METRICS和完整数组副本保存在 `legacy_offline_inprocess/`，同时记录原执行边界偏差。独立重放未重跑搜索，也未改变负结果或历史ledger。此次零物理处理独立计入 `E2_OFFLINE_REPLAY_JOB.json`；它不冒充原18.13秒物理stage的一部分。
