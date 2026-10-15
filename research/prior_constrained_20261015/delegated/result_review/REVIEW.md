# 独立有限证据审阅

审阅状态：**未发现会把本轮 NO_GO 改成 GO 的错误；保留有限搜索与实现偏差限定。** 最终科学判断由父 Codex 负责。

本次阅读任务指定的 `CODEX_NEXT_ROUND_ZH.md`、本轮 `PROTOCOL_FROZEN.json`、`E0_DECISION.md`、`METRICS.json`、`E1_DECISION.json`、`EXECUTION_DEVIATIONS.json`、`RESOURCE_SUMMARY.json`、`COST_LEDGER.jsonl`，以及 `src/prior_constrained/witness_search.py`、`prior.py`。按父线程追加授权，补读 `common.py` 的公式和结构门、其读取的 `STRUCTURE_CALIBRATION.json` 中三项阈值、A26 `src/a25/dataset.py` 的原 TRAIN 生成规则/`parent_material`/材料组装表达式及 `configs/a25/campaign.json`。读取的是生成规则源码与公开配置，没有打开 locked 材料、几何实例、真值或 clean-response 数据。只对已有 JSON/JSONL 做聚合和源代码静态核对；没有调用物理求解器、导入实验模块、修改协议或追加试验。

## 1. 现有结果只支持停止对照，不支持不存在竞争解

冻结要求至少两组 parent-disjoint、P0_EXACT、结构合格、`d_pred <= 2` 且两端各自 32 次固定噪声重放至少 16 次竞争的 pair 才能进入 E2。现有记录为零组；新 full 比较的最小 `d_pred=49.223622812870374`，缓存结构合格比较最小距离为 `44.826717731720585`。四个 anchor 的最近已验证 pair 在主门槛和 log(100) 敏感性下均无竞争重放行。

因此 `NO_GO_FOR_GENERATOR_COMPARISON` 与 `ONE_SUPPORTED_CANDIDATE_AMBIGUITY_UNRESOLVED` 的组合可保留，但后一状态只能解释为当前候选集/有限预算的证据状态。本轮没有完成充分搜索，没有排除其他参数邻域或未搜索 pair，不能报告全局唯一，也不能报告原 P0 下不存在竞争材料。

依据：`PROTOCOL_FROZEN.json:E1`；`METRICS.json:E0/E1`；`E1_DECISION.json`；`witness_search.py:150-176`。

## 2. 噪声、距离、ell 与重放公式一致

可见实现将 `sigma = RMS(F(anchor)-F(background))*0.1` 固定于 anchor，以对应 20 dB differential complex RMS；候选与另一 truth 端的重放均继续使用同一个 sigma，没有按候选幅度重估。`d_pred = sqrt(2)*||F(candidate)-F(anchor)||/sigma`，对应协议的 512 个实通道白化。重放生成 `y=F(truth)+sigma/sqrt(2)*(N_R+iN_I)`，每端使用同样预冻结的 32 个独立 seed；两端以及不同 pair 的重放相关性不能充当额外独立材料样本。

`common.py:20-21` 的 `pack` 按 source-major Re32 再 Im32 排列，`ell(pred,y,sigma)=vdot(pred-y,pred-y).real/sigma**2`。对于每复通道实/虚部各乘 `sqrt(2)/sigma`，这恰好等于 `0.5*||W residual||^2`；没有漏掉或重复乘 1/2。`TAU=sqrt(chi2.ppf(.95,512))`、`sqrt(2*ell)`、两端 noise gate 与 `|Delta ell| <= log(10)` 的调用关系正确。原 config 的 noise_rule 也规定 differential complex RMS、实虚部方差各为 `sigma^2/2`。

依据：`common.py:20-21`；`witness_search.py:12-14,34-48,112,123,131,140,150-153`；A26 `configs/a25/campaign.json:noise_rule`。

## 3. P0_EXACT 有明确参数域和原 rasterizer 重放入口，不应扩大到 PHASE 或先验密度

`prior.py` 明确把四类 family、连续参数范围、subtype `{0,1}`、ID 限制与振幅 `{0.3,1.2}` 分开核查；`verify_exact_endpoint` 在参数域检查后调用原 `parent_material` 重放并检验最大误差。`decode` 先取得公开背景再加入原生成器 pattern，固定相位检查作用于增量 `chi-chi_b`，没有把背景当作空气。局部 mutate 在原参数范围内产生新端点，`parent_seed=None` 不妨碍以完整合法参数重放；但不能把这种搜索后的参数分布称为原先验随机采样或赋予已知 posterior mass。

补读原源代码确认：原 `FAMILIES` 与本轮四类顺序一致，TRAIN 为 64 parents、split index 为 0，原 RNG 的 seed 结构、TRAIN 中 angle/center/scale/axis_ratio 和后续所有标量/数组范围及抽样顺序，与 `draw_train_row` 对应。原 `parent_material` 逐图除以非零 `value.max()`，再乘 `phase/abs(phase)`；config 的 phase 为 `(1,0.2)`、amplitudes 为 `[0.3,1.2]`，与本轮 `V` 和 `AMPLITUDES` 一致。原材料组装为 `chi=cache.chi_b+amplitude*parent_material(...)`，本轮 decode 相同。config 背景 weak=`0.1+0.04i`、ring=`2+0.1i`、环半径 `[0.2,0.35]` 也被本轮 background 直接读取。没有将 a1 字符串当成幅度 1。

P0_EXACT 路径比单纯 PHASE 更强，且 `phase_audit` 明确拒绝以相位检查建立 membership。此次源级核对补足了相位、振幅和参数支持的来源证据；E0 的 512 条精确重放、误差 0 仍是已有数值记录，本次没有重新运行重放。须区分两个方向：若某端点超出数值容差、违反已核对的原 P0 必要增量相位/非负/峰值归一化条件，可排除该端点属于原 P0；但通过这些必要条件仍不足以建立 P0 membership。单纯未找到生成器 retraction、且没有必要条件反证时，须保留 PRIOR_UNKNOWN/RETRACTION_NOT_FOUND，不能把检索失败写成不属于 P0 的证明。所有判断均需使用该 condition 实际幅度并区分 roundoff。

依据：`prior.py:15-18,35-90,93-159`；A26 `src/a25/dataset.py:14-15,50-89,127-186,394` 的 TRAIN/ID 路径与通用 rasterizer；A26 `configs/a25/campaign.json`；`E0_DECISION.md`；`METRICS.json:E0.TRAIN_generator_replay_rows/TRAIN_replay_max_error`。

## 4. 结构度量符合冻结定义；实际搜索有已披露的路径偏差

可见结构实现对归一化 rho 在 3 个尺度、5 个阈值上计算 occupancy、boundary 和 excursion-set 变化；continuous shape 项逐尺度去 DC，DC 单列。persistent 项为至少三个阈值的跨尺度平均 disagreement 达到 0.05，并单独报告组件/孔洞。四个保存的最佳 pair 的三项距离、persistent_levels 和 changed_cells 均越过冻结阈值；它们失败在预测距离与重复噪声竞争，而不是被错误地当作只有 DC 改动。

补读 `common.py:23-37` 确认实际 acceptance：descriptor 用 `rho_field(chi,bg)/amplitude`，逐尺度去均值；separation 对三个校准 minimum 全部使用 `>=`，同时要求 `persistent_levels>=3`、`changed_cells>=6`，changed-cell 差值阈值为归一化 rho 的 0.1；topology 和 DC 不参与单独放行。运行时读取的三项阈值逐值等于 `PROTOCOL_FROZEN.json` 内冻结副本。本次是公式与存档指标的静态核对，没有重新计算全部 39 次结构比较。实现从独立校准文件而非冻结文件读取阈值，若将来复用须保持该依赖不可变；此次未见阈值漂移。

协议原定 round 2 从 first-round local 端点继续；实现从 local+global 中选最小距离端点，三组 anchor 首轮没有合格 local。`EXECUTION_DEVIATIONS.json` 如实披露了这一偏差。应称“有偏差的有限自适应生成器参数搜索”，不能称完整执行预注册 PC-SP 或由此判断其性能上限；当前零 pair 下，该偏差不构成继续 E2 的依据。

依据：`prior.py:162-196`；`common.py:23-37`；`delegated/prior/STRUCTURE_CALIBRATION.json:recommended_minima_max_of_scales`；`PROTOCOL_FROZEN.json:structure/E1.round2`；`witness_search.py:114-148`；`EXECUTION_DEVIATIONS.json`。

## 5. 45 states / 680 RHS 与两次实现失败计费相符

独立聚合 `COST_LEDGER.jsonl` 的 STATE/RHS 事件，得到 45 个 state attempts 与 680 complex RHS，全部属于 E1；RHS 分解为 forward 360、receiver-transfer adjoint 320。ACTION 中共 45 次 factorization、10 次 receiver-transfer build，与资源汇总相符。不能把嵌套 ACTION RHS 再加到 RHS 事件总数上。

两次 JOB_FAILURE 分别是 scalar `.tolist()` 类型问题和重复 metadata `role` keyword，共收费 **3 states / 88 RHS**，而不是“两次失败等于两个 states”：第一次 1/40，第二次 2/48；完成启动为 42/592。总和 45/680，低于 E1 64/2048 cap。失败端点记录未被免费抹去。两次失败缺少独立进程计时，汇总使用 inclusive command upper bound；4.4133 秒应保留上界语义，不能作为精确纯物理用时。material-only 汇总为 1536/2048，未借用旧预算。

依据：`COST_LEDGER.jsonl` 第 20、53 行 JOB_FAILURE，第 21、54 行 JOB_BEGIN，第 598 行 JOB_END及全部 STATE/RHS 事件；`METRICS.json:E1`；`RESOURCE_SUMMARY.json`。

## 6. E2 停止是必要结果；新 prior-preserving N3 修正仍未运行

台账没有 E2 事件，`RESOURCE_SUMMARY.json` 的 E2 spend 为 0/0，`METRICS.json:E2.status=NOT_RUN`，`E1_DECISION.json:E2_authorized_by_evidence=false`。四臂 generator comparison、共享 C_P 的新效果、新 prior-preserving N3 修正及其 MAE/likelihood 权衡均没有本轮运行证据，不能写成通过、有效或优于旧方法。历史 E0 审计只支持其已记录的诊断范围。

可继续保存“本轮未建立重要竞争材料见证；停止 E2；是否另立研究任务由父线程决定”。不得在本轮 NO_GO 后追加训练、seed、physics 或降低先验/结构/noise 门槛。

## 7. 已执行旧版的正结果分支缺口；当前交付代码已阻止自动授权

**实际执行的旧版**曾用 `len(good)` 直接填充 `qualified_parent_disjoint_pairs` 并授权 E2，未显式检验跨 pair 的 material parent 集是否相交；global restart 使“来自不同 anchor”并不自动等于两端 parent-disjoint。局部 robustness 虽会收费运行并保存，但其结果未进入旧版 `E2_authorized_by_evidence` 判定。本轮 `good=[]`，这一旧分支没有触发，因此没有发生错误授权。

**当前交付代码**已由父线程修改并经本次静态复读确认：`witness_search.py:172-181` 固定 `qualified_parent_disjoint_pairs=0`、`E2_authorized_by_evidence=False`；如果未来出现原始合格 pair，只能得到 `PROVISIONAL_PAIRS_REQUIRE_REMAINING_PROTOCOL_CHECKS` 和 provisional count，不能自动打开 E2。它仍未实现 parent ancestry、prior-tail 和 robustness 的完整正结果资格审查，但已阻止未审查的计数冒充授权。该补丁没有重跑物理、没有重写旧数值证据，也不会把当前负结果改成正结果。

最终限定：`FULL_VALIDATED_NUMERICAL` 与很小 backward residual 只能表明当前离散计算结果；没有严格 forward-error bound、没有连续 Maxwell 证书，`posterior_mass=null` 必须保留。
