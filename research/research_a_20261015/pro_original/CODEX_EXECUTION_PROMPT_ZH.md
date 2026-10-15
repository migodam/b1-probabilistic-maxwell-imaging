# 给 Codex 的下一轮执行任务
# Research A — Maxwell Compensation-Driven Material Ensembles

你负责实现与执行一个有硬预算、真值隔离、可以得到 NO-GO 的研究试验。
先读同目录 `RESEARCH_REPORT_ZH.md`；它是本轮科学契约。
不要直接训练 Flow/MDN，不要启动新的大型 Maxwell bank，不要修改历史项目的 Gate。
本轮默认只运行 E0/E1；E2 为显式配置开启；E3/E4 默认关闭。

## 0. 项目定位与默认边界

目标不是完整 posterior：从一个合法基础材料重建，通过有限结构编辑与 nonlinear compensation，得到材料显著不同、full-Maxwell 数据相容的 candidate pool。
主指标是完整材料场的 Best-of-4 MAE、Selected-MAE、Coverage@4、物理不同解族数量和完整计算成本。
禁止用 event-probability MAE 代替材料 MAE；禁止用真值选方向、筛图或调算法。

默认：
- `max_stage: E1`
- `allow_gpu: false`
- `allow_training: false`
- `allow_new_response_bank: false`
- `enable_2d: false`
- `enable_3d: false`
- `enable_locked_campaign: false`
- P2 新增 full Mie `material x frequency attempt` cap=2304。
- P2 CPU wall hard stop=1800 s，不含人工阅读与报告，不承诺跑完。
- 任何预算用尽立即写 PARTIAL/STOP；失败、拒绝、cache miss 均登记。
- 旧Joint P2剩余541次不是本任务额度，不借用其他研究账本。
- 不得因为实验长、失败或候选不足而自动扩大预算、换noise或换support。

在现有工作区定位相关仓库。不猜本地绝对路径，不新建远端repository、不推送、不覆盖历史数据。
可在本地创建 `research/research_a_compensation/` 或独立本地工作目录，全部新增文件写在其中。
若依赖、权重或bank不在本机，列明具体缺失项；先完成可独立运行部分，不重造一套不一致后端冒充复现。

## 1. 必须读取的上游文件

A17 `migodam/a17-maxwell-certified-exchange`：
- closure_final/START_HERE.md
- closure_final/manuscript/CORE_THEORY.tex
- closure_final/reports/A17_CLOSURE_REPORT_ZH.md
- closure_final/manuscript/FINAL_REFERENCES.tex

A23 `migodam/a23-fast-material-imaging`：
- A23_EXECUTIVE_DECISION.md
- BAYESIAN_FEEDBACK_INVERSION.md
- OPM_REDUNDANCY_AND_MINIMALITY.md
- PHYSICS_NATIVE_PRIORS.md
- THEOREM_AND_COUNTEREXAMPLE_LEDGER.md

A24 `migodam/a24-source-anchored-om-feedback`：
- README.md
- A24_RESEARCH_REPORT_ZH.md（若存在）
- docs/SCIENTIFIC_REVIEW.md（若存在）

A30 `migodam/a30-maxwell-observable-material-manifold`：
- A30_E0_E1_DECISION.md
- 对应实验配置、noise packing、候选与材料场资产
- 不把Sparse-64的离线oracle改善当在线算法收益。

B1 `migodam/b1-probabilistic-maxwell-imaging`，branch `joint-p2-20261014`：
- research/joint_p2/RESEARCH_REPORT_ZH.md
- research/joint_p2/THEORY_HANDOFF.md
- research/joint_p2/METHOD_COMPARISON.md
- results/FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz、ACCEPTED_REFERENCE_AUDIT.json、MATERIAL_FIELDS.npz（以报告中实际相对路径为准）
- 只用accepted参考；不重新采用已被修正的旧积分结果。

C1 `migodam/c1-physics-state-generative-imaging`：
- p1/C1_P1_RESEARCH_REPORT_ZH.md
- P2对应的模型配置、冻结MDN权重和distillation target metadata。
- 若共同四维target的checkpoint缺失，不用P1两参数强度模型替代。

建立 `UPSTREAM_MANIFEST.json`，记录repo commit、文件hash、资产hash、历史暴露状态、硬件、数值库版本。
记录 `SOURCE_ACCESS_GAPS.md`；不得声称读过不存在/未取得文件。

## 2. 文件结构与接口

实现：
- configs/research_a_e1.yaml
- configs/research_a_e2.yaml
- src/target.py
- src/cache.py
- src/budget.py
- src/material_metrics.py
- src/structure.py
- src/edits.py
- src/compensate.py
- src/archive.py
- src/selector.py
- src/export_package.py
- src/evaluate_truth_offline.py
- tests/
- results/
- DECISION.md
- RESEARCH_REPORT_ZH.md

提供统一CLI：
`python -m research_a.run --stage E0 --config configs/research_a_e1.yaml`
`python -m research_a.run --stage E1 --config configs/research_a_e1.yaml`
`python -m research_a.run --stage E2 --config configs/research_a_e2.yaml`

路径和模块名可按原repo调整，但README必须给出可直接运行的实际命令。

核心接口：
```
Target.forward(material, nuisance, fidelity, budget) -> Prediction
Target.jvp(material, nuisance, direction, budget) -> DataTangent
Target.vjp(material, nuisance, data_vector, budget) -> MaterialGradient
Target.validate_material(material) -> MaterialValidity
Target.materialize(parameters) -> FullMaterialField
Target.physics_identity() -> immutable target manifest
```
所有full调用只能通过budget wrapper。cache key必须覆盖：
material/physical law、grid、自项、背景、source、receiver、frequency、polarization、packing、solver tolerance、nuisance policy及forward版本。
多源一个state和vector RHS分开计数；cache hit新调用为0，但历史setup成本仍写入cold ledger。

## 3. 真值隔离

generator/corrector/selector只能收到 observation、noise、acquisition、prior/support、基础解及预算。
真值、truth MAE、truth topology只由离线评分进程持有。

新truth可由独立数据生成模块产生，但向算法输出文件必须删除truth字段。
先冻结配置、选择器、候选集合及其hash，再运行离线评分。
在线模块禁止导入评分模块；增加自动测试防止任何truth字段通过函数参数进入在线路径。
检索旧oracle候选只用于独立离线capacity比较，不作为在线proposal训练标签或起点。
算法seed、noise seed、parent seed分开；noise/seed不得作为独立parent增加样本数。

## 4. E0：物理一致性与数学检查

复跑配套 `validation/verify_algebra.py`。
在实际backend中用最小可支付fixture检查：
1. j=A(chi)(Einc+Gj)，所有源共享同一chi。
2. finite deltaA使用a(chi1)-a(chi0)，不替代为a'(chi0)deltaChi。
3. relative identity:
   X=L0^-1 deltaA E0; T=L0^-1 deltaA G;
   deltaJ=(I-T)^-1 X。
4. 局部选择矩阵U：
   Q=L0^-1 U; K=U*GQ;
   deltaJ=Q(I-DK)^-1 D U*E0。
5. 标量材料一致性：
   a_i=sum(conj(E_i,s)*j_i,s)/sum(abs(E_i,s)^2)
   over所有向量分量和照明；
   同时检查正交缺陷，不是仅计算division。
6. 一个receiver-dark随机扰动不应自动通过material realizability。
7. dual correction:
   rp=b-Ljhat; rd=S*-L*Zhat;
   y-ycorr=rd*L^-1rp。
8. proper whitening与原数据一致，P2和A30约定不要相互替换。
9. 一参数global材料更新的matrix rank通常不是1。
10. 相同材料的inactive latent表示必须dedup。

本包独立八偶极检查不是实际repo backend重现。实际backend新增fixture调用必须计费。

## 5. 指标与候选集合

材料指标使用完整场和体积权重：
- real/imag MAE
- normalized volume MAE（固定动态范围，不依赖truth）
- complex-modulus MAE
- support IoU
- boundary distance
- cavity size/location/count、最小特征尺度与多阈值persistence
- full residual、材料合法性、多源本构缺陷
- Best-of-K，K=1/2/4，primary K=4
- Selected-MAE、selection regret
- coverage at frozen validation thresholds
- valid K、distinct K、cost per added valid world
- rejected/failed candidate counts、ROM false-rejection audit
- cold/warm/setup/RHS/wall/memory账本

P2径向误差用4*pi*r^2积分，不用latent MAE。
非激活空腔latent不能产生多个world。
各方法相同K和pool选择规则；不复制候选凑满K；合法基础解必须包含在嵌套pool。
最优真值候选只由离线评分选择。

研究合格阈值按noise/acquisition层的validation强基线中位MAE乘0.90冻结；10%为拟议最低收益，不是工程保证。
若数值波动已超过这个差异，标记无分辨能力，不能从test重新调门槛。
工程阈值未提供时明确NOT_SPECIFIED，不造一个固定数字。

## 6. 主算法：有限edit + nonlinear compensation

先做无NN版本。
基础解使用实际已有、同target且经验证的强常规重建；若不可行先修复并计费。

编辑族按固定次序轮询：
- cavity birth/death、位置/半径变化
- multiple-hole split/merge（仅在表示支持时）
- internal layers
- low-frequency boundary deformation
- localized contrast swap/high-contrast patch
- loss redistribution
- wavelet/multiscale completion

P2必须服从原径向decoder支持，不在该平台假装能移动偏心空腔或生成任意三维拓扑。
二维/三维才引入对应几何编辑。

每个edit：
1. materialize合法proposal；
2. 计算finite secant或有误差状态标签的代理；
3. 以编辑点的material/nuisance Jacobian评估允许补偿能力；
4. 用bounded LM/SQP/trust-region corrector求解完整数据相容；
5. 保护提出的结构区别，例如新空腔不能被corrector直接填掉；
6. 若未合格，记录真实失败/预算停止；
7. 通过full-Maxwell endpoint validation才进入archive；
8. 再做结构diversity增广，不允许diversity抵消数据不合格。

20% proposal预算给非贪心全域/随机edit，比例在validation冻结。
local continuation仅拓展已知分支；不宣称能发现所有disconnected families。
保留候选生成父节点、路径、分支未知状态。
临时结构保护属于搜索子域约束，不要把它说成完整prior。

在小支撑上优先试exact local Green/Schur engine，不先搭大型ROM。
full solver fallback始终保留。
必须承认PNGF/Woodbury已有先例；same engine也可给基线用。
A17的DeltaJ公式只在frozen material的模型变化下用，不能预测跨材料finite response。

## 7. E1：P2四维机制与材料MAE

冻结原zc,zs,epsilon_real,epsilon_imag支持与prior；保留inactive fibers。
gain为原shared real truncated-normal，不增设free phase/per-channel gain。
第一频kR=1.6；第二频2.4只做预定负控制，不为每个candidate自由选择是否加入。
现有fresh64parent*2noise集现为exposed regression。
新确认panel在算法冻结后生成32parent*1noise；可复用原bank进行proposal。

四臂：
A. finite-edit compensation
B. global multi-start material optimizer
C. fixed-budget prior/tempered population
D. frozen C1 MDN + 相同corrector（checkpoint缺失则NOT_RUN，保留缺失）

所有臂使用同一已付费proxy权限；最终Mie调用预算相同。
每parent每臂最多16次新material-frequency attempts：32*4*16=2048。
共用256次用于真值、固定前8parent的第二频率负控制及必要audit。
全阶段2304 cap；CPU1800s；GPU0；禁止新bank/训练。
连续gain的解析拟合不重复调用原始材料forward，但其每次CPU成本仍记录。

如果只用16次调用的SMC未混合，要称budget-limited sampler；不把它当“已收敛Bayes”失败。
候选全部保存full material与full predicted measurement；proxy-only的图不得算正式valid world。

输出：数学PASS/FAIL、MAE/coverage曲线、哪个structure edit成功/失败、selector-oracle gap、prior sensitivity、实际成本。
不报告“普遍r^3到r^5”；已有证据不支持。
第二频更易区分的既有趋势需独立核验，不能拿它当新生成器收益。

## 8. E2：可选、硬预算二维反证

只有 `enable_2d=true` 且E1未出现物理/接口失败时运行。
保留A30原24x24 TM、8源32接收器、complex128、原whitening、材料范围。
每个原四family按固定ID顺序选择第一个parent，weak/ring两条件，共8，标exposed regression。

三主臂：
A. finite-edit + secant-guided compensation
B. ordinary material multi-start
C. randomized material weak/null edits + same corrector

每condition每arm12个fullstate，合计288。
base最多32、拒绝点audit48、数值复查16。
全阶段cap384材料状态x频率，vector RHS9216，wall3600s，GPU0。
任何full J/adjoint/setup均计费。
支持允许多尺度活跃系数扩张，但不偷偷越过同一材料域、改变噪声或gain。
同一baseline可以用同一local engine。

机制消融：同edit无nonlinear correction；同搜索full vs cached engine。
判NO-GO：所有parent均无新增合法结构族；或者只有微小参数变化；或者同预算多起点的MAE-coverage-cost全面支配。
尚未收敛/预算耗尽标明，不能宣称材料唯一。

## 9. E3/E4保持关闭

本Prompt不自动触发三维campaign。
在E2支持机制后，只写出128-state的三维pilot配置、实际成本预估和依赖检查。
三维必须非同心、全矢量、多源共同材料、完整材料场、网格误差审计。
四parent演示不是总体成像优势；最终强论文需独立parent确认实验，预算另行冻结。
不因“论文需要三维”越过本任务边界。

## 10. MaterialHypothesisPackage输出

遵守同目录schema。
每个正式candidate必须有合法材料字段、full validated measurement和完整cost。
默认 `weight_reliability=candidate_only`、`weights=null`、`posterior_weight=null`。
残差softmax不允许填成posterior。
outside-support/unexplored uncertainty默认unquantified/unknown。
pairwise measurement distance必须声明是否profile nuisance和是否为observed-noise target。
候选clustering不等于证明mode独立；连接关系只有验证路径才能标verified。
给B提供完整material law/grid/solver identity和可调用的新action forward契约，不能用旧预测替代未计算的新measurement。

## 11. 最终交付文件

必须产出：
1. `RESEARCH_REPORT_ZH.md`：完整结果、数字、失败、解释，不只是状态摘要。
2. `DECISION.md`：GO/NO-GO/PARTIAL，区分物理存在、算法搜索、MAE和计算效率。
3. `METRICS.json`：逐parent、逐方法、逐K，含regret与失败。
4. `COST_LEDGER.jsonl`：包括拒绝、失败、setup、RHS、wall和历史成本。
5. `MaterialHypothesisPackage/`：manifest、材料、预测、候选对、hash。
6. `UPSTREAM_MANIFEST.json` 和 `TRUTH_ISOLATION_AUDIT.json`。
7. `TEST_RESULTS.json`、完整配置、实际执行命令、环境版本。
8. `FAILURE_CASES.md`：最差case、A30条件是否仍失败、ROM误拒绝。
9. `NEXT_STAGE_CONFIG_ONLY/`：仅条件性三维计划，不自动执行。

报告第一屏回答：
- 是否真找到更多结构不同的合法材料世界？
- 是否真的降低Best-of-4材料MAE？
- 相同信息下Selected-MAE改善多少，oracle差距多大？
- 与相同预算多起点相比是否值得？
- 下一阶段应停止哪项假设，而不是只说继续扩大规模？

请先做短plan并冻结参数，然后执行E0/E1，不停留在只输出计划。
遇到本地资源缺失，完成独立可做部分并准确报告阻断；不要猜文件、假成绩或悄悄更换target。


## 可选的零新增全波诊断：集合风险与选择风险
仅使用已存在的 P2 参考材料样本，可计算报告式(25)的 B(H|y)、S(H|y)，以及有限 bank 上的最佳已找到 K-medoids 风险。标为 reference-assisted diagnostic，不混入普通在线 selector 主表；有限 bank/有限优化值不是连续全局 D_K 的认证界。不得为此新建 response bank、训练模型或增加物理调用。此项不阻塞核心 E0/E1，预算不足则 NOT_RUN。
