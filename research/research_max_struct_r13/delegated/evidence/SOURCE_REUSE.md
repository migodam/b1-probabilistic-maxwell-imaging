# R13 来源证据与复用表

2026-11-12（新加坡）。工作只增加本目录的审计产物；没有重跑 Maxwell、访问新的封存真值、改写历史文件、提交或推送。父线程独立负责核心证明与科学裁决。本表中的“已核对”区分源码/回执一致性与独立物理实验。

## 交付包和可复核程度

Pro 包共 **86个文件**，逐项大小/角色见 `PACKAGE_INVENTORY.csv`。主要研究文字为 `report.md`、`theory_audit.md`、`evidence_audit.md`、`literature_review.md`、`paper_experiment_protocol.md`、`CODEX_TASK.md` 和书目；HTML/PDF/图是这些报告的呈现，不作为独立证据。检查了机制 probe/followup/评分与候选审查源码及已保存 JSON/CSV。NPZ只检查已有机制状态的数组头/维度，不重新读取评分真值。`SHA256SUMS` 等历史记录不重新验证。

包中包含三例 `states.npz`、全部分支 JSON/材料数组与导出评分，**没有**十个历史 AP/DR `CANDIDATES.npz`；`paired_selected_candidates/retrieval_receipt.json` 是历史提取回执，其 `/workspace/scratch/...` 路径不是当前工作树路径。因此五例候选池结论在本次是“原导出表+程序+历史回执”的交叉核对，尚未独立重载源候选数组。不得将这个级别写成完整独立复算。

`CLAIM_LEDGER.csv` 含25项逐主张证据、精确公式/统计、假设和限制。文献证据分 A全文/B摘要；`A-local` 是本地源码/回执证据，**不是**文献证据等级A。

## R4–R12 必须继承的结果与实现

阅读范围：逐轮 `CHAT_HANDOFF.md` 与 `RESEARCH_REPORT_ZH.md` 的结论/实验定义；对R6/R8/R12直接进入影响本轮归因的实际源码；不是重新审计历年所有运行。

|来源|本轮应复用的结构/实际证据|不得重复或扩大主张|
|---|---|---|
|R4|`src/cie_onehop.py` 的四cell残差传播；一个旧弯条正例；新六对象三臂均4/6；32网格普通global亦恢复|one-hop不是普遍更好；不能以更低辅助目标代替材料改善|
|R5|`THEORY_DERIVATIONS.md` 的有限残差完整乘积导数；只改材料块可复制旧例；阶段交接失败；同号cross在正反例均出现|不能把cross正号或冻结metric解释成充分成功条件|
|R6|`src/fft_tsom_ops.py`、`src/benchmark_solver.py::arm_config`；真实非周期FFT动作/伴随；新六对象共同4/6|`mf12`参数是追加阶段；不能误读为从第一步完整Fourier；已有FFT不再重写|
|R7|A+条件Möbius/MM实现、精确低秩有限响应及其失败；30/20两cohort分开|标准Woodbury/更精确有理动作不保证好材料；禁止把后选家族20与先30合并为事前50|
|R8|`src/cie_fft.py`、`NEURAL_WARM_STARTS.md`、`WARM_DEVELOPMENT_ANALYSIS.md`；硬DPCS/Fourier丢失分解与physical-complement anchor；42对象NN→direct强基线|“重新加入物理电流余量”已做；初始缺陷降低不保恢复；NN携带额外训练先验|
|R9|完整current松弛/低罚权旧壳正例及新六例无稳定增量；公开噪声选DPCS秩使各基线共同2/6→4/6；原尺寸3D和精度失败账|解锁current不是首次；60³→30³干净响应差约0.535，不能用matched图代替连续物理资格|
|R10|Green19完整材料分母、源/接收缺陷传播；新环改善与壳/桥反例；已存在低秩分母Woodbury代码|`19`是低秩公开Green投影，Full不等于完整G；不重新扫描rank挽救结论|
|R11|完整分母导数/消融；48对象Global/G19同23成功集合；配对网格环/壳条件增量；真3D r120均值未胜Dual|不能只讲48组无新增而抹掉真实条件正例；也不能用8组覆盖48组限制|
|R12|`supplement_paired/runtime/field_current.py::DenseFieldCurrentProjectors` 与 `run_projection`；共享实材料精确投影、全局noise ball、同投影AP/DR；同24条件救援与正确器损伤|DR是既有算法；11起点/44运行和3机制对象均已暴露；9/24等合取指标不是全部图像能力|
|Formula Atlas|`FORMULA_ATLAS.md` F-BASE-001–004、F-OP-004–007、F-BF-001–002、F-BF-006/009；物理resolvent、共同材料lift、Hadamard Gram、非根曲率与有限secant限制|固定G/S不等于固定材料Jacobian；current暗空间不是合法材料kernel；恒等式不作为R13新定理|

## 与新候选真正不同的因素

|因素|既有算法/协议|R13候选/控制|解释边界|
|---|---|---|---|
|状态|DR保留完整z；AP只投影交替；G19含current、stage、双beta、PR历史|事件改变z的法向坐标|同材料不等于同完整优化状态|
|精确跳步|普通DR每步同gap平移|精确AP对沿一条全局射线直接到首整数tie/crossing|完全相同后续轨迹，只能减少重复投影；不能提高恢复|
|近似proposal|普通DR从原AP辅助状态继续|从`c+t d`开始，整数`floor(t*)+1`或分数`1.001t*`|改变轨迹，无等价/GT保证；本包C收益属于此项|
|表示|legacy `J_D+(I-P_D)F_mf a`，最终完整APCS仍锁DPCS|full APCS从首步；NFFT `J_D+F_full a`；soft/elastic retained coordinates|full空间解锁已是NFFT先例；显式ball会另改目标，须分开|
|读出机会|G19阶段末；AP/DR每50步|共同付费审计日程、阶段内候选|B虽进入mf12，却未审计其终点；是预算下真实pipeline因素|
|corrector|原Phi下全局候选incumbent|保留raw/repair/prior-corrected/selected各阶段|可行时最小位移repair原样保留，只保证不动，不保证GT|
|源—接收诊断|已有G*/S*和投影multiplier|匹配同一个PA调用的lambda/residual，必要时小子集独立adjoint|`p=-d_X/(lambda s)`是重表达，不能当独立验证|

## 三例回执独立加总

`PRO_RECEIPT_COST_RECOUNT.csv/json` 由本目录脚本从九份JSON和三个NPZ头重新加总：A、B各12照明，C为6；每case为23个probe LU + 3×61个分支LU + 两个分支计时外的initial audit = **208 LU**。总计 **624 LU / 6240 complex RHS / 9投影H Cholesky**，CPU已记录小计 **241.448472642秒**，与Pro manifest一致。

这个加总没有运行物理求解器。CPU小计包括probe中的AP3000和两个500步试探，不能叫AP-only；六个后续进程的setup/import/initial audit/serialization等没有完整计时。端到端CPU与峰值内存仍为空，不据此宣称加速。

## 父线程需要核准的重点

1. 活动边界公式仅适用于实共享材料、固定分离度量、非退化有限区间；固定host cells、dark cell、精确ties必须单独处理。近似AP defect小不保事件附近的projection稳定。
2. KKT的lambda是半平方ball约束乘子，必须来自与残差相同的投影调用。旧a与`PA(c)=a′`的lambda混用是潜在实现错误。
3. 规定射线阈值、真实DR大box翻转、达到球内、材料恢复、全成本更快是五个主张。B准确翻转仍失败；C事件结果仍弱于原始NN入口DR raw。
4. 包内候选池数组缺失；本次不得声称重算十池。已有导出表足以提出归因限制，不足以替代今后的必要源数组审查。
