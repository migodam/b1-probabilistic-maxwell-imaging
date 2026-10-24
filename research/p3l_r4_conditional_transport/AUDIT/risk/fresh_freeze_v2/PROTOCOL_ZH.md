# 新锁定512与独立确认256：只准备，未生成/执行

主风险族固定为 Spline-fusion，依据三seed的**既有R2VAL256截断均值**0.0071267302458485，小于匹配Diag fusion的0.0071301649635037。差距很小；没有因DEV中Diag更好而切换主方法。两组新样本采用同一风险规则、45个既有动作、同三seed、原始相位/物理阈值与raw/矩形截断定义，L1以后不适应、不换方法。

## 唯一供主线程绑定的组件

`RISK_COMPONENT_FROZEN.json`（本目录fresh_freeze_v2）含45方法完整checkpoint身份、三个provider来源与K32/32/0、43份input-feature-state快照、TRAIN缩放/几何frame/源实现身份、验证选择源、所有比较与父样本/种子。模型参数和原闭合文件未改。元数据坏K只使用只读兼容器与已封存TRAIN/validation源gate恢复。未来观测字节尚不存在，不能虚称已知其身份；采集后、进入特征前另行封存。原始观测必须与context观测逐项相同，parent metadata必须对应。

原provisional fresh_freeze_v1从未绑定或生成；在任何新观测前加强了源观测/parent/有限anchor守卫，v2取代它。45模型、provider、验证选择与科学规则完全相同，v1保留。不要绑定v1。

当前两个TOP_LEVEL_*_TEMPLATE.json均是NOT_RELEASED模板，不在configs/FROZEN_*活跃位置。父线程完成VI选择后同时生成真正的`configs/FROZEN_LOCKED.json`与`configs/FROZEN_CONFIRMATION.json`，填入模板中的risk_component_identity/risk_rule/parents/generation_seed，methods_frozen=true、adapt_after_L1=false、rule_frozen_before_L1_generation=true，以及父线程自己的VI绑定，released=true。两个风险规则必须在L1观测生成前承诺。

父样本分别固定R4L1_00000..00511、R4L2_00000..00255；生成seed20261025020/20261025021。原世界rng=SeedSequence([seed,index,0])，原噪声rng=[seed,index,1]，不挑选残差/样本/父样本。物理生成与native96+probit64上下文由主线程队列单独授权，本风险队列不生成、不建上下文、不训练、不调Maxwell。

## 顺序与入口

1. 观测生成前调用`tools/r4_fresh_risk.py bind --cohort locked|confirmation --component 本组件 --top-freeze 主线程FROZEN文件 --out BOUND_RISK.json`。已有该cohort任何观测则拒绝补冻。
2. 父线程独立完成相应固定数量的原观测与在线context（所有66维，无GT网络特征）。
3. `tools/queue_r4_fresh_risk.py plan`只列出命令。`run`要求另一份父线程release（experiment=fresh_risk_panels_v1、cohort、component_identity、top_freeze_identity、allow_features_infer=true、released=true）。只有allow_offline_score=true才启动最后评分进程。队列不等待/补建缺失世界，不扩展样本数。
4. features只用原在线observation生成三provider固定特征；infer按CPU **64父样本批量**计算全部45动作。每批的raw/截断/失败标志以lossless stored NPZ即刻封存，aggregate只逐项拼接这些同一数组。实际分数只由这些保存数组定义；不以重跑近似相等的数组替换。全512与256均由完整64批次组成。
5. export在GT评分前封存raw Full诊断/代表元manifest。随后独立score校验模型/源/数组/ONLINE_COMPLETE/时间/PID后才打开隔离offline_only。按父样本先平均三seed损失，4000配对bootstrap，seed20261025099、CI98.333333%。全部45方法和全部父样本保留；预测失败明确anchor回退。没有可匹配原anchor则整批停在GT门外，不能删父样本后评分。

主要比较是Spline-fusion截断动作对冻结R2、同信息strong spectral、匹配Diag fusion、空间置换fusion；全部raw、截断、其他方法/相位修订也完整披露。全576复MSE/MAE、ROI与相位分层沿用既有独立评分器。

## 历史replay失败与新分数定义

原2e-6 batch1/64 replay检查为FAIL，原FAILED账本和部分输出保留。v2的1e-5检查只是一项float32技术启发式核对，**不把旧门槛改成PASS，不形成科学资格，也不是科学阈值放宽**。验证选择使用原来保存的64批量预测/分数。旧DEV闭合实际数组保持原样；未来两cohort只认本协议新保存的canonical batch64数组，不靠任何replay容差决定风险分数。

## 固定代表元与raw Full诊断（均未启动）

代表元只取L1前64父样本：Spline seed2026102501与冻结R2 seed33，各raw与clipped，共4目标；不平均seed图像。每父样本一个共享ordinary Full-Phi控制，独立于所有mu的提议与选择，和每代表元均使用原协议每臂16新states/256复RHS。共5臂，计划上限64×5×256=81920复RHS。ordinary成本记一次、不能免费；既有anchor/receiver获取成本单列原owner。原Full阈值282.873792、source_adjoint8、全66原decoder/无新prior均不变。

`tools/r4_fresh_representatives.py prepare`在ONLINE_COMPLETE后冻结实际mu及资格/skips；`run`需要独立匹配prepared_request_identity的物理release（experiment=fresh_locked64_representatives_v1），绝不默认运行。`score`是实际BookII COMPLETE之后的独立进程，校验时间/数组/arm结果后才读GT。共享ordinary只有一个实际字段文件，各目标比较引用它；全64分母、未资格/缺失/回退明确披露为NOT_EVALUATED或fallback，不丢弃。逐父NPZ/CSV/JSON报告Phi、prior、RMS、风险变化与父样本bootstrap，不冠以后验称谓。

raw材料Full前向诊断对Spline_s1、R2、strong_s1的raw+clipped六动作，覆盖全部512+256父样本，独立于anchor可行性筛选。在线export的RAW_FORWARD_FROZEN.json可交既有`tools/r4_raw_material_forward.py run`，L1计划24576复RHS、确认12288复RHS，共36864；另需主线程物理release/预算/排程。只报告非线性响应、被动性、ROI和矩形边界；任意NN图像的decoder可实现性为NOT_PROVEN，绝不当作P1合法世界或后验样本。

## 本次执行范围

只生成准备元数据和既有模型feature-state快照；未启动新观测/上下文/Maxwell/训练、未读LOCKED或确认GT。合成合同检查仅验证stored批64逐项保存、失败分母和未release模板拒绝；不构成科学证据。最终科学Go与物理队列均由主线程掌握。
