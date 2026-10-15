# 共同先验正确器 C_P：冻结候选合同

状态：**IMPLEMENTATION_NOT_ACTIVATED**。本文件是CODEX_NEXT_ROUND_ZH.md第6–7节的可审合同，不宣称已实现或验证。E0不足以激活E2。E1若未通过两组不共享parent的P0_EXACT竞争pair条件，则本合同保持NOT_RUN；不比较四臂、不追加训练。

## 输入、隔离与材料域

输入仅为在线y、其已知sigma、公开acquisition/background、N3、TRAIN/material-only prior、合法搜索历史、raw proposal及其合法prior seed。隐藏benchmark端点、clean means、geometry参数、truth identity均不得读入。原24×24 scalar TM/VIE、8 sources、32 receivers、complex128、gain1与S/G固定。

材料变量为非负实rho，chi=chi_b+v*rho。BOX、PHASE、P0_EXACT/P0_APPROX/P1_RELAXED/PRIOR_UNKNOWN分开保存。正式P0分支必须返回合法原生成器参数与可重放端点；PHASE不等同P0。连续field surrogate仅供提案，原rasterizer不得替换；retraction失败返回PRIOR_RETRACTION_FAILURE，不调用无限制box GN救援。

## 共同目标与约束

四臂M/T-soft-adapted/G/P共享同一实现、prior、reservoir、validator、thresholds、budget、stopping；每臂至少1/3 raw slots为独立global prior proposals。优先找到同时满足prior、likelihood、补偿和结构保留约束的可行端点，再最小化相对raw的改动，不无界追逐最低残差。global修正半径参照自身prior seed，不能统一施加N3小半径。

固定eta=.5：corrected-seed沿raw-seed结构方向的投影不低于raw的一半；另检查最终多阈值/多尺度离散与持久结构，投影不能代替结构验收。rho/a阈值.1/.25/.5/.75/.9、尺度0/1/2 cells、4前景/8补集、>=4cells组件诊断沿用材料校准。DC、低频、edit外补偿分列；local提案分别设外部及DC/低频预算。其数值必须从TRAIN/material-only自然变化尺度在激活前冻结，所有臂相同；本E0未给出已验证corrector补偿阈值，保留MISSING_PENDING_ACTIVATION，不能临场放宽。

最多两次接受的nonlinear correction updates；任何试探/拒绝/line search/full validation均计入本臂states/RHS。导数/adjoint/factorization/wall-time单列；已缓存可复用但另报standalone cold-start成本。预算不足即停。

## 似然reference与门槛

先执行与y无关的prior gate。ell=.5||W(F(x)-y)||²；512实通道，W每实虚通道sqrt(2)/同一sigma。主Delta ell<=log(10)，log(100)仅预声明敏感性；noise gate同时记录tau=23.785449004057032。

online reference为共同公开seed中prior合格的best-known likelihood；最终以四臂全部合法输出union reference统一重判。off-prior旧GN、hidden truth/witness不得进入reference。只称best-known-relative；没有全局证书不得称全局likelihood optimum。N3-anchor penalty不能作为独立proper prior。posterior_mass=null。

## 保存、验收与E2资源合同

保存raw→phase/prior retracted→corrected→full validated→archived及每步淘汰原因。无严格输出error bound的全模型数值验证只能FULL_VALIDATED_NUMERICAL；边界样本需tolerance refinement或BORDERLINE_NUMERICAL。普通full Jacobian含resolvent，不删feedback造消融。

若E1 gate通过，E2最多4观察×4臂×10 states=160 states；余32仅公共基础/必要数值复核/小型corrector attribution，总E2<=192states/8192complexRHS，整轮<=256states/10240RHS且任何更小有效上限优先。每臂最多3 raw proposals，K<=4含合格基础；不强制凑数。相同小组raw用旧corrector与C_P独立归因；旧off-prior输出不进正式archive。若PC-SP在同一C_P下不优于M/G，不能把corrector收益归于新生成策略。

激活前主线程必须冻结补偿/步长/alpha/正定prior metric/结构保留测量与可执行预算receipt，并审查online/offline隔离。以上是实现完成条件，现状没有已运行C_P结果或四臂结论。

E1构造与局部扰动reservoir均属offline-only。若E2条件满足，online臂另生成独立新seed的256-parent reservoir，不能使用E1隐藏见证的parent/参数/clean means；全轮材料生成仍受2048上限与父PROTOCOL_FROZEN.json约束。
