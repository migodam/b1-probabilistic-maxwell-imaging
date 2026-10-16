# L/Z/X 集成审计：已实现的核、不能混合的目标与可复用资产

2026-10-16。附件只读整合；本分工无新增物理、无网络、无训练、无脚本重跑。以下为可执行接口核对与冲突提示，核心数学有效性、实验方向及最终裁决留主线程。

## 1. 三包实际交付状态

| 包 | 正式报告 | 实际可运行内容 | 仍只有spec/假设的部分 |
|---|---|---|---|
| L | **存在完整正式MD与PDF**，另proof/experiment/literature附录 | `code/algebra_checks.py`：合成Schur缺陷、inverse precision、directional residual、finite-feedback二/三阶余项、prior横截/边界/不连通/体积反例 | L-A tangent atlas、L-B LIS-MCMC、L-C lifted HMC；没有真实N3 chart、合法anchor、路径或采样数组 |
| Z | 完整中文报告与理论/算法/协议 | 只有JSON Schema和文档；没有sampler、flow、decoder或训练脚本 | topology-open PZ、全尺度field latent、conditional flow、逐状态metric-MALA、SMC、DA、A17低秩metric应用 |
| X | 完整中文报告、四命题、28项文献表 | `code/sampler_core.py`通用Gaussian-prior冻结MH/pCN核；`run_synthetic_checks.py`合成/有限状态检查 | P2 likelihood/density adapter、标准tempered sampler、非Gaussian/混合P0、RJMCMC、Maxwell current refinement、物理误差证书、DA代码 |

X的JSON/stdout已保存合成结果；这不等于新的N3/Maxwell实验。本轮只读，未重跑原L代数或X检查。

### X通用核逐项

- `GaussianProposal`：给定mean/Cholesky的Gaussian sample/logpdf。
- `DefensiveMixture`：冻结可评价混合，包含标准Gaussian prior入口；独立MH用完整mixture密度。
- `FrozenPCN.from_jacobian`：冻结一个real prior-coordinate A，构造尺度有下限的pCN；其矩阵不以anchor为新prior均值。
- `FrozenPCN.validate`：检验R对称、RRT+BBT=I、B满维。
- `CachedOracle`：float64参数字节键缓存同一个negative-log-likelihood；计unique oracle evaluation；非有限返回直接报错，不当作零posterior悄悄拒绝。
- `run_frozen_sampler`：固定概率global/local调度；local锚点selector带反向/正向权重；保留拒绝重复；可选严格区间只跳过确定拒绝的full oracle。
- **边界缺口**：当前入口要求 `0<global_probability<1`，不能直接以0/1跑local-only/global-only；正式五臂driver需明确修改或另实现纯臂，而不能用ε调度冒称纯臂。
- **计费缺口**：该核计oracle cache miss，不自动给所有proposal尝试、reference密度、selector、Jacobian设置、温度链、failed evaluations和wall计费。集成需要外层共同账本。
- **数值合同缺口**：代码只实现标准Gaussian prior；不能直接把uniform x变量传入pCN并沿用简化接受率。`certified_phi_interval`若没有真包含界必须关闭。

## 2. 必须修正的包版本冲突

X的报告§1.1、README、NEXT和INPUT_MANIFEST把其当时收到的L记为四文件algebra/interface supplement，并写正式报告缺失。这是**X旧输入版本事实**。

本轮L是新完整包：`RESEARCH_L_FULL_REPORT_ZH.md` 65,730 bytes，并含正式PDF、证明/实验/文献附录。故本轮应记 `L_FORMAL_REPORT_PRESENT; NO_PHYSICAL_LOCAL_INSTANCE`。不要沿用 `L_REPORT_NOT_SUPPLIED`，也不要改写旧X文本使它似乎当时已经读过新L。

另有两项证据深度提醒：X `SOURCE_INDEX.md` 对C1 P2记 `first_150_source_lines_read`，`SOURCES.json`记关键定义/reference/results段；当前文献与项目事实须由主线程对本地正式文件核实。R10/DDMC、R28等仅摘要，不能据此宣称正文定理已被本轮验证。

## 3. local、proposal与prior不能互换

| 对象 | 正确合同 | 集成时避免的误用 |
|---|---|---|
| N3 | 确定性Top-1估计/合法online context | 不自动是prior draw、P0_EXACT anchor或posterior sampler |
| 原P0 | 四family原生成测度、fixed increment phase、真实背景/幅度、硬raster与归一化 | 1152>512不能推出P0内高维fiber；phase/box通过不能升级P0_EXACT |
| Z的PZ_GEOM/PZ_MATERIAL | 新先验/新target | 全尺度随机场不是P0无害换坐标；新支持不能覆盖P0历史负结果 |
| L iso-prediction fiber | 固定f(anchor)的等值集合 | 不等于含噪y的完整posterior；不合格anchor的等预测路径也不自动合格 |
| L noise tube / likelihood region | 有声明阈值的集合 | 都进noise ball不等于点likelihood相近，更不等于区域质量相同 |
| L local conditional posterior | 明确区域截断后的条件分布 | 不能免费给该区域在全局posterior中的质量 |
| Z/X global proposal | 冻结有密度q或明确prior-reversible核 | q频率、archive条数、softmax residual不是posterior mass |
| deterministic material corrector | 默认只用于discovery | 多对一GN/retraction后不能沿用raw q；DA也不能补未知推前密度 |
| X fixed-material current refinement | 所提chi不变，只提高同一likelihood求值精度 | 不是新的材料候选，也没有自动A17加速证据 |

原P0硬rasterizer局部常值的参数方向、decoder gauge与inactive coordinates，要按物理材料图和原测度分别处理。`ker Dg`不是材料新方向；去除冗余参数时也不能丢掉它携带的prior体积。

Joint P2中latent coefficient-uniform与cavity-radius-uniform不是同一prior。已有12个事件bins可沿同一连续ridge分布，不能预设为12个mode。`P2_PROXY_C`与`P2_FULL_MIE_C`需不同target ID；少量full端点误差不足以认证全球event mass。

## 4. 物理和A17边界

三包统一保留真实多源共享材料、固定白化/source order、passivity、背景和原gain约束。二维24×24 scalar TM证据不是一般三维vector-Maxwell；径向Mie也不是非同心开放3D支持。

A17指定两Galerkin空间差为 `J_C-J_A=Y Γ^-1 N`。对full误差仍需 `J_F-J_C`，不能删余项。普通full Jacobian已经含Maxwell resolvent feedback，不能删除它造“无A17”弱基线。

L协方差/weak-direction/volume推导需要一致prior metric、real packing、可逆性、同秩和必要谱隙；LM阻尼本身不是Bayesian prior。常秩结论要求合法光滑chart，不能直接跨hard raster、active boundary或rank transition。

X提出的严格区间安全拒绝需要真实包含同一target的likelihood区间；toy脚本区间直接用toy真值构造，**只证明接受逻辑不变，不证明计算省算**。无可靠界时可考虑合法DA，但X当前未实现DA，且本阶段不能因缺界追加物理。

本分工不裁决条件命题是否足以构成新理论；通用MH/SMC/pCN/coarea/Schur都不能仅靠改成Maxwell符号主张原创。

## 5. Schema并不直接兼容

### L原生

`L/Research_L/interface/LocalFiberPackage.schema.json`：`schema_version=research-l-0.1`，以physics_contract/prior_contract/anchor/local_region/geometry/candidates/probability_contract/cost/status组织局部sidecar。连续候选的posterior_point_mass强制null；global_region_mass未知可null。schema允许许多字段缺失/附加，故仍需真实语义核对：target身份、文件可读、full prediction与likelihood一致、路径与端点等级。

### Z候选

`Z/.../interfaces/PRO_L_CANDIDATE_INTERFACE.schema.json`：`schema_version=ZL-0.1`，逐candidate，顶层及子对象多处additionalProperties=false；要求material/generation/validation/cost，不同于L atlas聚合结构。L原生JSON不能只改文件名塞入Z。

`material.chi_checksum`是Z旧schema必需字符串；用户不要求SHA-256，本轮不能据此创建新的hash gate。需要不覆写原schema的版本化overlay或映射策略，以源commit/数组路径/实际材料语义核对替代“认证值”；不要伪造checksum或用空字符串冒称已核验。

L可给candidate-only点/路径；Z若将其作为统计proposal，仍需实际proposal law/正反核。没有真实instance时，只交schema与接口状态，不生产合成N3 geometry作为物理交付。

## 6. 可复用本地路径（已确认存在；本分工只枚举，不声称新数值验证）

工作树：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-research-lxz`，本次读取时HEAD `0b9024b11b1bbb338ae334a873e731fb46e4bb66`。

| 资产 | 相对上述工作树路径/外部绝对根 | 可复用范围 |
|---|---|---|
| 原P0复原/phase/provenance | `src/prior_constrained/prior.py`、`pairs.py`、`cached_audit.py` | 只读复用已封存prior审计，不重跑pair campaign |
| 原二维真实backend接口 | `src/research_t/target.py`、`handoff.py`、`corrector.py` | 合同/缓存审计；本阶段不能调用新physics |
| P2 target/decoder/likelihood | `src/b1_joint/{material,response,likelihood,reference,interface}.py` | 主线程核实原C proxy后适配已有oracle |
| P2银行/观测 | `research/joint_p2/results/COMMON_RESPONSE_BANK.npz`、`VALIDATION_OBSERVATIONS.npz`、`TRAIN_OBSERVATIONS.npz` | 已付费资产；fresh/locked读取受角色约束，不因为文件可见就打开 |
| P2参考与样本 | `research/joint_p2/results/*_REFERENCE_PARTICLES.npz`、`packages/*_samples.npz` | 数值target和accepted版本须由主线程核实，不能名字带reference就直接当gold |
| P2已有修正proposal | `research/joint_p2/results/CORRECTED_PROPOSAL_{diag,full,gaussian}_distill.npz` | 旧输出证据；不是自动可评价q的采样接口 |
| 实际C1 checkout | `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/C1/c1-joint-p2` | 本地存在；常见 `c1-physics-state-generative-imaging/p2` 路径在本次探查不存在，不应猜路径 |
| C1四维P2 density | 上一根下 `p2/c1p2/density.py`、`p2/C1_P2_RESEARCH_REPORT_ZH.md` | 需要P2接口，不要误用 `src/c1/density.py`或P1二维权重 |

## 7. 推荐父线程集成顺序（不另起授权）

1. 校正本轮L正式报告已存在的版本事实，保留“无真实local instance”。
2. 先完成target/reference/anchor本地证据审计；P0负结果只读引用，不重跑原pair/compensation。
3. 若已有P2 C proxy、accepted continuous reference、实际normalized proposal均可用，在同一个measure下做uniform↔Gaussian transport；Jacobian要进q。端点舍入需声明数值合同，不静默clip。
4. 使用X冻结通用核时补纯臂配置、外层完整proxy计数与wall、拒绝重复、四seed/checkpoint；所有输出标PROXY_TARGET_MECHANISM_ONLY。
5. 尚无可信多个分离区域时，比较event mass/材料函数而不虚构cross-mode成功；reference不够则明确停止该比较。
6. 当前零新physics约束内，A17只做已有矩阵/历史证据分析；没有真实误差界和同成本增量，不做认证筛除或新full求值。

最终实验状态与继续/停止由主线程基于实际结果决定。本审计没有把任何附件建议预算或网络训练量变成新授权。
