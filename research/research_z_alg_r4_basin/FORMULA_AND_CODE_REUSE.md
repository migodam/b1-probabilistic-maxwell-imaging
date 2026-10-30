# R4 公式与代码复用审计

状态：`READ_ONLY_SOURCE_REVIEW / NO_NUMERICAL_RUN / NO_GT_ACCESS`。本文件不证明实现、计算收益、成像成功或原创性。只拥有本表与 `theory/BRANCH_NAVIGATION.md`；未修改历史项目。

## 定向来源

- R2：本 repository 内 `research/research_z_alg_r2/FORMULA_COMPENDIUM.md` 的 F06、F09、F11、F12。
- R3：`research/research_z_alg_r3/MATHEMATICAL_THEORY.md`、`FORMULA_REUSE_MAP.md`，以及 `src/physics.py`、`src/model.py`、`src/solver.py`、`src/portfolio.py`。
- R9：本地 `Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/` 的 `PHYSICAL_VS_MATERIAL_SUBSPACE.md`、`src/stable_roi_operator.py`，以及公式集 03、06 章。随后按主线程指定路径补读 `Gaussian/B1/b1-p3l-r11-stable-digits/research/p3l_r9_stabilized_subspace/formula_collection/06_DERIVED_OPTIMIZATIONS.md` 的前置共同后端、端口成本和 varpro 段；不假定两份 checkout 的字节等同性。
- R10：本地 `Gaussian/B1/b1-p3l-r10-bilateral-structure/research/p3l_r10_bilateral_structure/` 的 `BILATERAL_OPERATOR_THEORY.md`、`src/transfer_campaign.py`、`src/reduced_response.py`。
- FP-T：初始定向文件名搜索没有命中，随后主线程提供实际路径，已读取 `Gaussian/B1/b1-p3l-r12-bilateral-fixedpoint/research/p3l_r12_bilateral_fixedpoint/TRANSVERSE_CONTINUATION.md` 与 `FIXED_POINT_SPECTRUM_AUDIT.md`，并定位 `agent_b/engine.py` 的接口。最终状态为 `ACCESSIBLE_READ`，不能继续记不可访问。

## 公式与代码逐项裁决

| 来源 | 可合法复用的内容 | R4 实际意义与限制 | 状态 |
|---|---|---|---|
| R2 F06 | 固定参考后的 \(Y=Y_0+L(I-\Delta K)^{-1}\Delta q\) | U=I、\(\Delta=\alpha D_a\) 得全 576 材料幅度族；全部活跃时仍是 576 阶，不是低秩 | `ADOPT_IN_THEORY` |
| R2 F09 | \(Wq=T, W=(L+WK)\Delta\)，测量先消元和共同实材料条件 | 可解释 R3 GlobalDual，但其辅助维度不为 576 材料自动降维；不是幅度射线公式的必要组件 | `REFERENCE_ONLY` |
| R2 F11 | 补集固定时精确 Schur 复用 | 若真的做条件小块可用；本轮全部材料幅度 ray 不虚构已知补集或 576→256 ROI 收益 | `NOT_REQUIRED_FOR_MAIN_ARM` |
| R2 F12 | 单坐标 rank-one 响应 Möbius、损失驻点二次式 | 可作强坐标对照；不能把二次式套到一般全秩 h。后者驻点多项式上界为 \(4N-2\) | `VALID_STRONG_CONTROL` |
| R3 `physics.Model` | 同一离散 G/S/Einc，原状态方程、实材料和 packing | 公共物理模型应共同使用；避免另建模型、换源接收顺序或偷偷换噪声 | `REUSE_BACKEND_IF_COMPATIBLE` |
| R3 `PublicModel.forward/jacobian` | 同一 LU 上 s 个源 RHS、min(N,r) 个 Jacobian RHS，逐状态计费 | N576、r通常小于N时 Full-J 为 s+r RHS，不是 N×s；耗时、完整分解和收缩也计 | `REUSE_ACCOUNTING_AND_DERIVATIVE` |
| R3 `ReducedKernel` | 独立材料选择列、参考 LU、固定核与切片限制 | 可作为全参考 setup 的参照；不把稠密 ray h 作为单个坐标 selector 绕过 rank 检查 | `ADAPT_FORMULA_NOT_FAKE_RANK` |
| R3 `Run.coordinate` | 二次驻点、端点及 incumbent 的有理坐标评分 | 已实现普通 exact coordinate control；当前只 correct 评分最好的 archive 点，不是多低谷逐个 correct 的实现 | `EXISTING_CONTROL` |
| R3 `Run.globaldual` | 接收对偶硬数据约束、共同材料 variable projection、端点 full forward 与 full-J polish | 每次内层重新 QR；其目标是材料相关度量，零点等价不表示局部极小等价，近极点可虚降 | `CSI_VARPRO_EXISTING_CONTROL` |
| R3 Guarded portfolio | 保留公共普通求解器候选，最终按原全部观测残差排序 | 只能保证不劣于已保留的候选；不能保证好于相同总预算普通法或材料 GT 更准 | `REUSE_INCUMBENT_RULE` |
| R3 theory §2 | 较高原始损失的新根先全变量正确再比较 | 本轮多低谷操作已有明确项目内思想前身；不能将操作本身包装成首次提出 | `DIRECT_CONCEPTUAL_PRECEDENT` |
| R9 PHYSICAL / stable ROI | 真实材料兼容、宿主精确消元、s+r 端口 Jacobian | R9 原先固定宿主和 ROI=256；本轮全576未知不能复制其降维收益 | `BOUNDARY_AND_BACKEND_REFERENCE` |
| R9 F-OPT-006/007/008 | 完整稳定 QR 基、仿射弱坐标约束、强坐标 profile | 全基换坐标仍是同一个 GN；冻结弱坐标改变可行域；profile 不是积分 | `STANDARD_CORRECTOR_TOOLS` |
| R10 B2/B3/B7 | 端点相依精确双边差、真实材料 Jacobian 与双边 Gram | 可解释强/弱方向，但固定 Born Gram 不能代替原非线性资格 | `OPTIONAL_DIRECTION_DIAGNOSTIC` |
| R10 B8/B9 | 固定形状的幅度有理族、Schur 复用、一二阶幅度导数 | R4 实幅度 b+αa 是直接同结构应用；任意非 GT shape 仍须重做其分解 | `DIRECT_FORMULA_PRECEDENT` |
| R10 `transfer_campaign.Shape` | 构造 B/L/R、复 Schur、缓存 LQ/QR、三角求解幅度响应及导数 | 可参考其核心数值结构；必须去掉 Gaussian66 decoder、固定256、复幅度 prior、历史σ/数据文件等专属假设 | `PORT_NUMERICAL_PATTERN_ONLY` |
| R10 `reduced_response.py` | 短基/低阶响应与原方程缺陷、输出误差、完整回退 | 若主线只需精确全阶扫描可先不引入；不能把过去64阶/256格资格迁移到576实材料 | `OPTIONAL_APPROXIMATION_NOT_QUALIFIED` |
| FP-T audit | 零残差根处的阻尼导数、非零残差 Hessian 项、算法单位模态 | 只能支持正则局部横向行为；不推出跨 basin、GT、全局纤维或 posterior | `ADOPT_LIMITING_THEOREM` |
| FP-T `transverse_correct` | 当前精确 J、固定弱绝对坐标、边界与回溯、付费失败 | 若需要保留远分支可做备选；本轮先用全变量 correct 避免新增约束可行性问题 | `OPTIONAL_STANDARD_CORRECTOR` |

## 与普通线搜索、CSI 以及 R10 的差别

| 方法 | 搜索对象 / 在线目标 | 与 R4 多低谷版本的关系 |
|---|---|---|
| 回溯线搜索 | 沿下降方向不断缩步，只接受局部下降 | 不会主动越过原目标上坡区；R4 扫整个合法区间且保留远低谷 |
| 普通 global line minimization | 在同一 ray 只选择原目标全局 argmin | 使用同一个精确有理 evaluator 时，数学查询完全相同；R4 的增量是多个低谷分别 correct 后才排序 |
| R10 固定形状幅度 profile | 相同类型的精确幅度 resolvent；历史主要为固定形状重复响应/条件任务 | R4 公式直接继承；候选保留、不同低谷校正和全576实材料预算需要新实现/证据 |
| CSI / R3 GlobalDual | 辅助电流/接收传递与材料共同兼容；varpro 改变非零残差度量 | R4 始终查询原有理响应的原数据损失，沿有限材料 ray 组织候选；它没有提出新的 CSI 目标或新的状态方程 |
| basin hopping / 多起点局部校正 | 产生远种子、逐个局部优化、按原目标选 | R4 属于这个标准机制家族；它的可检验特征是用固定形状的谱复用产生多个低谷种子，而不是随机扰动种子 |

## 实现交接：最小函数合同

建议独立函数 `propose_amplitude_valleys(public, m_current, b, bounds, budget, config)`：

1. 只读取公共 G/S/Einc/Yobs、材料向量、bounds 和原观测 metric；形成 \(a=m_{current}-b\)。拒绝零方向、空区间、非有限输入并计状态。
2. 共享背景 setup；每个 a 形成全阶 \(T_a,B_a\)，默认 complex Schur，缓存 LQ/QR。setup 全部进同一 ledger。
3. 在固定合法区间内作冻结有限节点查询；细化检出的低谷。返回 incumbent、端点和所有未被候选上限截断的分离低谷，不按原始损失提前筛掉较高低谷。
4. 对每个低谷在同等冻结子预算内调用共同完整正确器。原始高损失只是诊断；box、数值失败与实际总预算仍是停止条件。
5. 逐个保留原全部观测评分；全部正确后再排序、材料去重。保存 `PRE_CORRECTION`、`POST_CORRECTION`、失败和截断记录。

建议 ledger 字段至少包含 `full_LU_576`、`full_RHS_576`、`Schur_576`、`triangular_RHS_576`、`matrix_builds`、`ray_evaluations`、`derivative_evaluations`、`corrector_trials`、`wall_seconds`、`failed_paid`。RHS 与分解分别报告，不能把 LU、Schur、QR 折成一个含糊的 RHS 数字。

有限节点扫描的状态名应为 `DETECTED_VALLEYS`；只有可靠根隔离覆盖才可声明 `ALL_STATIONARY_POINTS`。最终材料候选可标 `PHYSICS_VALIDATED`，仍不能自动标 `GT_RECOVERED`、`DIFFERENT_BASIN`、`POSTERIOR_SAMPLE` 或 `NOVEL_METHOD`。

最有判别力的结果是：在同预算下，多低谷版本产生一个原始损失较高、correct 后优于只取线 argmin 的合法端点；若实际 Maxwell 未出现此现象，解析反例仅保留为算法区别解释，收益结论 `HOLD`。

## 最终执行补记（主代理）

V3实际复用了冻结R3-DIGITS的Model、packing和Regularizer，新的forward/Jacobian没有第二实现。Full-J、data-adjoint与Gram调用都共享其成本计数。V3所有L-BFGS用(m−0.5)/19.5，保持物理可行集和原正则；其尺度处理是标准仿射链式法则。主代理另实际定向读取A18 Route A当前`Gaussian/Research_W/publish-w-p0/research/research_w_gsda_r4/src/gsda/fwi.py`（stage_scale、coord_g和L-BFGS段），以及Route B R11的`shared/invert_digits.py`（trust cap与alpha回溯段）。只移植尺度模式，不移植已知ROI、复损耗、多频数据、GT支撑或旧成功率。A17X的receiver/source端口与材料合法性通过R9和R3物理实现、已冻结跨session审计作来源参照，没有声称完整移植其反馈算法。

V1/V2 source snapshot、最终V3 snapshot和全部开发失败都保留。最终独立fresh没有GT提升，因此所有ADOPT_IN_THEORY行只表示数学/实现采用，不表示论文收益已经成立。FP-T/R12本地资产已访问，不能记作不可访问。
