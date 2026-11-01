# N1 独立实现审阅

审阅者：独立 N2 worker；2026-11-01。范围为 `joint_current_model.py`、`profile_schur.py`、`controlled_release.py`，并只读相关 `current_normal.py`、`current_mode_selection.py`、父线程物理 wrapper 与 E1 配置接线。没有修改这些实现，没有追加数值测试、读取 GT 或作理论优先权判断。父线程单次 E0 运行另承担导数、profile 和改基数值证据。

**初始审阅状态：HOLD，需关闭下面的 P1 same-span 数值秩问题；其余核心代数静态检查通过。** 本文反映审阅时源码，后续修复需单列回执，不能把静态检查当运行验收。

## 确认正确的接线

1. **复数据对实变量的实化一致。** Joint model 对数据、状态与各 Jacobian 作用都使用同一个 `model.pack`；父 wrapper 的 pack 为 source→receiver/cell→real/imag。`current_normal.mode_matrix/real_unpack` 使用另一套自洽的 concatenate(real,imag) 仅做内部几何投影，没有把该顺序混入 Joint data residual。空间复 atom 经 `spatial_modes` 明确展开为 2p 个实系数，没有将复自由度记成一个实变量。

2. **非物理电流下的材料导数保留 DP·e。** `linearize` 的逐行材料注入为 `-w Q + dw e`，其中 CSI 为 w=1,dw=0；CIE 为 w=beta/(1+beta m),dw=-beta²/(1+beta m)²；normalized CIE 为 w=1/(1+beta m),dw=-beta/(1+beta m)²。current 导数使用同一 w 乘 L V。除 fixed state_weight 必须为声明的固定量外，未发现把 B(m) 或 CIE 材料相关权重当常量的问题。

3. **Schur/profile 代数与两个通道一致。** H=Bᵀ(I-P_A)B+Λ+LM，g=Bᵀ(I-P_A)b+ell。新增 residual-QR range 给 D=QᵀB、v=Qᵀb；软释放 alpha=1/(1+gamma) 后同时 downdate H 与 g，不能只改 Hessian。恢复 `u=-alpha(v+Dp)`、`c_C=lstsq(R,u)` 后重新求旧系数，符合剖面消元。

4. **软释放惩罚在实际残差 QR 坐标。** `thin_basis` 返回 R=QᵀA 的原系数顺序，因此 penalty 是 gamma||R c_C||²/2。它没有把非零 atom 缩放当逐渐释放，也没有在任意原始系数中裸加 gamma I。同容量 replacement helper 以共同 retained core 分别计算两个增量，并检查实际实系数数量。

5. **材料任务分数保留残差与曲率项。** `physical_task_gain` 使用 -(r+Jp)ᵀJdp-||Jdp||²/2 和共同 prior 的线性/二次项，没有把 LM 计成真实物理信息。Born 预评分明确标 `Born_no_error_bound`；完整 finalist 通过付费 `paid_physical_jvp` 重评。该量仍是局部模型预测，不是有限物理收益证书。

6. **有限 seed 先固定实际物理基，再 hold/release。** 由 u=±tau w 恢复候选与 core 系数，当前材料用实际 box QP 重算，不能把无约束位移的分数归给 clip 后材料。hold 阶段只开放 core，因此已插入的 candidate 物理动作保持固定；free 阶段再释放该固定候选基。路径中辅助状态的成功没有直接成为最终物理解，checkpoint 必须调用完整物理正演，原 incumbent 被保留。

7. **最新 bounded-scale 修复显式改变尺度合同。** 审阅过程中 `finite_normal_seed` 从固定半径版本更新为 bounded feasible scale。当前程序先记录请求尺度，再联合测试±两支，必要时只做最多八次有限 bisection、每次支付 QP 与 Green 作用；最后保留 actual_material_rms、tau_requested、scale_status、normal activation、材料最小有效位移。它没有证明区间可行性单调，也没有声称找到全局最大可行 tau。若正文仍说“严格达到固定请求 RMS”，必须随源码更新为“请求值为上限、实际值单列”。

## 必须处理的 P1：same-span 候选可被舍入噪声伪增广

位置：`current_normal.py::thin_basis`，初读行 46–55；调用方 `profile_schur.py::release_candidate`，初读行 71–75；`quotient_normal_signature` 存在同样投影后判秩入口。

`thin_basis(A)` 的判秩阈值只有 `rtol * max(abs(diag(R)))`。当 candidate 本来位于旧 A 的 range 内时，Z=A_C-Q0(Q0ᵀA_C) 理论上为零，但浮点 Z 通常只有约 machine-epsilon×||A_C|| 的残差。对 Z 再以它自身的 tiny 最大对角作相对阈值，会把这段 roundoff 识别为非零秩。随后 Q 将噪声归一化为单位方向，D、v 和材料步就可能获得没有物理新增自由度的假变化。

**影响：** 这是 same-span/no-op 与 dressing 负控制的必要数值边界，不能只用普通 full-rank 新候选的改基测试替代。一般 C→CT+UK 且 C 确实带新方向时常能通过改基测试，而纯 C∈range(U) 的 no-op 仍可能失败。

**建议关闭方式：** 判定投影残差秩时保留原 A_C 的尺度基准，或先用相对于原候选／旧基的绝对 residual-norm 门排除纯数值零空间；同时保留正常小尺度候选的缩放合同。该项需由实现 owner 修复，并纳入父线程已经安排的一次 E0 same-span/no-op 证据。独立审阅者未追加测试或代替 owner 修改实现。

## 其他需要保留的范围与实现限制

- **P2：QP 停止状态没有返回。** `box_quadratic_step` 在 L-BFGS-B 分支直接返回 result.x，没有保存 `success`、message、投影梯度或达到迭代上限的事件。该步始终符合 box，后续真实非线性下降检查提供保护，但不能把每个数值 QP 都称为精确最优。因此 profile 差界是理想模型公式；实际被截断 QP 应记录数值误差／状态。`unconstrained_step_bound` 已正确命名，不应在 box 活动时宣称它一定界住实际约束步差。

- **P2：父 wrapper 的 prior 曲率是 IRLS surrogate。** `physics.py::regularizer` 初读行57–60返回 Dᵀdiag(2smooth+tv/q)D。这是共同正曲率模型；不是 exact smooth-TV Hessian（TV 系数应为 epsilon²/q³），也不是 `regularizer_residual` 的精确 GN Hessian。它可合法充当本文局部二次模型 Λ，但正文必须称共同 surrogate Λ。真实 objective 与 gradient 使用的是正确原 prior，finite acceptance 也重算原 prior，因此这项不必否定算法，只限制“精确 Hessian”表述。

- **配置显式字段名称需对齐。** 初读 E1 runner 给 `feasible_seed_search/seed_backtracking_trials/minimum_seed_material_rms`，最新实现读 `bounded_seed_scale/seed_scale_trials/seed_min_material_rms`。当前默认分别同为 True/8/.02，所以当前默认运行不受数值影响；但改 runner 字段会被静默忽略，冻结配置需要改成实际被读取的名字。

- **白化合同。** callable data_whitener 被同样用于 residual 与导数，正确性要求它是固定实线性算子；当前 scalar data_scale 满足。不能把任意非线性 callable 解释成已正确求导。

- **可选 callback 的预算范围。** ReleaseCost 在 callback 前后检查总预算，但可选内部 corrector callback 未逐次受 ReleaseCost 预付 gate 约束，可能先超支再发现。当前 E1 使用外部共同 corrector 的独立剩余 wall 分配，因此这不是当前接线路径的阻断项；不要把 API 可选 callback 的任意配置称为严格原子级预算保证。

- **A17 与 prior-art。** 当前 common pool 包含 Fourier、receiver、Green-Krylov、公共 random 模式；没有仅凭命名就建立独立 A17 理论增量。方法可能验证共享材料 task selection，但不能用本次静态实现审阅宣称理论首创或已经超越成熟 CIE。

## 等待的验收证据

本审阅不增加公式误差地板测试。需要父线程已经安排的单次 E0 将 derivative、DP e、直接联合消元对 profile、QR 软惩罚、actual box step、pure-span no-op／dressing 各项逐一标注 PASS／FAIL／NOT_RUN。P1 未关闭前，不应把“改基无增量”负控制标为完成；最终科学判断与 E1 是否支持新方法由父线程保留。

## 结案补记：实现修复与现有回执，未重新运行测试

2026-11-01，父线程要求仅依据现有代码与回执关闭审阅。本审阅者未重跑 N1、N2、公式或成像实验。

**原 P1 实现问题已修复。** 当前 `thin_basis` 接受 `reference_scale`，实际阈值使用 residual 自身 QR 对角尺度与原候选尺度中的较大者。`release_candidate` 传入 `norm(A_C)`，`quotient_normal_signature` 传入 `norm(AC)`。因此纯 same-span 舍入残差不会再仅因相对于自身非零而归一化为“新增方向”。父线程明确转交 owner 的实际 pure-span 回执为 residual rank=0、step=0；本次文件检索尚未找到其指定的独立 `n1_implementation_checks.json` 路径，故该单项运行数值先标为 **PARENT_REPORTED，独立 receipt 文件待定位**，不冒称已经读到文件。

已独立读取 `audit/E0_PRODUCTION_WIRING.json`，状态为 PASS。其中新生产接线的 joint-CIE 方向导数相对误差 2.2798e-10；物理切向法向相对量 3.6998e-16；hard profile 的 H/g/step 差约 1.70e-17／2.64e-19／8.41e-18；dressing 材料步差 2.30e-17；soft-coordinate 步差 2.72e-17。该现有回执包含 9 次完整 LU、114 个实际复 RHS，相关成本没有因本次独立审阅重复发生。它支持本文已检查的生产接线，并不构成成像收益或理论优先权认证。

其他关闭与保留项如下：

- E1 runner 现已使用 `bounded_seed_scale/seed_scale_trials/seed_min_material_rms`，与实现读取字段一致；原来仅靠默认值相同的接线缺口已关闭。
- prior Λ 按父线程最终解释明确为共同 IRLS 正曲率 surrogate；不称 exact smooth-TV Hessian。该解释关闭表述问题，不改变其有限物理收益预测仍需验证的范围。
- 新 observable backbone 规则只依赖公共 S 的奇异谱与固定相对 cutoff=1e-3，保留 4×4 的 12 个、6×6 的 16 个接收模式，和额外公共 Fourier atom。原固定小 backbone 的数据表示容量失败仍单列，不能说它已经完整检验或否定 normal-direction 机制。已读取 `audit/N1_INITIAL_REPRESENTABILITY.json`：其中24行保留先前配置的 raw receiver-tail 审计，文件实际全局范围约 **0.55390–0.90068**；应引用精确行值，不能把口头近似 .55–.85 当覆盖所有记录的范围。文件中的 rank=12/16 字段不表示旧 core 的 receiver-data rank，也不代表换 backbone 后重算的剩余 tail。
- 父线程保留原264次新 N1 运行中的两次 QP active-box 舍入导致 TRF 初值越界失败。当前共同 corrector 在进入 TRF 前统一将初值投影回同一原 box，属于 bounded QP→bounded TRF 的数值交接修复。没有为这两次历史失败重新运行，因此历史失败仍为 FAIL，不能回填成 PASS。
- QP 分支仍直接返回有限 L-BFGS-B 迭代的 result.x；未返回停止状态的 P2 可观测性限制继续保留。无约束 profile 步差界不自动覆盖 active-box、有限求解误差或远处非线性路径。
- 本次关闭的是已指出的 same-span 数值实现阻断项及配置／术语接线，不是声称整个 N1 搜索算法已证明全局收敛、成像胜过成熟 CIE，或所有实现限制都不存在。原始失败、容量变更、费用与实验版本须在最终报告分开保留。N2 的封存结果不受本次关闭影响，也没有新增 N2 运行。

**结案判断：代码层 P1 已关闭；生产接线 E0 PASS 已独立核对。pure-span 单独数值回执暂按父线程报告保留来源限定；剩余 P2 与科学有效性边界继续保留。**

父线程随后澄清：pure-span 的 rank=0、step=0 原本仅由 implementation owner 在消息中报告，并不存在遗漏的原生测试日志。父线程将该消息级结果保存为 `audit/N1_PURE_SPAN_OWNER_RECEIPT.json`，明确标记 `OWNER_REPORTED_RESULT`。本审阅将它仅作**间接 owner-reported 控制证据**引用，不作为本审阅者独立复跑、原生终端日志或已读取的新数值运行；代码修复本身已静态核对，生产 E0 的独立文件证据仍为上面的 `E0_PRODUCTION_WIRING.json`。这一来源限定不阻止关闭代码层 P1，也不升级其科学证明力度。
