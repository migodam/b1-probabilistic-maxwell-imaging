# R13 基线源码与资格审计

本次只审查源码、冻结配置，以及旧的已暴露 A r10_twod_0004 receipt 元数据。没有读取任何 fresh 材料、fresh 材料质量/评分或 GT，没有启动运行，也没有修改活动源码或在线 C2。最终科学解释、运行完成性验收与材料结果验收仍由 parent 决定。

## 1. 先区分两个完全不同的 GlobalDual

| 稳定名称建议 | 当前入口 | 实际未知量与目标 |
|---|---|---|
| **Dual-CIE / global PR** | `global_dual` | DPCS/Fourier 辅助 current + 实材料交替；使用一个跨 source 的 PR 系数/步长；不是物理 total-field source-Gram 目标 |
| **physical GlobalDual (β=0)** | `physical_globaldual` | 只有 N 个实材料变量；每次评价解原 Maxwell，使用 total-field source Gram 及其完整导数 |
| feedback (β=.2)，若另行运行 | 历史 `method='feedback'` | 使用 `Q+.2J` 的物理 feedback Gram；不能标成历史 GlobalDual |

历史 `research/research_fbg_r3_experiments/src/full_maxwell_corrector.py:69-79` 对 `method='globaldual'` 明确强制 `beta=0`，而 feedback 默认 `.2`。当前 `src/run_experiment.py:129-151` 明确调用同一个 `MaxwellModel.data_value_gradient(m,Y,beta=0,reference=initial)`，并保持未缩放的 L-BFGS-B `ftol=1e-12, gtol=1e-8, maxls=30`。proposal 的 `maxiter=3000,maxfun=12000` 是本轮冻结预算合同，不应写成历史每次运行都用了这两个上限。

### physical GlobalDual 具体目标

原物理实现位于 `research/research_fbg_r3_experiments/src/physics.py:79-102`。令原 Maxwell 物理状态的总场为 Q，接收残差为 R，则 β=0 时

`F_GD(m)= c/2 * tr[R (Q*Q)^(-1) R*] + R_prior(m)`，

`c=||Q_ref||_F² / (ns ||Y||_F²)`。

这里的 ns×ns Gram 随 m 改变；原解析伴随包括 Gram 导数项，不能替换为固定权重的普通物理 LS。当前 runner 原样复用它，未删去该项。其每个目标试探点已经支付物理 LU/adjoint/Gram，随后 `record` 单独计算原物理 Phi；**主选择器始终按原 Phi，而不是 F_GD 选点**，与历史 corrector 的 physical incumbent 习惯一致。

历史 reference 是均匀材料 `.1` 的物理状态 (`full_maxwell_corrector.py:72`)；当前 reference 固定为同一个公开初始物理状态的 Q (`run_experiment.py:133-148`)。本轮 `confirmation_local_reference16.json` 明确声明均匀 `.1` 初值，因此该 reference 约定相同。若未来改为 NN/非均匀初值，当前 reference 会改变 c 相对于 prior 的权重，应明确称为适配，不能自动声称与历史 reference 完全相同。

reference 显式传入也避免 `physics.py:89` 默认另取 reference 导致 cache 指向 reference 状态的问题。当前传入的 initial row 确实含所需 Q；该目标无需 reference 的 Y 字段。

未发现当前 physical GlobalDual 的关键接线错误。但 source-Gram 使用原 Cholesky 路径；失败时应保留实际失败与已付费候选/成本，不应臆造收敛或自动伪逆恢复。该物理 GlobalDual 不是一个原作者 NFFT/CIE 算法复现声明。

## 2. 24×24 标量合成基线的实际原生设置

以下表描述普通 scalar synthetic 合同。共同材料变量始终 N=576，材料盒 [.02,20]，原 prior `smooth=1e-4,tv=1e-4,epsilon=.01`。令 rD 为由公开 observations/noise 冻结的 receiver DPCS rank，F_mf 为 native unitary Fourier mask，P=I−VdVd*，JD 为原 DPCS offset。

| 方法 | current / material 更新 | native stages 与 beta | current 初始化 | current family 与 DOF |
|---|---|---|---|---|
| CIE (`cie_single`) | sourcewise PR+ current；原共享实材料 scalar-CIE readout，含同 prior | mf3/6/10 各500次，再 mf12 1500次；βJ=βm=1/.2/.2/.2 | cold：系数0，J=JD；m=m0 | `JD+P F_mf A`；每 stage 的独立实 current DOF 为 `2*rank(PF_mf)*ns` |
| Dual-CIE (`cie_dual`) | 同上 sourcewise PR+，material readout 不变 | 同 masks/次数；βJ始终1，βm=1/.2/.2/.2 | 同 CIE cold | 同 CIE |
| global-PR Dual-CIE (`global_dual`) | 全 source 一个 global PR+/line step；scalar material readout | 与 Dual-CIE 相同 | cold | 同 CIE；改变 PR 聚合方式，不改变空间 |
| A+ (`a_plus`) | global PR scalar-current；原 onehop 材料损失的 Möbius/four-color MM readout，默认2 sweeps | 与 Dual-CIE 相同 | cold，保留原 R7 约定 | 同 CIE；材料仍576自由实变量 |
| G19 (`g19`) | 同 global PR scalar-current；固定公开 Green-strong Galerkin denominator，完整 moving-denominator 导数，材料 L-BFGS 内迭代上限8 | 与 Dual-CIE 相同 | warm：原物理 j(m0) 以同一 LS 投入原 affine DPCS/Fourier family；不是直接保留 j | `JD+P F_mf A`；同 masks 的 current family 与 CIE 相同 |
| G19 full APCS | 与 G19 完全相同的 PR/material/beta/次数 | 每个 stage 从第一步就用 mf12；beta 与各 stage 次数保留 | 物理 j 的原 LS warm；DPCS 被 JD 固定，因此仍可有投影损失 | `JD+P F_full A`；独立实 current DOF=`2(N-rD)ns` |
| G19 full current | 与 G19 完全相同的 PR/material/beta/次数 | 每个 stage 从第一步就 mf12 | 完整 Fourier correction 可将初值保留至 roundoff | `JD+F_full A`；独立实 current DOF=`2Nns`；JD 是可被纠正的 offset |
| NFFT-equations / joint L-BFGS-B | 无 CIE contraction；全 current correction 与材料 jointly 优化普通 LSIE Eq3–6 + 同 prior；明确的 L-BFGS-B 适配 | 无 beta/mf continuation；runner `maxiter=6000`，adapter `ftol=1e-12,gtol=1e-8,maxls=20` | Z=0，所以 J=JD；m=m0；不是 J=0，也不是物理 j warm | complete Fourier 的 unitary 等价空间 Z；`2Nns+N` 个联合实变量 |
| physical GlobalDual β=0 | 原物理 source-Gram objective + 完整解析伴随，L-BFGS-B | 无 Fourier/current continuation；proposal cap3000/12000与实际 CPU deadline | m=m0；J 始终为该 m 的物理解 | 仅 N 个独立实材料变量 |
| direct physical L-BFGS-B | 原 Phi/adjoint，无辅助 proposal | `maxiter=6000,maxfun=24000,ftol=1e-12`；返回梯度/目标乘1e6、scaled gtol .01、maxls30 | m=m0 | 仅 N 个独立实材料变量 |
| physical Full-J TRF (`physical_trf`，后续固定16 reference) | 原物理 residual + exact analytic Full-J + 同 prior residual；trust-region reflective least squares | whole-process cap；`max_nfev=3000,x_scale='jac',ftol=xtol=1e-10,gtol=1e-8` | 已声明盒内的 m0；没有 CIE/current 初始化 | 仅 N 个独立实材料变量；Full-J 是导数矩阵，不是额外材料/current未知量 |

source 依据：

- 原 mask、beta、归一化、inner60 和 outer coefficient_tol0：`research/research_max_struct_r6/src/benchmark_solver.py:34-80`。
- A+ 的 arm 绑定与 global/sourcewise 区分：`research/research_max_struct_r8/src/benchmark_solver.py:13-33`；实际 MM 默认两 sweep：`research/research_max_struct_r7/src/aplus_readout.py:103-128`。
- G19 的 fixed public Green basis/inner8：`research/research_max_struct_r10/src/lowrank_denominator.py:86-123`，经 R11 `full_material_readout` 复用。
- 三种 current family、cold/warm 政策、actual rank/DOF：`src/full_current.py:102-118,195-220`；runner cold/warm 显式设置：`src/run_experiment.py:188-190`。
- NFFT 的无 halves LSIE objective、JD source normalization与初值：`src/nfft_control.py:23-35,85-100,168-176`。`wd=1/||Y_s||²`，`ws=1/||JD_s||²`；普通 CIE 使用 `ws=1/||Einc_s||²`，两者不能混称相同目标。

普通 CIE 的 one-step PR+ 使用精确固定材料 quadratic line step，material 更新后继续原 PR 历史，stage 边界 reset；没有为审计分块重启优化。标量 dual-beta G19 并没有被证明共享一个整体共同目标或 Lyapunov 函数。

## 3. Fourier 容量、实际 rank 与旧 metadata

mf3/6/10/12 在24²上生成36/144/400/576个 native Fourier coefficients/source，但 DPCS 投影后独立 rank 必须按实际空间计算。**旧 mf12 已经是 full APCS**，所以 full APCS-from-start 主要消除早期截断/continuation，而不是增加旧方案的最终 APCS 容量；full current 另外解除固定 DPCS。

仅提取旧的已暴露 A r10_twod_0004 `representation_v1` receipt 元数据，rD=13，独立 spatial ranks 为36/144/400/563；对应 real source-current DOF 为864/3456/9600/13512。full current 是13824；另加共同576个实材料变量。提取记录为 `BASELINE_EXPOSED_METADATA.json`，未打开 HISTORY、材料数组或 GT。

旧原始 receipt 的 `state_real_DOF=27648` 是两份复数场的 **stored ambient scalars** (`4*N*ns`)，不是 CIE 独立系数 DOF，更不是 direct 变量数。没有修改旧 receipts。当前 runner 使用 `stored_field_real_scalars` 与明确独立计数、`auxiliary_receipt.current_spaces` (`run_experiment.py:317-319`)。DR 两场 iterator 的 ambient scalars 也不能当作传播/材料约束流形的内禀维数；AP 内禀维数未证明则保持 unknown。

## 4. 付费物理审计机会与选择

新 confirmation 使用共同的绝对 **process CPU 2秒 bucket**。runner 在一次 proposal 物理记录后更新下一个 bucket (`run_experiment.py:70-98`)；辅助 CIE/G19 在完整 material update 后查询，NFFT 在 accepted joint update 后查询。旧 development/representation 缺少该 callable 时仍按原 every50 fallback，不能倒写成相同审计协议。

共同 bucket 不意味着审计次数、精确时刻或 wall time 相同：一次 inner solve/line search/LU 可以跨过 bucket；CIE 原生初始/stage audits、NFFT 初始/terminal audits仍保留且支付成本。物理 GlobalDual/direct/corrector 每个 objective trial 本来已支付物理 forward，因此保留每个合法、及时的物理候选，而不是丢掉已有候选来制造相同数量。

物理 GlobalDual proposal 的 `record` 记录原 Phi 与另一个字段中的 Gram auxiliary value (`run_experiment.py:136-148`)。所有主选择器都是最小原物理 Phi；secondary 为原 discrepancy 内最小 prior。共同 corrector 从已付费原 Phi incumbent 启动，使用 beta=None 的原物理目标，而不是继续 source-Gram 或 auxiliary-CIE 目标。低辅助值、低物理残差或更多审计机会本身均不是材料正确性证据。

## 5. 有限预算与原生 schedule 完成性

`confirmation_core48` / `confirmation_mechanism16` 的总 process cap120秒，proposal deadline为绝对72秒，包含 imports/setup；corrector至多追加48秒且受总 cap约束。并非 proposal 从 setup 后另获72秒。`confirmation_local_reference16` 的 cap360秒、proposal216秒、corrector最多追加144秒，native masks/betas/辅助步骤数不变，且配置明确禁止把 Apple local CPU 时间与 i9 pooled。

较大 cap 的目的可表述为降低工作截断；**不能从配置宣称一定完成原生 schedule**。必须查实际 termination、完成 PR steps、reached stage 和已付费审计记录。receipt `status='COMPLETE'` 是子进程正常写出产物，也可能同时有 `termination.proposal='PROPOSAL_PROCESS_CPU'`，不是3000次 schedule 完成或数值收敛声明。

旧暴露 A4 的 metadata 给出一个不涉及材料质量的实例：CIE/Dual/global-PR Dual 完成3000步；A+仅2071、G19仅2007、full APCS仅1873、full current仅2159步后到 proposal CPU cap。这些都是同一120秒总预算下的实际已付费 prefix，不能改写成各方法都完成全 schedule 的比较。它也不能证明任何 fresh/local360 运行的完成状态。

需要同时报告：整进程 CPU（external exit优先，internal receipt采样不含最终 exit）、physical LU/RHS、projector G/S/H 操作、setup与失败成本、实际已审计 prefix。原 aux physical delta 已包含在最终 physical model counters 中，不能重复相加。有限预算截断证据与科学材料结论应分开。

## 6. 模型/来源资格边界

- 标量 CIE/Dual 是已有 source-audited real-material/common-prior CIE 适配，不是未改动的原作者复杂材料算法逐位复现。A+/G19是项目自己的明确后续材料 readout/metric，不能当作原论文原法。
- Austria classic cold CIE/Dual 使用原 `paper_cie_fft`、公开 metadata 的 mf/beta/stopping、rD15、[0,20]与不缩小的 equal-area radius。当前 factory 为 mf5/8/10、βm8/2/.5；Dual/global PR βJ8，CIE βJ=βm；paper loop 支持每-stage coefficient tolerance。G19 为保留物理 warm 而使用 R8 warm loop+equal-area Green/OneHop factories，并保留自己的 global coefficient stopping，不能称为原论文 G19。
- 3D vector 和 measured source-specific Fresnel 使用各自兼容的原生适配，不能套本表的 scalar current rank。已有接口不支持的 source-specific G19/A+ 或 scalar-NFFT vector transplant不能伪造为已比较基线。
- 新 physical GlobalDual 入口当前审计范围也是普通 scalar/common-port MaxwellModel；vector 的同名接口不接受 beta/reference，measured adapter 明确排除 feedback 参数，不能把这一 scalar 接线自动扩展到它们。本轮 local-reference16 manifest 声明的正是 scalar TM 2D。
- NFFT-equations control 明确为 **NOT_ORIGINAL_NFFT_SOM**：完整 Fourier/JD family与 Eq3–6规范化的已声明实材料/box/prior/joint-L-BFGS适配。2018短论文未完整给定原 CG/材料步/初始化/停止，当前没有声称 mature author reproduction。详见 `NFFT_PRIMARY_AUDIT.md` 和 `NFFT_EQUATIONS_CONTROL.md`。
- 当前 `physical_globaldual` 与旧 `global_dual` 标签必须始终分开；不能通过改图例把 global-PR Dual-CIE 补写成已跑的 physical GlobalDual。

本次未发现新 physical GlobalDual 接线的关键技术错误。审计未评价任何 fresh 材料结果，也未宣称任一方法完成了尚未验收的运行、优于另一方法或满足论文级复现资格。

## 7. 后续 physical Full-J TRF 源码补充审计

本节对应 parent 声明的 TRF source `2cd5516cb` 与 ready `32bfdfa4c`；查阅时工作树 HEAD 为后者。只读 source/config，没有读取该活动队列或其任何 fresh 材料结果，没有运行测试反演或修改源码。

`src/run_experiment.py:120-146` 新增 `physical_trf`，走 `phase='direct'` 与 whole cap，而非60% proposal + L-BFGS corrector。residual 是

`[pack(Y(m)-Y)/||Y||F ; r_prior(m)]`，

Jacobian 是

`[FullJ(m)/||Y||F ; J_prior(m)]`。

这与旧 `research/research_fbg_r3_experiments/src/full_maxwell_corrector.py:47-65` 的 native/TRF 分支一致；`trf`、`x_scale='jac'`、三个容差也相同。`max_nfev=3000` 与绝对 process CPU guard 是本次预声明运行预算，并非历史默认120评价预算。旧 wrapper 会先 clip m0到盒内 (`:10-13`)，新分支直接使用公开 m0；本次声明初值已在盒内。若以后换成盒外/QP舍入种子，不能忽略这项 wrapper guard 差异。

### Prior residual 与原 prior 是否相同

已直接核对 `research/research_fbg_r3_experiments/src/physics.py:62-71` 和原 `Regularizer.evaluate` (`research/research_z_alg_r4_diversity/src/physics.py:215-222`)，而不是从“least_squares”名称推断。令 `d=Dm`、`q=sqrt(d²+epsilon²)`、边数为 M：

`r_prior=[sqrt(2*smooth/M)*d ; sqrt(2*tv/M)*d/sqrt(q+epsilon)]`。

因 `d²/(q+epsilon)=q-epsilon`，其半平方范数恰为

`.5||r_prior||² = smooth*mean(d²) + tv*mean(q-epsilon) = R(m)`。

源码的 prior Jacobian 使用第二段对 d 的实际导数

`sqrt(2*tv/M) * [1/sqrt(q+epsilon) - d²/(2*q*(q+epsilon)^(3/2))]`，

再乘同一 D；不是把 TV residual 当常权差分。原数据 `pack` 保持复数残差 Frobenius 范数，Full-J 的 source/receiver/real-imag 排列与 pack 一致 (`research/research_z_alg_r4_diversity/src/physics.py:98-107,151-154`)。因此该 scalar 合同下的 **least-squares residual objective** 等于原 Phi，prior residual/Jacobian表达一致；这不意味着 Gauss–Newton/TRF 的局部曲率是完整非线性目标的精确 Hessian，也不意味着全局收敛或材料可辨识。

Full-J 沿用保留 LU 的 receiver-transpose多 RHS解析实现，计数为 `FULL_ANALYTIC_J` / `FullJ` / `receiver_transpose_RHS`；不能把 `njev` 当作额外 N 个材料未知量，也不能只数物理 forward LU 而漏掉 Jacobian RHS与线性代数成本。

### 选择、预算和覆盖

每个 residual trial 的已付费物理状态按原 Phi记录，保留合法及时候选；raw是初始物理状态，corrected与primary为已付费最小Phi。`audit_predictions.py` 当前虽然仅把字符串 `direct` 分类为 direct，但 TRF所有运行候选的 phase也为 `direct`：其 raw pool仍只含initial，corrected无corrector pool时fallback到primary。因此现有索引推断在此路径一致，无须修改活动源码。尚未据此验收任何活动 TRF receipt。

正常返回时保存 message/nfev/njev/optimality摘要；预算异常用 whole-process终止标志，不能宣称达到求解器收敛。配置的360秒 whole cap用于固定16个已声明 fresh对象的同机 reference，另一个已暴露B诊断为120秒；它们不能与120秒 i9 core结果直接合并为同硬件同预算榜。配置文本中继承的“corrector最多144秒”不适用于 TRF direct分支，实际分支使用整段 cap且没有另一次corrector。

此补充配置的 task列表实际是16 fresh条件运行 +1 exposed诊断；当时复制的 `expected_tasks=96` 是陈旧摘要，审计/交付分母应由实际task列表与parent总账纠正，不能沿用它。Parent本轮明确的固定总账为828个fresh条件运行+89个exposed方法运行=917个固定方法运行，rescue的 `5×实际触发数` 与其它诊断另列；这些数目不是917个新材料对象或900个新对象。

未发现新增 TRF residual/Jacobian接线的关键错误。它填补的是原 Full-J TRF 强基线，不是未经验证的新优化原理或原作者所有条件逐位复现。完整 optimizer/internal-Jacobian轨迹仍未保存，见 `AUXILIARY_STATE_COVERAGE.md`；不能从材料快照承诺精确 TRF continuation。
