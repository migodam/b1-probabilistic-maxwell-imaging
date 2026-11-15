# R14 跨物理复用接口审计

审计日期：2026-11-15。范围：只读现有代码、公开输入打包清单、历史物理资格和费用回执；没有打开材料 GT 数组、没有打开 Fresnel heldout 数组、没有运行新求解或实验，也没有修改任何历史文件。本文是工程复用清单；新公式、科学可行性和是否启动扩展由主线程判断。

**结论：存在真实的 3D vector Maxwell 复用路径。** 最小实质性扩展应独立包装已有 `VectorPR3D.current_pr` 和 `VectorCIE3D.material`，保持每个 voxel 一个实 χ、三个耦合电流分量、原生 2×2×2 空间块、原 prior/权重、全局 PR 和所有四组合相同的 L-BFGS-B material inner allowance。它是 **3D 几何对应实验**，不是 R7 四色 Möbius/Arc 的逐字基线。Austria 的 equal-area scalar seam 更小；真实 Fresnel 还必须保留逐源 receive/DPCS/APCS，不能只删除 `physics_family` 检查。

## 1. 实际材料量、电流和 Green

| 项目 | 已有实际接口和约定 | 代码位置 |
|---|---|---|
| 通用 tensor forward | `TensorDDA(points, volume, k=2π)`；`solve(epsilon,directions,polarizations)` 接受 `(N,3,3)` 复对称 ε；它不是已有逆问题的材料参数化 | [vector_backend.py](../../research_max_struct_r1/src/vector_backend.py#L96-L151)，特别是 L125–144 |
| 当前 3D inverse | `VectorInverse` 和 R13 `DenseVectorModel`/`VectorPublicModel` 使用 **N 个实 isotropic χ=εr−1**；R13 明确 `lossy=False`；不优化 polarizability | [vector_inverse.py](../../research_max_struct_r7/src/vector_inverse.py#L1-L20)；[vector_adapter.py](../../research_max_struct_r13/src/vector_adapter.py#L18-L41) |
| 物理自项 | `s = −1/3 + i k³ V/(6π)`；`G = V*Goff + s I`；`Goff` 是无自项的 dyadic kernel，三个分量不能拆开 | [vector_backend.py](../../research_max_struct_r1/src/vector_backend.py#L96-L111) |
| χ 与 CM+RR | R7 原文已明确 `alpha = V*chi/(1-s*chi)`，恰是 CM+RR；优化 χ 与优化 alpha 不同 | [vector_inverse.py](../../research_max_struct_r7/src/vector_inverse.py#L1-L5) |
| 电流约定 | inverse 的 `J` 为 susceptibility density current `j=χQ`，不是含体积的偶极矩 `p=Vj`；原 TensorDDA `VectorState.J` 又表示 `−ik*j`，名称必须按 API 区分 | [vector_backend.py](../../research_max_struct_r1/src/vector_backend.py#L135-L151)；[PHYSICS_EVIDENCE_BOUNDARIES.md](../../research_max_struct_r13/delegated/physics/PHYSICS_EVIDENCE_BOUNDARIES.md#L7-L11) |
| 行序 | `J,Q,E` 为 `(3N,P)`，cell-major x/y/z；`d=np.repeat(chi,3)`；一个材料梯度必须把同 cell 的三个分量求和 | [vector_adapter.py](../../research_max_struct_r13/src/vector_adapter.py#L51-L74)、[L93–101](../../research_max_struct_r13/src/vector_adapter.py#L93-L101) |
| 全物理状态 | `L=I−diag(repeat(chi,3))*G`；一份 retained dense LU 解全部 P 个 source RHS；`Q=E+GJ`；`Y=S J`；状态不合格残差 `>1e−9` 报错 | [vector_adapter.py](../../research_max_struct_r13/src/vector_adapter.py#L51-L75) |

原 `TensorDDA.solve` 在 total-field 坐标使用 `I−Gχ`，R13 current-state 在 J 坐标使用 `I−χG`；它们是相同构成关系的不同未知量。不能据矩阵左右顺序把二者当不同物理，也不能把 tensor forward 的材料自由度冒充已存在的 tensor inverse。

### 原生 local block 的确切含义

`VectorCIE3D.__init__` 使用 `block=ijk//2`，只有同一 **2×2×2 voxel block** 内的相互作用保留；三个 component 维用 `repeat(...,3)` 同时展开。`H=(G-sI)*same_block_mask`，`Hop=csr_matrix(H)`。side=6 时有 27 个空间块；每块含 8 个 voxel、24 个 complex current components，是 **24×24 dyadic block**。同 voxel 的完整 3×3 self block 被 `sI` 扣除；不同 voxel 之间的交叉分量保留。[vector_cie3d.py L25–29](../../research_max_struct_r8/src/vector_cie3d.py#L25-L29)

这对应 2D 中“在每维两格的不重叠 local-cell block 内做一次 off-diagonal 交换”的几何规则；它并非固定四个 scalar rows。把 3D 排成 2D stack、把 xyz 当三个独立材料变量、只保留同分量、把 2×2 scalar 块随意贴到 3D，都不保持这个已有接口。

## 2. 两个 metric family 已经存在的精确入口

以下只转录现有已实现关系，不提出新推导。设 `d=repeat(χ,3)`，`r=J−dQ`。

| family | current 入口 | material 入口 | 已有行为 |
|---|---|---|---|
| P0 | `b.current_pr(z,Y,'CIE_TSOM',beta0,mf,True)` | `b.material(z,J,'CIE_TSOM',beta0)` | `P0=diag(beta0/(1+beta0*d))` |
| P1 | `b.current_pr(z,Y,'A_TSOM',beta1,mf,True)` | `b.material(z,J,'A_TSOM',beta1)` | `b1=beta1−s`，`den=1+b1*d`，`u=beta1*r/den`，`a=d/den`，`P1 r=u+a*(Hop@u)` |

Current 实现：[vector_pr3d.py L11–34](../../research_max_struct_r8/src/vector_pr3d.py#L11-L34)。它先构建 `K=I−diag(d)G` 与 `rhs=diag(d)E`，再用同一个 P 变换两者；`SB=S@B, KB=P K B`。PR axis 为 `None` 时是全局实 PR+；保留旧 `gam,w,direction,mf`，使用 source-normalized exact quadratic line search。无需新 current optimizer。

Material 的全导数已经包括 **dP/dχ**：[vector_cie3d.py L68–77](../../research_max_struct_r8/src/vector_cie3d.py#L68-L77)。P1 分支含 `du=-beta*(Q+b*J)/den²`、`v=R+Hop* (conj(a)*R)`，以及 `conj(R)*(Hop@u)/den²` 项；随后在每 voxel 内汇总 xyz，再加入原 prior。P0 分支也保留自己的完整材料相关分母导数。不能用固定 P surrogate、不能漏掉 xyz 汇总。

R7 的 `VectorInverse.material(...,'A_material_joint')` 已有同类完整导数，且非根残差检查明确激活了 dP 项：[vector_inverse.py L61–76](../../research_max_struct_r7/src/vector_inverse.py#L61-L76)；[NONROOT_DERIVATIVES.json L2–35](../../research_max_struct_r7/RESULTS/vector_3d/NONROOT_DERIVATIVES.json#L2-L35)。它的 acquisition/current solve 与 R8 不同，只能作为已存在的导数实现和历史检查，不能改称新 3D factorial 已验证。

### 绝不能直接移植的 Arc

`VectorPR3D.aplus` 在 L44 把候选端点硬编码为 `[0,4]`；`optimize_pr` 在 L51 也固定这个 box。R13 的原公开 vector box 是 `[0,8]`。此外该 coordinate routine 是逐 voxel 顺序更新，既不是 R7 四色 2D MM，也不是 R13 adapter 已支持的算法。[vector_pr3d.py L35–62](../../research_max_struct_r8/src/vector_pr3d.py#L35-L62)

因此本次最小方案 **不调用 `aplus`、不调用 `optimize_pr`、不调用 R7 `mobius_material_readout`**。四组合都调用同一 inherited `material` value/gradient，再用同一 L-BFGS-B inner allowance=12 和 `model.bounds`。这保持真实 3D 变量/全导数/原 prior，同时不把它包装为 R7 Arc 复现。

### real 与 complex readout 边界

R8 `VectorCIE3D(lossy=True)` 可以用 `z=[Reχ,Imχ]`，其 `decode/encode/grad` 和 `material` 实现了分量导数；原 `optimize` 另有 Imχ∈[0,1.2]。但 **R13 public adapter 只构造 `lossy=False`**，公开材料、prior 和 bounds 均为实 χ。[vector_cie3d.py L32–37、L68–77、L91–104](../../research_max_struct_r8/src/vector_cie3d.py#L32-L104)；[vector_adapter.py L18–31](../../research_max_struct_r13/src/vector_adapter.py#L18-L31)

拟议四组合只覆盖这个已有实 isotropic 接口；不顺带增加 lossy、tensor 或 polarizability 逆问题。

## 3. prior、source metric、basis 与原始状态规则

| 项目 | 必须复用的原值/规则 |
|---|---|
| 原物理选择目标 | `Phi=mean_s ||S J_s−Y_s||²/||Y_s||² + prior`，**没有 1/2**；不套用 scalar 全 Frobenius 归一化 |
| 原 prior | `2e−4 * mean_edges |χ_i−χ_j|² + 1e−6 * mean_cells χ_i²`；3D 最近邻，不是 scalar smooth/TV prior |
| auxiliary weights | data 为每源 `1/||Y_s||²`；state 为每源 `1/||E_s||²`；`rho=10`；material 目标的 state 部分含 `rho/P`；这些因子保持原实现 |
| DPCS | `b.prepare(Y)`，S 的一个 SVD，rank 阈值 `.1*s[0]`；`Jd=V ((U*Y)/s)` |
| APCS | 原生 3D IFFT，每个 Fourier mode 同时生成 xyz；投影 DPCS 后默认 `whiten_apcs=True`，经 SVD 取 orthonormal range |
| mf | 原 `mf=1→2`；side6 时实际 basis dimensions 为 24→192；这不是 R14 2D 原生非正交 Fourier coefficient metric |
| PR state | `b.prstate={gam,w,direction,mf}`；同 mf 保留 memory；**mf 改变时 gam 置零并重置 PR**，没有 2D core 的 nested-mask zero-padding transport |
| 物理传播 | dense dyadic G/S、retained LU；3D IFFT 只负责 basis construction，不能声称 FFT Maxwell complexity |

目标/梯度 API：[vector_adapter.py L80–116](../../research_max_struct_r13/src/vector_adapter.py#L80-L116)。prior/source/basis：[vector_cie3d.py L25–37、L43–77](../../research_max_struct_r8/src/vector_cie3d.py#L25-L77)。PR stage state：[vector_pr3d.py L17–34](../../research_max_struct_r8/src/vector_pr3d.py#L17-L34)。R13现有 wrapper schedule/inner：[vector_adapter.py L145–197](../../research_max_struct_r13/src/vector_adapter.py#L145-L197)。

`run_vector_auxiliary` 当前只接收成熟 CIE/Dual/Global 名称，固定使用 `CIE_TSOM` 两个 block，明确拒绝 scalar A+/G19。其限制不能通过 relabel 绕开；下面的独立 wrapper 才是有明确变化范围的新入口。

## 4. 给主线程审阅的最小 wrapper API（未实现、未运行）

```python
run_vector_metric_pair(
    model, Y, m0, pair, config, check, audit_callback,
    *, audit_due=None, state_callback=None, cost_callback=None,
)
```

参数合同：

- `model` 必须为实际 `VectorPublicModel(public)`；使用 `b=model.legacy`，不重新实现 Green、acquisition、prior 或 material derivative。
- `pair` 仅 `00/01/10/11`。第一位选 current；第二位选 material。0 一律使用该阶段 beta0，1 一律使用该阶段 beta1。
- `config` 显式保存 `stages`、`material_inner_maxiter=12`、`global_pr=True`、`CPU_cap=120`、`proposal_CPU_cap=72` 和 source/prior/bounds 摘要。四组合相同 config；没有 Arc 参数。
- 若依照 R13 **global-Dual current** 阶段：前 60 cycles 为 `mf=1,beta0=6,beta1=6`，后 60 为 `mf=2,beta0=3,beta1=1`。若主线程选择原 CIE-global reference，beta0 后半应为1，必须另行明确冻结；不从方法名称暗中推断。
- 其 dispatch 仅为 `current_pr(z,Y, kind[pair[0]], beta[pair[0]], mf, True)` 与 `minimize(lambda q: b.material(q,J,kind[pair[1]],beta[pair[1]]), z, jac=True, method='L-BFGS-B', bounds=model.bounds, options=...)`。其中 `kind={'0':'CIE_TSOM','1':'A_TSOM'}`。
- Material 四组合统一继承 R13 options：`maxiter=12,ftol=1e−12,gtol=1e−8,maxls=25`。`maxiter` 是上限；保留实际 `nit/nfev/success/message`，不强迫恰好12 accepted steps。不要改变 `.rho/.E2/.prior/.bounds`。
- `check()` 使用原整进程 CPU；proposal 至72秒，随后共同 original-Phi corrector 至120秒。72与120是绝对进程时刻，不是在 setup 后重新计时；任何超限完成操作仍收费，候选不纳入预算内选择。若120个原 cycles先完成，按 parent 固定合同进入共同 corrector；若希望填满72秒继续 proposal，必须另列新 iteration allowance，不能默默改成另一算法。
- `audit_due()` 是原每2整进程CPU秒首次完成 material update 的 paid audit；保留初始/阶段/结束 audits。`audit_callback(m,paid_state,Phi,metadata)` 不再求解同一 material。选择依据调用 `model.data_terms`/`physical_score`，不复用 scalar runner 的 D 归一化。
- common corrector 使用原 `model.data_value_gradient` 与 `[0,8]^N`；原 full physical prior、数据目标、noise radius 均不变。新120秒 cap与历史180秒协议不同，时间只在本次同 host/同合同内比较。
- `state_callback` 在完成 update 后接收 `m,J,Q,b.prstate,b.bases,b.V,b.Jd,E2,rho,stage/cycle,source weights,bounds,flags,cost` 的 detached 数组；若实现 checkpoint，必须保留 vector 原生 stage reset 规则。2D `StatefulDualMetric` 不能充当这个 state container。
- setup/current/material/failed attempts/physical LU/source RHS/adjoint RHS/SVD/basis matrices/Hop actions/序列化都收费。P1 current 增加 `Hop@K` 和 `Hop@rhs`，material 增加 `Hop@u` 与 adjoint action，必须新增明确 counters；不能沿用旧 CIE counts 假称它们免费。
- 报告名称应为“same local-cell geometry 的3D two-metric analogue，继承 full analytic material readout”，不能叫“R7 Arc baseline”或字面 CIE2019 author reproduction。该 native-beta0 合同下 `00` 也不能叫原 global-Dual：其 scalar material beta 现同为 beta0，而原 Dual material 后半是1。

此 API 规格已足够把扩展变化限定为独立 metric dispatch 和完整账本；没有要求新 solver/backend。主线程仍应审查 actual matrix/current/material dispatch、box 和候选选择，再决定是否启动16条新轨迹。

## 5. 已有可用 2–4 case 的真实路径与费用

实际旧 R13 worktree 根：

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r13_worktree/research/research_max_struct_r13`

本次目录检查确认以下 **四个 public NPZ 实体均存在**，各约1.9MB；没有打开其数组。当前 R14 git worktree 中只有 tracked packing receipt，不能把相对路径当成已复制的数据。

| 固定 case | 实际 public 相对旧 R13 根路径 | 可支持的范围 |
|---|---|---|
| slotted_high_frame | `delegated/physics/new3d/runner_public/slotted_high_frame.npz` | 已曝光 matched6³ 实材料 inverse；216实 χ，648复 current rows，60源、180接收分量 |
| unequal_coupled_steps | `delegated/physics/new3d/runner_public/unequal_coupled_steps.npz` | 第二独立旧几何，同模型/目标 |
| mixed_crossed_bracket | `delegated/physics/new3d/runner_public/mixed_crossed_bracket.npz` | 第三独立旧几何，同模型/目标 |
| slotted_high_frame_crossgrid8 | `delegated/physics/new3d/runner_public/slotted_high_frame_crossgrid8.npz` | 原frame的8³ generation→6³ inverse压力对；**不是第四独立物体** |

公开打包声明无truth、G/E/S形状见 [PACKAGING_RECEIPT.json L1–74](../../research_max_struct_r13/delegated/physics/new3d/runner_public/PACKAGING_RECEIPT.json#L1-L74)。冻结的 acquisition、box、prior、noise与对象数见 [FROZEN_PROTOCOL.json L4–37](../../research_max_struct_r13/delegated/physics/new3d/FROZEN_PROTOCOL.json#L4-L37)。

最小2-case可固定取前两项，不按旧恢复分数筛选；3-case覆盖三独立旧几何；4-case再保留crossgrid反例压力。四组合对应8/12/16 trajectories。它们都是已经曝光的机制后续，不是新盲确认。整 cube 是未知域，没有 GT support mask。

### 实际旧运行与资格测量（不是新预算预测）

- `QUALIFICATION.json` 的 adapter/legacy prediction、J、Q、objective、gradient differences 均为0；原 physical residual 为 `2.1878e−15`；physical objective directional derivative relative difference `1.0432e−8`。该回执共有6 LU、360 source RHS、120 adjoint RHS，0.115705 CPU秒/0.126420 wall秒；这是有限接口/解析参考检查，不是3D恢复效果证据。[QUALIFICATION.json L2–16、L135–151](../../research_max_struct_r13/delegated/physics/QUALIFICATION.json#L2-L151)
- 原2-cycle interface回执为0.073229 CPU秒，2 paid LU/120 source RHS，2 PR steps、4 material calls、1 receiver SVD；含另一次原physical-gradient检查的总数是3 LU/180 source RHS/60 adjoint RHS。**只可称 interface exercise**，不能把它当 substantive factorial结果或据此估计完整 inversions。[INTEGRATION_CHECK.json L432–505](../../research_max_struct_r13/delegated/physics/INTEGRATION_CHECK.json#L432-L505)
- R7 非根 dP 检查使用 constitutive relative residual .204808；A-material directional derivative relative error `8.6841e−9`，1 physical LU/6 source RHS/9 material evaluations，wall .017500秒；旧6-plane-wave acquisition与本次60-dipole acquisition有别。[NONROOT_DERIVATIVES.json](../../research_max_struct_r7/RESULTS/vector_3d/NONROOT_DERIVATIVES.json#L2-L35)
- 历史4-case×3成熟vector基线的12个 Windows children 均 returncode0/COMPLETE；queue外部整进程CPU总数1007.359375秒，单条58.75–97.515625秒，包含 imports和最终serialization。来源为 [queues/physics_v1/RECEIPT.json](../../research_max_struct_r13/queues/physics_v1/RECEIPT.json) 中这四case的 `cie/dual_cie/global_dual` 字段，未读取offline材料分数。它们不是本次120秒合同的速度承诺。

旧 `RESULTS/physics_v1/<case>__global_dual/RECEIPT.json` 四条实际记录如下；这些 receipt 在上述原 R13 worktree 中实存，当前R14 tree并未复制。physicalcounts是全旧run（含corrector），proposal列来自`auxiliary_receipt`：

| 旧 case | whole_process CPU s | proposal CPU s | proposal cycles / material calls | proposal LU / source RHS | full-run LU / source RHS / adjoint RHS |
|---|---:|---:|---:|---:|---:|
| slotted_high_frame | 92.5625 | 24.1875 | 120 / 1155 | 14 / 840 | 660 / 39600 / 38760 |
| unequal_coupled_steps | 96.6875 | 23.421875 | 120 / 1221 | 13 / 780 | 732 / 43920 / 43140 |
| mixed_crossed_bracket | 86.75 | 23.015625 | 120 / 1218 | 12 / 720 | 621 / 37260 / 36540 |
| slotted_high_frame_crossgrid8 | 60.21875 | 22.8125 | 120 / 1138 | 13 / 780 | 357 / 21420 / 20580 |

每条proposal还有1 receiver SVD、216 basis columns、120 current matrix builds、12960 KB columns、12960 SB columns，以及7200 KJd/SJd source columns。原 corrector 有的报告projected-gradient收敛，有的报告process-CPU停止；不能从status COMPLETE推定全部优化收敛。回执的旧protocol cap为180秒，不能拿这些总CPU直接声称新72/120合同表现相同。

### continuum 与 mismatch 边界

原 CM+RR dyadic backend 有独立解析 Mie 对照，但没有一般 continuum certificate。旧高contrast eps9/radius.09 的 staircase h=.02/.015 far errors为9.015%/7.064%；volume-averaged χ对应19.76%/16.93%。h=.03等粗格负结果仍保留；近机器精度方程残差不消除边界/材料偏差。[PHYSICS_EVIDENCE_BOUNDARIES.md L32–34](../../research_max_struct_r13/delegated/physics/PHYSICS_EVIDENCE_BOUNDARIES.md#L32-L34)、[MIE_REFINEMENT_DIAGNOSIS.json L111–144](../../research_max_struct_r13/delegated/physics/MIE_REFINEMENT_DIAGNOSIS.json#L111-L144)

8³→6³ frame的旧clean mismatch为8.006%，大于所加1%噪声；保留原noise-only99%球，只能报告coarse-model失配压力，不能抬radius把它变成noise认证。[PHYSICS_EVIDENCE_BOUNDARIES.md L24–26](../../research_max_struct_r13/delegated/physics/PHYSICS_EVIDENCE_BOUNDARIES.md#L24-L26)

## 6. Austria equal-area：可合法复用的较小 seam

现有公开路径在原R13根的 `delegated/physics/classic_public/austria_150to32.npz` 与 `austria_matched32.npz`，本次只核对文件存在。二者同一物体、两个离散化条件，不是两独立材料。

`paper_physics.acquisition(side,frequency)` 明确使用 h=2/side、a=h/√π、20 plane waves、40个3.75m接收点；内部/外部Green和analytic self一致。[paper_physics.py L32–45](../../research_max_struct_r8/src/paper_physics.py#L32-L45) `EqualAreaGreen` 复用相同有限Toeplitz FFT action，仅允许正确的equal-area quadrature radius；它不是把.45h disk物理半径扩大后仍叫不重叠圆盘。[L19–29](../../research_max_struct_r8/src/paper_physics.py#L19-L29)

`PaperOneHopAbsorption` 对dense G **直接调用原 `OneHopAbsorption.__init__`**；故2×2四cell P1及完整材料导数完全相同。对`ShapeOnlyGreen`，它只从相同equal-area公式生成实际local blocks/self；不是另一个P1。[paper_physics.py L80–91](../../research_max_struct_r8/src/paper_physics.py#L80-L91)

已有R13 `run_full_current` 在`physics_family='equal_area_disk'`时：cold mature loop用`paper_cie_fft`；warm G19 loop用原R8 PR code加private `EqualAreaGreen/PaperOneHopAbsorption` seam。[full_current.py L85–101](../../research_max_struct_r13/src/full_current.py#L85-L101) 因而扩展R14 2D core的family选择有具体历史依据，**不是删guard后仍调用原不重叠disk构造器**。

最小合法改动范围应仅是：显式family→已有Green/block factory，保留densepublicG、真实radius、原materialprior、box[0,20]、originalaudits；checkpoint记录family并按family恢复已有cachedspectrum/block对象。现有R14 `_bind` 的family拒绝在 [dualmetric_core.py L217–228](../src/dualmetric_core.py#L217-L228)，尚未修改。`EqualAreaGreen`另有`diagonal`属性；若未来调用其物理iterative接口，checkpoint也必须保存该属性。只使用其继承apply/adjoint时原shape/n/pad/spectrum足够，但不要遗失family身份。

配置不能用R14默认mf12/rank-noise选择替代。原 [classic_auxiliary_config L183–208](../../research_max_struct_r13/src/physics_inputs.py#L183-L208) 固定DPCS15、mf5/8/10、materialbeta8/2/.5、Dual/globalcurrentbeta8、zeroχ/APIC、box[0,20]、原prior；stage allowances来自public。cold paperloop使用逐stage coefficient thresholds，[paper_cie_fft.py L309](../../research_max_struct_r8/src/paper_cie_fft.py#L309) 与当前R14 global threshold不同。若保留R14 global stopping，必须标成新factorial合同，不冒称历史cold baseline逐轨迹相同。

经典原文合同实际cache在 [delegated/paper_baselines/PAPER_CONTRACTS.json L3–23](../../research_max_struct_r8/delegated/paper_baselines/PAPER_CONTRACTS.json#L3-L23)：NIE2016 forward150²/inverse64²、rank15、无附加prior；当前inverse32²/commonprior/corrector是公开已声明adaptation，不是author-code reproduction。历史noise是10%fixed-norm方向，ball为12%discrepancy，不是CN99。[PHYSICS_EVIDENCE_BOUNDARIES.md L28–30、L38–44](../../research_max_struct_r13/delegated/physics/PHYSICS_EVIDENCE_BOUNDARIES.md#L28-L44)

## 7. Fresnel：P1本身可复用，common-port state不能复用

可用原R13 publictraining包：`delegated/physics/fresnel_public/dielTM_dec8f.npz`、`twodielTM_8f.npz`，本次只确认文件存在。旁边`*_HELDOUT_OFFLINE_ONLY.npz`没有打开。

实际R8模型是equal-area scalar-TM，2GHz、24²/.16m、order2 empty-target multipole continuation，逐源实际receive rows；R7早期commonsector/.45h/cylindricalpilot是另一历史合同，不能拿其sharedS替代当前qualified路径。R8 `load_file` 在L26共轭total−incident与emptyincident到exp(−iωt)，`operators`在L30–37用h/√π；`SourceReceiverModel.receive/receive_adjoint`严格用 `Ss`。[fresnel_adapter.py L16–37、L90–108](../../research_max_struct_r8/src/fresnel_adapter.py#L16-L108)

R13可携模型为 `SourceSpecificPublicModel(public)`，输入`S_ports`或原S+`source_receiver_indices`，不打开observations/GT。[physics_inputs.py L64–91](../../research_max_struct_r13/src/physics_inputs.py#L64-L91)

已有合格的source-specific算子接口已经存在于 [fresnel_cie_fft.py L18–36](../../research_max_struct_r8/src/fresnel_cie_fft.py#L18-L36)：

- `SourcePorts(model.Ss)` 实现逐源 `S_s @ J_s` 与匹配adjoint；
- 每源一个DPCS basis `Vd[s]`，每源`NativeFourierAPCS` action/adjoint；
- 相同denseG和`PaperOneHopAbsorption`，所以local四cell P1不需改公式；
- 既有 [fresnel_solver.run L19–31、L67–70](../../research_max_struct_r8/src/fresnel_solver.py#L19-L70) 已保留diagonalsourceweights，`a_plus`分支可调用R7完整Möbiusreadout。

但R13已qualification的production wrapper `run_fresnel_auxiliary`只接受CIE/Dual/global，明确拒绝未资格的A+/G19 transplant：[physics_inputs.py L94–115](../../research_max_struct_r13/src/physics_inputs.py#L94-L115)。有source-specificAreadout代码，不等于R14四组合/statefulresume已经实验合格。

R14 core目前 `np.asarray(model.S)`、单个Vd、公共SU/adjoint不符合实际`Ss`；模型的commonS全receiver数也不匹配训练39rows。因此**仅启用equal_area_disk仍不足以支持Fresnel**。若做新statefuladapter，必须复用上述SourcePorts/逐源Fourier factory，checkpoint保存完整Ss/masks、逐源SVD/DPCS列表、旧PR及频率状态；所有四组合保留原diagonalmetrics，不宣称source-mixing covariance。物理audit仍用原source-specificreceiveadjoint和modelstate。

现有物理证据范围：1053trainactualchannels与711heldoutactualchannels、empty-fieldvalidationmedians2.515%/2.408%、15%training-normdiscrepancyball。噪声covariance与receivehorncalibration未知，不作CN99。[PHYSICS_EVIDENCE_BOUNDARIES.md L16–28](../../research_max_struct_r13/delegated/physics/PHYSICS_EVIDENCE_BOUNDARIES.md#L16-L28)

旧authorized heldout评估有27unique LU/972source RHS、wholeevaluator0.822801CPU秒/0.947316wall秒（本地Apple）；它的约9–10%heldoutfit和原corrector作用已写入旧报告，不是materialtruth，也不是新A+优势证据。[PHYSICAL_VALIDATION_REPORT.md L35–47](../../research_max_struct_r13/delegated/physics/PHYSICAL_VALIDATION_REPORT.md#L35-L47) 另有50-updateperiodicaudit接口回执：27逐源SVD，2physicalLU/54sourceRHS，0.366945CPU秒；它是callback付费检查而非算法比较。[FRESNEL_PERIODIC_AUDIT_CHECK.json L39–77](../../research_max_struct_r13/delegated/physics/FRESNEL_PERIODIC_AUDIT_CHECK.json#L39-L77)

## 8. 交付边界

工程上可实施且最小的下一步是上面的**真实3D current/material dispatch wrapper**，配合同一[0,8]实χbox、已有2×2×2dyadicblock、全导数materialL-BFGS-B、全局PR及父线程明确的72/120CPU合同。四个旧case公共输入实存，足够冻结16条conditional expansion；没有必要新建backend或假2Dstack。

尚未实施的缺口是：独立pairdispatch wrapper及其P1额外计数、完整vectorstate保存、newcontractmanifest和parentreview；尚无R14跨3D的实测factorial结果。Austria需要小的equal-areafactory/stopping合同适配；Fresnel还需要逐源stateful接口。本文不授权任何运行，不得把已有物理资格、旧方法完成回执或低残差改称新的metric恢复结论。
