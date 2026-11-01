# N1：受控辅助电流释放实现合同

本实现把新增电流方向作为有限的搜索自由度，最终材料仍交由完整 Maxwell 正演及共同材料 corrector 验收。代码没有 GT、标签、已知另一根或隐藏 seed 参数；材料是完整实向量。科学结论、E0 独立审计和 E1 冻结由父 Codex 负责。

## 文件及入口

| 文件 | 作用 |
|---|---|
| `src/current_normal.py` | 实电流坐标、多源 lift、共享实材料法向、current 商空间诊断 |
| `src/joint_current_model.py` | CSI/CIE 真非线性辅助目标、完整实导数、成本账本 |
| `src/profile_schur.py` | joint-GN 剖面化、硬/软释放、同秩替换和物理任务收益 |
| `src/current_mode_selection.py` | 公共 Fourier/receiver/Krylov/random 池、预筛、付费物理 JVP |
| `src/controlled_release.py` | 有限正负 seed、固定坐标阶段、软释放 rollout、完整物理回读 |

主入口：

```python
from controlled_release import run_release
result = run_release(model, Y, m0, config)
# result: m, candidates, history, cost, stop, physical_score, config
```

`model` 使用父实现的 `G,E,S,n,shape,pack,state,tangent,regularizer`。`state(m)` 返回 `j,Q,Y,lu`。所有数据与状态 residual/derivative 都通过同一个 `model.pack`，因此保留冻结 R4 的 source/receiver 排列。若没有 `tangent`，允许复用 `state['lu']` 的直接 tangent solve，并显式计费；本研究已有的包装器提供 `tangent`。

`result['m']` 只更新为已完整物理回读且原目标更好的材料。`candidates` 同时保留完整物理 score、原数据相对 residual、辅助目标/数据/state loss、辅助 defect、辅助材料和 current、duplicate/qualification 标签。辅助 residual 对应 `m_auxiliary,j_auxiliary`；如果配置了 corrector，它与校正后物理材料不同。

`config['corrector']` 可接受父提供的 `corrector(model,Y,m,config)` 回调；也可以由父对输出 proposal 统一分配 corrector 配额。后者是 E1 默认集成方式。未提供阈值时，`PHYSICAL_VALIDATED` 只表示完成原 Maxwell 回读，不宣称是精确数据根、合格新 basin 或 GT 恢复。

## 实坐标和正确导数

一个实电流基张量为 `(n,p,d)`，系数严格为 `d` 个实数。一个复空间 atom 在每个 source 上具有独立复系数，所以占 `2*p` 个实自由度；不是一个实自由度。`spatial_modes(U,p)` 完成这个 lift。替换要求新旧实际实自由度相等。

辅助目标为

\[
\Psi(m,j)=\tfrac12\|W_y\mathcal R(Sj-Y)\|^2+
\tfrac\rho2\|\mathcal R(P(m)[j-m(E+Gj)])\|^2+\mathcal R_m(m).
\]

父的共同 `model.regularizer(m)` 返回 value、gradient 和共同 IRLS 局部先验曲率；后者不是非二次先验的精确 Hessian。物理 gain 使用该明确声明的先验二次模型。LM 是联合求解器阻尼，不加入物理收益或测量信息。

取 `Q=E+Gj,e=j-mQ,L=I-diag(m)G`。状态导数是

\[
D(Pe)[h,\delta j]=PL\delta j-P\operatorname{diag}(h)Q+(DP[h])e.
\]

`chart='csi'` 时 `P=I`；`chart='cie'` 时 `P_i=beta/(1+beta*m_i)`；另提供明确命名的 `cie_normalized`，其分子为一。二者都完整使用 `P'_i e_i`。`state_weight` 是共同固定的额外归一化。非物理 rollout 不删除这个导数。

## 剖面 Schur 和连续软释放

`profile_schur` 消去 retained current，形成

\[
H=B^T(I-P_A)B+\Lambda+\lambda I,\quad
g=B^T(I-P_A)b+\ell.
\]

对 `Z=(I-P_A)A_C=QR` 使用带 pivot 的秩揭示 QR，`R=Q.T@Z` 保留原 candidate 系数坐标。秩阈值参考原 `A_C` 的尺度；纯 retained-span 候选的 roundoff 不会被归一化成假方向。current 商空间也使用原 candidate 尺度。`D=Q.T@B,v=Q.T@b`。硬释放是 `H-D.T@D,g-D.T@v`。软释放在 `u=R c_C` 上加 `gamma*||u||^2/2`，因此 `alpha=1/(1+gamma)`，而非缩放基列：

\[
H_\gamma=H-\alpha D^TD,\quad g_\gamma=g-\alpha D^Tv,
\quad\gamma=\max(0,\eta_0/\eta_{target}-1).
\]

所有材料 proposal 都先解实际 box 二次问题，再对实际受限位移计算物理 gain。不会将未受限步的 score 归给 clip 后的步。`operation='replace'` 的旧组和新组都从同一 retained core 计算增量；`operation='augment'` 增加一组。`operation='dressing'` 只做 `C -> C+U K`，不产生新范围；它是 raw/dressed 同空间的 no-op 对照，不能算机制收益。

在非线性 rollout 中，QR 正则约束每一步的新增 current 增量。接受条件使用固定阶段的真辅助目标。`rho_schedule` 若存在，会明确改变阶段目标；每步重建导数，没有跨阶段的 BFGS 历史。

## 公共池及任务排序

所有 selector 都获取、实化和评分同一个 8 组公共池。默认 family 按 Fourier、receiver、Green-Krylov、random 循环，每组一个空间 atom。API保留的原型默认 backbone 是一个 receiver atom 加一个 Fourier atom，用于复现首轮容量失败；正式修正实验明确使用 `backbone_receiver_atoms='observable'`、`backbone_receiver_rtol=1e-3`，补入全部可靠接收方向。池仅依赖公开算子、观测、共同材料锚点和预先固定随机数。实际4/6网格的active实DOF为168/216，见正式ALGORITHM和冻结run manifest。

`selector` 可为 `material_step,energy,normal,random,fourier,krylov,receiver`。后面三个只在同一池中按 family 优先和固定原顺序排列，仍支付完整池的获取和评分费用。`normal` 是按实际候选动作的法向能量排序的明确消融，不是主方法。

`material_step` 先用冻结锚点 Born map 生成 proposal 分数；它没有尾误差界，记录为 `Born_no_error_bound`。若父明确提供 `physical_J`，则可使用该共同已付费参考。默认只有固定数量 finalist 获得真实物理 tangent 核验：旧步 JVP 缓存一次，每个 finalist 再付一个完整 source block。原物理 gain 包含共同先验的梯度、曲率两侧作用，不含 LM。finite 探索不要求此 gain 为正。

normal trace 先在固定 `||L V||` current 度量中取 retained 商空间、整组白化。source-complete 复空间组的 trace 在本配置下可恒为 `1-1/(2p)`，所以它仅作为角色诊断。另记录新旧联合解的实际 current 动作及其法向比例；零动作单列。

## 有限 seed 和一次明确修正

在 residual QR 坐标中，先用 `v+D p` 的归一化方向；该向量为零时使用 `D H^-1 D.T` 最大特征向量。`seed_rule='random'` 或 `'fourier'` 提供同池对照，并使用同一可行尺度规则。局部响应是 `h(u)=p-H^-1 D.T u`，实际带 box 的材料点仍用 QP 求解。

最初固定半径版本在公开 DEV 低材料锚点失败。父在 E1 前明确授权一次修正：`seed_material_rms=2` 改为声明的**最大**有限材料 RMS，current/defect cap 不变。默认 `bounded_seed_scale=True`：

1. 同时测试 `+tau_requested` 与 `-tau_requested`。
2. 若失败，测试零尺度；若零点本身不满足 cap，只报告 `NO_FEASIBLE_TESTED_SCALE`。
3. 在 `[0,tau_requested]` 上测试最多八个共同二分尺度，每次两个符号都解 box QP 并真实计算 defect。
4. 返回一个共同绝对 `tau` 的已测试可行 pair，并标记 `CAP_LIMITED`，记录 requested/achieved RMS 和原请求失败 receipt。由于非线性及活动面变化，不声称可行性单调或找到了全局最大尺度。
5. 任一 achieved RMS 小于 `.02` 都记 `NO_OP_MATERIAL_SCALE`；无有效法向动作也拒绝。

`bounded_seed_scale=False` 保留最初固定尺度负控制。有限探索中的非正局部物理 gain 不会自行删除候选。固定阶段只允许 retained current 更新，因此固定的 candidate 物理分量不随新 QR 改变；随后在同一个原始 current basis 中解除约束。所有失败、重复、拒绝与尺度搜索试验均计费。

## 父冻结的 E1 配置键

```python
scale_E = np.linalg.norm(model.E)
config = dict(
    chart='csi', rho=1., data_scale=np.linalg.norm(Y),
    state_weight=1/scale_E, eta_target=.5, lm=1e-4,
    backbone_receiver_atoms='observable', backbone_receiver_rtol=1e-3, backbone_fourier_atoms=1,
    group_atoms=1, pool_groups=8, finalist_jvps=2,
    episodes=1, hold_steps=3, joint_steps=8,
    finite_seed=True, execute_local=False,
    seed_material_rms=2., seed_min_material_rms=.02,
    bounded_seed_scale=True, seed_scale_trials=8,
    seed_current_scale=scale_E, seed_defect_scale=scale_E,
    seed_current_cap=3., seed_defect_cap=3.,
    rollout_current_scale=scale_E, rollout_defect_scale=scale_E,
)
```

一个 episode、finite pair、禁用 local slot，最多产生两个物理 endpoint，供父分配相同的两个 60-evaluation TRF corrector 配额。无 finite seed 的对照使用 `finite_seed=False,execute_local=True,joint_steps=8`。若增加 episode 或开启额外 local slot，会增加 proposal 机会，必须另行冻结；默认代码的可配能力不是扩展 E1 预算的授权。

## 成本和停止

`cost` 顶层保留父冻结后端 Counter 的实际差值，包括 `LU,actual_complex_RHS,source_RHS,tangent_RHS`；缓存命中不会伪造 LU。`auxiliary` 单列 current/Green/receiver 作用、dictionary 获取、QR、Schur/QP、局部 objective、尺度搜索、seed 插入和失败分支。`wall_seconds` 覆盖全部本次调用和可选 corrector。`expected_lu/expected_complex_rhs` 仅是无后端计数时的显式回退诊断，不与实测 Counter 相加。

支持 `max_wall_seconds,max_lu,max_complex_rhs,max_full_states,max_green_vector_actions`。`max_full_states` 限制 API 调用次数，实际 LU/RHS 以 measured backend 为准。预算中断保留已经验收的 incumbent。父在外部执行 corrector 时，须将其成本加回同一总账；本模块不会把它称为免费一步。

## 已完成的 DEV 实现验证

以下只用于可执行性和导数；不支持成像优势或新物理解的主张。DEV 公共包存在披露材料字段，但算法未引用该字段或材料误差评分，未使用另一端点生成 seed。

| 检查 | 结果 |
|---|---:|
| 五个模块 Python 编译 | PASS |
| 原 candidate 完全包含在 retained span：硬/软释放 | rank=0，材料步差=0，新增系数=0 |
| 原 current candidate 完全包含在 retained current span | real DOF=0，duplicate=True |
| 非物理 CIE 点的完整 `(h,delta j)` 导数中心差分相对误差 | `7.627445097894848e-11` |
| 上述点的非零本构 defect norm | `.3941309346966354` |
| 固定尺度负 receipt：DEV 4x4、6 source、m=.1，rho=1、共同归一化 | 1 LU，24 复 RHS，约 `.0123 s` BLAS1 |
| 该次真实 finalist gain（selected receiver_5） | `.008641730262715122` |
| 该次 candidate 实 DOF / retained 实 DOF | `12 / 24` |
| source-complete normal trace | `.9166666666666666` |
| 固定 RMS2 的实际正/负材料 RMS | `1.3871923929646286 / 1.3646533686188342` |
| 固定 RMS2 的 current/E 正/负 | `4.513913365732972 / 4.5116255184647684` |
| 固定 RMS2 的 defect/E 正/负 | `10.172651197955652 / 9.500225941665677` |

最初固定半径 pair 均被共同 box/预算拒绝，没有执行 rollout 或 corrector。这一失败 receipt 保留，不能被修正后的可行 pair 覆盖。父的独立 E0 和 E1 才判断新实现的数学条件、真实物理收益及科学范围。

修正后以同一个冻结 DEV 配置仅执行一次完整实现验证：3 次 hold、8 次 release，未改变参数，未运行 corrector。公共 `|tau|=3.9946358505923047`，请求值 `8.246990143158307`，状态 `CAP_LIMITED`；20 次带符号可行性评估包含请求点、零点及八次共同二分试验。正/负实际材料 RMS 为 `.6741048757473197/.662266383054832`，current/E 为 `2.194405357286386/2.1936017270807646`，defect/E 为 `2.9935934838802467/2.7844395743096673`。得到两个完整物理 endpoint，耗时约 `.10213 s`，实测 3 LU、36 复 RHS。物理相对 residual 约 `.9851401`，辅助相对 residual 约 `.906917`；第二个材料为 near-duplicate。此验证只确认有限流程可运行，明确不支持盆地或成像质量主张。实现随后冻结，后续只修复集成缺陷。

## 独立扩展：native_tsom 仿射变体

这一节是父 Codex 在早期 N1 的表示能力与 APCS 容量诊断后单独授权的方法变体。它不是对前述负结果的静默改写；旧 flag、旧材料锚点方法及已保存结果继续保留。先前 observable-backbone 虽有较大总 current DOF，但可见方向占据了大部分自由度，不能据此宣称已在合理 native Fourier-TSOM APCS 容量上检验过 N1。

新选项 `affine_mode='native_tsom'` 只允许 6×6 网格。它采用公开观测确定的固定 DPCS：

\[
j_D=V_d\Sigma_d^{-1}U_d^*Y,\qquad
j=j_D+U_A a,\qquad P_D=I-V_dV_d^*.
\]

receiver 相对 cutoff 固定为 `1e-3`，6×6 预期 DPCS rank=16。重用 `cie_tsom_baseline.fourier_apcs(shape,mf=2,Vd)` 的原 native 四角 Fourier 投影，其 span rank=16。以该 span 的正交基前 15 列作共同 retained APCS，原第 16 列作 removable old group。每个候选只替换这一列，因此始终是 16 个复空间 APCS atom、每 source 独立复系数；6 source 时自由实 DOF=`16*2*6=192`。DPCS 是固定向量，没有自由变量 DOF。完整 dark 空间有 20 个空间方向；本变体匹配 mf=2 的 16 维 span，不声称已实现 20 维完整 APCS。

公共八组池的第一个候选 `fourier_native_old` 精确恢复原 mf=2 span。其余候选按 Fourier/Krylov/random 循环，全部先投影到 DPCS 正交补和共同 15 维 APCS 的商空间。新增 Fourier 候选来自原 mf=2 span 之外，避免把同一旧方向多次包装成替换。固定 `selector='fourier'` 保持池中原顺序，选择第一个 old candidate；所有 selector 获取与评分同一八组池、相同实 DOF。未默认加入额外付费精确 tangent 候选。

在 4×4 上，dark rank=4，若使用 core3+1 则替换商空间只剩一维，任何满容量选择都是同 span 的 no-op。因此代码明确拒绝 native 4×4，未新增该网格实验。

### 物理状态与 joint 状态严格分开

字典返回 `joint_anchor_current=jD`。`evaluate_modes(...,joint_current=jD)` 和非线性 rollout 都从这个非物理 joint 点开始；它不等于 `state(m)['j']`，代码不假装真实材料物理流形属于 native 仿射表示。一般材料导数完整保留 `DP[h]e`。APCS 更新不改变 DPCS 分量，endpoint 记录 `fixed_dpcs_change_norm`。

完整 Maxwell 的 `state(m)` 单独用于原物理 `rF`、所有物理 JVP、完整 nonlinear readback 和 incumbent。`whiten_physical` 与 joint 的 `whiten_data` 是不同入口。显式 `physical_data_scale=||Y||` 时，物理残差严格为原 `pack(F-Y)/||Y||`，不会继承辅助目标的 source weights。原有分支未设置 physical override 时继续原 whitening 行为，可重放旧配置。

joint 目标与 constant-beta CIE 的 source 归一化一致：

\[
\tfrac12\sum_s\frac{\|Sj_s-Y_s\|^2}{\|Y_s\|^2}
+\tfrac12\sum_s\frac{\|\beta j_s-r(m)(Q_s+\beta j_s)\|^2}{\|E_s\|^2}
+\mathcal R_m(m),\quad
r_i(m)=\frac{\beta m_i}{1+\beta m_i}.
\]

配置使用 `chart='cie',beta=.2,rho=1,state_weight=1`，并设置 `state_source_whitener=1/||E_s||`；它乘入 residual、全部 material B 和 current A，包含非物理点完整 `DP e`。数据 whitening 为 `1/||Y_s||` 在同一 `model.pack` 坐标中的重复权重；默认由 `pack((1+1j)*ones_like(Y)/norm(Y,axis=0))` 生成，调用者明确传入的数组不会被覆盖。

该 joint-GN 变体与 native CIE 的目标和 mf=2 span 可以相同，但其正交坐标、GN/LM 与 nonlinear acceptance 不是原论文 PR-CG 的逐步算法复刻。raw normal 指标仍参照另付费求得的精确物理 current，明确标记 `raw_unweighted_exact_physical_anchor_diagnostic_only`。这是非物理 joint 状态处的参考诊断，未证明 source-weighted 目标下的法向度量不变性，也不把物理锚点的 0/2/4 阶结论直接归给 jD。

固定 Fourier 的 old/new 材料步差严格为零，不应因此阻止其实际辅助目标下降。此变体对所有控制统一使用 `operation='replace',local_gain_min=-inf`。仍只允许完成原物理回读后的材料更新 incumbent，corrector 由父统一分配。`finite_seed` 在 native 默认关闭；若父明确启用，seed 从同一个 jD 开始，normal reference 仍单独保留真实物理状态。

### 单次公开工作验证配置及 receipt

父指定公开 case：`obs_5a04b0ce1306e592`，即 `data/e1_case_ids.json` 中 side6/layoutuniform 的字典序最小项。只读取其公开算子及 Y；起点为公开常数 `m=.1`，没有读取 GT、truth 目录或材料评分。只执行下列一次优化，没有配置搜索、内部 corrector 或新 4×4 运行：

```python
config = dict(
    affine_mode='native_tsom', operation='replace',
    chart='cie', beta=.2, rho=1., state_weight=1.,
    state_source_whitener=1/np.linalg.norm(E, axis=0),
    data_whitener=model.pack((1+1j)*np.ones_like(Y)
                            /np.linalg.norm(Y,axis=0)[None,:]),
    physical_data_scale=np.linalg.norm(Y),
    selector='fourier', release='hard', local_gain_min=-np.inf,
    finite_seed=False, execute_local=True, joint_steps=80,
    max_wall_seconds=2.6, episodes=1, pool_groups=8,
    group_atoms=1, finalist_jvps=2,
    backbone_receiver_rtol=1e-3, lm=1e-4,
)
```

机器可读的实际执行摘要如下；wall 包含本次全部 setup、selection、rollout 和原物理 readback。未包含父随后另行分配的 corrector：

```json
{
  "kind": "N1_NATIVE_AFFINE_PUBLIC_IMPLEMENTATION_RECEIPT",
  "case_id": "obs_5a04b0ce1306e592",
  "grid": [6, 6], "sources": 6,
  "initial_material": 0.1,
  "status": "EPISODE_LIMIT",
  "selected": "fourier_native_old",
  "fixed_DPCS_rank": 16,
  "core_APCS_rank": 15,
  "active_APCS_rank": 16,
  "active_real_current_DOF": 192,
  "physical_step_difference_from_old": 0.0,
  "physical_gain_from_old": 0.0,
  "joint_anchor_is_physical": false,
  "joint_anchor_defect_norm": 1.0729281172651546,
  "initial_joint_data_loss": 4.7759950828823846e-30,
  "initial_joint_state_loss": 0.13341595373061915,
  "joint_steps_declared": 80,
  "joint_steps_attempted": 65,
  "joint_steps_accepted": 50,
  "joint_steps_rejected": 15,
  "physical_endpoints": 1,
  "physical_score": 0.0680850780934225,
  "physical_relative_residual": 0.36898724257709303,
  "auxiliary_objective": 0.014103638291296025,
  "auxiliary_relative_data_residual": 1.4459913576188847e-13,
  "auxiliary_state_loss": 0.014094352790196768,
  "auxiliary_defect_norm": 0.5493694457762865,
  "fixed_DPCS_change_norm": 1.6010772509218584e-13,
  "wall_seconds": 0.8156615840271115,
  "actual_LU": 2,
  "actual_complex_RHS": 30,
  "source_RHS": 12,
  "tangent_RHS": 18,
  "paid_physical_JVPs": 3,
  "green_vector_actions": 86948,
  "receiver_vector_actions": 77783,
  "joint_linearizations": 66,
  "joint_objective_evaluations": 178,
  "rank_revealing_QR": 156,
  "material_QP_solves": 140,
  "small_Schur_solves": 140,
  "receiver_SVD": 1,
  "native_APCS_SVD": 2,
  "receiver_tail_SVD": 9,
  "Fourier_atoms_generated": 52,
  "DPCS_projection_columns": 16,
  "internal_corrector_calls": 0
}
```

source-complete real lift 在当前实现中使用稠密 current tensor，因此辅助 Green/receiver 计数包含实际数组乘法中零填充的 source 列，不把这些操作隐藏成免费传播。wall 和后端 LU/RHS 都单独保留。

同一公开 case 只额外执行一次必要的非物理 source-weighted 联合方向 wiring check，未作完整 E0 重放：`(h,delta j)` 导数中心差分相对误差 `7.992940407649957e-10`；直接 paper-CIE 目标与 `JointCurrentModel` 的 value 差为 `0`；原 mf=2 Fourier span 的投影相对误差 `1.2907944888141269e-15`；`whiten_physical` 与 `raw/||Y||` 差为 `0`。native APCS receiver 数值秩使用明确的 `1e-11*||S||` 阈值，避免把理论零动作的 roundoff 解释成接收可见方向，诊断 rank=0；该阈值不改变基、目标或 optimizer。

这个独立 wiring check 在常数 `m=.4` 另构造了一次新的 36 阶物理参考：额外 1 LU、6 source RHS，须在全局实现验证成本中另列，不计作上面 2 LU/30 RHS 方法 trial 的免费准备。

这次验证确认了 192 实 DOF、固定 DPCS、非物理导数及物理/辅助目标分离的实现。辅助数据 residual 接近零而原物理 residual 仍为 `.369`，所以不能据此宣称新的物理解或 Maxwell 恢复成功。实现至此冻结；父负责后续匹配 constant-beta/mf2 CIE、完整成熟 CIE、公开五起点、总成本配额和科学判断。


## 最终执行裁决（5748条主批次唯一回执中的N1证据）

1. 旧physical-anchor原型的容量失败、observable修正及CSI/CIE状态权重实验均保留。旧分支新确认未显示超越强基线组合的覆盖；不能据此否定原生APCS空间。
2. 最后native-affine分支保持固定DPCS、16个free APCS、6源，即192实current DOF。660次开发运行中所有ordinary selectors及matched CIE为6/12联合恢复；finite task/Fourier/random/Krylov均7/12，共享同一个新增高对比双岛。这是generic finite探索的有限正面证据，不是task专属。
3. 冻结后另生成独立六对象，648预测封存后评分。普通Fourier、task、成熟CIE、直接材料seed的clean/noisy K3联合均6/12；finite task为3/12，Fourier/random finite为6/12。开发finite增量没有转为新的确认覆盖。
4. **新鲜次级正面证据**：task local在一个clean、幅度16的弯条上恢复真实材料，两个公共初值终点NRMSE≈.000310、r≈8.9e-5；普通Fourier/强基线没有恢复该观测。IoU=.625与GT自身Otsu相同，按固定联合门仍FAIL。不能把这项材料恢复抹掉，也不能改门后声称总体联合提升。
5. 已付费日志显示task选择随机字典的random_native_3，192→192DOF、各一个公共corrector seed。两起点一正一负的原物理局部预测收益最后都成功，所以简单正收益门不是全局恢复证书。
6. 针对该已披露对象另做24次posthoc组件诊断，普通normal和soft也在各一个clean起点恢复同一材料；random/energy未恢复，所有noisy控制未恢复。它不是新的fresh统计样本，更不能证明A17分数唯一必要。

最终为 **有限机制正面／通用优越性HOLD／N1_576优越性NOT_RUN**。源版本、输入、原物理端点、实际成本与失败见主报告及独立审阅。当前方法是A17任务相关性启发的auxiliary-current研究，不是A17原生Galerkin exchange后端或新全局收敛算法。

