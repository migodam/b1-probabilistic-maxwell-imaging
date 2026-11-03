# R5 公式与轨迹诊断独立审查

2026-11-03。只读审查；没有读取 GT/真实材料数组，没有执行 forward、Jacobian、HVP、公式恒等式测试或参数扫描，没有修改诊断实现。审查读取了 R5 推导、Stage I 协议、诊断/分析/轨迹记录源码，以及 R4 one-hop、R2 ordinary readout、FBG 物理包装器与冻结物理/regularizer API。另只读汇总已经生成的 168 个 checkpoint JSON。核心定理、原创性、实验支持何种科学结论均由父线程裁决。

## 结论与需要修正的输出口径

未发现使本批 finite-r 梯度、共同实材料 LS 或完整 Maxwell Jacobian 失效的实现错误。两个输出口径需要修正或显式限定；均不需要新物理计算。

1. **`exact_minus_GN_data` 混入 prior 曲率。** `src/trajectory_diagnostics.py:108–114` 的 `exact=hᵀH_Phi h` 来自 `inherited_hvp`，而 `gn=||J_real h||²/||Y||²` 只有 data GN。继承函数 `research_max_struct_r4/src/cooperative_geometry.py:43–55` 最后一行明确加入精确 prior Hessian。因此当前差值是 `hᵀH_residual-feedback h+hᵀH_prior h`，不能称纯 data residual feedback。建议保留 full-Phi 曲率，同时后处理分列 prior directional curvature、exact data curvature 与 residual-feedback curvature，或把旧差值明确称 full-Phi-minus-data-GN。已即时通知父线程和 diagnostics worker；对方确认安排后处理。
2. **ordinary/self 的零 transport 没有可解释的方向角。** 两种传播均为逐 cell 非零复标量，对共同实材料无正则自由 LS 的权重逐 cell 抵消，故解析 `d=d0`、`transport=0`。`frozen_ls` 用浮点 `ds-d0`，`cosine` 只检查范数是否大于零，因而对舍入误差也输出有限夹角。168 点已存 JSON 中 ordinary 最大 transport 范数 `5.796060485725511e-14`，self 最大 `5.656377614321901e-14`；例如开发 global 的 s0/i500 self transport 与实际步夹角余弦为 `-0.1717287`，这不是方向证据。应按解析结构把这两类 transport 方向角标为 `N/A / analytically zero`，原始浮差可以保留；无需引入人为数值 gate。已即时通知父线程和 diagnostics worker。

## 实化、source metric 与 P 的伴随

固定合同采用默认 `Omega_s=diag(ws)`，`ws_s=1/sum_i |Einc_is|²`，不同源共享同一实材料向量。该固定合同没有自定义 `state_metric/state_weights`。`analyze_checkpoint:134–135` 与原 `run_onehop_cie:266–269` 的该默认 metric 一致；诊断函数不是任意 full source metric 的通用入口，不能在未来传入 full Omega 后仍沿用此硬编码。

`frozen_ls:45–68` 中 `K[b,o,s,i]=P[b,o,i] Q[b,i,s] sqrt(ws_s)`，Gram 与 rhs 分别取 `Re(KᴴK)`、`Re(KᴴPr)`。这是把复 residual 乘 source 平方根后实化的共同实材料正规方程，没有把每源或实/虚分量变成独立材料。`d0_i=Re sum_s conj(Q_is) r_is ws_s / sum_s |Q_is|²ws_s`；`r-Q*d0` 是这一实化 C_mat 的正交补缺陷，顺序正确。

对 one-hop，代码构造 `P=(I+s diag(a)Goff)diag(beta/d)`，其中乘在右边的 `beta/d` 对输入 cell 缩放，和原 one-hop 一致。`analyze_checkpoint:170–174` 先算 `R Omega_s`，再算 `(I+s Goffᴴdiag(conj(a)))R Omega_s`，最后左乘 `diag(conj(beta/d))`；这是正确的 Pᴴ 作用。`-rho Re sum conj(Q) PᴴP r Omega_s` 因而是正确 frozen-metric 梯度。`phase_blocks` 只替换副本的 goff，正反相及 self 不修改原真实物理 G。

## 非根梯度与 connection

继承的 `onehop_material_value_gradient:101–120` 在固定 J/Q/Omega 时使用

```text
du_i/dm_i = -beta (Q_i+b_i J_i)/d_i²,
da_i/dm_i = 1/d_i²,
adj_u = R Omega_s + Goffᴴ diag(conj(a)) R Omega_s.
```

这同时保留 `du` 和 `da` 项。诊断先减去同一点 prior，再减 frozen 梯度，得到的 `g_connection` 正是 `rho <Pr,(dP)r>` 的实材料梯度；没有把 prior 混入 connection。正反相的半差/半和正确给出 `g_cross` 与 `g_quad`，含相同 rho。

ordinary 精确梯度和 frozen 梯度 `trajectory_diagnostics.py:185–191` 同原 `_readout:83–94` 相符，并且可以直接化简为

```text
g_connection,ordinary,i
  = -rho beta³ sum_s |r_is|² ws_s / (1+beta m_i)³.
```

本批 beta>0、box 为 [.02,20]，所以该项非正是结构事实；负 connection 推动的无约束方向是材料增加方向。它自身没有接收数据方向，不是 physical Phi/GT 收益证书；prior、box 与其他梯度仍影响真实 readout。self 项同样可由 `Re(b_i/d_i)` 控制符号，本批负值观察不能据此推广到任意复 Green/材料域。

`cross_value=2 Re<u,v>_Omega` 是残差平方范数展开的 cross，而 `g_cross` 是 `0.5*rho||u+s v||²` 中 `rho Re<u,v>_Omega` 的导数。因此 `g_cross=(rho/2) grad(cross_value)`，不是 cross_value 的未缩放导数。现有协议显式采用前一种 cross_value，数值计算没有因子错误；报告应保持这一口径。

## 冻结 C_mat、四 cell Gram 与实际 readout 的边界

四 cell Gram 的方向次序和全材料回填正确；每个 cell 恰在一个不重叠 2×2 block 内。对可逆 Gram，`ds-d0` 等于理论中的 `G_M^-1 C_matᵀ M(I-Pi_mat)r`。已存 168×4 类 LS JSON 全部 `fallback_blocks=[]`，本批没有触发伪逆语义。

潜在回退边界：若将来 `la.lstsq(gram,rhs)` 丢弃数值方向，`ds-d0` 除 transport 外还含 `-(I-G†G)d0`；不能继续无条件解释为精确 transport 公式。并且对 normal matrix 做 SVD 已将原 K 的条件数平方。现有回退/rank 字段可揭示该分支；本批未触发，所以不构成本批结果污染，也不需要补跑实验。

R5 推导第15、17、47行和 Stage I 第16–20行已正确声明 frozen LS 与真实材料算法的差别。代码另在相同 `m_before,J_after_current,Q` 上调用原 ordinary/one-hop readout，保留 nonlinear P、prior、box、确定性 seed 与原60步最大额度，随后用原完整 forward 评价。它是算法反事实，不能称已达到材料全局最优或精确驻点；应同时保留 `material_inner_success/message/nit`。self 是 goff=0 的 one-hop inner solve，已经明确不是历史 blockself 的逐字节重放。实际 `delta_m` 是单次 readout 步：recorder 每次 callback 都更新 previous，不是两个稀疏检查点间的累积步。

因此：根处或冻结 LS 的相关性可以作为局部诊断，不能被报告改写成实际带 prior/box 非根 readout 的公式、路径因果或 material recovery 保证。当前文字没有作此错误等同。

## 物理 Jacobian、曲率、评价位置与 LU

FBG `MaxwellModel.data_value_gradient` 返回的是 full Phi 值/梯度，含共同 prior。诊断正确另取 `g_physical_data=gphysical-gprior`，实际 Phi 差与线性预测均使用 full Phi。投影梯度使用同一 box。

冻结物理 `Model.jacobian:98–107` 使用保留 LU 的 receiver-transpose 多 RHS，构造 source-major 的完整 complex J。其 real J 将各 source 的 receiver real rows 与 imaginary rows 分组；`pack` 使用每 receiver 实/虚交错，二者只差行置换。`J_realᵀJ_real=Re(J_complexᴴJ_complex)`，所以当前 Gram、特征值与 `||J_real h||` 正确；不能把这个未重排行序的 J 直接乘 `pack(residual)` 做 residual-J contraction。本代码没有这样做。

complex Gram 的谱允许 complex coefficient directions，而实际材料是实数；物理实材料谱应以 real Gram 为主，complex Gram 只能作另列的代数比较，不是第二套实际实材料信息量。微小负 Gram 特征值及浮点 rank 已用 numerical semantics 标注，不应解释成严格核或负信息。

`inherited_hvp` 使用同一保留 LU 的 base adjoint、tangent、incremental adjoint，并使用精确 smooth-TV Hessian；没有把 `Regularizer.evaluate` 返回的 IRLS weights 当精确 TV Hessian。全 Phi HVP 本身正确，错误仅在上述差值字段的口径。

**评价位置需写清：** 当 `iteration==500` 时，`endpoint_geometry` 得到的是最后一次材料 readout **之前**的 `m_before`，方向是该次 readout 的 `m_after-m_before`；它不是保存的 stage_end `m_after` 的谱。此选择适合预测最后一次有限步，但“阶段端点谱”字样需要加此限定。counterfactual 也从同一 m_before 出发，比较口径一致。

LU 复用链条正确：当前 m 的 value/gradient、Jacobian、HVP 共享同一 state/factor；之后才计算 m_after 的完整 forward 和四种 counterfactual 的真实 forward。这些不同材料的 LU 有单独收费，没有以缓存命中冒充新 forward。小 Gram、complex metric RHS 与物理 LU/RHS 分开计数。`small4_metric_complex_RHS` 只覆盖 Pr 和 P-unexplained 两次 block action；K 构造另随 Gram-build 字段计，不能把该 RHS 字段当全部辅助算术的 FLOP 总数。

## 分析 runner 的范围

runner 不改采样 checkpoint；已有 completed checkpoint JSON 被复用，并在 receipt 中保留 historic checkpoint 成本与本进程成本的区别。其 `COMPLETED` 表示文件分析流程完成，不表示 readout inner 收敛、机制被证实、GT 改善或科学接受。CSV 未展平 ordinary 梯度和所有 endpoint 字段，但 JSON/NPZ 保留这些量，不属于数值漏算。没有因本审查新开实验。
