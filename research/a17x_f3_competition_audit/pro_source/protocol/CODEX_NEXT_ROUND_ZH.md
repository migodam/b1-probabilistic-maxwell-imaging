# ARCHIVED PRE-F3 PROTOCOL -- NOT AN AUTHORIZATION TO RE-RUN P0

New F3 evidence at pinned commit a323ef223731870b390f875c55f041aa8cc855b4 has completed the quotient prior, legal anchors, four-arm transport, short MH and path diagnostics. Start from `protocol/F3_FIRST_GATE.md`. The proposed 128 full states / 4096 RHS below is a potential ceiling only, not a permission to re-open F3/P0; P1 readiness and new budget authorization are separately required. F3 high-order arm D is NOT A17 Schur exchange, which remains NOT_RUN.

---

# Codex：A17-X 最小决定性实验

你负责执行 **A17-X / finite-source-transport pilot**。先读随附 `report/RESEARCH_A17X_FULL_REPORT_ZH.md`，重点是第 5–7、11、15–16 节。不要训练 NN，不要扩大 P1，不要修改 A17 原论文或覆盖任何旧决定。

目标只有两个，而且必须分别回答：

1. 在原先验支持内，有限 Maxwell 等效激发是否能产生有意义、全部照明相容的材料变化？
2. 对同一机制，A17 是否比精确局部 T-matrix / Woodbury、普通 full-J corrector 和 standard goal-oriented ROM 更可靠或更便宜？

不能只比较 GN+SVD，然后把所有 finite-feedback 收益记给 A17。

## 0. 工作区与已有资料

优先使用当前已有 B1 克隆，在独立目录 `research/a17x_finite_source_transport/` 保存新代码/结果。没有可用本地克隆时，先完成纯代数检查并报告实际缺少的资产；不要猜测 `/Volumes/...` 路径，不覆盖正在运行的 P1。

主要仓库：

```text
https://github.com/migodam/a17-maxwell-certified-exchange
https://github.com/migodam/b1-probabilistic-maxwell-imaging
https://github.com/migodam/a26-om-feedback-prior-completion
https://github.com/migodam/a30-maxwell-observable-material-manifold
https://github.com/migodam/c1-physics-state-generative-imaging
```

必须读取 A17 的 `closure_final/manuscript/CORE_THEORY.tex` 和 `CONDITIONAL_EXCHANGE_THEORY.tex`，以及 B1 commit `8b620b315b2b4b212d0f1e4c143ced49f22f02e0` 下 `research/research_f_intrinsic_transport/RESEARCH_REPORT_ZH.md`。最新 P1 若仍运行，只记录 UNKNOWN，不读取半成品来推定成功。

报告中的主要公式须按实际本构、复数实化、mass、source packing 实现，不能直接把三维 DDA 的极化率当 χ。

## 1. 先做零新增真实物理审计

运行随附：

```bash
python checks/check_finite_source_transport.py
```

需要 NumPy/SciPy。这里的线性方程是小型合成矩阵，不计为真实 Maxwell states，但单独记录 CPU 用时。不要把合成 G/H/O 正例写成真实天线采集成功。

建立 manifest：target_id、prior_id/law、noise_id、material decoder、gauge、coordinates/mass、source packing、active constraints、数据 split、缓存来源与预算。

区分：known TRAIN 几何控制；原 P0；条件 smooth-P0；radial proxy/full-Mie；新的 P1。禁止混合这些 target 的成功率与概率。

原 F2 最优性门槛未通过的历史状态仍为 NOT_RUN。理论上探索不要求 MAP 驻点；若本轮另做 feasible-anchor pilot，必须在新 protocol 中事前明确，仅要求真实 prior 与数据可行性，不声称这些点已通过旧 F2 或已到 MAP。P0 七维 smooth 仍只作控制实验。

## 2. 新预算与停止规则

本轮建议新硬上限：128 full material states，4096 complex RHS，GPU 训练 0，一般三维新 campaign 0。

E0/E3 至多 24 states/768 RHS；E1/E2 至多 80 states/2560 RHS；失败/独立审计预留 24 states/768 RHS。每个 source/adjoint RHS、重试、失败、basis solve、数值导数、full endpoint audit 都计入；Green matvec、分解、内存、wall time 另列。

若第一阶段发现真实后端成本不适合这个矩阵，只能在不读取质量结果的 microbenchmark 后整体缩小配对设计，不按方法表现单独增减配额。不得借用 P1 或旧 F 的余额。

没有可信 prior/noise law 或没有非平凡合法位移，E4 posterior 与 E5 learning 均保持 NOT_RUN。发布、push、提交远端和大训练不是默认任务。

## 3. 必须实现并独立核验的接口

对固定完整 anchor：

```text
DeltaX = U D U*
P = L0^{-1} U
H = U* G P
O = S P
B = U* E0_total
Q(D) = solve(I-D H, D)
T = Q(D) B
delta_current = P T
DeltaY = O T
```

低维 latent 不意味着 DeltaX 低秩。先确认 U 的真实物理支持；若原 prior 不允许局部材料块变化，不能把它放宽成 box 代替。

A17 增广必须用：

```text
P_V = R_V U
barC = (I-R_V L0) C
Gamma = C* L0 barC
N_U = C* (U-L0 P_V)
P_new-P_V = barC solve(Gamma,N_U)
```

对同一个有限 D，核验：

```text
DeltaY_new-DeltaY_V
 = (S + O_V Q_V U*G) barC solve(Gamma, N_U Q_new B)

rho_D = U D B - L_D (P_V Q_V B)
      = (U-L0 P_V) Q_V B
```

保存 full−new 的未解析余项。不要把已解析 exchange 当 full−old。

有可靠 eps_H、eps_O 时核验：

```text
theta = ||Q_V|| eps_H
bound = (eps_O + ||O_V Q_V|| eps_H) ||Q_V B|| / (1-theta)
```

只有 theta<1 且常数确实有效才能使用该充分界。浮点 SVD/经验 calibration 默认标 NUMERICAL_BOUND_ESTIMATE；没有严谨下界来源，不标 deterministic CERTIFIED。theta>=1 标 bound unavailable，不宣称 full Maxwell 必奇异。

## 4. E0：冻结 weak geometry

用相同材料坐标对照 full J、旧 ROM J、A17-exchanged J。Full J 是隔离 evaluator，不能免费泄露给新方法 selector。各方法必要的 full J 构造允许高效 receiver-adjoint/JVP 缓存复用，并按实际成本计费。

保存 weak-projector error、false-weak directions、方向响应区间、Fisher defect、finite secant error、剩余 defect 和谱隙。没有谱隙时用 cluster/soft filter，不比较单个向量的任意符号。

专门检查：GN 步相同却 weak geometry 不同；一次 exchange 很小但未解析 tail 很大；full 稳定而 Galerkin core 不稳定。

## 5. E3：材料可实现性，先于学习

实现单源标量 source chart：

```text
e = b + H t
d = t/e
Dd_Dt = diag(1/e) (I-DH)
logabsdet_real = 2 logabsdet_complex(I-DH) - 2 sum(log(abs(e)))
```

场为零、I-DH 奇异或本构逆不合法必须拒绝或单列退化情况，不能用隐藏 eps/伪逆/clip 把非法材料变成功。

多源共享材料检查：每一行 T_i 必须与 (B+HT)_i 共线；三维各向同性按 3×sources 整块共线，而不是各分量独立读材料。

实现反馈闭合联合核：

```text
T1(a)=diag(a)B
O_B(a)=vec(O T1(a))
gamma_i(a)=(H T1(a))_i B_i*/||B_i||^2
C_B(a)=(H T1(a))_i - gamma_i(a) B_i
E_B=ker(O_B) intersect ker(C_B)
d_i(a)=a_i/(1+gamma_i(a))
```

E_B={0} 仅说明该 straight-source 子类为空，不能证明没有其他弯曲 fiber。即使 a 是 prior tangent，d(a) 也可能离开 fixed-phase/shape prior，必须重新完整验证。

对一般 a，保存 R=HT1-diag(gamma)B，并核验精确关系：

```text
d(t)=t*a/(1+t*gamma)
DeltaY(t)=t O T1 + t O Q(d(t)) R
```

必须包含：单源成功但其他源失败；两独立源人工相容正控制；固定相位失败；实际小型 TM 多源测试。人工 H/O 例不算真实材料 world。

## 6. E1：同材料距离的五臂比较

只有明确原 prior 内的合格 anchor 才启动。优先两个独立 parent；只有同 parent 的多条件时明确标 n_parent=1，不能当独立统计样本。

预注册 2 anchors×1 common direction×2 material radii，共四个比较块。方向由共同的 prior/粗物理规则产生，不读取 full evaluator 或 GT；方向构造成本按冻结共享规则列账。若想研究各方法自己的方向选择，另作 end-to-end 组，不混入这个机制隔离实验。

五臂：
1. frozen material-J weak/null step；
2. refreshed full-J standard corrector；
3. ordinary prior/material-native predictor-corrector；
4. exact local T-matrix/Woodbury + 同一个有限材料 proposal；
5. A17 finite-source transport。

每臂每块最多三个中间新状态及一次独立 endpoint full audit。匹配材料 RMS 距离，不匹配 latent norm。所有方法用相同 prior、Phi_cap、progress、termination tolerances；无进度/no-op 不算成功。若 requested displacement 不可行，明确失败，不偷偷降低成功尺度。

若原 P0 保留既有有意义尺度，就不能把小于该尺度的记录升级为新结构；新 P1 的尺度须由 TRAIN/物理任务事前确定。

没有 A17 对 exact-port 的增量，报告 NO_A17_INCREMENT；不得只强调它胜过 frozen J。

## 7. E2：current action 的付费比较

在相同 E1 提案上比较 forward、GN、weak classification、finite-response/likelihood-decision selectors。候选 current pool 的生成规则相同，所有试算和未选中的候选也要计费。

优先使用 rho_U Q_V B 与有限双侧 return，而不是旧 r 的 GN gain。评分中的 upper bound 没有可靠来源时标 heuristic。记录同总成本下的 full finite error、可靠判定率、合法材料进度及拒绝次数。

不要假设 Galerkin enrichment 单调改善 Fisher、残差或 likelihood；真实反转必须保留。

## 8. Probability 和 NN 的硬边界

普通 physics corrector 只标 candidate discovery。没有正反向 proposal law / Jacobian / 可逆 augmented map，不叫 posterior。

明确 prior law 即使没有 density evaluator，也可检查是否有合法 prior-reversible kernel；不能直接把 implicit prior 判死刑。但一个无概率定义的候选库不是 prior。

P1 若已有可信目标且 E1/E3 非平凡成功，只输出下一阶段 E4 的独立预算与 reference 需求，不在当前额度中强行做高维 posterior reference。E5 不启动。

## 9. 必须交付

输出：

```text
START_HERE.md
PROTOCOL_FROZEN.md
TARGET_MANIFEST.json
E0_WEAK_GEOMETRY.csv
E1_MATCHED_TRANSPORT.csv
E2_ACTION_COST.csv
E3_REALIZABILITY.csv
CANDIDATES.jsonl + endpoint arrays
RESOURCE_LEDGER.jsonl
RESEARCH_REPORT_ZH.md
DECISION.md
REPRODUCE.md
```

每个未运行的部分写 NOT_RUN 和原因。报告 physics / exploration / probability / risk / cost 五张独立结果表。GT 只允许进入隔离 evaluator，不进入 anchor、方向、chart、prior、basis 或在线候选选择。

最后逐项回答：真实合法有限位移有没有；row-realizability 是否阻碍；exact-port 是否已足够；A17 是否修复 false-weak；新 selector 是否在相同总成本下更好；最关键剩余 tail 是否可控；是否值得继续，或应停止什么。

本轮不是论文包装任务。允许结论是：有限传输有用但 A17 无增量；多源/先验使机制不可用；或者简单 full-J/原 prior 方法更好。