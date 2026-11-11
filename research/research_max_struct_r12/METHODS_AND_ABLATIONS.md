# R12 方法、消融及证据边界

本文件整理已完成的双场投影核心、四对象 NN 接口开发与高对比度 true3D 补充。独立24确认已完成192次执行与整组封存，并已完成离线评分；科学解释见总报告与审查；本文件没有运行新实验。数学总稿见 [THEORY_DERIVATIONS.md](THEORY_DERIVATIONS.md)，最终复核问题见 [audit/FINAL_THEORY_REVIEW_NOTES.md](audit/FINAL_THEORY_REVIEW_NOTES.md)。

## 1. 状态、尺度与材料共享

继承有限离散物理关系

\[
J=D_mQ,\qquad Q=E^{inc}+GJ,\qquad F(m)=SJ(m),
\]

其中材料是实各向同性标量。标量 TM 中每个 cell 的材料共享所有照明；true3D 中还共享该 cell 的三个矢量分量。核心要求 cell-major 行序：`cell0-x,y,z,cell1-x,y,z,...`，不对每个分量或每个照明独立拟合材料。

取固定尺度 \(s>0\)，令 \(X=J/s\)、\(\alpha=m/s\)，使用度量

\[
\|\Delta X\|_F^2+\|\Delta Q\|_F^2
=s^{-2}\|\Delta J\|_F^2+\|\Delta Q\|_F^2.
\]

改变 \(s\) 会改变投影和轨迹。完成的 NN/high3D 开发均冻结 \(s=6\)，没有按测试材料反调尺度。

两个集合为

\[
\mathcal C=\{(X,Q):X_i=\alpha_iQ_i,\ m_i\in[\ell_i,u_i]\},
\]

\[
\mathcal A_\delta=\{(X,Q):Q=E^{inc}+sGX,
\ \sum_j\|W_j(sS_jX_j-Y_j)\|_2^2\le\delta^2\}.
\]

\(S_j\) 可包含合法的按源 receiver mask，\(W_j\) 是按源白化块。当前核心不实现跨源相关白化。所有源共享一个全球数据球和一个乘子，不能把它改成逐源独立球。

## 2. 精确有限维共享实材料投影

对一个 cell，将所有源和矢量分量排列为复向量 \(x,q\)。固定实斜率 \(\alpha\) 后，最近场对为

\[
\widehat q=\frac{q+\alpha x}{1+\alpha^2},\qquad
\widehat x=\alpha\widehat q.
\]

令 \(a=\|x\|^2\)、\(b=\|q\|^2\)、\(c=\operatorname{Re}(x^*q)\)。消去场后的距离为

\[
\varphi(\alpha)=\frac{a-2c\alpha+b\alpha^2}{1+\alpha^2},
\qquad
c\alpha^2+(b-a)\alpha-c=0.
\]

枚举 box 内全部合法实根和两端点，比较真实平方距离，即得到每个 cell 的全局最近场对。代码使用稳定二次根、cell 内同尺度归一和平方范数比较，避免根消减及 Gram 差值消减；零交叉项、暗 cell、端点最优与并列最优均保留。并列时确定性选择第一最小值，不宣称材料标签唯一。

这是实共享斜率 TLS 的有限维构造；不是独立源拟合，也不是复介电或各向异性张量投影。原非局部材料 prior 没有包含在 \(P_\mathcal C\) 中。代价为 \(O(NcP)\)，其中 \(c=1\) 或3，\(P\) 为照明数。

实现：[src/field_current.py](src/field_current.py) 的 `material`。原始 Pro 第4节及 toy 保存在 [sources/pro/MAX_STRUCT_RESEARCH/](sources/pro/MAX_STRUCT_RESEARCH/)，未改写。

## 3. 固定背景传播—数据投影

设 \(B=sG\)、\(H=I+B^*B\)。无数据球约束时

\[
\overline X=H^{-1}\{X+B^*(Q-E^{inc})\}.
\]

逐源定义

\[
T_j=sW_jS_j,\quad y_j=W_jY_j,\quad
V_j=H^{-1}T_j^*,\quad K_j=T_jV_j.
\]

若 \(\sum_j\|T_j\overline X_j-y_j\|^2\le\delta^2\)，返回无约束投影。否则使用共同 \(\lambda\ge0\)：

\[
X_{\lambda,j}=\overline X_j-
V_j\lambda(I+\lambda K_j)^{-1}(T_j\overline X_j-y_j).
\]

利用缓存的 \(K_j\) Hermitian 特征分解，对全球残差平方求一个单调标量根。\(H\) 与材料迭代无关，只分解一次。相同的 \(T_j\) 复用 \(V_j\) 和 port 特征系统。

零特征空间残差产生可达下界；球小于该下界时抛出 `UnreachableDataBall`，落在边界时使用伪逆极限，不伪称总存在有限正乘子。代码的数值 port rank 阈值为 `32*eps*max(dofs,receivers)*max_abs_eigenvalue`，边界能量容差为 `128*eps*energy_scale`。因此“精确投影”描述有限维数学问题；浮点实现带有明示的 rank/边界容差，不能等同无误差的符号计算。

合成CN噪声合同使用 raw iid circular complex 噪声；Fresnel与Austria的模型差异/固定范数噪声球另见BENCHMARKS，不能套用以下概率资格。若 \(E|n|^2=\sigma^2\)，全局99%半径为

\[
\delta=\sigma\sqrt{\chi^2_{2n_{obs}}(.99)/2}.
\]

任意白化 API 的调用者必须提供与其噪声协方差和 mask 相容的 \(\delta\)；核心不自动把 raw 半径迁移到任意 \(W_j\)。

## 4. AP、DR 与物理候选验收

令 \(c_k=P_\mathcal C z_k\)。AP 为

\[
a_k=P_\mathcal A c_k,\qquad z_{k+1}=a_k.
\]

DR 为

\[
a_k=P_\mathcal A(2c_k-z_k),\qquad
z_{k+1}=z_k+\omega(a_k-c_k),\quad0<\omega<2.
\]

完成研究默认 \(\omega=1\)。同一完整物理初态满足 \(P_\mathcal Cz_0=z_0\)，故 AP/DR 的第一次更新同为 \(P_\mathcal Az_0\)。这是场状态恒等式；暗 cell 的材料读出标签可不唯一。之后才可能出现反射记忆差异。球内错误物理解同样是固定场状态；没有全局逃逸或材料恢复保证。

每50步将 C shadow 的材料交给原物理 backend，重新求完整物理场和原 \(\Phi=D+R\)。候选只按原 \(\Phi\) 选择，初态保留为 incumbent。核心不读取 truth。投影审计保存轨迹前/后状态、两个 shadows、物理候选/最佳场、材料、计数和 CPU/wall 时间。iteration0另作baseline记录：其 `gap=0` 和 A字段是初始状态占位，尚未做 A投影，不能当作consensus或应用shadow-gap bound。预算退出不强加额外免费正演。

AP/DR 的 auxiliary 数据球不意味着每个材料 shadow 的完整物理预测都在球内；gap 与材料 resolvent 会影响该差别。随后的原物理 corrector 优化 soft \(D+R\)，也不保证维持球内可行性。详见 theory 中的 conditional gap bound 和本轮 high3D 保存状态分解。

## 5. 固定背景成本及 grouped v2

令 field 维度 \(d=Nc\)，相同 port 块去重后的接收维度为 \(p_r\)：

- setup：形成 \(G^*G\)、一个 \(d\times d\) H Cholesky、每个不同 port 的 \(p_r\) 个 H RHS、port Gram 与 eigensystem；dense 主 setup 为立方级，存储平方级。
- 每次 A：\(P\) 个固定 H RHS、G/G-adjoint 的多 RHS action、接收及缓存的 H-inverse 接收伴随 action、标量乘子计算。
- 每次 C：跨源/分量的逐 cell TLS，不做材料相关 Maxwell LU。
- 物理审计/corrector：原 backend 的材料相关 LU、forward/adjoint RHS、原 prior 和优化费用，全部另计。
- 完整 audit archive 的压缩/写出是实际 CPU/wall 费用，不能只报61个物理 LU 来代表3000步投影成本。

`H_*_rhs` 以 RHS 列计数，一个 Cholesky solve 内含两次 triangular solve；不是“每列一个新分解”。`receiver_adjoint_action_*` 表示缓存的 \(H^{-1}T^*\) 应用，不是额外 H 求解。`ProjectionRun.costs` 是 runner 增量；audit 的 `costs` 和 projector 的累计 ledger 包含 setup。不要相加嵌套 timer/counter。

`grouped_ports=False` 保留初始 v1。可选 v2 (`True`) 只将完全相同 port 的源列打包为矩阵 RHS，并聚合每个特征空间的 residual power；求解的投影问题、尺度、数据球和目标不变。聚合顺序会改变浮点舍入，不保证长非凸轨迹 bitwise 一致。最初已启动的 basin/低对比 vector v1 结果保留；NN/high3D 和确认冻结使用 v2。

一个已完成的 Pro 100-dof、6源、24接收、200次固定 A 投影 workload：v1 **0.085336CPU s**，v2 **0.035540CPU s**，该 workload 比值 **2.40×**；两场输出相对差异低于 \(10^{-15}\)，H RHS均1200。它不是重建速度、CUDA速度或材料恢复结论。原始计数见 [src/projection_grouping_benchmark.json](src/projection_grouping_benchmark.json)。

## 6. 已完成的比较及其隔离范围

|比较|保持相同的条件|可回答的问题|不能单独推出的结论|
|---|---|---|---|
|同完整物理初态 AP vs DR|A/C、尺度、审计、物理目标及预算口径|第一次相同更新之后的反射轨迹差异|全局收敛、所有局部驻点消失|
|NN-direct vs NN-complete-physical→DR|冻结 NN、原数据、box、prior、共同 corrector|四个暴露 R8 对象的 NN 接口保留/损伤|没有 NN→AP 臂，不能将全部差异唯一归于反射|
|NN-CIE vs NN-DR|同 NN 材料、原物理选择及共同 corrector|不同 current/field 接口对信息的影响|current 表示、辅助目标与 proposal 成本也不同，不是单因素反射消融|
|high3D CIE/Dual/AP/DR|相同真3D forward、初态、原 prior、box、共同 corrector|两个 matched 对象及一个配对 cross-grid 上的有限表现|配对不是第三独立对象；不同 proposal 不是等支出|
|raw→同 \(\Phi\) corrector|同一已付费 raw 候选、原 data/prior|材料变差来自 data/prior 权衡还是两项同时改善|更低 \(\Phi\) 不构成 GT 或唯一性证书|
|grouped v1 vs v2 固定 A workload|同一数学投影问题和输入|矩阵多 RHS/聚合的工程吞吐|独立科学机制或普适速度优势|

NN 已完成 **4对象×4臂=16预测**，是公开/已暴露开发回归。high3D 已完成 **2独立几何+1配对×4臂=12预测**。各 cohort 的 seal 均先于该 cohort 的离线 truth 读取。数值及成本见 [NEURAL_AND_SELECTION_REPORT_ZH.md](NEURAL_AND_SELECTION_REPORT_ZH.md)。

## 7. Primary 与独立报告的 secondary estimator

Primary 仍是已付费完整物理候选池中的最低原 \(\Phi\)。新增 secondary 为

\[
\arg\min_{m\in\mathcal P_{paid}}(R(m),\Phi(m))
\quad\text{s.t. legal box and }\|F(m)-Y\|\le\delta_{CN99}.
\]

这里是词典序 tie break；池无合格候选时明确 unavailable，不用 primary 顶替。它不向迭代反馈、不增加 forward、不用 GT，只在已有完整物理审计池内筛选，附加检查/复制/写出费用仍计入预算。不同方法的审计池覆盖率不同，必须分别披露。

来源是 high3D 暴露开发负结果，不是独立证据。实际注册文件 [delegated/confirmation/SECONDARY_SELECTION_FREEZE.json](delegated/confirmation/SECONDARY_SELECTION_FREEZE.json) 记录 **24观察已生成、确认逆解启动0、GT评分0** 时冻结此 secondary；不能写成它在数据生成前就已冻结。该文件规定八方法共用规则、whole-child240CPU、proposal144/corrector96及30/60/120/240前缀。

独立24的 primary/core 在 [ROOT_CORE_FREEZE.json](delegated/confirmation/ROOT_CORE_FREEZE.json) 中更早冻结。新增 secondary 的实现/注册不等于确认成功；在完整 cohort 封存、独立评分、失败计入和最终评审前，不能宣称已验证，也不能把潜在 secondary 收益记为纯 DR 反射机制收益。
