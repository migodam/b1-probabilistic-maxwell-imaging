# R3-DIGITS 计算规模与精确性独立审查

第一版：2026-10-29 19:38 UTC。审查者实际运行上下文为 `gpt-6.1-sol / high`，从本审查子会话的 `turn_context` 核验；没有启动 Astra 或其它子代理。只读公开算子、公开 DEV 观测、源码和相邻 R3 理论；未读取 `private_truth`、种子或标签。此文不是最终科学裁决。

## 1. 数学审查：条件下通过

固定有限 scalar-TM VIE，令 `D=diag(m)`、`M=I-DG`、`B=I-GD`。`M` 可逆时，`Q=B^-1 E0=E0+GJ`，`J=M^-1 D E0=DQ`。这里没有交换 D 与 G，也没有 Born 近似。

`MATHEMATICAL_SCALING.md` 的 D1 继承相邻 R3 `theory/ALTERNATIVE_STRUCTURE.md` 的 A8–A9。因为 `E Q=Y-F(m)`，硬数据约束下的最小本构残差范数恰为

\[
\psi=\operatorname{tr}[R(Q^*Q)^{-1}R^*],\quad R=F(m)-Y.
\]

需要 `rank(E0)=s` 和物理解可逆，因而 `rank(Q)=s`。这些是可逆域内的有限矩阵恒等式，不是唯一性或全局优化定理。

设 `A=R H^-1`、`H=Q*Q`、`C=A*A`。微分得到

\[
d\psi=2\Re\langle A,S\,dJ\rangle_F-2\Re\langle QC,G\,dJ\rangle_F,
\quad M\,dJ=\operatorname{diag}(dm)Q.
\]

所以 `M* Λ=S* A-G*(QC)`，`g_i=2 Re Σ_a conj(Λ_ia) Q_ia`。源码 `src/physics.py::Model.smallgram` 包含第二项，没有把变化的 Gram 冻结；LU 伴随用 `trans=2`，完整 Jacobian 的接收传播用普通转置 `trans=1`。这两个选择都正确。

## 2. 薄 QR 放置与实际实现

旧 R3 每次分解 `A.T`，尺寸 N×(N−s)。本轮主路径先解完整物理系统，再形成 `H=Q*Q`，对 s×s Gram 求特征值诊断、Cholesky 与三角回代。出现 `cond(H)>1e10`、非有限条件或相对解残差>1e-8 时，源码回退到 **Q 的 N×s 薄 QR**，无对角 jitter。回退费用另记 `SMALL_GRAM_THIN_QR_FALLBACK`。

若 `Q=UT`，`H=T*T`，右解 `W=R T^-1`、`A=W T^-*` 与源码两次三角求解一致。薄 QR 放在总场矩阵上，而不是材料空间或原 Green 矩阵上。它不免除 N 阶 LU，也不把物理状态降到 s 维。Gram 会平方条件数，回退仍可能因严重源失秩而失败；失败必须原样保留。

审查建议：回退目标可直接用 `||W||²` 避免 `Re <R,A>` 的数值抵消；当前实际公开点两者相符，未据此发现阻断错误。小的线性方程残差本身也不能证明 resolvent 条件良好。

## 3. 实际成本，不能只比较 RHS

|一次新材料处的评价|完整 LU|源 RHS|额外传播 RHS|其它实际费用|
|---|---:|---:|---:|---|
|普通数据目标+伴随梯度（后增归因对照）|1|12|12 伴随|原数据收缩与共同正则项|
|GlobalDual 目标+精确梯度|1|12|12 伴随|576×12 总场、12×12 Gram/诊断/Cholesky；必要时薄 QR|
|Full-J/GN|1|12|48 接收转置|1152×576 实 Jacobian、与阶段基相乘、法方程、64/144/576 维 LM 解、试探 forward|

主耗时包含 O(N³) LU 与 O(N²s) 传播；Gram 形成约 O(Ns²)。小 Gram 相对旧宽 QR 的改变是真实的，但不构成优于普通伴随梯度的成本定理。Full-J 的 48 个接收 RHS 已复用同一 LU，没有故意使用 576×12 条逐材料 tangent 或有限差分。全变量 GN 的材料线性代数必须与物理传播一起记账。

`Run.started` 在 `Model` 构造前设置，包含模型准备、多尺度基、优化、失败和归档的数值 wall；`load_public`、文件读取、threadpool 管理和最终 JSON 写盘在此计时外。最终报告必须把它称为 **solver 数值预算**，另列外层进程 wall（若记录），DEV的75秒和主实验授权的90秒都不能声称为完整端到端时限。操作计时未覆盖全部 Python 开销，求和不应冒充总 wall。一次原子操作可越过截止，overshoot 必须报告。

BLAS 通过 `threadpool_limits(1)` 约束，CLI 也设置线程环境；结果保存 threadpool 信息与进程 CPU。两并发 case 不等于免费资源：需保留并发数、硬件、每方法 CPU 与完整 wall。

## 4. 多尺度、倒数坐标与归一化

三方法都在同一个576单元物理系统上工作。正双线性基 `P_h` 的行和为1，阶段为8²、12²、24²；24阶段明确检查 `P_24=I_576`。native 为 `m=P_h c`；倒数方法优化 `x=1/c`，仍是 `m=P_h(1/x)`，所以两个方法在每阶段允许的物理材料图像集合完全相同。它没有采用会改变粗阶段图像集合的 `m=1/(P_h x)`。

物理图像正则项是

\[
\mathcal R(m)=\lambda_s\operatorname{mean}[(Dm)^2]+
\lambda_{TV}\operatorname{mean}\{\sqrt{(Dm)^2+\epsilon^2}-\epsilon\}.
\]

其中差分除以19.5，三臂使用同一权重、epsilon、box 和常数初值1。GN 使用共同原数据目标 `.5||F-Y||²/||Y||²+R`。GlobalDual 提案采用 `a*psi/(2b)+R`，`b=||Y||²`，`a=||Q_ref||²/s`。V1封存源码只在 `dual_factor is None` 时设置该因子，此后跨阶段固定；V2明确从第一次已付费的共同粗阶段均匀材料1前向取得Q_ref，且fine阶段不重设。故完整梯度乘同一个固定常数，没有缺失归一化导数。它不是原likelihood：平均总场能量缩放不能消除各源方向的Gram度量。最终候选始终回到共同原目标排名，故不能把提案目标下降单独算成成像收益。

`all_unknowns_released` 当前只检查24阶段存在已评价候选。最终审查还必须核对该阶段实际 Full-J/伴随调用及迭代记录；仅评价一个粗阶段延拓初值不足以证明实际进行过576变量搜索。返回 top1 可以来自早期阶段，但必须明确其来源，不能把变量维数与已恢复细节混为一谈。

## 5. 本审查新增的数值检查费用（不可并入算法成绩）

只对公开 DEV 文件 `dev_aeab4b38a05d.npz`、公开算子和常数材料1进行实现检查。未优化材料，未读取 GT，未使用或生成随机种子。当前 DEV 噪声身份的更正见 `EXPERIMENT_NOISE_CORRECTION.json` 及后续 generation receipt，不能把本次检查称为新鲜留出实验。

- `cond(Q)=25.0955944447`，`cond(H)=629.788860532`。
- Cholesky 目标0.20810266224889082；独立薄 QR 目标0.20810266224889076；相对差2.6675e-16。
- 一个确定性方向，步长1e-3、1e-4、1e-5的中心差分相对误差分别为4.60e-9、1.70e-8、6.86e-7；解析方向导数−3.7950631199e-5。有限差分只用于审查，求解器不使用。
- 新增7次576阶 LU，84源 RHS、12伴随 RHS，合计96 complex RHS；7次12阶 Gram 评价，1次576×12薄 QR比较，`np.linalg.cond` 另做一次 Q SVD 与一次 H SVD。
- 工具记录该 Python 命令 envelope wall=0.756445625秒；逐操作 wall 与 CPU 未输出，状态为 `NOT_RECORDED`，不得补推。另一次默认 Python 因缺 NumPy 在导入处失败，没有物理求解。

这是同一个离散模型的代数/实现检查，不是独立物理验证或性能样本。父线程已要求停止额外数值身份检查，审查不再增加该类调用。

## 6. 当前结论

`MATH_PASS_WITH_CONDITIONS`：D1–D5与当前常用路径正确，源 Gram 条件保护缺口已在本审查期间补入。`PERFORMANCE_HOLD`：尚未审阅冻结九数字结果、外层预算收据及完整576阶段记录。强新数学、全局性、连续PDE、唯一性与普通梯度普遍加速均无证据。


## 7. 19:48 UTC 静态复核与审查更正

此前第一版把开发过程早期未封存草稿的 `psi/psi_initial` 描述误套到已评分V1，并写成阶段重置；这不符合封存V1源码，现已在本文件、最终科学审查和公式表全部更正。依据是 `SOLVER_DEV_V1_SOURCE/src/solver.py` 的 `if self.dual_factor is None` 与当前 `Run.score/dual`；实际V1为一次固定场能量因子，V2为明确固定均匀材料1参考。此处承认审查文字错误，不改变历史预测或实验。

父线程批准的V2开发调整把8、12粗阶段改为三臂同一native Full-J，各阶段最多320次迭代、按共同wall安排停止，各臂重新运行并自行支付相同路径费用，24阶段才按方法分流。直接读取 `experiments/frozen_predictions/dev_v2` 六份JSON而不读取GT：

|公开DEV ID|三臂fine初值逐元素完全相同|共同粗阶段已付费RHS / Full-J|fine阶段完成GN迭代 native / reciprocal / dual-polish|GlobalDual fine伴随评价|
|---|---|---|---|---:|
|dev_2a5c221387b1|是，576元素|43524 / 640|514 / 525 / 334|979|
|dev_703ae1b0abb1|是，576元素|43536 / 640|661 / 620 / 263|849|

六个top1均来自24阶段。上述证明DEVv2确实进行了576变量的导数/迭代，不能替代主九例的同样检查。DEV是已知开发数据，二次DEV不算新鲜盲确认。

对 `SOLVER_DEV_V2_SOURCE` 与当前源码逐行差分的静态审阅：最终恒等基绕过乘I，由 `run` 内P24=I检查保证资格；精确稀疏差分组装 `Dcsr.T @ diag(weights) @ Dcsr` 与原IRLS二次模型一致，保留每次新权重及倒数链式因子；Gram_RHS记48个小矩阵右端，不能混入物理RHS；L-BFGS接受callback只增事件与计数，不求物理。它们不改变实数代数公式，但浮点加法顺序、计账开销和wall限下的迭代吞吐会改变，所以不是DEVv2数值预测的逐位重放。没有新增数值验证，按父线程要求只作静态审查。

主实验freeze已记录于 `HOLDOUT_FREEZE.json`，时间为19:47:38.116883 UTC；`configs/main_methods.json` 明确三臂各90秒、两worker、每worker单线程BLAS、相同先验/box/初始化及320次粗阶段cap。逐个比较 `HOLDOUT_SOURCE_SNAPSHOT` 的16文件与当前原文件字节，16/16完全相同，未计算哈希。父线程确认freeze时holdout公开输入数为0，随后授权生成9例；当前公开generation receipt确认9例、零噪声、9LU/108源RHS。相同粗算法只在未因wall截断产生分歧时保证相同fine初值；必须在全部主例检查实际数组，不能只凭源码宣称完全相同。

最终静态结论：`SOURCE_FREEZE_STATIC_PASS`。已审读冻结版runner、batch和scorer：原始top1与oracle分离，scorer从嵌套top1取残差；独立验算用同模型总场LU，ZGECON估计器另记时间且内部RHS标记opaque，不能冒称零额外费用。主实验结果仍为 `PERFORMANCE_HOLD`，无新增物理运行。


## 8. 主实验后、全部GT评分前的两个补充：静态审查

`budget_fill.py`/`configs/BUDGET_FILL_CONTRACT.json` 于19:58:09.572410 UTC另行冻结；`data_adjoint_control.py`/`configs/DATA_ADJOINT_CONTROL.json` 于20:02:06.447952 UTC另行冻结。它们都在原九例Y生成后设计，父线程要求在全部63份预测snapshot前不开主GT评分。它们不是生成前预注册的新鲜方法确认，也不是新数据。各补充六个快照文件与当前源码6/6逐字节相同；原主实验16文件仍16/16不变。只做源码审查，没有新增PDE/梯度数值检查。

**Budget fill：** 每臂从原在线top1继续最后576维同一目标的GN，倒数臂仍为倒数GN，其余为nativeGN。新预算 `max(0,90-original_wall)`，不是再给90秒；加载、Run准备和新LU都包含在补充计时内，保留原材料，没有复用不存在的LU。1e9迭代cap只是避免原1000次cap提前结束，正常收敛、无接受步或wall仍可停止。原计数加一次增量生成root累计计数；总研究费用只能加 `incremental_accounting`，不能把累计root再加一次。合并原候选后按共同原物理目标去重/排序，故保留的原incumbent保证top1目标不增；新计算中即使随后异常，之前已完成forward、合法且更好的付费候选仍可能入档。需要同时报告 `budget_fill_status/error`，根status继承原值不能单独代表补充成功。

**普通data-adjoint对照：** 复用同一共同粗路径、fine初值权限、45% L-BFGS-B时间安排、选项、box和正则，其提案目标就是原评分。因为数据项为半平方范数，梯度确为

\[
M^*\Lambda=S^*R,\qquad
g_i=\Re\sum_s\overline{\Lambda_{is}}Q_{is}/\|Y\|_F^2+\nabla_i\mathcal R.
\]

实现没有多余系数2。每次新状态1LU+12源RHS，梯度复用LU加12伴随RHS，Gram调用应为0。`DataAdjointRun.dual` 重写提案而继承公共调度；收尾只改输出身份，内部历史 `dual_fullj_polish` 名称已明确是通用尾段标签。

该对照有助分辨优化器组合与Gram目标效应，但不是只有目标一项的逐位匹配：ordinary对照GN cap为1e9，原主实验为1000；应优先与budget-filled GlobalDual比较，而fill仍有从incumbent重启/再次LU的不同。对照90秒包括load，原90秒数值预算不包括load，存在轻微不对称；wall受并发、回调开销和迭代吞吐影响。必须报告这些差别，不能将单次胜负当严格因果识别。实际fine初值/计数和全九例成绩尚待最终审阅。


## 9. 父线程提供的576规模表示内核计时

`SCALING_KERNEL_BENCHMARK.json`记录一个公开均匀材料1状态、每种表示一次调用，机器同时运行两路成像worker。本审查只读该收据，没有重跑。旧表示固定setup wall0.0488524秒/CPU0.033611秒，逐目标setup wall0.00440671秒/CPU0.004271秒，均单列；旧576×564 QR目标/梯度内核wall0.0668536秒/CPU0.050002秒，48个LS小矩阵RHS。新完整物理求解＋小Gram精确梯度wall0.0247559秒/CPU0.016秒，费用明确为1次576阶LU、12源+12伴随=24物理RHS，另有12阶Gram及48小矩阵RHS。目标分别为0.1770092634534417和0.17700926345344176。

该次wall比为2.7005，直接支持“本公开状态下，新表示单次内核比旧宽QR更便宜”的有限观察。它不是重复计时统计、成像比较或相对普通数据伴随/GN的普遍加速证据；旧表示setup也不能不加说明地摊进或摊出所有目标。此额外调用费用应在总研究台账单列，不能塞进任何盲方法预算。


## 10. 核心63结果后的计算审查更新

原27、补齐27、ordinary9已依序在20:04:41、20:08:25、20:13:14 UTC封存并随后评分。原三臂逐case完全相同的fine576初值已由coverage核实；ordinary对native九例也全部精确相同。原主27例均有实际576维导数和接受步；ordinary9例有2–1005次data-adjoint，Gram=0，六例GN无接受步如实保留。63端点全有限/box合法，验算差≤3.33e-16；不是仅凭 `all_unknowns_released` 字段作结论。

读冻结预测累计计数：补齐native wall739.0504秒/CPU645.5355秒，25,787 LU/1,054,932物理RHS；reciprocal739.0432/649.7174秒，26,255 LU/1,075,716 RHS；Gram761.2932/673.9434秒，28,578 LU/1,041,180 RHS及7,539 Gram；ordinary448.8433/344.2859秒，15,098 LU/533,220 RHS。fill新增wall分别138.7868、164.0241、50.1672秒，不能重复加累计root。所有fill error=null；最大累计overshoot0.16662秒。

所以计算规模正确、旧宽QR到小Gram的单次内核收益存在，但实际成像并没有普遍wall优势；Gram的LU更多、总wall更高。核心视觉/材料收益集中于数字9，并非计算少就必定恢复。完整科学结论见 `FINAL_SCIENTIFIC_REVIEW.md` 第9节。噪声扩展及完整总成本仍待封账。

最终封账状态更新：`SCIENCE_FREEZE.json`于20:41:43 UTC关闭新增物理。噪声仅14/27完成、13NOT_RUN，禁止部分均值；三TRF、单方向Jv和两低幅度控制的最后解释见最终科学审查§12。它们不改变本文件核验的代数与预算事实，也不建立成熟普通法支配或原高对比任务全恢复。全部成本由最终cost台账更新，本审查没有新增物理。
