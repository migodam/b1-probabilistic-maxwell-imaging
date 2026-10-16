# F3 独立概率与驻点审查

审查日期：2026-10-17。初始审查对象为既有 F1 端点与 `src/research_f`，以及父线程提出的 F3 设计。新运行源码与结果的复核在末尾单列。此次独立审查新增 Maxwell states、LU、RHS、Mie、GT 读取均为 **0**。没有导入项目物理解算模块；数值交叉核对只使用独立 NumPy/SciPy 代数与一维积分。旧文件、旧 gate、旧账本均未更改。

## 1. 共同尺度积分后的正确先验

令两权重独立服从 U[a,b]，a=0.4，b=0.9，L=log(b/a)。定义

- rho=log(w0/w1)，范围 [-L,L]；
- c=(log(w0)+log(w1))/2；
- w0=exp(c+rho/2)，w1=exp(c-rho/2)。

从 (c,rho) 到 (w0,w1) 的 Jacobian 绝对值为 exp(2c)。可行纤维为

`log(a)+|rho|/2 <= c <= log(b)-|rho|/2`。

因此对共同尺度纤维积分后，相对于 d rho 的密度恰为

`p_rho(rho) = [b² exp(-|rho|) - a² exp(|rho|)] / [2(b-a)²]`，当 |rho|<L；支撑外与端点密度为 0。

此密度非均匀。它在 rho=0 有尖点，在两个端点线性趋零。数值稳定形式为

`b² exp(-|rho|) * [-expm1(-2(L-|rho|))] / [2(b-a)²]`。

独立积分为 1.0000000000000002，误差估计 1.11e-14；p_rho(0)=1.3，两个端点恰为 0。见 `INDEPENDENT_PROBABILITY_ALGEBRA.json`。这些是概率代数核对，不是物理求解。

若七维计算坐标为 q=(p0,...,p5,t)，其中前六个 p 为原始均匀形状坐标，rho=(2t-1)L，则相对于 d q 的先验密度为

`p_q(q) = p_rho((2t-1)L) * 2L`。

`log(2L)` 只加一次。前六个均匀坐标在单位盒内密度为 1；若直接改用有物理单位的形状参数，则要使用相应常数 Jacobian。常数虽在单一 MH 比率中抵消，记录归一化密度时仍需保持明确。

共同尺度的积分不等同于在规范代表上代入原始二维权重密度，也不等同于把七维商坐标设成均匀。忽略纤维长度会改变先验。

## 2. 规范代表、材料体积与离散重复

规范代表可取

`w0*=sqrt(ab) exp(rho/2)`，`w1*=sqrt(ab) exp(-rho/2)`。

两权重均属于原闭支撑 [a,b]；rho 端点分别给出 (a,b) 或 (b,a)。原图像先把两个加权 Gaussian 相加，再除以该和的最大值。对两个权重同时乘以任一正数不会改变此归一化图像，因此规范代表保存材料图像。独立正数组代数控制的图像最大误差为 1.11e-16。

此结论依赖实际生成器正是共同最大值归一化。它不容许改变 Gaussian 中心、宽度或相对权重，也不容许在存在额外未同尺度缩放的 Gaussian 项时直接套用。

在 q 坐标中运行 MH，目标应为 `exp(-Phi(q))*p_q(q)`，其中 `Phi=||r||²/2`。**不得再乘或除材料 Gram 行列式**。若另行报告相对于七维材料 Hausdorff 体积的密度，在局部浸入且单射的片上才有 `p_material=p_q/sqrt(det(Dm.T@Dm))`；存在多个前像时还需求和。把此材料密度直接当成 q 坐标 MH 目标会重复计算体积。

除去连续共同尺度 gauge，并未自动除去两个 Gaussian 的标签交换。交换两组中心、宽度并令 rho→-rho 仍表示同一材料。保留有标签坐标、采用上述原始推前密度是合法的；若强行排序去掉离散重复，则须推导折叠后的密度。不能宣称当前七维坐标已是全球唯一的材料坐标。中心重合、形状退化或最大值单元并列时，局部秩与光滑性也要另行判断。

## 3. Likelihood 驻点与 posterior 密度是不同问题

既有 F1 的两个 frozen gate 失败仍保留：ring optimality=0.0452733958，weak=0.00200098872，门槛为 0.001。新 q 坐标与新判据只属于 F3；即使通过，也不能回写旧 gate 为通过。

对于 likelihood 目标 Phi=||r||²/2，g=A.T@r 为 q 坐标梯度。在单位闭盒的精确边界，KKT 符号为：

- q_i=0 时 g_i>=0；
- q_i=1 时 g_i<=0；
- 0<q_i<1 时 g_i=0。

正确符号的边界梯度可从投影梯度中置零；其余保留。**七个坐标均须检查**，特别不能只检查六个形状坐标而跳过 rho 的可行符号。仅当 rho 边界梯度方向正确时，较大原始 rho 梯度才不阻止可行驻点。

使用 1e-8 数值 active-set 时，必须明确：严格内部点的真实切锥仍然无约束，近边界置零属于预先声明的数值约定，不是精确 KKT 定理。建议同时报告原始梯度、精确边界投影梯度、该数值投影梯度、各坐标到边界的距离及 TRF bound-scaled optimality。

可用有界线性最小二乘生成候选：

`min_delta ||r + A delta||²/2`，约束 `-q <= delta <= 1-q`。

该步骤尊重原先验支撑；结合真实非线性 likelihood 回查和可行回溯，可作为继续优化的方法。若将距离边界小于预设阈值的候选直接吸附到精确边界，吸附点是新的材料状态，必须占用新 state 并重新得到真实 likelihood。不能把旧端点预测当作吸附后预测。

记录 `||Dm delta||`、线性预测 Phi 改善、实际 Phi 改善及 projected gradient，有助于区分剩余可行改善与不可行的权重比外推。有限盒内 GN 改善微小本身并不等价于切锥一阶驻点，也不证明全局最优或非线性邻域无改善。

新判据可作为预先冻结的有限精度 likelihood anchor qualification；不能因其与旧判据不同便声称修复了旧 F1 的科学失败。最终阈值与实际通过状态由父线程冻结协议和运行证据确定。

posterior 相对于 d q 的负 log 密度为 `Phi-log(p_rho)-log(2L)`。其驻点还包含先验梯度，与 likelihood 驻点不同。在 rho≠0 内部，

`d log(p_rho)/d rho = -sign(rho)*coth(L-|rho|)`。

rho=0 为尖点，端点 prior 密度为零、log prior=-infinity；不能直接用普通光滑梯度条件讨论这些点。Likelihood 最优点完全可能位于 posterior 密度为零的支撑端点。因此需要 likelihood anchor 的方法可受新驻点 gate 约束，但 finite-noise MH 的定义与有效性不依赖 MAP 或 likelihood 驻点成功。

## 4. MH 目标、边界与初始化

有限噪声 Gaussian likelihood 在支撑内不是硬 noise-ball 截断。正确 MH 接受 log 比率为

`min(0, -Phi_new+Phi_old + logp_q_new-logp_q_old + logK(old|new)-logK(new|old))`。

只有已经证明提议对称时，最后两项才能消去。noise-ball 可以作为输出诊断，不能据此提前拒绝高残差点并仍宣称目标保持原 posterior。

固定、对角高斯随机游走加逐坐标模 2 反射到 [0,1]，得到对称边界核。反射必须处理任意次数越界；简单的一次镜像或 clipping 不足以定义该核。clipping 会产生边界原子，不能继续使用纯连续目标的简单密度比。

**固定但相关的高斯协方差加逐坐标反射，一般不对称。** 反射混合核包含 `phi_Sigma(S*y+2k-x)`，其中 S 为各坐标独立符号翻转；只有增量密度也对这些翻转不变时，才能直接交换 x、y。独立二维反例给出 K(y|x)=15.9809646392、K(x|y)=20.2794357321，比值 0.7880379341；同方差对角控制的双向密度完全一致。若使用相关固定预条件提议，可对盒外提议直接拒绝形成 self-loop，或计算完整反射提议比率；不能省略修正。状态相关协方差即便无边界通常也需完整 proposal correction。固定于起点、以后不更新的矩阵不是状态相关矩阵，但相关反射问题依然存在。

旧 F1 实际端点严格位于原支撑内部，rho prior 虽很小仍为正，故可直接以有限 log target 开始 MH。其后进入内部可能使先验密度大幅上升；即使 likelihood 略变差，仍可合法接受。不能为保持旧 likelihood 最优值而增加非 posterior 的单调下降限制。

若另一个 arm 确实从精确 rho 端点开始，则 log target=-infinity。必须显式定义从零密度起点到有限密度候选的接受概率为 1，并避免 -infinity - (-infinity) 的 NaN。起点为零密度不使其成为 posterior 样本；需保留初始化与有限步非平稳性说明。更简单的独立 MH arm 可以使用旧有限密度 F1 起点，完全不依赖新 anchor gate。

被拒绝时必须保留重复链状态；accepted proposals、候选库、deterministic continuation endpoints 和完整 MH chain 是不同产物。种子、提议尺度、白化 sigma、目标 prior、接受次数、总提议次数以及越界处理应写入运行证据。

## 5. 允许与不允许的结论

本审查支持上述商先验、规范代表与 MH 公式；其实现是否符合公式须复核新源码。后续只运行极短链时，合理结论是“此核的目标为所声明 posterior，并实际产生这些有限步接受/拒绝”，而不是已达到平稳分布、混合成功、可信区间校准或发现了具有某 posterior 质量的新材料模式。

这些病例仍是同一公开暴露 parent 的两个条件背景；no-GT 运行路径不使其成为新盲测。保留原 sigma、原在线观测、原 prior 支撑、严格 state/RHS/wall 预算与 no-GT 路径是必要条件。数值响应满秩、noise-ball 相容、likelihood 驻点和 posterior 覆盖分别报告，不能互相替代。

## 6. 新源码与实际结果复核

### 6.1 运行前源码审查

已逐段阅读 `src/research_f3/prior.py`、`common.py`、`run.py`，未执行物理入口。

- `logprior` 使用正确的稳定 rho 密度，并加一次 log(2L)；规范 decoder 精确设置两个闭端点，避免浮点权重低于 0.4 一个 ulp 导致支撑验证失败。
- q 导数的形状链式系数 0.26、0.045 与 ratio 列的 ±L 正确；共同最大值归一化导数的商法则正确。
- `kkt` 检查全部七维梯度，保留数值 active-set 与精确边界两种诊断；线性 GN 使用原盒约束，继续优化仅接受真实 Phi 不增的候选。这个单调限制只用于 likelihood 优化，没有应用于 posterior MH。
- MH 的固定对角尺度仅在每个 arm 开始时计算一次；逐步不更新，反射实现允许多次越界。接受比使用完整 log prior 与 likelihood，没有 noise-ball 截断。实际链文件保留拒绝时的重复状态，并明确八步、没有去 burn-in、精度未评估。
- 初始 -infinity log target 到有限 log target 的接受分支显式处理；没有 inf-inf 运算。接受后正常使用有限比率。
- ledger 在付费动作前预留；新增阶段与旧 F 账本路径分离。此处为静态读取，不独立验证物理解算器和账本运行结果。

运行前已向父线程提出三项修正：anchor qualification 还须包含 noise-ball、材料/响应秩与物理残差；D 的差分 stencil 必须整段在合法盒内，不能以 clipping 代替固定方向二阶导数；MH 提议记录不得把被拒绝的候选标为 posterior sample。这些修正已在实际执行/恢复源码中复核：最终 qualification 包含原 noise-ball、材料秩7、响应秩7及物理 backward residual；D 在付费 Hessian 前检查合法差分 stencil，必要时缩小步长；提议记录明确 `posterior_samples=false`、`posterior_chain_proposal=true`。

### 6.2 实际运行

已读取最终 `RUN_FINAL_CLOSED.json`、追加式 ledger、48 条候选记录、两份八步 MH chain 及 F1.5 保存数组。独立重算见 `INDEPENDENT_F15_STATIONARITY.json` 和 `INDEPENDENT_MH_AND_LEDGER_AUDIT.json`。这些审查不调用物理解算，也不读取 LU/pivots、GT 或离线评价。

**F1.5 新判据的数值通过得到独立支持。** 两个最终点 t 恰为 0，所以数值 active-set 与精确边界的投影梯度一致；rho 梯度分别为正的 0.0119838821 / 0.1106716306，满足下边界可行符号。

| 新 quotient likelihood 诊断 | ring | weak |
|---|---:|---:|
| 新 states | 3 | 1 |
| 白化残差 | 22.2395547577 | 21.5121122553 |
| 精确边界 projected-gradient inf | 0.000579260410 | 0.000535627363 |
| 有界线性 GN 材料 RMS | 5.29818e-8 | 6.76656e-8 |
| 线性预测 Phi 改善 | 9.30783e-10 | 4.56623e-10 |
| 材料 / 响应 rank | 7 / 7 | 7 / 7 |
| backward residual | 1.65076e-15 | 1.52181e-15 |
| 该点 prior density | 0 | 0 |

重建梯度与报告逐分量完全一致。新 likelihood qualification 不改变旧 F1 Gate2 FAIL，也不构成 posterior MAP。因为两个点均在精确边界，此处的实际通过无需依赖 1e-8 近边界约定。

**MH 公式与实际有限链一致。** 独立使用 `1.44*sinh(2L*min(t,1-t))*2L` 重算先验，重放种子 26101771 的全部高斯增量和接受均匀数。16 次提议的 q 与保存记录完全一致；接受/拒绝、重复保留状态、Phi、log target、初始零密度扩展标志全部一致。log prior 最大差 8.88e-16，链 log target 最大差 2.84e-14。每例接受 3/8 次，五次拒绝均正确保留旧链状态。源码和运行均没有将 noise-ball 作为 posterior 硬截断。

| 八步有限 MH | ring | weak |
|---|---:|---:|
| 接受数 | 3 / 8 | 3 / 8 |
| 达任务 RMS 的提议数 | 1 / 8 | 2 / 8 |
| 达任务 RMS 的保留链状态数 | 0 / 8 | 0 / 8 |
| 保留链最大 RMS | 0.00468347 | 0.00475189 |

达到任务位移的提议全部被拒绝，不能当作 MH 链输出的任务级 teacher。两个 chain 都很短，从零密度 likelihood 边界开始且未去 burn-in；目标核正确不等于已平稳或已校准覆盖。chain 文件中逐步字段 `initial_boundary_zero_density` 实际记录更新后当前状态是否仍为零密度，故第一步为 false；真正的初始扩展事件由对应 proposal 的 `boundary_initialization_extension=true` 记录。解释时应以这两个时刻的区别为准。

**四臂比较没有支持 D 优于普通 native。** 本段 teacher 仅指保存的筛选条件同时满足：P0、正 prior density、原 noise-ball、相对冻结原 F1 的 deltaPhi<=2、任务位移及结构检查；不表示 posterior 质量已估计。

- ring 四臂均为 0/2 teacher；C_native 两点的 prior density 均为 0，D 第二点也为 0。
- weak 的 C_native 与 D 都只有第一个半径通过（各 1/2）：C 的 deltaPhi=1.140604，D=1.632730；第二半径为 3.098031 / 4.230149，均超过 2。因而首半径 C 更好，不能将 D 的存在性结果说成算法优势。
- weak mixed-native 池有 2/8 teacher，ring 为 0/8。池里包括两个 C_native 确定方向端点的复用及六个随机方向，不能称为八个独立同分布随机提议。
- B_updated 的第二点在两例都未达到指定第二半径（ring RMS 0.00540745、weak 0.00653120，目标 0.01055111）；ring A_frozen 第二点 RMS 0.01012286，也超出冻结 2% 匹配容差。这些点可以报告，但不能作为严格同位移的二级比较。C_native 与 D 的四个目标距离均匹配。

**RHS 对比须区分真实增量与 cold 等价。** ledger 每例 D 为 2 states / 56 RHS（16 forward +8 Maxwell Hessian current +32 一次 Jacobian refresh）。每例 C_native 为 2/16，A 为 2/16；B 因第一点与 A 共享缓存，实际增量为 1/40。mixed-native 池仅新增 6/48，另两个点复用 C_native，因此其八点 cold forward 成本为 64 RHS；前七点 cold 等价为 56 RHS，与 D 的 56 相同。实际 incremental 48 与 D56 不相等，不能混写。此比较同 anchor、同噪声，公共 anchor 几何对各方同等视为共享设置；如果计入公共设置，须向各方一致添加。前七 mixed-native 点 weak 仍有两个 teacher，D 有一个；同预算结论不支持 D 胜出。MH 前七次提议 cold 为 56 RHS，但 teacher 比较须用实际保留链，不能用被拒绝的 proposals。

**恢复与成本保留完整。** 首次运行在 ring C_native 第二半径因规范权重浮点值 0.39999999999999997 低于严格下界而失败；失败记录及 8 states /224 RHS 均保留。恢复只修复规范权重的浮点边界表达，并用 `gap=2L*min(t,1-t)` 防止极近边界 log prior 的减法消去；没有修改真实支持、q 提议、先验测度或已冻结筛选条件。恢复载入已付费状态/缓存，未重跑 F1.5，并将新增消耗追加同一 ledger。两个 JOB_END 的受控时间合计 16.303935709 秒。

最终累计 **46 states /640 RHS /16.303935709 秒**，F1.5 共 4 states，均在 64/2048/1200 秒与 F1.5<=16 的授权上限内。最终 closure error=null；初次失败没有被擦除。独立审查本身新增 states/LU/RHS/GT 均为 0。

### 6.3 最终审查结论

商先验推导、固定对角反射 MH 的目标核及保存有限链操作均通过本次独立审查；新 likelihood 边界驻点满足冻结数值判据。没有发现支持额外付费运行的必要性。实验证据不支持 D 的优势，不支持任务级 MH 链 teacher，不支持平稳 posterior 或覆盖结论。最终科学方向、novelty 与后续投入由父线程裁决。
