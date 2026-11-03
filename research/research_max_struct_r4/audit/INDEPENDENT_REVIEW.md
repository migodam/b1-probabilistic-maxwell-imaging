# R4 独立批判审查

日期：2026-11-03。性质：静态数学／实现审查，供父线程最终科学判断。

本审查读取指定理论、A/B 实现、开发审计、prior-art 摘录和冻结 confirmation 合同，并为核对实际调用而读取既有 `cooperative_paths.refine`、R2 global CIE 的源指标／读出调用、R3 native readout 和 `MaxwellModel` 薄适配器。没有启动求解器、梯度检查、恒等式检查、计时或评分；没有读取 confirmation GT、私有映射或 confirmation 结果。下述开发数字仅复述 `DEVELOPMENT_AUDIT.md`，没有重新验证原始预测。

## 结论

未发现 A 的 exact-block／onehop 实材料梯度或 B 的固定 eta profile 梯度存在致命代数或复伴随接线错误。A 的同零集和 Neumann 界是带明确条件的有限矩阵结论；B 的 profile metric 是精确消元结果。二者都不能直接提供成像优势、普遍盆地扩大、收敛或新物理信息。

当前最需收紧的是机制归因：`THEORY_DERIVATIONS.md:128–133` 的残差展开式正确，但“actual mechanism”若理解为实验已经单独确认有符号交叉项，则证据过强。self 对照同时删除 `2 Re<u,Lu>` 和 `||Lu||²`，不是仅删除交叉项。现有证据支持“该公开开发对象上，加入块内非对角缺陷输运的完整算法成功，而指定 self／global／native 端点失败”；尚未辨明线性交叉项、正二次项及其优化轨迹各自是否必要。

## A：可逆界与零集

1. 令 `D=I+X(beta I-diag G)`、`L=D^-1 X Goff`、`u=beta D^-1 r`，则 `K=D(I-L)`，故 `R_exact=(I-L)^-1 u`、`R_one=(I+L)u`。只要所有 `d_i!=0` 且在所声明的诱导范数中 `||L||<1`，`I-L` 和 `I+L` 都可逆。于是两种状态罚函数的零点均恰为原离散状态方程 `r=0`；这不意味着非零残差处的目标、驻点或轨迹相同。

2. `R_exact-R_one=(I-L)^-1 L²u` 及 `||error|| <= q²/(1-q)||u||` 正确。代码 `cie_onehop.py:70–98` 对正实盒内 `m/|1+bm|` 的逐行极值计算正确：导数符号由 `1+Re(b)m` 决定，检查端点及合法临界点。`min |d_i|` 的投影二次极值也正确。该函数不读取 GT。

3. 文中开发网格的 `0.16810 / 0.03784` 是现有记录的 infinity-norm 数值，本审查未重算。对应结论只对产生这些界的实际 G、分块、盒和 beta 生效，不能把 infinity-norm 误差因子作为 Euclidean／源指标加权 Frobenius／H(curl) 同因子。零集可逆性本身不依赖这种范数转换。

4. `whole_box_bound` 是一个显式诊断函数，`run_onehop_cie` 没有自动调用或拒绝不满足界的通用输入。冻结 confirmation 沿用指定公共算子与 beta=.2/1，因而可引用同一算子的既有界；不能把该 Python API 描述成对任意 G 自动提供可逆保证。换网格／材料盒／块／beta 后原数值界不自动沿用。

5. exact-block 的另一个充分条件 `Im_H Gin > 0` 推理正确：`K=X(X^-1+beta I-Gin)`，对正实可逆 X，后一个矩阵的虚 Hermitian 部严格负定，不能有零向量核。半正定情形不足以得到该严格论证。该条件在文中是条件性定理，不是本审查已经测得每块严格正定的事实。

## A：复杂伴随与实材料梯度

`cie_block_absorption.py:137–150` 与推导一致。对固定 J、Q、固定 Hermitian 源指标 Omega，有

`dR=-K^-1 diag(dm)(B R+beta Q)`，

所以以 `K^-H(R Omega)` 作伴随，取逐行复内积实部，得到实材料梯度。代码确实包括 K 对材料的导数，没有把左预处理器错误冻结。

`cie_onehop.py:101–120` 同样正确：

- `du_i=-beta(Q_i+b_i J_i)/d_i²`，`da_i=1/d_i²`；
- `W=R Omega`；对 u 的伴随为 `W+Goff^H diag(conj(a))W`；
- 对 a 的额外项为 `Re sum_s conj(W_is)(Goff u)_is/d_i²`。

实现中的 `off_action(...,adjoint=True)` 同时共轭并交换块矩阵两轴，正确使用 Hermitian 伴随。材料导数保留复平方 `d_i²`，不是误用 `|d_i|²`。这些梯度仅是固定 J、Q 的材料内层梯度；它们不是消去／沿着实际有限次 current 更新轨迹后的总导数。

current 更新对固定 m 使用同一个线性残差变换作用于 H/T，随后使用 `H^H(R Omega)`；与所声明的二次 current 子问题一致。两 beta 阶段中 current 使用 beta=1，材料使用 beta=.2，因此整体是既有 transported 两 beta 调度，不是一个共同辅助目标的精确交替最小化。有限次 PR／L-BFGS 亦不提供辅助目标全程单调或原 Phi 单调保证。

## 同 prior、单一物理 Phi 与比较公平性

静态接线支持以下有限结论：

- confirmation 的 `metric` 是 `0.5||F(m)-Y||²/||Y||²+model.reg.evaluate(m)[0]`；各方法都调用相同 `MaxwellModel` 和同一个 prior callable。A 改变辅助状态罚函数，没有另加 material ridge 或形状先验。冻结 payload 参数检查与 `Regularizer` 的既有默认参数相匹配。
- `physical_sigma=||Y||` 使辅助阶段的物理候选排序也是该 Phi。confirmation 再显式重评 host 和真实 stage endpoints，按原 Phi 选一个 seed，用同一个 `refine(maximum=1000)`，最终保留 host／seed／refine 中的原 Phi incumbent。
- `refine` 的目标和梯度均为同一个物理目标，乘 `1e6` 的优化缩放对各方法相同。它保留所有已评估点的最好物理 Phi，不必等同于优化器最后返回点。没有看到按 GT 选择 seed、候选或根的调用。
- 这些是模块静态边界，不替代对外层 runner 是否只传 ALLOW payload、封存顺序及私有文件访问的执行审计。此模块的 `truth_access='NONE'` 本身不是操作系统级隔离证明。

公平性仍有三个必须保留的限定：

1. **同 prior 函数不等于同辅助 landscape／同 prior 相对作用。** A 对状态残差施加材料相关左变换，B 对数据残差施加 M(m)；这会改变辅助阶段与未缩放 R(m) 的相对权衡。可以说最终同 Phi、prior 未改，不能说所有辅助优化问题相同，或从相同 prior 名称排除全部隐式偏好。

2. **不是等算力结果。** 开发审计记录 onehop 57,669 次材料评估、global 7,500 次；native_extended 还有更大 mf12 容量和 8,000 次辅助更新。onehop／self／block 的材料内层框架、上限与容差相近，self 本身也有 70,051 次评估，因而 self 失败是有用的局部对照，不能简单归咎于 onehop 独有 L-BFGS；但它仍不能把跨 global 的全部差异归结为物理输运，或给出效率优越性。开发 wall 为单次测量，保留实际局部成本即可，无需新增任意预算门。

3. **同 corrector 上限排除了“给 A 换 corrector”这一直接混杂，但不是 basin 证明。** 不同 seed 会产生不同线搜索／评估次数；开发端点均在 1,000 次上限停止，未认证驻点、局部最小或其吸引域体积。允许报告失败端点到成功端点的经验性改善，不能称全盒盆地扩大。native_extended 失败也不代表所有普通 CIE 日程必败。

冻结确认集只有六个同算子、同 acquisition、同对比度等级、指定三类形状的 fresh instances；其结论应保留这个分布范围。`block4_cached` 是 compatibility arm，三主方法的候选池成本按三方法之和记账；不能把候选池的原 Phi 最小值视作单一 onehop 的成绩。冻结合同已经明确禁止 equal-compute 声称，这一点应继续保留。

## A：恒等式、机制和既有工作的边界

exact block absorption 是可逆状态残差的左变换；同物理零集来自代数，并不自动是新算法贡献。onehop 是该局部残差变换的一阶截断；它没有截断最终 forward 的多重散射，因为物理排序／corrector 仍调用全 Maxwell resolvent。

开发数据中 self 失败，而 exact／cached／onehop 均通过原开发联合门，支持“完整块无限循环不是该案例成功的必要条件”和“指定 self-only 算法不足以复制该案例收益”。它不证明一般材料中一次交换足够，也不证明 exact 与 onehop 的优化器／material 结果等价；即使残差近似误差有小界，非凸优化轨迹也可能分离。

prior-art 摘录已明确涵盖空间 beta 的 NIE/CIE、局部 Green 抽取、离散 self 重整化及一般矩阵 T 重整化。未核实完全相同四格罚函数／算法仅是检索未决，不能推成 first。此审查使用该摘录作风险边界，没有额外检索或独立核验其全文。若保留新颖性候选，只能针对具体 off-shell 目标、一次交换实现和已限定的经验表现继续判断；不宜声称发明了 Green 分解、Neumann 级数、可逆重整化或辅助非物理电流。

## B：profile 梯度正确，但范围有限

`q=A(j-j*)`、`C=SA^-1` 的消元给出

`q*=-C^H[eta I+(1+lambda)CC^H]^-1 e`，

`M=(eta I+lambda CC^H)[eta I+(1+lambda)CC^H]^-1`。

`port_profile.py:49–63` 正确形成 receiver Gram 和正定小系统，不需要求密集全 current 逆。`receiver_resolvent` 以 `A^T` 的转置求解再转置得到 `S A^-1`，这里使用 `trans=1` 而非共轭伴随是正确的。

材料梯度可在原 lifted current j 坐标使用 envelope theorem：`Cq=Sj-F`，所以固定 j 时 `d(Cq)=-dF=-C diag(dm)Qphysical`，而 `dq=-diag(dm)Qhat`。这给出文中的

`-lambda Re sum conj(z)Qphysical - eta Re sum conj(q)Qhat`，

再统一除以 `||Y||²` 并加 prior 梯度。`port_profile.py:67–78` 使用 `trans=2` 求真正的 `A^-H S^H port`，正确；没有把 `M(m)` 当常数而漏导数。

该结论要求：原 forward A 可逆；全维、无约束 current 被精确消元；eta 为固定正实标量，lambda 为固定非负实标量；源／接收权为当前 Euclidean 范数；m 是实对角材料；Y 与归一化固定；G/E/S 不随 m 改变。若 eta 随材料调节、profile current 受子空间约束／近似求解、换成向量张量材料或自适应 metric，需要相应新推导，不能沿用此梯度公式。

对 `CC^H` 的特征值 h，M 的特征值为 `(eta+lambda h)/(eta+(1+lambda)h)`。这说明有限 eta 确实改变材料相关权重；没有加入观测。eta→0 时仅在 range(C) 上趋于 `lambda/(1+lambda)`，其正交补仍为 1；零噪声、有限 eta 或加上未同步缩放 prior 时，都不能随意说与原目标完全相同。lambda=0 是 state／wavefield relaxation 对照，lambda=1 的特定开发结果更低残差却更差材料，不能据此给 B 成像优势。B 不在当前冻结三主方法中，不能用 A 的确认结果替它获得支持。

## 唯一建议的后续验证（本轮不执行）

若下一步重点是区分真实 Green 相位方向与一般同强度 coupling，最有信息量的是 **onehop 的单个 pi 相位反转控制**。在另一个明确的开发研究中，将 onehop 状态项

`|| (u+Lu) Omega^(1/2) ||F²`

仅改为

`|| (u-Lu) Omega^(1/2) ||F²`。

仅在辅助残差变换中置 `Goff -> -Goff`，真实 forward 的 G 不变。保留同一分块、d/a、源指标、prior、初始点、beta 日程、current／材料上限和相同物理 corrector，完整包含该控制目标自身的材料导数；current 和 material 接线都要一致反号，不能冻结 L 的导数。它保留 `|Goff|`、块拓扑、逐行可逆界 q、`||u||²` 和全部 `||Lu||²`，只将 `2 Re tr(u^H Lu Omega)` 的符号反转。由于 diag(L)=0，两个 metric 的对角项 `diag[(I±L)^H(I±L)]` 也相同；局部矩阵应用／标量除法操作结构相同。`q<1` 同时保证 I+L 和 I-L 可逆，因此两控制仍有同一物理状态零集。

用既有物理 Phi 排序，预算终止状态照实保留，不新增成功门或 GT 选择。若 plus 成功而 pi 失败，支持“相位相对于自残差的方向在这个案例／日程中有作用”，并排除“仅幅值与该同强度反号 coupling 就足以复制”的具体解释；不能排除所有其他任意 coupling，亦不能证明真实 Green 相位普遍最优。若两者均成功，则正确符号不是该案例成功的必要条件。该单一控制不需要随机种子或扫参；不得把它插入已冻结的 confirmation 或事后改写其分母。

## 简短问题清单

1. 将“actual mechanism”限制为正确的结构解释；在相位控制前，不称交叉项／真实 Green 相位已有独立因果确认。
2. 保留同 prior、最终同 Phi，但明确辅助状态／数据 metric 和相对权衡改变；双 beta 不是单一辅助目标下降定理。
3. 不作等算力、收敛、local-minimum 或普遍 basin enlargement 声称；corrector MAXITER 和不同内层成本必须随结论保留。
4. A 的可逆界只覆盖已用算子／盒／beta／范数；B 梯度只覆盖固定 eta／lambda 和精确全 current profile。
5. 既有重整化／CIE／variable projection 有直接先例；未找到同式实现不构成新颖性证明。
6. 本审查不读取 confirmation 结果，不判断 fresh-confirmation 是否支持成像优势；2D scalar-TM 证据不能外推为未运行的 3D vector-Maxwell 结果。

## 补充审查：第 7 节 material readout 的 metric 不变子空间条件

本补充仅静态读取新增 `THEORY_DERIVATIONS.md:187–230`，不实验、不读取任何 GT 或新增结果。第 5 节现已将先前过强的“actual mechanism”表述收紧，明确同时存在交叉项和二次项，且不将双 beta 扫描称为单一辅助目标下降；原审查针对这些措辞的问题已得到回应。

**代数正确，条件明确。** 对实化后满列秩 `C in R^(k×n)` 及固定实对称正定 M，有 `D_M C=I`，从而

`D_M Pi=D_M C(C^T C)^-1 C^T=(C^T C)^-1 C^T=D_0`。

因此 `D_M-D_0=D_M(I-Pi)=(C^T M C)^-1 C^T M(I-Pi)`。由于 `C^T M C` 可逆，`D_M=D_0` 等价于 `C^T M(I-Pi)=0`；M 对称使之等价于 `(I-Pi)M C=0`，即 range(C) 对 M 不变。实对称 M 对该空间不变时，其正交补也不变，故又等价于 `[M,Pi]=0`。这里的 iff 只针对“对所有右端 b 的整个线性读出算子相同”，不存在遗漏的单向条件。

**零缺陷处冻结 P 的一阶使用正确。** 对 `P(m)(b-C dm)`，在基准缺陷为零的一阶线性化中，`dP` 乘基准零缺陷，故不贡献一阶残差；读出的局部二次问题使用基准 `M=P^T P`。这并不授权在原先非零缺陷的实际 material inner objective 中冻结 P：那里 `dP` 项必须保留，原 A 代码也确实保留。满列秩要求每 cell 的实化照明向量非零；全灭场 cell 的非唯一读出不在该逆矩阵公式范围内。

**局部标量消去及 square-C 边界均正确。** 固定源指标被吸收到字段后，cell 内所有实／虚／source 分量共享的正标量 metric 是各 cell 上的 `w_i I`，保持由 C 的分离支撑列张成的空间，因而无正则、内点的一阶读出不变。特别地，所有 sources 共享的 cell-wise 复标量 P 在实化后是缩放旋转，产生 `P^T P=|p_i|²I`，同样消去；不能把它与 sources 间不同或任意块 metric 混同。若 C 方阵且可逆，`Pi=I`，确实所有正定 weighting 都消去。每 cell 一实源、一实方程、非零照明的 toy 属于这种边界。

**解释需保留的两个限定：**

1. 不变性被打破只保证存在某个 b 使两个线性读出不同。对具体实际缺陷，必要且充分的差异条件是 `C^T M(I-Pi)b != 0`。若 `b in range(C)`，两种 metric 对它永远给出同一步；即使有补空间且 M 非局部，也可能保持 range(C)，或实际 b 的补空间分量落在 `C^T M` 的零作用子空间。因此“非局部”及“补空间非空”本身都不是观察到收益或实际一步改变的充分条件。当前“generally does not”可作泛型结构说明，不能作为当前 Maxwell J/Q 已检验违背不变性的数值证据；实际 Fourier/current 子空间还能限制可达 b。

2. 一个复 source 配合共享实材料时，每 cell 已有两个实方程，一般便有非零本构残差补空间；多个复照明不是该补空间存在的必要条件。原文所说多复照明确实保留补空间，并未明确声称其必要性，但建议交付解释时将边界表述为“恰一实方程的 toy 会丢失机制；单个复照明已可能保留”。第 212 行的 sign/phase “matter”宜理解为 “can matter”，而非从代数推出必然影响。

`+L/-L` 的 metric 可以因交叉项符号不同而对 range(C) 产生不同耦合，因此第 213 行使用 “possible difference” 合适；公式不预测哪一个更好，亦不推出两者必然不同。父线程告知了 pi 物理拟合失败、plus 恢复；本补充未读取或独立核验该结果。即使该对照成立，也只是与条件性结构解释一致，不能据此把零缺陷一阶 LS 公式提升为实际非零缺陷、带 prior/box、有限内层、双 beta 全迭代的因果或收敛定理。

将本节标为既有加权 LS／R3-D1 消去现象的条件精化而非原创投影定理，范围适当。现有第 224–230 行对 prior、box、非零缺陷和不精确求解的限制充分，建议再明确上述“存在某个 b”与“实际 b”的区别即可；无需追加实验、门或机器精度恒等式检查。
