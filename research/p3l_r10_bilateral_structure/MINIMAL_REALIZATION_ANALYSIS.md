# R10 最小输入输出实现、全材料消去与响应等价

日期：2026-10-30。标准有限维实现理论在固定形状 Maxwell 传递族上的应用；不是新定理声明。本文不运行物理，不以阈值数值 rank 代替严格 rank。

## 1. 正确的压缩对象是传递族

固定形状 A 时，令 \(P^2=A\)、\(B=PCP\)、\(L=OP\)、\(R=PE\)。完整响应为
\[
G(\delta)=F(\delta)-F_h=\delta L(I-\delta B)^{-1}R.
\tag{M1}
\]
把 \(z=1/\delta\)，有 \(G(1/z)=L(zI-B)^{-1}R\)。因此是标准严格真有理传递函数，B 为内部状态，R 为输入，L 为输出。这里“时间系统”只是数学实现解释；材料幅度 \(\delta\) 不是物理时间。

在零附近
\[
G(\delta)=\sum_{k=0}^{\infty}\delta^{k+1}M_k,\qquad M_k=LB^kR.
\tag{M2}
\]
级数在例如 \(|\delta|\|B\|<1\) 下收敛。M1 则在全部无极点点成立，不依赖 Born/Neumann 收敛。

## 2. 可达、不可观测与最小阶

定义
\[
\mathcal R=\operatorname{span}\{R,BR,\ldots,B^{q-1}R\},\qquad
\mathcal N=\bigcap_{k=0}^{q-1}\ker(LB^k).
\tag{M3}
\]
由 Cayley–Hamilton，二者都对 B 不变；\(\mathcal R\) 包含全部会被源激发的状态，\(\mathcal N\) 包含在任意传递阶都不被接收看见的状态。精确最小实现状态空间是
\[
\mathcal R/(\mathcal R\cap\mathcal N),\qquad
r_{min}=\dim\mathcal R-\dim(\mathcal R\cap\mathcal N).
\tag{M4}
\]

只删 \(\ker L\) 不合法：\(x\in\ker L\) 可能在下一阶有 \(LBx\ne0\)。只看 rank(R) 或 rank(L) 也不够；单输入单输出系统仍可有 q 个可达可观测模态。

令块 Hankel 矩阵
\[
H_{a,b}=\begin{bmatrix}
M_0&M_1&\cdots&M_{b-1}\\
M_1&M_2&\cdots&M_b\\
\vdots&\vdots&&\vdots\\
M_{a-1}&M_a&\cdots&M_{a+b-2}
\end{bmatrix}
=\underbrace{\begin{bmatrix}L\\LB\\\vdots\\LB^{a-1}\end{bmatrix}}_{\mathcal O_a}
\underbrace{[R,BR,\ldots,B^{b-1}R]}_{\mathcal C_b}.
\tag{M5}
\]
当 a,b 足够大（a,b≥q 是充分上界），\(\operatorname{rank}H_{a,b}=r_{min}\)。这是经典 Ho–Kalman/实现理论。短 Hankel 的 rank 仅是最小阶下界；小奇异值截断给经验有效阶，不是证明状态可精确删除。当前 q=256 时，显式构造巨大 Hankel 未必必要，双侧 Krylov/可达可观测基可更经济，但其费用也须计入。

## 3. 有限矩匹配不等于全传递等价

两个实现若全部 M_k 相同，则在零附近解析相同，并由有理函数延拓在共同无极点区域相同。若两个已知有限状态实现维数分别为 q1、q2，令差实现为 \(\widehat B=\operatorname{diag}(B_1,B_2)\)、\(\widehat L=[L_1,-L_2]\)、\(\widehat R=[R_1^T,R_2^T]^T\)。Cayley–Hamilton 表明：若
\[
M_k^{(1)}=M_k^{(2)},\qquad k=0,\ldots,q_1+q_2-1,
\tag{M6}
\]
则后续全部矩也相同。这是精确代数的充分条件，不是建议以数值近零矩作廉价认证。

反之，匹配少数 M_k 仅是局部幅度近似。对截断 p 阶的模型，若 \(|\delta|\|B\|<1\)，直接尾界为
\[
\left\|G(\delta)-\sum_{k=0}^{p-1}\delta^{k+1}M_k\right\|
\le\frac{|\delta|^{p+1}\|L\|\|B\|^p\|R\|}{1-|\delta|\|B\|}.
\tag{M7}
\]
当条件不满足，此界失效，不代表原 resolvent 发散或真实模型无解。强散射时应使用目标幅度域上的有理近似和实际端口误差，不能把低阶 Born 匹配当全域性质。

有限幅度节点的完全一致也通常不足以认证解析相同；若已知严格分子分母阶并有足够精确独立节点，可以用有理插值唯一性，但浮点噪声、未知阶和极点邻近使其变成条件数问题。R10 有限 sketch 只作为选择/诊断，不作 gauge 认证。

## 4. 相似变换是内部坐标等价，不是物理材料变换

两个最小实现产生同一传递函数时，存在可逆 S 满足
\[
B_2=S^{-1}B_1S,\qquad R_2=S^{-1}R_1,\qquad L_2=L_1S.
\tag{M8}
\]
这表示同一输入输出行为的状态坐标可以不同。一般 S 不保留 \(B=PCP\)、\(L=OP\)、\(R=PE\) 的共同物理因子结构，也不保留对角材料、固定几何或 decoder 像；因此不能解释成两幅物理材料之间的 gauge。

三层必须分别命名：

1. **有限观测相容：** 一次 y 或有限节点上预测差在容差内；R10 主候选资格属于此层。
2. **条件传递等价：** 固定采集、某个形状幅度族在全部允许幅度上相同；需要强得多的证明。
3. **材料 gauge：** 原材料族内存在明确变换，使声明的完整观测族不变，且通常要求可连续连接/可逆作用等结构。前两层不自动推出这一层。

有限数据无法区分的材料可以是离散、多分支或近似相容，不必形成光滑纤维。用词“possibility set/候选集合”比“已证 exact material fiber”更符合目前证据。

## 5. 全 D 固定状态投影为何可能不降维

要区分三种不同问题：固定形状只变标量幅度；原受限 decoder 的 D(x)；全部独立对角 D。后者更强，不能把其不可能性未经证据搬给 Full66。

考虑一个**固定线性状态嵌入** \(h=Vz\)，要求对全部允许 D 与所有源激励精确，并在约化空间内闭合。充分且标准的条件是
\[
\operatorname{range}(DE)\subseteq\operatorname{range}V,
\qquad DC\operatorname{range}V\subseteq\operatorname{range}V
\quad\text{对所有 D}.
\tag{M9}
\]
若 D 可独立变化且 E 每一行至少一个源分量非零，则 \(\operatorname{span}_{D}\operatorname{range}(DE)=\mathbb C^q\)：取 D 为单格坐标投影即可得到每个坐标轴。因此 M9 的固定状态空间必须为全 q 维，甚至不必使用 C 连通性。这是对“所有内部状态都由一个固定小基表示”的直接障碍。

相同结论可由 \(D=t\operatorname{diag}(d)\) 在 t=0 的导数得到；只要允许族含独立对角变化的开邻域，所有首阶 \(DE\) 已激发完整空间。若 decoder 只有受限形状变化，这一论证的前提不成立。

## 6. 输出商空间的更精确条件

内部状态全维不排除某些共同不可观测状态可以商掉。对要求共同线性闭合的 reduction，定义最大的共同不可观测空间
\[
\mathcal N_{all}=\{v: O A_{D_k}\cdots A_{D_1}v=0
\text{ 对所有有限字及所有 }D_i\},\qquad A_D=DC,
\tag{M10}
\]
其中空字条件为 Ov=0。它是包含于 ker O 且对所有 DC 不变的最大空间。一般这是**同时约化整个算子族**的条件；独立字比冻结同一个 D 的物理轨迹更强，因此不可据此宣称排除所有可能的静态输入输出公式。

给出一个透明的充分无压缩条件：若 D 是全部对角矩阵，O 的每一列非零，且 C 可逆，则 \(\mathcal N_{all}=\{0\}\)。证明：对 \(v\in\mathcal N_{all}\)，一字条件要求 \(ODCv=0\) 对所有 D。取单格 D 得 \(O_{:,i}(Cv)_i=0\)，故 Cv=0，进而 v=0。更一般在 O 无零列时恰有 \(\mathcal N_{all}=\ker O\cap\ker C\)。

结合 E 无零行、C 可逆、O 无零列，可判定标准固定线性共同状态嵌入/商 reduction 没有严格降维空间。它不是计算复杂性下界，也不排除：材料相关 Schur 补；直接有理公式；低秩近似加误差校正；某个受限 decoder 的特殊结构；统计任务压缩。

若从“所有对角 D 与 C 的公共不变子空间”出发，所有 D 不变迫使空间为坐标子集，C 的跨块耦合再排除非平凡子集。这个常见论证同样只针对对应的联合不变性类别。不能把 \([D,C]\ne0\) 这一条现象扩大成一切双边压缩都无效。

## 7. 非对易性、可观测性与端口误差互不替代

非对易性：\([D,C]_{ij}=(d_i-d_j)C_{ij}\)，通常阻止复用一个共同对角化/Schur 三角结构。它不直接给出模型误差数值。

可观测性：由 \(LB^k\) 的全阶闭包决定，不能只看 L 或 O 的零空间。

端口误差：若近似状态满足真实残差 r，则 \(\Delta F=O(I-DC)^{-1}r\)。它包含实际材料 resolvent 和最终端口；\(\|[D,C]\|\)、\(\sigma_{k+1}(C)\) 和压缩 rank 都不能直接替代。

极简结构反例不需新增数值实验：令 L=[1,0]、R=[0,1]^T、B 的 (1,2) 元素非零。则 LR=0（Born 矩为零），但 LBR 非零（第二散射阶可见）。这说明“输入子空间与输出子空间不相交”并不意味着没有响应；内部传播连接两者。相反，C 可以高秩而某个给定形状的 reachable/observable quotient 很小。

## 8. 对 R10 实验与论文的具体建议

- 如研究固定形状幅度扫描，报告达到所需端口误差的有效阶、每形状 setup、每幅度评估、原始模型回退及实际摊销点；不要只展示 B 或 Hankel 的奇异值曲线。
- 如研究 Full66 材料覆盖，最小实现只提供结构解释。需证据证明 shape/amplitude 分块在原 decoder 内合法，不能用条件族成功替换全材料结果。
- 若用 transfer sketch 选择候选，应同时保留材料距离。相同 transfer 行为的不同材料是反演不确定性的一部分，不应全部去掉。
- 主张应是结构引导的机制和有限证据；标准最小实现、相似等价、块消元和 resolvent 恒等式不作为新颖性主体。

经典实现来源：[Ho and Kalman, 1966, NASA bibliographic record](https://ntrs.nasa.gov/citations/19670049337)。其 metadata 用于身份；本文件公式给出独立标准推导，未声称逐行核对不可获取的原始全文。相关实质文献证据与检索限制见 [PRIOR_ART_AND_NOVELTY.md](PRIOR_ART_AND_NOVELTY.md)。

## 9. 已完成压缩结果的准确定位

R10 实测了标准 block Krylov/Galerkin 在固定形状幅度族上的效用。八个 DEV 形状、每形状三十节点中，32 维数值残差/势函数双门为 0/240；64、96、128 维均 240/240。随后冻结 64 维，在十个相对 DEV heldout 的锚点形状上新取三十二个幅度节点，共 320/320 通过、零精确回退。实际 M32 冷核比精确全阶 Schur 节省 63.76%，比直接 LU 节省 51.82%；[完整记录](REDUCED_RESPONSE_FRESH10.md) 区分实际 M32 与其他 M 的外推。

因此，**64 是所测基族、节点、数值门和费用定义下的经验有效阶**。这不等于：M3 的严格 Hankel rank；任意基族的最小准确 rank；所有幅度连续域的认证阶；所有形状的共同压缩阶；或 Full66 后验/材料空间维数。32 维门失败不能证明所有 32 维近似无效，因为该充分误差界可能保守。三百二十个通过节点也不能证明严格可达—可观测商只有 64 维。

该另组形状核未使用真值，但材料 GT 在其他支线已全局开放，故仅称 heldout-shape kernel probe，不称新的盲法材料测试。严格最小实现计算、向外舍入证书及全幅度统一误差认证仍为 `NOT_RUN`。该证据支持“固定形状查询有可用低阶近似”，并不与第 5–6 节针对更强 all-D 共同状态约化的限制矛盾。
