# R10 相关工作与可检验的区别

R10 的主要结果是用有限照明所得的复散射数据，严格排除整片材料盒。设 \(d=k-\operatorname{rank}B\)。把 \(d+1\) 个试探材料列投影到 \(B\) 的左零空间，取带符号余子式 \(\lambda_j\)，得到 \(K=\sum_j U_{i_j}(q_j)\lambda_j=B\alpha\)。因此 \(\alpha\) 可以直接组合已采集的物理读数。只要被选中的 \(d+1\) 个材料值正确，精确模型便满足该组合给出的恒等式，其余 \(k-d-1\) 个材料可在原盒任意变化。将整组复加性误差与一个共同复增益的联合椭球传入恒等式，得到对每个试探变量分别至多二次的必要不等式。对最后一维的二次式检验全区间正性，Bernstein 细分只遍历剩下 \(d\) 个坐标。这个外包覆盖声明误差及已存离散系数下的所有可行材料。 [R10 命题一、二和正性证明](MATHEMATICAL_THEORY.md)

## 行列式公式与区间证书的直接先例

- **Boyer、Garzella 与 Guevara Vasquez，2016。** 离散 Schrödinger 反问题的式 (4.12) 已把含节点势的边界响应写成 \(L_{BB}-L_{BI}(L_{II}+\operatorname{diag}q)^{-1}L_{IB}\)。从该式用伴随矩阵清除分母，每个节点势分别至多出现一次，这是行列式的多线性。式 (6.1) 另将两个完整边界响应之差写成势差与两组内部解乘积的配对，Theorem 6.2 用该乘积张成空间证明条件性的几乎处处两两可辨识。R10 的材料列来自单接收响应行，余子式取在有限照明的缺失空间，作用是把试探材料的检验压回同一组已采集读数。 [Boyer 等原文 §4、§6](https://arxiv.org/pdf/1510.02848) · [期刊 DOI](https://doi.org/10.1137/15M1043479)
- **Garloff 与 Graf，1999。** §3.1 式 (8) 以张量 Bernstein 系数的最小值和最大值夹逼整个参数盒上的多项式值，§3.2 用 de Casteljau 型递推二分。§4–5 已把正性判据用于严格多项式不等式，按可证满足、可证排除和待定分类。R10 复用这套盒上正性逻辑。其低维搜索来自物理约束留下的逐变量二次结构，每个枢轴控制项的被消去变量是一元二次，可以在整段上解析检验。 [Garloff–Graf 原文 §3–5](https://www-home.htwg-konstanz.de/~garloff/SolvPolyIneq.pdf)

余子式消元、逐变量多线性与 Bernstein 判定各有成熟先例。R10 的直接六维对照也能排除同一批已证排除区域，它从完整正演出发，先精确消去共同复增益，再在所有材料坐标上做 Bernstein 证明。两种方法在同一离散模型和误差界下的实际差异是维数、保留的联合误差信息，以及证完同一片区域所需的工作。针对「有限照明缺失维数 → 可观测余子式 → 带联合误差的全盒外包」这一整条构造，所核读的原文仍不足以判定优先权。 [六维判据与比较](MATHEMATICAL_THEORY.md#8-六维直接对照也可以严格证明相同的事) · [逐区域记录](results/direct_box_unselected.json)

## 已知根消去和低维提升

- **Farrell、Birkisson 与 Funke，SIAM J. Sci. Comput.，2015。** §2 的式 (2.6)、(2.8) 对残差乘以 \(\|u-u_\star\|^{-p}+\alpha\)。Fig. 2.1 演示无移位会给出远处的虚假小残差。§4.2 的 Lane–Emden 方程测试用同一零初值与 \(p=1,\alpha=1\) 找到两个解。R10 沿 \(z=z_\star+tn\) 的精确锚点展开为 \(F=a\,\Re t+b\,\Im t+c|t|^2\)。令 \(s=t/|t|^2=1/\bar t\)，则 \(F/|t|^2=a\,\Re s+b\,\Im s+c\)。按 \(\|n\|^2\) 归一后，左侧就是 \(p=2,\alpha=0\) 的范数消去残差。可检验的区别是这里能用一次线性最小二乘求出复幅的两个实分量，方向仍需搜索。\(s=0\) 无有限逆像，带噪排序和收敛判断应回到原始兼容残差及独立正演。若锚点只近似满足兼容性，原始残差还会产生 \(F(z_\star)|s|^2\)，线性子问题实际求解的是扣除锚点残差后的方程。 [Farrell 等原文](https://arxiv.org/pdf/1410.5620)
- **Golub 与 Pereyra，SIAM J. Numer. Anal.，1973。** 可分离最小二乘以伪逆消去线性参数，优化余下的非线性变量。R10 的 \(s=-A(n)^+c(n)\) 正是这一类方法。比较方向和复幅联合优化可辨别收益是否来自变量投影，比较同一物理残差上的多初值与移位消去可区分换元和一般消去的收益。 [Golub–Pereyra 原文](https://doi.org/10.1137/0710036)
- **Gao 等和 Huang–Xu 的 affine phase retrieval。** Gao 等 Lemma 3.1 已用复平面三非共线中心的平方距离差恢复单复变量。Huang–Xu §4.3 Lemma 4.7 已用增广坐标 \((z,1)\) 的 Hermitian 提升刻画碰撞，两个解之差属于测量映射的核，秩至多为二，右下元素为零，最后一列的非末项非零。R10 右下元素固定为一，因此自由坐标数 \(d^2+2d\) 直接来自 Hermitian 块和复列的维数。满列秩是线性提升唯一的充分条件。一旦与模型或数据不符，线性最小二乘的矩阵可离开秩一曲面，应检查 \(Z=zz^\mathrm H\)。\(d=1\) 时，可直接求测量仿射子空间与 \((x,y,x^2+y^2)\) 的交。有效线性秩为二时至多两个孤立解，秩更低可以有连续族，非退化与秩都应由实际算子判断。Huang–Xu Theorem 4.6 的 \(4d-1\) 泛型阈值假定平方仿射范数的特定矩阵族，无法直接给 Maxwell 行约束计数。PhaseLift 的 Thms. 1.1–1.2 则对随机强度读数的 SDP 给出精确性和噪声界。R10 直接使用复数散射样本，二次等式来自材料兼容性。 [Gao 等原文](https://arxiv.org/pdf/1608.06117) · [Huang–Xu 原文](https://arxiv.org/pdf/1805.07899) · [Candès 等 PhaseLift 原文](https://arxiv.org/pdf/1109.4499)
- **代数全解参照。** Tisseur–Meerbergen 的二次特征值问题用伴随矩阵铅笔线性化 \(Q(\lambda)v=0\)，假定一个多项式谱参数。R10 的 \(t,\bar t\) 原样不具有这种形式。实际的低维代数参照是把实部和虚部分开，直接做消元或数值同伦，验回所有原方程，只保留实且有限、材料可行的候选。一般零空间维数下，复直线方向搜索没有提供全根枚举保证。 [Tisseur–Meerbergen §1、§3](https://doi.org/10.1137/S0036144500381988) · [T.Y. Li 数值多项式同伦](https://doi.org/10.1017/S0962492900002749)

## 电磁反演里已经做到哪一步

- **SOM 是最贴近的全波反演对照。** Chen、Zhong 与 Agarwal 的原文 §3.1、§3.2 用数据算子的 SVD 得到确定电流，再在余下子空间优化电流与对比度，目标含散射数据及状态失配。§4.2 的 3D 包覆立方体用 30 个电偶极子和 30 个接收器。早期 twofold SOM 对弱观测子空间再约束。R10 的物理未知量是单接收数据约束下的电流系数，因此在同一测量矩阵下，以原兼容系统多初值求解和 SOM/CSI 式双残差求解，才有可比较的分支命中率。 [Chen 等原文 §3–4](https://archive.intlpress.com/site/pub/files/_fulltext/journals/maa/2010/0017/0004/MAA-2010-0017-0004-a006.pdf) · [Zhong–Chen twofold SOM](https://doi.org/10.1088/0266-5611/25/8/085003) · [van den Berg–Kleinman CSI](https://doi.org/10.1088/0266-5611/13/6/013)
- **最接近材料歧义的非线性算例。** Chen、Lu、Karniadakis 与 Dal Negro 的 §IV Fig. 7 用单平面波照射不对称的两圆柱。PINN 从模拟总场求一个等效介电剖面，该剖面和原目标明显不同，重新 FEM 正演的总电场在计算窗内相对 \(L^2\) 差为 0.30%。作者指出该等效剖面可能不唯一。这个例子支持近似场匹配下的材料歧义。R10 的同数据分支还需要报告材料差异、复散射数据误差和真实兼容残差。 [Chen 等 Optics Express 作者稿 §IV Fig. 7](https://arxiv.org/pdf/1912.01085)
- **生成多候选后用正演筛选已有实验。** Chehelgami、LoVetri 与 Khoshdel 从同一散射数据条件化的潜空间扩散产生多个介电图，用 2D 标量矩量法正演，输出原数据失配最低的一幅。作者稿 §II.B.4 和 §III.A–B 用未见合成目标及校准实测数据分别评价，Table I 分别报告数据和介电图 MSE。R10 的代数兼容可用于检验候选是否确实对应同一数据方程，追加的独立激励用于区分已发现的候选。 [Chehelgami 等作者稿](https://arxiv.org/pdf/2510.25729) · [IEEE Open J. Antennas Propag. 论文](https://doi.org/10.1109/OJAP.2026.3727556)
- **局部二次电磁约束已有直接先例。** Molesky 等从 Maxwell VIE 消去总场，以局部见证算子 \(Q_j\) 和极化电流写出二次等式，再做对偶或 SDP 松弛，继而从可行或近似可行极化推回材料作为光子器件设计起点。其 §I.A 式 (2) 的等式以给定填充材料生成，灰度情况可放宽为不等式。全局光学定理只约束整体收支，局部版本保留各空间投影的约束。R10 使用每个未知实对比度的局部条件 \(\Im(w_i\bar v_i)=0\)，与复接收数据一同求解。恢复材料时还取决于 \(v_i\ne0\)、材料范围与正演成立。 [Molesky 等式 (2) 与 §I.B](https://arxiv.org/pdf/2504.14083) · [Molesky–Chao–Rodriguez 局域化 T 算子原文](https://arxiv.org/pdf/2008.08168)
- **两种相邻的不可辨识性。** Li–Wang §3.4 和 §4 给出固定频率 Maxwell 非辐射电流的核空间刻画及显式例子，所求是源。真实材料分支的电流还须来自该材料在指定入射下的自洽极化。Kilgore–Moskow–Schotland §3 和 Thm. 2 从正演 Born 项的系数递推构造逆级数，引入正则化的一阶伪逆并给局部收敛充分条件。它提供弱散射和局部反演控制，级数收敛结论本身不比较远离锚点的多支。 [Li–Wang 原文](https://arxiv.org/pdf/2402.10407) · [Kilgore 等原文](https://www.math.drexel.edu/~moskow/papers/MaxwellBorn.pdf)

## 多解如何用实际读数区分

- **二候选最优照明已有原样的判别算子。** Bouchet 等 Phys. Rev. Lett. 2021 式 (2)–(3) 直接最大化单位入射光子数的 \(\|(S_2-S_1)u\|^2\)，最优输入是 \((S_2-S_1)^\dagger(S_2-S_1)\) 的最大本征态。Figs. 1–3 用 532 nm 光、DMD 和双正交分量探测，在两扩散片之间分别有、无微球时测两张透射矩阵，并以未优化照明作对照。发表实验针对目标存在性报告入射光子用量下降超过两个数量级。R10 在功率 Gram \(G\) 和噪声协方差 \(\Sigma\) 下若使用 \(u^\mathrm H\Delta S^\mathrm H\Sigma^{-1}\Delta S u/(u^\mathrm HGu)\)，是同一准则在发射与接收度量下的版本。关键证据是候选先由反演生成、再预测并获取新的物理读数，以留出的读数判别。 [Bouchet 等原文式 (2)–(3)](https://dbouchet.github.io/assets/publications/PhysRevLett.127.253902.pdf)
- **全波实验的参照尺度。** Liu 等 SOM-Net §IV.E 在 Fresnel FoamDielExt 分别做 3 和 4 GHz，8 个发射源、241 个接收器，与 SOM、U-Net 比较，Table V 分列材料 SSIM/RMSE，Figs. 13–14 另画电流和散射数据。Yago Ruiz 等 Inverse Problems 2024 §5 把合成训练后的模型应用于独立实测 Fresnel 数据，实测有限视角，Fig. 7 和 Table 2 报告泡沫与塑料的定量估计差异，塑料偏差更大。两种报道都支持并列呈现数据拟合和材料误差。Eom 等 2026 在线发表的 TE Fresnel 原文摘要直接比较单发射和多发射，报告单发射表现随频率和源位置变化，增加发射可改善小目标定位。只有单点照明也未必意味着单频单接收，例如 Le 等的 Carleman 算例用一个点源、面上接收和 \([\pi,2\pi]\) 内 121 个波数，输入另含 10% 合成乘性噪声。 [Liu 等原文](https://arxiv.org/pdf/2209.03567) · [Yago Ruiz 等原文](https://iris.uniroma1.it/retrieve/b76805aa-90bf-4c70-9e9a-a170bd68b59c/Yago_Ruiz_Deep_2024.pdf) · [Eom 等发表页与摘要](https://doi.org/10.4208/cicp.OA-2026-0059) · [Le 等 §5 原文](https://arxiv.org/pdf/2506.21699)

## 从点分支到有界误差材料集合

R7 的虚拟实验已经在指定接收精度下计算材料可行区间，并在一个算例里通过额外接收选择两段中的一段。Jaulin–Walter 的 SIVIA §5 给出通用的非线性按误差界反求参数集合方法，同时获得材料集合的内、外夹逼，他们也将 set inversion 用于分子光散射数据。R10 若取 \(d=1\)，可沿用这套集合解释，但可以先借 Maxwell 条件降低约束的维数。设 \(nB=0\)、\(U_i(q)=e_i-qH_i\)、\(b_i(q)=nU_i(q)\)。则 \(K=U_i(q)b_j(p)-U_j(p)b_i(q)\) 属于 \(\operatorname{ran}B\)，所以 \(\alpha=B^+K\) 给出可由旧源合成的读数。若 \(q,p\) 是对应真材料，理想数据满足 \(y\alpha=qR_i b_j(p)-pR_j b_i(q)\)。加入共同接收增益和有界误差后可用低阶不等式对 \((q,p)\) 作可行投影，新的物理照明可以继续收缩这个集合。能证明的区间外排除依赖正确传播噪声和增益界，以及消元时保留退化分支。 [R7 虚拟实验](../research_z_alg_r7/index.html) · [Jaulin–Walter 1993 原文](https://www.ensta-bretagne.fr/jaulin/paper_automatica93.pdf) · [Jaulin 等光散射应用](https://doi.org/10.1088/0305-4470/30/22/012)

检索覆盖已知根消去、affine phase retrieval、Maxwell 局部功率、全波子空间反演、非辐射源和最佳照明。Semantic Scholar 返回 HTTP 429，以上对照依据开放原文或已注明的出版页摘要。

R10 对共享误差的分析还使用了白化后各投影子空间的几何。Casazza 与 Kutyniok 定义 3.15、命题 3.16 给出加权投影之和及其谱界，本轮用于控制多个短检验共同遗漏的误差。本文的尖锐最坏放大率来自单位球上的最大投影再取最小，有限真子空间不覆盖整个空间是基础线性代数。与物理模型相接的是从实际入射和材料得到这些特定投影，并在保存的一组全波读数上严格构造能通过所有两项检验的错误材料。[Casazza 与 Kutyniok 原文](https://arxiv.org/pdf/math/0311384)
