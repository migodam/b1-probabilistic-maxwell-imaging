# 能量、被动性与互易：能复用什么，不能凭空增加什么信息

主线程推导；不运行光学定理机器精度测试。本章的标量连续式与离散代数式分开。所有“功率”解释需要正确场归一化、单元质量与外边界；32个有限采样receiver不是全角功率计。

## F-BI-E01：二维TM的局部能流

采用e^{-iωt}、实正磁导率μ、E_z=u。Maxwell方程给出

\[
H=(\partial_yu,-\partial_xu)/(i\omega\mu),\qquad
\langle S\rangle=\frac{1}{2\omega\mu}\operatorname{Im}(\bar u\nabla u).
\]

若Δu+ω²μεu=−iωμJ_z，则乘以\bar u并取虚部直接得到

\[
\nabla\cdot\langle S\rangle
=-\frac\omega2\operatorname{Im}\epsilon\,|u|^2
-\frac12\operatorname{Re}(\bar uJ_z).
\]

在无外加源的区域积分，外向净通量等于负的吸收功率。若背景吸收，应包含背景的吸收，不能只计ROI增量。当前P1使用固定有限网格宿主与原Green合同；把宿主消去并没有创造无损参考背景。

**性质：标准条件物理定理。** 未测量闭边界场/法向导数，不宣称在本轮实际验证全能量平衡。

## F-BI-E02：离散VIE的功率型恒等式

以均匀单元质量w为例，原方程j=D_mE、E=e+Gj。左乘j*并取虚部得

\[
\operatorname{Im}(e^*j)
=\sum_i\operatorname{Im}(m_i)|E_i|^2+j^*\operatorname{Im}_H(G)j,
\quad \operatorname{Im}_H(G)=\frac{G-G^*}{2i}.
\]

左右同乘单元质量与相同物理归一化即可对应消光、材料吸收和向背景耦合。非均匀网格应在同一质量内积下重写，不能直接沿用欧氏矩阵。推导没有除以m_i，零材料格仍成立。

该等式作为离散代数对任意G成立；将最后一项称为非负辐射/背景耗散，还需要离散\operatorname{Im}_H(G)的相应半正定性质。点配点/自项近似并不自动继承连续算子的所有能量结构。

精确宿主消元后，同样对h=D(E_R+Ch)有

\[
\operatorname{Im}(E_R^*h)
=\sum_i\operatorname{Im}(d_i)|E_{ROI,i}|^2+h^*\operatorname{Im}_H(C)h.
\]

这是**增量系统的代数关系**，不是完整材料相对宿主的总功率差；后者含宿主场变化及干涉项。

## F-BI-E03：T算子的耗散结构

在非零材料子域D可逆时，T=(D^{-1}−C)^{-1}。由逆矩阵差恒等式

\[
\operatorname{Im}_H T
=T^*\left[-\operatorname{Im}_H D^{-1}+\operatorname{Im}_H C\right]T.
\]

若D为对角，−Im(d_i^{-1})=Im(d_i)/|d_i|²。该式说明被动材料与背景耗散共同限制响应。D有零格时使用E02或在活跃支撑上限制，不能用伪逆无条件替换。

这是经典光学定理/散射算子约束的离散形式。Le Ru等2013研究T/K矩阵的辐射修正；Molesky等2020利用功率约束构造散射上界，并非本轮新定律。[原始论文](https://doi.org/10.1103/PhysRevA.87.012504)、[T算子界](https://doi.org/10.1103/PhysRevResearch.2.033172)。本章自身代数推导不借这些摘要替代条件证明。

## F-BI-E04：固定形状的非正规稳定性下界

令P=√A为实非负对角、B=PCP，K_δ=δ^{-1}I−B。非Hermitian B无需正交本征分解。定义

\[
\beta=\lambda_{\min}(\operatorname{Im}_H B),\qquad
b_\delta=\frac{\operatorname{Im}\delta}{|\delta|^2}+\beta.
\]

若b_δ>0，对任意v，−Im(v*K_δv)≥b_δ||v||²；由Cauchy–Schwarz得到

\[
\sigma_{\min}(K_\delta)\ge b_\delta,
\qquad \|K_\delta^{-1}\|_2\le b_\delta^{-1}.
\]

它不要求B正规，故比单看极点距离更稳妥。合法幅度矩形Reδ∈[.8,1.8]、Imδ∈[.03,.18]内，因为b<a，b/(a²+b²)随a下降、随b上升，故

\[
b_\delta\ge \frac{.03}{1.8^2+.03^2}+\beta.
\]

对数值β必须说明特征值误差。计算出的正下界是浮点条件诊断，不自动成为严格区间证书。实际8个形状的β略负，因此本报告**没有**将离散ImB称为PSD；正的材料损耗项仍足以给出数值正稳定余量。数据见 `results/transfer/COERCIVITY_ZERO_RHS.json`。

## F-BI-E05：可用于压缩的残差—端口误差关系

精确状态u=K_δ^{-1}R，近似\tilde u的完整残差r=R−K_δ\tilde u，端口响应为Lu。于是

\[
F-\tilde F=LK_\delta^{-1}r,
\qquad \|F-\tilde F\|_F\le\frac{\|L\|_2\|r\|_F}{b_\delta}.
\]

这是离散条件误差界；固定背景/形状只需准备β和||L||。若计算残差比精确三角解更贵，或者界过松，则没有算法收益。可用后验端口伴随进一步收紧，但其费用也必须计入。本轮是否真正实现压缩必须以实验文件为准，不能由该式直接宣称加速。

若白化预测误差上界为e，近似白化残差\tilde r_y，则高斯势差满足

\[
|\Phi-\tilde\Phi|\le\|\tilde r_y\|e+\tfrac12e^2.
\]

这将物理误差连接到候选noise gate或积分likelihood误差；不能连接到缺失材料分支质量而无需额外假设。

## F-BI-E06：互易不等于Hermitian，也不等于新增数据

若同一离散质量下C^T=C，D对角，则

\[
T(D)^T=D(I-CD)^{-1}=(I-DC)^{-1}D=T(D).
\]

形状B=PCP也为复对称，但一般B*≠B。其复Schur分解合法，Hermitian实谱理论不合法。观测O T E_R是32×8矩形；源/接收不共址、标签不同，不能把该矩阵随意转置当作同一测量。互易只在匹配源接收、极化与归一化时对应交换实验。

被动性原本就在P1 prior里，重加一次不能提供新的信息。完整功率需要边界/角积分；||32receiver data||²只是该采样协议的范数。没有新的功率观测，本轮不增加能量罚项，也不将其当作新的候选筛选优势。

## 裁决

保留E04/E05作为物理压缩的可计算稳定性门；这是标准耗散算子理论在当前形状传递中的有用组合。全功率实验、能量增强成像和新光学定理均为NOT_RUN。固定频率改变δ与改变频率不是一回事；频率因果性、Kramers–Kronig或矩阵值振子理论不能直接套到材料幅度参数上。[相关频率理论](https://www.nature.com/articles/s41467-023-43221-2)。
