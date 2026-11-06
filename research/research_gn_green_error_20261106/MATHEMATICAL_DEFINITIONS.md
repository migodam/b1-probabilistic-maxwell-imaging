# 误差、能量与 Two-fold 的计算定义

## 1. 实际算子和误差对象

本轮固定二维scalar-TM/VIE、24×24网格，源为列。沿用冻结R4数据中的 G、S、E0：

$$L(m)=I-D_mG,\qquad j(m)=L(m)^{-1}D_mE_0,\qquad Q(m)=E_0+Gj(m),\qquad F(m)=Sj(m).$$

这里 G_in=G 为576×576，G_ob=S 为48×576；均从576维复电流域出发。没有另乘单元面积或k²，存储矩阵已包含pulse积分。不是用左奇异向量把576维材料塞入48维观测空间。

分析两个主要误差：

$$e_m=\widehat m-m_*,\qquad E_j=j(\widehat m)-j(m_*).$$

材料误差投向原G的右奇异坐标，是同一网格上的坐标诊断，不是材料Jacobian可辨识性证书。复投影的单独材料分量一般为复数，不能当作合法实材料候选。另算 [Re G;Im G] 的实右奇异基作保实控制，结果标为material_real，不替换主表中的原G定义。

电流误差则有直接物理意义：

$$Q(\widehat m)-Q(m_*)=G_{\rm in}E_j,\qquad F(\widehat m)-F(m_*)=G_{\rm ob}E_j.$$

进一步由两状态方程相减，精确得到

$$L(\widehat m)E_j=D_{e_m}Q(m_*),\qquad
E_j=L(\widehat m)^{-1}D_{e_m}Q(m_*).$$

这保留全部多次散射，不是Born。总场乘法及resolvent会混合固定谱方向，所以材料误差的G投影不能代替实际电流误差。

## 2. 单算子的强／弱／补划分

令 A=UΣV*，在右域补足V至576维。使用正交投影：

$$P_b=V_bV_b^*,\qquad X_b=P_bX,\qquad
\eta_b(X)=\frac{\|X_b\|_F^2}{\|X\|_F^2}.$$

主定义：σ/σ_max≥.1为strong；1e-12<σ/σ_max<.1为weak；其余为complement。最后一组包含宽接收矩阵的结构性右零空间和数值近零方向，不是“前几个强方向以外全部都叫严格零空间”。

主数据上Gin秩576，strong19/weak557/complement0；Gob秩48，strong9/weak39/complement528。补空间这里指range(A*)的右域正交补，不是输出空间range(A)的补。

各部分正交，故能量比例之和为1。普通范数比例为sqrt(η)，不能加成100%。逐案例表同时给分量RMS、平方能量和norm_pct。

传播后的分量能量有

$$\|AX_b\|_F^2=\sum_{i\in b}\sigma_i^2\|v_i^*X\|_2^2.$$

因此“误差幅度主要在补/弱空间”和“传播后的场误差主要由强方向承担”可以同时成立。所谓能量是离散平方Euclidean/Frobenius范数；未计算完整E/H电磁储能或Poynting通量。

## 3. 两种串联顺序：都能作100%能量表

对完整的Gin和Gob三分组，主表使用内部后接收：

$$X_{ab}^{I\to O}=P_{{\rm ob},b}P_{{\rm in},a}X.$$

对每个a，Σ_b||P_ob,b P_in,a X||²=||P_in,a X||²；再求和a，得到

$$\sum_{a,b}\|X_{ab}^{I\to O}\|_F^2=\|X\|_F^2.$$

反序 X_ab^(O→I)=P_in,a P_ob,b X也满足该式，但各格通常不同。所有向量分量之和仍等于X，**不同a,b的分量不必两两正交**。因此这是可复算的串联能量分配，不是证明存在一个共同强弱本征分类。

若把Gin或Gob再作用于这些分量，传播后不一定仍正交：

$$\|AX\|_F^2=\sum_{a,b}\|AX_{ab}\|_F^2+
2\sum_{(a,b)<(c,d)}\operatorname{Re}\langle AX_{ab},AX_{cd}\rangle.$$

代码保留有符号cross energy。不能将各格的孤立传播能量归一化后冒充实际场误差的唯一百分比分摊。

## 4. 真正交集：必须保留剩余项

令 A_a、B_b分别为两套正交子空间，Q_ab为A_a∩B_b的正交基。数值计算在A_a内求解

$$V_{{\rm ob},b^\perp}^*V_{{\rm in},a}c=0,\qquad Q_{ab}=V_{{\rm in},a}N_{ab}.$$

使用固定sin-angle容差1e-9，并保存1e-10/1e-8下的秩敏感性。不以cos≈1的舍入误差直接声称解析精确交集；这里是数值交集，无符号认证。

这些交集相互正交，但不必覆盖整个域。令

$$P_\cap=\sum_{a,b}Q_{ab}Q_{ab}^*,\qquad P_R=I-P_\cap.$$

则交集能量加剩余能量为100%。若省掉P_R，就会把本来无法用共同强弱标签分类的误差消失掉。本轮复域主定义中，交集秩分别为 weak-in/weak-ob20、weak-in/complement-ob509，其余0；共同分类域529维，剩余47维。

## 5. 旧TSOM的投影空间：先投影，再正交化

项目已有表达为

$$Z_{ab}=P_{{\rm ob},b}V_{{\rm in},a},\qquad
P_{Z_{ab}}=Z_{ab}Z_{ab}^{\dagger}.$$

这与上一节交集不同，也不同于把 P_ob,b P_in,a X 当作正交投影。计算先对 V_ob,b* V_in,a作SVD，再运输回Gob子空间，用相对原单位基的绝对阈值1e-10取rank；没有把极小投影残差相对自身归一化成虚假方向。

不同Z_ab的空间可以重叠，故||P_Zab X||²/||X||²之和允许超过100%，不能强行归一化。项目来源：R9 formula_collection/02_FIXED_OPERATORS_AND_PROJECTIONS.md 的 F-OP-005，以及本地a3_rom/basis.py的twofold_chart；未把它称为新TSOM实现或重新验证原2009论文。

旧实现的实际两fold chart还可直接写为

$$Z=(I-P_{{\rm ob},S})V_{{\rm in},S},\qquad U_T=[V_{{\rm ob},S},\operatorname{orth}(Z)].$$

两个块正交，其余为I−U_TU_T*。本例为9+19维，剩余548维。Z投到联合的Gob weak⊕complement，保持内部系数共享；分别对weak和complement投影后允许独立系数，一般会扩到更大的空间，不能混为同一个chart。

## 6. 为什么本例TSOM总百分比恰等于100%加交集剩余？

这是两个互补内部空间的标准投影谱结论，本轮直接用于解释结果，不主张原创。

本例Gin补空间为空，所以 A_weak=A_strong^⊥。固定某个Gob组B，在B内考察自伴算子 C=P_B P_A P_B，谱位于[0,1]。

- λ=1对应A∩B；λ=0对应A^⊥∩B。
- 0<λ<1的方向不属于任一交集，却同时属于range(P_B A)及range(P_B A^⊥)。

因此两个投影后空间的正交投影之和，加上两个交集投影之和，等于2P_B。再对Gob分组求和：

$$\sum_{a\in\{S,W\},b}P_{Z_{ab}}=I+P_R.$$

从而

$$\sum_{a,b}\frac{\|P_{Z_{ab}}X\|_F^2}{\|X\|_F^2}
=1+\frac{\|P_RX\|_F^2}{\|X\|_F^2}.$$

高对比主材料表的59.29%共同分类剩余，因而对应159.29%的重叠投影总量；电流为32.79%与132.79%。不是多出物理误差，也不是概率超过1，而是同一方向被两个投影后空间重复覆盖。

如果内部也有第三个非零组，该二互补空间恒等式不能直接照搬。

## 7. 聚合、阈值与解释边界

主表每个method/amplitude/noise有12个形状，使用

$$\eta_{\rm pooled}=\frac{\sum_t\|P e_t\|^2}{\sum_t\|e_t\|^2}.$$

不是先求各对象百分比再无说明地平均。CSV另有逐例、均值和既有PASS/FAIL条件表。不同初值/噪声/参数化结果不是独立材料抽样，本轮不报人群统计显著性。

.1是明确的计算定义，不是自然常数。预先声明.03/.3敏感性；另外从已保存模态系数给出.001/.01等完整阈值响应，不重跑物理，也不按GT挑最佳门槛。Gin在.001门下全部576方向都会被称strong，说明不能脱离阈值谈某个唯一“强方向百分比”。

G和S在所有材料与迭代中固定。本轮看到的是误差向量在这些固定空间中的位置变化；没有测量或声称G随对比度旋转。材料可辨识性、局部Jacobian几何、全局吸引域和GT唯一性仍是不同问题。

