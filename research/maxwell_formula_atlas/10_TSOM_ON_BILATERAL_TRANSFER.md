# 10 · TSOM 如何作用于双边 Maxwell 传递

2026-10-30，回应用户关于 $O(I-DC)^{-1}DE_R$ 的追问。本章为解析推导与算法定位，新增 Maxwell RHS=0；不改写 R10 已冻结实验。结论是：**可以使用 SOM/TSOM，材料与 Green 不对易并不构成否决理由。** 要区分绕开内层完整正演的电流辅助反演，与降低每次完整正演成本的近似算子求解。

实质文献依据为作者公开的 [Xu、Zhong、Song、Chen、Ran，MR-FFT-TSOM，DOI 10.1109/TGRS.2014.2329032](https://renchengsong.github.io/papers/TGRS2329032.pdf)，§II，式(3)—(11)。该文明确重述接收 SVD、内部传播基的投影及数据/状态双残差。以下符号吸收了原文电流中的频率常数；原论文算法、这里的噪声修正、variable projection 和低秩正演必须分开归属。没有声称已逐式读到2009原篇。

## F-BI-019 · 从双边逆算子还原电流方程

令 $N=256$，$n_s=8$，$h\in\mathbb C^{N\times n_s}$ 为固定宿主消元后的 ROI 对比电流，$y=Y-F_h$。定义 $K_D=I-DC$。当 $K_D$ 可逆时：

$$h=K_D^{-1}DE_R\quad\Longleftrightarrow\quad h=D(E_R+Ch),\qquad F-F_h=Oh.$$

这是精确等价式。原完整正演只需解8个右端，并不显式构造256×256逆矩阵或整个 $\mathcal T(D)$。SOM/TSOM 的入口是将 $h$ 作为辅助未知量，以固定投影组织它，同时满足上述状态方程与数据方程。它不要求 $D$ 与 $C$ 对易，也不要求真实材料解集线性。

## F-BI-020 · 两次筛选是投影像，不是交集

设 $S=WO$，其中 $W$ 为已固定的逐源噪声白化；$S=U\Sigma V^*$。保留噪声下可靠的前 $\ell$ 个方向：

$$h_d=V_\ell\Sigma_\ell^{-1}U_\ell^*Wy,\qquad P_\ell=V_\ell V_\ell^*.$$

令 $V_{C,r}$ 为固定内部算子 $C$ 的所选右奇异向量，或使用明确声明的 Fourier 基 $V_{F,r}$。TSOM式表示为

$$Z=(I-P_\ell)V_{C,r},\qquad Q=\operatorname{orth}(Z),\qquad h=h_d+Q\alpha.$$

这里 $\operatorname{orth}$ 删除相关/零列，保留相同列空间。$\operatorname{range}Z$ 是投影后的像，一般不等于两空间交集；$ZZ^\dagger$ 才是该像的正交投影。截断接收强方向后，$OQ$ 通常非零，不能省略数据残差。原作者的 twofold 内部基/Fourier构造支持这个位置；它不是从材料 Jacobian 的 weak8 方向构造电流。

**固定缓存：** $Q,CQ,OQ,V_\ell,CV_\ell$；同一观测的 $h_d,Ch_d$。宿主消元后的 $C$ 一般不再具有完整均匀网格的卷积形式，Fourier作为基仍合法，但不能未经推导声称 $C$ 乘法自动具有FFT卷积加速。

## F-BI-021 · 共享材料与含噪数据的联合目标

原 P1 材料保持 $D=D(x)=\operatorname{diag}(\delta(x)a(x))$，$x\in\mathbb R^{66}$。不得逐源恢复不同材料，或用自由像素代替原 decoder。为不把含噪接收强系数硬锁在 $h_d$ 上，可加入其修正：

$$R=[V_\ell,Q],\qquad h=h_d+R\beta,\qquad \beta\in\mathbb C^{q\times n_s}.$$

采用明确的候选生成目标（归一化权重需单独冻结）：

$$\mathcal J_\lambda(x,\beta)=\|S(h_d+R\beta)-Wy\|_F^2+\lambda\|K_{D(x)}(h_d+R\beta)-D(x)E_R\|_F^2+\mathcal R(x),\quad\lambda>0.$$

对全部源共享同一个 $D(x)$；$\mathcal R$ 可为原先验惩罚或事先固定的随机化目标。没有 $N$ 阶逆出现在目标计算中：$K_DR=R-D(CR)$，因此固定的 $CR$ 可预计算，材料更新主要是对角行缩放。

有限 $\lambda$ 允许状态残差，所以这是优化辅助目标，不等于原 Maxwell likelihood。作为材料可能性生成器完全可以使用；候选需由原完整正演验证。不要求先解决精确 posterior 才能研究这个生成器。若把接收强系数固定不动，则应将其标为传统截断数据近似，而非原含噪模型的精确等价。

## F-BI-022 · 电流系数的小规模消元

对固定 $D$，将上一目标的两项堆叠：

$$T_\lambda(D)=\begin{bmatrix}SR\\\sqrt\lambda K_DR\end{bmatrix},\qquad B_\lambda(D)=\begin{bmatrix}Wy-Sh_d\\\sqrt\lambda(DE_R-K_Dh_d)\end{bmatrix}.$$

若 $R$ 满列秩且 $K_D$ 可逆，则 $T_\lambda$ 满列秩：其下块已经满列秩。于是

$$\beta_*(D)=T_\lambda(D)^\dagger B_\lambda(D),\qquad \min_\beta\mathcal J_\lambda=\|(I-T_\lambda T_\lambda^\dagger)B_\lambda\|_F^2+\mathcal R(x).$$

采用薄 QR/SVD 求解，不应为了书写方便在实现中显式形成 $(T^*T)^{-1}$。若 $q\ll N$，主要内层费用约为 $O((N+n_r)q^2+q^3)$，多个源共享同一QR；不再需要每一步 $N\times N$ LU。若 $q$ 很大或需要大量恢复迭代，收益可能消失。

令 $r=T_\lambda\beta_*-B_\lambda$，正交性 $T_\lambda^*r=0$ 给出

$$d\|r\|_F^2=2\operatorname{Re}\operatorname{tr}\{r^*[(dT_\lambda)\beta_*-dB_\lambda]\}.$$

因此约化目标求导不用显式求 $d\beta_*$，但它仍是上述有限罚目标的导数，不是原 Maxwell posterior 的梯度。这是标准 variable projection；是否比精确256阶Full-J更好需要真正的同质量成本比较。

## F-BI-023 · 若目的是直接压缩中间逆：低秩 Green 与 Woodbury

另一条路线固定近似 $C_r=AB^*$，$A,B\in\mathbb C^{N\times r}$（例如截断SVD因子）。令 $g=DE_R$。在 $I-B^*DA$ 可逆时：

$$\widetilde h=(I-DC_r)^{-1}g=g+DA(I-B^*DA)^{-1}B^*g.$$

这是标准 Woodbury 恒等式，对近似算子 $C_r$ 精确。对每个新材料仍需更新 $B^*DA$ 并解 $r$ 阶系统，但不用重新分解 $N$ 阶矩阵。若 $C=C_r$ 是实际精确低秩，则也是原离散物理的精确表达；一般截断情况下只是近似。

这与 TSOM 有共同的固定传播压缩哲学，但算法对象不同：TSOM主要减少反演中未知电流系数；此式直接压缩正演求解器。不能用“把 $C$ 分解一次”推出“所有 $D$ 下的逆都已分解一次”。同样不要求 $D$ 与 $C$ 对易。

## F-BI-024 · 完整物理缺陷给出投影的失效位置

对任意近似电流 $\widetilde h$ 定义真实缺陷

$$r_D=DE_R-K_D\widetilde h.$$

当 $K_D$ 可逆，精确误差恒为

$$h-\widetilde h=K_D^{-1}r_D,\qquad F-\widetilde F=OK_D^{-1}r_D.$$

在 F-BI-023 的低秩算子求解中，进一步有

$$r_D=D(C-C_r)\widetilde h.$$

这说明被删传播的实际影响同时受材料、内部场和resolvent放大控制；只看固定Green的奇异值尾不够。若有稳定性界 $\|K_D^{-1}\|\le b_D^{-1}$，则 $\|W(F-\widetilde F)\|_F\le\|WO\|_2\|r_D\|_F/b_D$；无可靠界时用完整校正/验证。

对纯TSOM，缺陷可能来自电流子空间不足、有限罚优化未收敛或材料与电流失配，不能全部叫“流形不稳定”。误差超门时先考虑扩大固定电流基或精确校正，不必默认更新材料弱方向。

## 本次判断与未执行边界

1. **成立：** 双边式天然允许TSOM式固定电流表示；完整非线性仍在状态方程中，没有使用Born替代。
2. **成立：** 固定基缓存+小系数QR有实际避开内层完整LU的路径；这比只说“可能做SVD”更具体。
3. **尚未实测：** 原66维decoder内的上述TSOM/variable-projection候选生成器是否改善材料覆盖及总费用。本次没有执行，不将R10双边Gram随机目标的失败当作它的失败。
4. **已有但不同：** R10固定形状的rank64传递压缩与条件幅度积分得到有范围的正面计算结果；它不是这里完整的TSOM反演实验。
5. **下一项有辨别力的实验：** 共享同一精确ROI验证器和材料decoder，对照普通Full-J随机目标与固定TSOM电流系数+材料优化，按合格材料覆盖、遗漏/拒绝及完整费用比较。先不扩大弱8射线，不把准确posterior作为启动门。

本章数学结论属于已知状态分裂、投影、最小二乘消元和Woodbury的具体组合。潜在论文增量只能来自原材料支持下更好的合法可能性生成或总费用；不能将这些标准公式本身写成新定理。
