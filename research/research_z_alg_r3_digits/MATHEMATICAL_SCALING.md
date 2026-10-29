# 从 12 到 576 材料：精确小 Gram 目标与伴随梯度

本文件保留决定实现的推导。模型为固定几何、固定频率的有限二维 scalar-TM disk-pulse VIE；材料为实数对比度，576 个物理单元全部未知。它不是 Born、连续 PDE 认证或三维矢量 Maxwell。下述矩阵求解均包含完整多次散射。

## 1. 统一物理变量

采用源为列的约定，$G\in\mathbb C^{N\times N}$、$E_0\in\mathbb C^{N\times s}$、$S\in\mathbb C^{r\times N}$ 固定。令

$$D=\operatorname{diag}(m),\quad M=I-DG,\quad B=I-GD,$$
$$J=M^{-1}DE_0,\quad Q=E_0+GJ=B^{-1}E_0,\quad F=SJ=SDQ.$$

假设 $M$ 可逆；由 $\det M=\det B$，$B$ 同时可逆。材料可以高度非均匀，不假设 $D$ 与 $G$ 对易。$Q$ 是完整总场，$J=DQ$ 是对比电流。

## 2. 继承 R3 的精确消元，改变计算方式

R3 用接收对偶矩阵 $W$ 满足 $WE_0=Y$，并最小化共同材料兼容性。其消元目标已严格推导为

$$\psi(m)=\operatorname{tr}\{R H^{-1}R^*\},\qquad R=F(m)-Y,\quad H=Q^*Q.\tag{D1}$$

适用条件是 $E_0$ 满列秩、物理系统可逆。因此 $Q$ 满列秩，$H$ 正定。该式与 R3 硬数据 variable projection **严格相同**，不是额外的物理降阶。两种形式的标量归一化也必须一致后再比较。

旧实现对 $A^T=[Z(I-GD)]^T\in\mathbb C^{N\times(N-s)}$ 做 QR，$ZE_0=0$。本轮 $N=576,s=12$，每次是 $576\times564$ 的分解。D1 改为一次完整物理解及 $576\times12$ 总场的薄 QR，避免该宽辅助矩阵。它**没有免除 576 阶 LU**。

若 $Q=UT$ 是满列秩薄 QR，则

$$H=T^*T,\qquad \psi=\|RT^{-1}\|_F^2.\tag{D2}$$

实现应使用三角求解，不显式求 $H^{-1}$。当源列条件尚可时也可对 $12\times12$ Gram 做 Cholesky；形成 Gram 会平方 $Q$ 的条件数，因此须报告病态状态与回退。增加未声明的对角 jitter 会改变目标，不能再声称严格等价。

## 3. 不构建 Full-J 的精确梯度

令

$$A=RH^{-1}\in\mathbb C^{r\times s},\qquad C=A^*A\in\mathbb C^{s\times s}.$$

对实材料微分，$dH=dQ^*Q+Q^*dQ$，于是

$$d\psi=2\operatorname{Re}\langle A,dF\rangle_F
 -2\operatorname{Re}\langle QC,dQ\rangle_F.\tag{D3}$$

这里利用 $dH^{-1}=-H^{-1}(dH)H^{-1}$；第二项是变化的材料度量的导数，**不能冻结 Gram 而省略**。

从电流方程微分得到

$$M\,dJ=\operatorname{diag}(dm)Q,\quad dF=S\,dJ,\quad dQ=G\,dJ.$$

用已有 LU 解伴随系统

$$M^*\Lambda=S^*A-G^*(QC).\tag{D4}$$

代入 D3 即得

$$\boxed{\frac{\partial\psi}{\partial m_i}
=2\operatorname{Re}\sum_{a=1}^{s}\overline{\Lambda_{ia}}Q_{ia}.}\tag{D5}$$

对实材料这是实梯度；复共轭转置不能与普通转置混用。若目标乘一个固定归一化常数，梯度同乘该常数；若常数也依赖材料，必须再求其导数。本轮不采用随迭代变化而未求导的归一化。

**算法后果：** 目标加梯度一次只需要一份 $N$ 阶 LU、$s$ 个源 RHS、$s$ 个伴随 RHS、小型 Gram/薄 QR 和矩阵乘法。没有逐材料有限差分，也不需要 $N-s$ 个辅助投影变量的每步 QR。

## 4. 与普通方法的真实比较

原数据损失 $\ell=\|F-Y\|_F^2$ 的伴随梯度是

$$M^*\Lambda_d=S^*R,\qquad
\partial_{m_i}\ell=2\operatorname{Re}\sum_a\overline{\Lambda_{d,ia}}Q_{ia}.\tag{D6}$$

所以最优普通伴随梯度同样只需 $s$ 个伴随 RHS。**小 Gram 不构成对一切普通反演的梯度成本优势定理。** GlobalDual 的差异首先是优化的非零残差几何，不是魔法减少物理维数。

普通 Full-J/GN 则通过

$$T_s=SM^{-1},\qquad
\frac{\partial F_{ra}}{\partial m_i}=(T_s)_{ri}Q_{ia}\tag{D7}$$

构造解析 Jacobian。$T_s$ 可用一个转置求解的 $r$ 个 RHS 得到；本合同 $r=48<N$。实化后 Jacobian 为 $2rs\times N=1152\times576$。GN 还需要形成或作用于法方程，以及材料空间步的线性求解；这些时间必须记账。若改用 matrix-free GN 或 L-BFGS，该费用会变化，不能只凭 RHS 数断言快慢。

倒数参数 $\eta_i=1/m_i$（本合同 $m_i>0$）满足

$$\frac{\partial F}{\partial\eta_i}=-m_i^2\frac{\partial F}{\partial m_i}.\tag{D8}$$

这改变优化坐标，不改变材料范围、物理前向或允许的图像。

## 5. 多尺度是参数化，不是 Maxwell 秩降阶

每阶段用已声明的确定性 prolongation $P_h$ 将 $h\times h$ 材料参数映射到同一个 24×24 物理网格：$m=P_hz$。Native 坐标的链式导数为

$$\nabla_z f=P_h^T\nabla_m f,\qquad J_z=J_mP_h.$$

若使用正权重行和为 1 的插值，则参数 box 的像仍位于同一物理材料 box。最终 $h=24$ 时 $P_h=I$，576 个物理材料自由度全部释放。粗阶段少的是优化未知量；$\operatorname{diag}(P_h\,dz)$ 一般仍可满秩，因此完整 Maxwell 系统始终是 576 阶。

正则项在物理图像上计算，并对所有方法使用同一公式、权重安排、材料范围和初值权限。若 coarse reciprocal 使用 $m=1/(P_h\eta)$，其粗材料空间与 $m=P_hz$ 一般不同，必须改成相同物理图像映射或明确申报差异，不能偷偷把两者当作相同先验。

## 6. 条件、极点与最终验收

D1 与原数据目标在可逆域中有相同精确零点，但通常有不同非零极小点。且

$$\frac{\|R\|_F^2}{\sigma_{\max}(Q)^2}
\le\psi\le
\frac{\|R\|_F^2}{\sigma_{\min}(Q)^2}.$$

接近极点时，$\|Q\|$ 增大可令兼容性目标虚降；R3 已给出严格反例。因此候选在线选择与最终验收必须回到所有源接收的原始数据损失及共同正则项，并检查材料合法与物理解的数值稳定性。小 $\psi$ 本身不是恢复材料的证据。

独立端点重求解只验证同一有限离散 VIE 的实现和保留预测，不是跨网格或独立连续物理模型验证。数值相符不证明唯一性；观测数 1152 大于材料数 576 也不证明 Jacobian 满秩。

### 本轮实现的归一化与共同先验

对每幅观测，$b=\|Y\|_F^2$。原数据评分为 $\ell_d=\|R\|_F^2/(2b)$。兼容性候选使用固定参考 $m_{ref}=\mathbf1$ 的总场，令 $a=\|Q_{ref}\|_F^2/s$，优化 $a\psi/(2b)+\mathcal R(m)$。$a$ 在本次求解内固定；即使首次兼容性阶段出现在24×24，它也不由变化的当前材料重设。这个尺度让平均总场能量可比，不把各个源方向的度量变成原数据目标。

在物理网格边上定义 $d_e=(m_j-m_i)/19.5$，实现共同正则项

$$\mathcal R(m)=\lambda_2\,\operatorname{mean}_e d_e^2+
\lambda_{TV}\,\operatorname{mean}_e\left(\sqrt{d_e^2+\epsilon^2}-\epsilon\right).$$

它是各向分边的平滑 TV 与二次平滑，不是数字模板。Native、reciprocal、GlobalDual 的在线排序全部使用 $\ell_d+\mathcal R$。兼容性方向和该排序的下降方向可以不同；每个候选仍保留原物理评分。原数据校正采用有界阻尼 GN 和标准滞后扩散/IRLS 二次模型，不把其正则化 Hessian 近似说成精确二阶导数。

在原图差分矩阵 $D_e$ 下，这个每步更新的正则曲率是

$$H_{\mathcal R}^{IRLS}=D_e^T\operatorname{diag}\!\left[
\frac{2\lambda_2+\lambda_{TV}/\sqrt{d_e^2+\epsilon^2}}{n_{edge}}
\right]D_e.$$

它保留相邻网格的固定稀疏结构。最终 $P=I$ 时，可直接稀疏组装再加到数据法方程；无须用稠密单位基反复做矩阵乘法。权重仍随材料改变，既没有冻结非线性曲率，也没有假设材料与传播算子交换。这是标准 IRLS 与稀疏代数的工程复用，不能称为新的物理近似或新定理。

## 7. 来源和原创性边界

D1 直接继承 R3 的全观测消元；QR/Gram、隐式微分、伴随和多尺度是标准工具。R3 已明确其与 CSI/SOM、variable projection 和材料相关 PDE 残差度量的直接关系。本轮首先是一项 576 材料盲成像的规模验证。是否出现新的算法证据，必须由所有九数字的冻结比较决定；不能仅凭把宽 QR 换成薄总场 QR 宣称新的 Maxwell 代数定理。

材料相关 Gram 残差度量的直接理论先例可参见 van Leeuwen 与 Yang 的[作者全文](https://arxiv.org/html/2403.15292v2)。该文处理 PDE 约束松弛并推导加权非线性最小二乘及解的 Gram 权重。本轮固定离散电流、硬数据消元和实材料梯度须按上述条件逐项推导；记号和约束不同不构成原创性的充分理由。

## 8. 这个目标怎样改变图像搜索：材料相关的照明归一化

这条解释不需要增加任何物理实验。对当前合法材料的总场作薄 SVD，$Q=U\Sigma V^*$，则

$$\boxed{\psi(m)=\|R V\Sigma^{-1}\|_F^2
=\sum_{a=1}^s\frac{\|R v_a\|_2^2}{\sigma_a(Q)^2}.}\tag{D9}$$

固定材料时，Maxwell 方程对照明线性，因此源组合 $E_0V\Sigma^{-1}$ 的总场恰为 $U$，内部场列两两正交且范数为 1。GlobalDual 相当于比较这些**按当前内部总场归一化的照明组合**的观测误差。普通数据目标按原始照明幅度计量；两者非零残差的几何自然不同。

这里的 $V,\Sigma$ 随材料变化，不是一次固定 SVD 解决所有材料，更不是额外测量提供新信息。D3–D5 已包含这种变化的完整导数。该重加权可以改变优化吸引域，但不保证改变方向有利，也不保证跨越所有局部极小点。

若观测噪声每个接收行在源坐标中是独立、等方差复噪声，方差为 $\sigma_\eta^2$，则固定材料下

$$\mathbb E\|\eta V\Sigma^{-1}\|_F^2
=r\sigma_\eta^2\operatorname{tr}(H^{-1}).\tag{D10}$$

所以小场奇异值方向的噪声会被放大；不能把该目标直接当作原噪声 likelihood。这个解析结果说明为什么无噪声成功不能自动推广为噪声鲁棒性。它不是本轮已完成的带噪实验证据。

最后，GlobalDual 提案由 L-BFGS-B 求解，普通主基线是 GN/LM；因此“GlobalDual 结果更好”本身不能区分**目标变化**与**优化器更换**的作用。本轮在任何主 GT 揭晓前另冻结普通数据伴随 L-BFGS-B＋相同 Full-J 校正的匹配对照，保留该补充在数据生成之后才冻结的边界。只有相对于这个对照的额外优势，才支持 Gram 目标具有独立搜索增量；仍不能据此证明普遍因果规律或全局收敛。
