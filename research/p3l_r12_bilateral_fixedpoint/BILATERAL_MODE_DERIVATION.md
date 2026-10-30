# R12 双边材料模态：当前场、端口对与原先验坐标

本文件先解释物理，再给可执行的坐标定义。所有公式针对冻结的256阶ROI离散模型；不宣称连续PDE定理或全局材料纤维。新增Maxwell RHS=0、LU=0。标准公式的成立与R12探索是否更有用必须分开判断。

## 1. 两端怎样看见同一次材料改变

在当前材料 $D=\operatorname{diag}(d)$ 下，内部对比电流满足

$$h=D(E_R+Ch),\qquad F=F_h+Oh.$$

记 $M=I-DC$、$N=I-CD$，假设二者可逆。当前总场与接收传递分别为

$$R_D=E_R+Ch=N^{-1}E_R,\qquad L_D=OM^{-1}.$$

这里 $R_D$ 每一行告诉我们各个源怎样照亮该材料格点，$L_D$ 每一列告诉我们该格点的电流怎样传到全部接收端。两个端口都已经包含当前材料导致的多次散射。对电流方程求导：

$$M\,dh=(dD)(E_R+Ch)=(dD)R_D,$$

因此

$$\boxed{dF=L_D(dD)R_D.}$$

不需要 $C,D$ 对易，不需要 $D$ 可逆。对比幅度 $d_i$ 是复材料增量；若参数是占比 $a_i$ 且 $d_i=\delta a_i$，则相应列还应乘 $\delta$。把已经含 $\delta$ 的占比Jacobian再与材料导数相乘会重复计入该因子。

宿主Born式只在 $D=0$ 时取 $L_0=O,R_0=E_R$。以 $O,E_R$ 生成固定方向可以作为明确的消融臂，但不等于当前Full-J。反之，若已付费取得当前 $L_D,R_D$，后续Gram或SVD不再产生新的Maxwell右端。R12 `shared/physics.py` 在同一付费 `need_jac=True` 状态缓存二者。

## 2. Khatri–Rao列与双边Hadamard Gram

逐源、每源逐接收端向量化。复像素材料Jacobian $J_m$ 的第 $i$ 列为

$$j_i=R_{D,i,:}^{T}\otimes L_{D,:,i},\qquad J_m=R_D^T\odot L_D.$$

$\odot$ 在此专指逐列Kronecker积。于是

$$\boxed{G_m=J_m^*J_m=(L_D^*L_D)\circ\overline{R_DR_D^*}.}$$

逐项验证只需对接收索引、源索引分别求和。右侧第二项为 $\sum_s\overline{R_{is}}R_{js}$，因此必须是 $\overline{RR^*}=\overline R R^T$，不能误用 $RR^*$。每个因子均为Hermitian半正定，Schur乘积仍半正定。接收单边 $L^*L$ 忽略照明；源单边 $\overline{RR^*}$ 忽略可观测性；任一个单独使用都不是原数据Hessian。

## 3. 真正的端口对模态

取薄SVD

$$L=U_L\Sigma_LV_L^*,\qquad R=U_R\Sigma_RV_R^*.$$

对第 $a$ 个接收模态、第 $b$ 个源模态，定义材料行

$$K_{(a,b),i}=\sigma_{L,a}\sigma_{R,b}\overline{(V_L)_{ia}}(U_R)_{ib}.$$

由

$$L\operatorname{diag}(v)R=U_L[\Sigma_LV_L^*\operatorname{diag}(v)U_R\Sigma_R]V_R^*$$

可知完整 $K$ 与原数据导数仅差输出端的等距酉变换：$K^*K=G_m$。在32接收、8源时最多有256个对，而不是把两个单边弱子空间取集合交集。保留乘积 $\sigma_{L,a}\sigma_{R,b}$ 最大的32个对，是R12声明的截断启发式；乘积强度没有包含当前decoder的可达性，因此不保证它是66维先验上最优的32项截断。

`agent_a/modes.py:pair_modes` 实现这些加权行。`fullGram` 保留全部双边信息；`pair` 默认只保留top32对。两者是计算表示/截断的区别，不是新增物理信息。

## 4. 原Gaussianized66先验的正确链式法则

原native先验坐标为 $(z,u)\in\mathbb R^{64}\times[0,1]^2$。R12实际优化使用 $x=(z,\xi)\in\mathbb R^{66}$，其中 $u_j=\Phi_N(\xi_j)$。该变量变换使声明的参考概率为标准Gaussian66；它与native Uniform2坐标下的MAP点目标不是同一个点估计。

设 $s=(\mu+A z)/\tau$，$a=S(s)$ 为冻结的quintic saturation，

$$\delta=0.8+\Phi_N(\xi_1)+i[0.03+0.15\Phi_N(\xi_2)],\quad d=\delta a.$$

真实复材料导数 $D_x\in\mathbb C^{256\times66}$ 为

$$D_{x,:,1:64}=\delta\operatorname{diag}(S'(s)/\tau)A_{ROI},$$
$$D_{x,:,65}=a\,\varphi_N(\xi_1),\qquad D_{x,:,66}=i\,0.15a\,\varphi_N(\xi_2).$$

$\Phi_N$ 是CDF，$\varphi_N$ 是PDF。两列属性不能漏掉CDF导数。平台上的 $S'$ 可为零，所以材料像的局部秩可低于66。

冻结白化为每个源依次排列32个实部、32个虚部，整体乘 $\sqrt2/\sigma$。对实参数 $x$：

$$J_w=W_{real}(J_mD_x),\qquad\boxed{H_{data}=J_w^TJ_w={2\over\sigma^2}\operatorname{Re}(D_x^*G_mD_x).}$$

此处不能把复矩阵 $D_x^*G_mD_x$ 直接作为实Hessian。原Gaussian先验再给 $I_{66}$，但本轮方向消融的“数据可见性”不把 $I$ 混入；预测/校正的正则化应由主线程协议单独声明。

## 5. 材料度量、局部decoder零空间与有限弱坐标

原材料RMS是全576格点的复差范数除24。ROI外材料固定，故即使只保存256行导数，也必须使用

$$G_x={1\over576}\operatorname{Re}(D_x^*D_x).$$

令 $B=[\operatorname{Re}D_x;\operatorname{Im}D_x]/24=U\Sigma V^T$。删除相对奇异值阈值 $10^{-10}$ 以下的局部零方向，取

$$T=V_+\Sigma_+^{-1},\quad T^TG_xT=I,$$

然后分解 $T^TH_{data}T$，保留最小特征值的 $k$ 个方向 $W=TU_{weak}$。于是 $W^TG_xW=I$。它们是单位材料RMS的低数据代价方向；不能把decoder不改变材料的latent零方向当作新材料世界。阈值定义的是有限精度局部商，不证明全局gauge或全局可积纤维。对于latent范数极大的方向，实际原prior代价也可能很高，必须在真实候选中记录。

原生自由占比加未知幅度模型存在 $a\mapsto a/c,\delta\mapsto c\delta$ 的同材料尺度gauge；原Gaussian64 decoder没有因此自动获得全局重缩放gauge。R12主臂使用原66维decoder，不能借用原生258模型的gauge解释。

对一个预测/校正步骤，把锚点 $G_0,W_0$ 固定，并定义

$$C_w=W_0^TG_0,\quad c(x)=C_w(x-x_0).$$

`coordinate_preserving_complement` 返回Euclidean正交 $Z$ 张成 $\ker C_w$。校正 $dx=Z\alpha$ 精确保留上述有限仿射坐标。若每次校正都更新 $G,W$，只能得到移动局部规则，不能声称保留同一个有限坐标。若需只在当前材料像内校正，还应与其明确声明的切空间参数化联合处理；$Z$ 本身包含局部decoder零方向。

## 6. 实验可以支持什么

完整当前双边Gram与直接Full-J具有相同局部数据信息；二者不应作为两种独立科学发现。固定宿主双边、接收单边、源单边与top32 pair截断是可比较的提案几何消融。只有按真实原decoder、真实非线性校正、真实观测噪声门、有限材料差异与完整费用评估，才能判断某个方向更有用。小数点精度的恒等式重放不构成材料覆盖、纤维存在或新颖性证据。

## 7. 主线程Addendum2：原先验感知方向

在前述纯数据审查发现高latent代价后，主线程明确选择原Gaussian先验感知几何。数据因子记为 $A_{data}=J_wT$（各消融用自己声明的因子），则

$$H_{post,local}=H_{data}+I,\qquad A_{aug}=\begin{bmatrix}A_{data}\\T\end{bmatrix},\qquad A_{aug}^TA_{aug}=T^T(H_{data}+I)T.$$

`full_prior,pair_prior,receiver_prior,source_prior,host_prior` 由这个增广因子的SVD选择weak16。它加入的是原Gaussianized66先验曲率，**没有新增观测信息**，也不是把先验改为材料Lebesgue度量。主线程冻结H1/H2/HA使用这一声明的几何；原纯数据结果仍保留为有限步可行性负证据。

为避免正规方程在最弱方向上消减精度，实现直接分解增广实因子，而不先形成latent Hessian再乘病态canonical lift。完整双边和Full-J的等价仍成立，先验项在二者中相同。单边或pair截断加入相同先验，也不能追回已丢失的观测信息。

`agent_a/DECODER_WEAK_FEASIBILITY.json` 对8个暴露R9端点、每个几何16列及最节省latent范数的弱空间组合、两种符号，在latent位移4处直接调用原decoder。纯数据272条被测射线中12条达到材料RMS0.02，原先验感知各臂均为272/272；这些结果没有使用Maxwell。少数纯数据组合经曲率可达到较大有限材料差，因此不能据此证明整个纯数据weak16信赖球不可行。这里只证明原先验感知几何在声明的有限射线集合中更容易产生足够材料移动，后续仍需真正的非线性观测校正。
