# B1：材料事件概率与 Maxwell 成像

2026-10-13。本文解释本项目实际使用的数学对象。已有 A17/A31/A32 结论保持原归属；以下 Bayes、importance、链式条件概率和插值稳定性不是新颖性声明。

## 1. 先说明图像所代表的世界

散射测量观察的是物体对入射电磁波的响应。一个实心球、一个有中心空腔的球和一个具有中心空腔及封闭球壳空隙的球，可以产生相近响应。算法应比较这些完整世界的概率质量，而不能只比较三张图的拟合残差。

设 H∈{0,1,2}，分别对应上述三种固定同心几何；外半径 R=1。H1的空腔为0≤r≤.50；H2的空腔为0≤r≤.24与.46≤r≤.60。第二个空腔是一个三维连通的球壳形真空区，不是两个任意位置的小球孔洞。二维图只是三维材料的中心切面。所有介质区域共享 ε=ε′+iε″；真空 ε=1。

每分支的 ε′∼Uniform[2.2,2.8]，ε″∼Uniform[.04,.16]，相互独立；初始 πH=1/3。连续材料精确单点的概率为零。本文给某张代表图标注的概率是整个H分支的质量，而不是那张图本身的质量。CLI图像代表是分支条件最大likelihood的近似材料；主要科学图/面板用条件均值ε作为合法代表，整体后验均值与可信区间另外报告。二者均不改变分支事件质量。

这是一项封闭的三模板模型。它给模板外材料零先验支持，并不能证明外面的真实物体不可能存在。更大材料世界中的posterior mass尚未定义；JSON为其保留null/未量化状态。Top2展示保留第三模板真实质量，绝不把两项重新归一化。

## 2. 三维矢量物理与概率目标

采用exp(−iωt)，μ=1，外部真空，kR=1.6。分层球展开同时保留电/TM系数a_n与磁/TE系数b_n，默认nmax=10。角度20°至160°，8个角度，两种正交入射偏振/散射平面投影，得到16个复远场振幅。已知接收距离相位被除去，A=i(S1,S2)/k；这不是二维scalar TM。

主实现逐多极匹配切向E/H。以适当归一化的径向Riccati函数u，TM跨界面连续的是u及u′/m，TE连续的是u及m u′，m=√ε；外部u=ψ_n−a_nξ_n或ψ_n−b_nξ_n。内核只有regular ψ，壳层同时允许ψ与第二类函数。每个polarization解2L阶边界系统，列/行缩放是可逆换元，不是物理jitter。独立admittance递推、未改写A13均匀球系数、A15已有treams多层系数及能量/极限核验分别约束实现。

复数据按全部Re后全部Im实化为y∈R32。噪声σ是实验固定的每实分量标准差：σc=√2σ为proper complex标准差。对F_H(θ)，

\[
L_H(θ;y)=(2πσ^2)^{-16}\exp[-\|F_H(θ)-y\|^2/(2σ^2)].
\]

分支prior density为1/(.6×.12)。证据与后验：

\[
Z_H=\int L_H(θ;y)\,\mu_H(dθ),\quad
p_H=\frac{π_HZ_H}{\sum_Gπ_GZ_G},\quad
\pi(dθ\mid H,y)=L_H(θ;y)\mu_H(dθ)/Z_H.
\]

Gauss–Legendre节点权重是在连续proper prior上的积分权重，不是一个被偷换成有限粒子prior的概率模型；3/5/7阶是数值近似。绝对log evidence保留Gaussian归一化，Bayes factor在相同噪声模型中约去共同常数。

带校准增益时y=gF_H(θ)+η，g∼N(1,s_g²)。边缘化得到C_H(θ)=σ²I+s_g²F_HF_Hᵀ。likelihood包含−½logdet C_H；不能只保留加权residual。协方差依赖候选，A32的共同固定白化签名不能不经修订用于该边缘目标。保存/重放观测必须携带s_g；本项目已经补齐NPZ元数据并保留遗漏导致的错误重放。

## 3. “材料token”的严格联系

假设一个proper prior树，节点s是已经生成的离散前缀；后续结构、几何、材料和nuisance具有明确条件prior。令Z(s;y)是该前缀的全部连续/离散后续likelihood积分，则

\[
Z(s;y)=\sum_t p_0(t\mid s)Z(st;y),\qquad
p(t\mid s,y)=\frac{p_0(t\mid s)Z(st;y)}{Z(s;y)}.
\]

沿路径相乘会telescoping，得到叶节点prior mass×叶证据/根证据。这给条件材料概率与完整Bayesian posterior的联系；网络softmax只有在估计了这些条件量并经过合法修正时才能承载这种含义。

计数、连接关系、封闭空腔出生等适合作为离散节点。半径、边界坐标和复本构变量应优先连续积分/采样；把它们tokenize只是另一种离散模型，需要明确精度和测度。Beam ranking不会计算被裁剪的子树证据，也不能把保留下来的路径质量自动归一化为1。本原型只实现一层三分支树，不声称已完成任意可变维结构生成。

## 4. 两条实际算法及其修正

主线是branch-wise quadrature。先枚举所有注册H，再积分θ，所以列出的拓扑不会因proposal未访问而消失。它仍可能遗漏模板外拓扑及分支内未解析的窄峰；枚举三个事件不是全世界mode discovery。

对照是每分支16个独立样本的防御importance。把θ线性变换到z∈[0,1]²；条件prior在z上为1。proposal为

\[
q_H(z)=.5+.5\,Beta(z_1;a_1,b_1)Beta(z_2;a_2,b_2).
\]

Beta参数只从order3后验矩估计；q≥.5，保留完整注册支持。每分支证据估计为N⁻¹Σ L/q，再乘π_H并归一化。θ与z的常数Jacobian在prior/proposal ratio中正确约掉。证据均值在此独立采样条件下无偏；三个随机证据的归一化ratio一般有有限N偏差。Bootstrap不认证未采到的窄峰；高ESS也不证明准确。

48次样本全波调用与27次proposal准备均计费。当前单案例最大概率差.02565，超过.02；因此没有用它宣称比积分更准。不存在NN、SMC链或RJ接受率被冒称已执行。

## 5. 场近似与概率积分是两个误差

若F_A=F+e，则势能差含residual–error干涉项。对共同σ，

\[
|\Phi_A-\Phi_F|\le\|F-y\|\|e\|/σ^2+\|e\|^2/(2σ^2).
\]

场光滑不意味着likelihood也容易积分；当σ降低，likelihood在材料坐标上可以迅速变窄。第一轮强信号ΔP虽仅约1e−9，低质量分支的logBF差仍.366，不能把概率表面稳定当完整证据门槛通过。

第二轮把已支付5/7阶场网格拟合为tensor-Legendre多项式，再对其likelihood做17/25阶积分。这些增密是廉价代数evaluation，逐项另记，绝不是新增全波状态。它定义了一个可计算的近似forward target；通过24个新full点和4个已有独立treams评估核验近似。经验检查不能提供整个prior支持上的uniform remainder，故标EMPIRICALLY_REFINED_REFERENCE。原生full-grid的REFERENCE_UNRESOLVED保持原样。

若未来能在全部支持上界定likelihood tilt δ的oscillation ω，已有bounded-tilt理论给TV≤tanh(ω/4)，从而界定每个事件概率误差。当前没有取得该全域界，不能将24点最大误差代入作为证书。

## 6. A17/A30/A32的实际连接与限制

A17在冻结Galerkin current改动下给ΔJ=YΓ⁻¹N，并严格保留shared real material/realification。A32共同proper prior与残差r下，两个affine模型的势能差为

\[
h=(Y^Tr+Y^TJ_0s)^TΓ^{-1}Ns+\tfrac12\|YΓ^{-1}Ns\|^2.
\]

所以dπ₁/dπ₀∝exp(−h)，并且给定双侧签名U=(Ns,Y^TJ₀s)的条件分布不变。事件变化等于Covπ₀(1_H,exp(−h))/Eπ₀exp(−h)。这是一种具体物理因子与一般Bayes tilt的连接，不是完整posterior天然低维。

本项目比较实际Mie的central-FD affine Jacobian和dipole截断的affine Jacobian，以SVD构造同形因子。它们不是A17 native Galerkin交换。三分支defect均rank2，签名最多4维，而材料只有2维，没有压缩收益。跨H使用不同anchor时必须把H和各branch normalization保留；此时H在签名中可知，不能宣称τ_H小。

非线性full-versus-affine剩余在错配案例的后验加权RMS可达3.62，不能用affine恒等式直接替换full likelihood。竞争例full和calibrated-L最大概率差.01232，信息不足例仅.000232；这是冻结模型诊断，不是部署加速资格。

A30 material observability在此用于prior尺度归一化的JᵀJ诊断，不把局部曲率逆当非高斯/多拓扑posterior。固定材料两点的full finite secant可计算白化距离与等prior二点误判率Φ(−d/2)；复合拓扑仍必须积分，不能把两张MAP的距离冒充复合Bayes error。

## 7. 不确定性与停止条件

空腔occupancy p_void(r)=Σp_H1_void,H(r)，反映这些固定模板的事件质量；中心已知，不能称位置估计。复材料SD通过全方差分解合并分支内Re/Im方差及分支均值差，不是独立像素noise或simultaneous可信带。

按真实H条件的材料coverage与H-marginal coverage不同。原来一项校准标签混用了名称，已修正并保留初版。24个连续材料parent各30个噪声，仅24个独立材料簇；bootstrap必须在parent内保留整组replicates。90%离散拓扑credible set可以因离散性包含97%以上质量，故它的coverage应同时比较实际报告mass，不能机械要求恰好.90。

Brier/log/coverage是模型内诊断；若没有充分参考、独立材料覆盖和SBC rank检查，不能宣称真实世界校准。部分分箱cluster区间不含0的负证据也必须报告。当前预算384/384，后续真实物理、未知位置/边界、复杂形状、严格认证和跨实例优势均停止在未授权的新阶段边界。
