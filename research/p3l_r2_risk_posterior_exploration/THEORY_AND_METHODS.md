# 风险目标、概率目标与几何机制

本文件区分数学性质与待检验算法。统计目标始终是原 P1；所有 FULL 指原 24×24 二维 TM/VIE 的完整离散模型。

## 材料 MSE 监督不需要 MCMC teacher

设材料在复数 Hilbert 空间中，内积为全网格平均实内积。令 μ(y)=E[χ|y]。对任意只依赖观测的估计 T(y)，条件期望的正交性给出

`R(T) = E||χ−μ(y)||² + E||T(y)−μ(y)||²`。

因此 joint simulation 下的实际材料 MSE 是合法的条件均值学习目标；需要独立 parent 的泛化评价，不要求先求出每份观测的 posterior。有限训练数据、模型容量和优化不保证得到 μ，也不保证可信区间或完整概率密度。原 P3-L 的失败 Gaussian 平均不构成此恒等式的反例。

若 q0、J 和残差都是 y 及已知 solver 的确定性函数，它们不增加新的观测信息。它们可能改进归纳偏好和有限样本/有限算力下的学习效率；独立增量必须通过同输入信息、容量与完整成本的消融检验。不能把多做一次 Jacobian 称作取得了新传感信息。

## 材料风险与弱方向不是同一排序

`D=Dg(q0)` 给出材料速度，`J` 是正确白化的 Maxwell 导数。任务尺度下可考虑 `(||Jv||²+||v_z||²)/||Dv||²_material`。分母核或极小奇异值必须单独处理，不能把 latent 饱和当真实材料自由度。该广义谱是局部候选生成的工具，不是 posterior 的新参考测度；不插入任意 Gram determinant 改变原 prior。

在局部线性 Gaussian 近似中，材料风险贡献涉及 `D C D*`，而非单独 J 的最小奇异值。真实问题的 C 未知；本轮分别比较观测几何、监督剩余材料风险与普通不分割估计，不能把一个 GN covariance 自动认作真实 C。

有限移动还包含 decoder secant 余项和 Maxwell response 余项。`Dg δq` 很小不排除有限 δq 跨越饱和界后产生大材料变化。轨迹记录实际 RMS `.02/.04`、signed progress、原 prior cost、完整 Φ 与谱旋转，不用参数距离代替材料运动。

## 概率合同与新 non-Gaussian 提案

原 q=(z64,u2)，z 标准 Gaussian，u 在单位方框均匀；`Φ=||F−y||²/σ²`，σ固定 .00119211798514519，实虚噪声方差各 σ²/2。在 `u=NormalCDF(η)` 下 x=(z,η) 的原 prior 恰为标准 Gaussian66。

本轮多中心 Student-t 非线性三角 proposal 采用 `η=cη+s tη` 与 `z=cz+L tz+B(NormalCDF(η)−NormalCDF(cη))`。其 Jacobian 为块三角，绝对 determinant 是 `det(L) prod(s)`；非线性 shear 不改变 determinant。每个完整支撑的 mixture density 可明确求值。独立 MH 使用真实 FULL likelihood 及两端 proposal ratio；数学不变性不保证有限时间 mixing。

另一条路线按八个 illumination 分块渐进同化；每个阶段保留累计 likelihood、增量权重、祖先、重采样和合法 mutation。终点必须真实到达全部八源原目标，未完成时只报告当时的 tempered/partial target。

固定 z 后对原 u 方框积分属于 **SLICE**，其 du 测度与 FULL marginal 不同。它可以独立验证 non-Gaussian kernel、normalization 与材料矩，不会使 66维参考自动合格。

## 轨迹与风险策略

零宽轨迹没有连续多维 posterior 的普通体积。几何轨迹本轮先明确作为 candidate-only 或监督风险策略的动作环境；其访问频率不作 posterior mass。上下文 bandit 直接拟合允许的 TRAIN GT 风险改善；有限时域 model-based policy 同时拟合下一步在线特征变化，规划未来风险回报。二者均在未读取测试 GT 时冻结，并须实际运行独立在线 rollout。稀疏行为覆盖和未来局部几何的近似不构成离策略价值证明。

## 先例及证据范围

- [Cui 等，2014，LIS](https://doi.org/10.1088/0266-5611/30/11/114015)：prior-relative 子空间和 Rao–Blackwell 已有成熟先例，不能宣称首次强弱分解。
- [Adler 等，Task adapted reconstruction](https://doi.org/10.1088/1361-6420/ac28ec)：监督风险与条件均值成像已有理论基础。Semantic Scholar 的 arXiv 年份为2018，期刊发表信息另按 DOI 记录，不混作同一日期。
- [Greenberg 等，2019，APT](https://arxiv.org/abs/1905.07488)：联合模拟下的归一化 density learning 与 proposal correction 有明确先例。弱坐标预测正确不等于 FULL posterior 正确。

三条记录已通过 ScholarQA 批量核验；本次发现查询部分429，保留原始返回。上述定位不是穷尽新颖性审查。独立 Maxwell 算法贡献取决于真实风险/统计精度及总成本的结果。

## 实际风险动作的数学范围

本轮选定的网络输出先经过一个包含全部 P1 材料的凸矩形投影。令 C 为 ROI 内实部/虚部合法范围的直积，ROI 外固定为原背景。对任意合法真值 χ∈C，Hilbert 投影性质给出

`||P_C(T)−χ||² ≤ ||T−χ||² − ||T−P_C(T)||²`。

所以这个投影不会增加任何合法 P1 真值下的材料平方误差。这是一个逐材料的数学性质，不是网络泛化保证。C 严格大于整个非线性 decoder 的材料 image；逐像素落在 C 内不能证明存在单个合法 q 产生它，也不能证明该图像与观测相容。网络输出以及若干轨迹的平均图统一称为 **材料风险动作**，不称为发现的物理材料世界。

对非线性 Maxwell 算子，`F(E[χ|y])` 通常不等于 `E[F(χ)|y]`。因此本轮对风险动作另做真实 forward 检查，其 residual 是部署诊断；它既不等于 posterior predictive check，也不能单独否定该动作的材料风险价值。真正 posterior predictive 需要一份有资格的材料 posterior 及其 pushforward，本轮 FULL 未合格时不作此声明。

## 为什么风险排序与信息排序不同

若仅为分析而取一个局部 Gaussian covariance C_q，线性化材料的不确定性为 `D C_q D*`。在 prior-whitened 参数中，能最大限度保留这份近似材料方差的 d 维子空间由 `C_q^(1/2) D*D C_q^(1/2)` 的最大特征方向给出；这与 J 的最小奇异向量一般不同。两个排序只有在材料度量、先验和 likelihood 几何同时具有相应共同特征向量等条件下才可能一致。

这里的 C_q 只是已声明的局部近似或一个有资格的条件 covariance；用 GN 逆矩阵代入不会使它自动成为真实 posterior covariance。实际监督风险以新 parent 的 GT 损失评价，完整 posterior 的资格另行检查。这个线性代数结论是研究工具，不能作为独立原创定理。

## 有限响应的两类余项与抵消

令 Δχ=g(q+h)−g(q)，D=Dg(q)，C_χ=DF(g(q))，J=C_χD，所有数据量沿用原 whitening。则有精确记账分解

`F(g(q+h))−F(g(q)) = Jh + C_χ(Δχ−Dh) + [F(g(q)+Δχ)−F(g(q))−C_χΔχ]`。

第二项是 decoder secant 对数据的影响，第三项是 Maxwell 有限响应余项。这两项通常既不正交，也不具有相同符号。总平方误差包含交叉项 `2 Re⟨e_decoder,e_Maxwell⟩`；不能把各自范数占比称作可相加的因果百分比。本轮保存的拒绝提案和离线 GT secant 采用此分解，不用任意一个余项的大小证明另一个机制不存在。

## 动态弱方向与可积性

新增局部诊断用谱 projector P(q) 定义向量场，并在材料 RMS 下归一化，避免任意特征向量符号进入比较。在光滑、秩固定及谱分离条件下，两个场的 Lie bracket 是否仍属于所选分布，是 Frobenius 可积性的局部必要检查。decoder Hessian 的对称项在 pushforward bracket 中抵消，因此可将 `Dg(q)[X,Y]` 投影到材料弱切空间的正交余部，检查 normal leakage。

本轮只对两个已暴露 anchor、两个预定差分尺度进行真实物理检查。有限差分误差、谱隙变化、材料导数数值核及回路离散误差均保留；一个不收敛的 leakage 数字不能证明不可积，一个小数字也不能证明存在完整流形。这项诊断不改变轨迹的 candidate-only 或监督风险含义。

## 定向确认：灵敏度图中的 Maxwell 信息究竟在哪里

原 P1 在 ROI 内具有 `g_i(z,u)=host+δ(u) h((A_i z)/τ)` 的形式。对任意固定正交 frame V，有逐像素恒等式

`||D_z g_i V||_2 = |δ(u)| h'((A_i z)/τ)/τ · ||A_i V||_2`。

因此本轮局部 CNN 的 weak/strong sensitivity 图，由当前 decoder 导数乘两个固定空间包络构成。V 是由 TRAIN 的平均 `J_z^T J_z` 冻结的 frame；它不是每个测试对象上动态更新的弱 projector。对这两个输入图，Maxwell 的作用只通过这份训练 frame 进入。两个包络之间的比例在导数非零处与当前材料无关，还可能为共享局部卷积提供空间位置信息。

这是该 P1 decoder 和本轮具体特征的代数性质，不是否定一切物理特征。它要求更强的反证对照：保持 decoder 导数、通道数、参数量、训练规则和输入边界不变，仅把 V 替换成预先固定、独立于 J/GT 的随机正交 frame，同时分别去掉 sensitivity 与 residual-backprojection 通道。

第二批新对象确认据此冻结七个模型家族、三项主比较。即使真实训练谱胜过这个随机 frame，也只说明该固定包络对这份有限网络有帮助，不构成在线反馈必要性、完整 posterior 正确性或所有随机空间上的最优性证明。若随机包络同样有效，应将机制归于 decoder/空间风险表征，而非强行归于 Maxwell 谱。

## 输入摘要下的 Bayes 风险上限

若模型实际上只接收固定摘要 c=h(y)，平方损失下该输入类的最优动作是 E[χ|c]，未必是完整观测下的 μ(y)=E[χ|y]。条件期望塔式律给出

`R(E[χ|c]) = R(μ) + E||μ(y)−E[μ(y)|c]||²`。

因此物理基础图、固定谱灵敏度和残差摘要可以改变有限样本学习的效率，却不自动保留完整观测的全部信息。本轮的监督模型没有获得充分统计量证明，不能以 joint-simulation MSE 训练这一事实宣称已经求得 FULL Bayes mean。这个经典风险分解只澄清目标，不是新定理。

训练/validation的上下文反演上限为native48＋probit48，新测试为更强的native96＋probit64，存在明确的context covariate shift。测试上的配对风险结果仍是实际测试部署规则的合法比较；它并不证明训练context与部署context的条件分布相同，更不支持从训练收敛推导完整后验均值一致性。
