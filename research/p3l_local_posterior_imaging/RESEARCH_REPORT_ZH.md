# 从合格材料重建到可信 posterior：P3-L 实验记录

本报告是正在运行的独立 campaign 检查点。历史 P0/P1/P2、F3、A17-X 的结果和账本保持原样。这里的新结果均来自原 P1 二维 scalar TM/VIE；没有三维 vector-Maxwell、Mie 或新大型网络的证据。

## 为什么基础图已经不错，概率成像仍然困难

本轮 32 个新对象都已得到满足原数据资格的物理基础图，其中 27 个达到预注册的驻点标准。这建立了一个较强的单点基线。可是，找到一份能拟合观测的材料，与知道其周围材料的概率质量，是不同的任务。

原 prior 是 64 个 Gaussian field 系数和两个独立 uniform 属性。它们经过饱和的非线性 decoder 形成被动复材料。一个 Gaussian 参数扰动可能在当前位置看起来很弱，有限移动后却激活原来饱和的像素，使完整 Maxwell likelihood 急剧变化。于是局部线性化得到的椭球，可以给很多实际低 likelihood 的材料过多权重。

## 已执行的完整 posterior 参考尝试

已执行固定度量 MALA、条件初始化后的 MALA、变化 decoder 度量的 MALA、正确 probit 初始化后的 MALA，以及 canonical pCN 对照。六组各有四条链，每条保留 1,024 warmup 和 4,096 production。接受率提高没有自动带来混合：变化 decoder 度量的配对实验多用了约 15.9% job wall，材料函数 Rhat 仍显著不合格。

这不是把 uniform prior 改成 Gaussian。我们使用精确坐标变换 `u=normal_CDF(eta)`，并把变换 Jacobian 纳入目标。在这些坐标中，原 uniform 的推送恰为标准 Gaussian。生产阶段不再按历史自适应；变化度量提案使用两端的 Gaussian determinant、drift 与二次型，likelihood 和梯度始终来自真实 full-Maxwell。

独立 tempered SMC 也已跑到 beta=1。128 个粒子在退火中只剩一个初始祖先，并未形成可信独立 reference。祖先数本身不是错误的定理证明，但结合有限 mutation、局部链分歧和缺失交叉验证，不能把它的权重 ESS 当成全局后验质量证明。

另完成四条 HMC pilot，每条 512 warmup 和 1,024 production。它使用完整 66 维目标、正确质量矩阵约定、leapfrog 和 Hamiltonian MH；有限能量误差不会被任意截断成“成功”。实际材料函数最大 Rhat 为 1.26968，最小 bulk ESS 为 11.47，材料均值 RMS MCSE 为 0.01015。性质坐标的最大 Rhat 仍为 1.51447。因此它改善了局部混合，却仍不足以回答完整 posterior 的材料风险。

下一项修订不扩大网络，而是使用非线性、全局可逆的三角坐标：在性质坐标改变时，预先拟合空间场的条件补偿中心，保留中心周围全部 64 个残余坐标。该变换行列式为 1，但先验必须在变换后的真实材料坐标中求值，不能另用标准 Gaussian 代替。它是采样预条件器，不是把零厚度中心面当作 posterior。

## 强方向积分的实际失败

四个事前按观测特征选定的对象，分别取弱维数 2、4、8、16、24，在零位移及固定 seed 位移处计算强条件 GN 参考，再做 128 次真实 likelihood 重要性校正。40 组中有 38 组 ESS 小于 2，最大 ESS 也只有约 5.27。它们都不具备 teacher 资格。

有限优化未完全达到 mode 标准时，单纯 `profile - logdet/2` 还缺少局部线性梯度项。我们保留了这个问题，而没有改写旧结果。进一步的零物理分析发现，10 组的该项 log-normalizer 修正上界小于 0.01，但权重仍发生明显退化；因此失败不能仅归因于优化没有再多走几步。更关键的是目标在所用 GN Gaussian 范围内强烈偏离二次近似。这个偏差包含 decoder、非 GN 曲率和更高阶响应，不能全部归为 Maxwell feedback。

解析饱和 decoder 矩已经通过独立高精度积分和 Monte Carlo 验证。它准确回答的是“给定这个 Gaussian 参考，材料均值是什么”。它没有把一个错误的 Gaussian 参考变成正确 posterior。

## 对材料图像有何后果

我们先冻结全部 32 个 development 对象的估计器、保存材料图，再在独立离线步骤读取 GT。优化基础图的平均 complex MSE 为 0.0138212；固定属性 GN Gaussian 的解析材料均值为 0.0206962，差为 +0.00687494，parent paired bootstrap 的 95% 区间为 `[0.0029612, 0.0120149]`。因此这一近似在本 development 总体上明显恶化成像。Gaussian 参考下的几何中位数接近基础图，平均 MSE 约 0.0138410。

这是一项对具体近似的负结果。它没有比较已知正确的 FULL posterior mean，不能证明 Bayes 成像本身无益。数据仍是 development，不能称为 locked 泛化结果。

随后完成全部 32 对象的联合属性修订。我们先在正确 probit 目标中重新优化，再以完整 66 维 GN Gaussian 为参考，通过两属性条件化及解析场材料矩，计算联合材料均值。GH 阶数 5、9、17 的材料均值已稳定，最大后两阶 RMS 差约 2.90e-6。这排除了本次观察由粗 quadrature 引起的解释。

新的坐标 mode 材料图平均 MSE 为 0.0120131，较原 native-q MAP 低约 13.1%；可是联合 Gaussian 材料均值仍为 0.0185700，较同一新 mode 高约 54.6%。四份固定对象上的 64 点 Full importance correction 都退化至约一个有效样本。因此，仅加入 property-field 相关性并未修复 Gaussian 参考；不能把坐标 mode 改善写成 posterior mean 改善。

这里还存在重要的数学区别：posterior 概率律在 probit 变换下保持不变，但 density mode 不是坐标不变量。在原 u 坐标中，该坐标 mode 等效于额外加入 `0.5*||normal_CDF_inverse(u)||²` 的点估计惩罚；它并非原 uniform-measure MAP，也不是 Bayes action。若将其作为实际估计器，需要独立新对象验证，不能用这组 development 数字宣称概率方法获胜。

[四份事前固定对象的实际材料对照](B/joint_gaussian_figures/index.html) 采用共同实部、损耗色标。方差图仅属于声明的 Gaussian 参考，未通过 posterior 校准。

## 一个有正确含义的低维条件参考

为区分“物理接口有错”和“高维参考困难”，在第一份冻结对象上构建了两弱一强的三维条件问题。其他 61 个 Gaussian 坐标及两个属性固定；有限 box、原条件 Gaussian prior、Jacobian 和噪声在求值前登记。所有节点都用真实 full-Maxwell 计算。

三阶 tensor quadrature 的最后两阶，材料均值 RMS 差约 8.46e-10，log-normalizer 差约 6.64e-8。这里的数值收敛说明受限条件问题可被准确计算。它不回答 box 外的质量，也不代表完整 posterior 的弱 marginal。没有用节点分箱展宽伪造连续 quantile 或 credible interval。

## 几何诊断与后续方法

在两个固定 anchor 上，以三个有限差分步长计算真实 full-J 导出的 rank-8 投影场。所测方向对的坐标法向 Lie bracket 收敛到约 0.0913 和 0.1295；冻结 Maxwell 材料响应、只更新 decoder 的对照已复现大部分数值。它是具体投影场的局部数值诊断，不是任意材料 posterior 不存在、exact fiber 存在或一般不可积性的证明。材料度量解释和谱隙都需要一并保留。

归一化 CNN/MLP Gaussian、MDN 和 Flow 的实现已完成。512 train 与 128 validation 的 joint-simulation 数据及在线 context 已准备好。这条监督路线与 posterior distillation 分开；尚无合格 posterior teacher。新模型仍需经过实际训练、强基线比较及 full-likelihood 校正检验，不能仅凭架构归一化宣布概率正确。

实际的 5 个模型族、3 个 seed 共 15 次训练已完成。每次均运行 1,000 updates 后触发预注册 validation patience，验证集都选择 update 0。独立复核恢复了全部 checkpoints：选中权重与各自随机初始化完全相同；训练 minibatch NLL 下降而验证 NLL 迅速变差。未发现 batching、参考测度 Jacobian、frame/sign 或 checkpoint 选择错误，但这不足以证明观测没有信息。我们将其记为本规模下的泛化失败，而不是训练成功。

进一步的简单密度对照使用同一 512/128 数据、train-only PCA、ridge 均值和收缩 covariance。验证 NLL 为：选中 PCA Gaussian 14.8138，无条件 Gaussian 15.1007，最佳神经初始化 15.2620。它只是归一化 proposal 的验证分数，不是 FULL posterior 校准。后续真实物理校正优先用于这两个已拟合简单模型，停止为验证集已否决的神经训练额外支付查询。

## 求值加速的独立实测

对同一 anchor 的 16 份预先冻结提案，已运行 optimized Full、完整 T-matrix、预条件 GMRES、rank-8/16 普通 Galerkin 以及同空间 A17 表达。所有比较使用同一观测、相同 proposal 和接受随机数，近似决定先于 Full 参考落盘。

rank-16 的条件误差区间提前判定了 5/16 次决定，其余 11 次需要 Full fallback；rank-8 没有提前判定。同样的随机数下没有观察到决定差异。这个结果基于明确的浮点误差模型，不是定向 interval arithmetic 的形式认证；失去 Neumann 充分条件时只记界不可用，不认定 Maxwell 算子奇异。

新付费成本为 960 RHS：576 RHS 的保存逆 setup、128 T-matrix、128 GMRES、128 Full 验证。后者作为共享 fallback 已计费，没有让 oracle 免费替近似挑选决定。GMRES 的 16 份提案均收敛。相同 current space 内 A17 与普通 Galerkin 的等价性不构成独立算法增量；是否节省总 wall 尚需包含 setup 和界计算后判断。

## 中断、费用与证据边界

本轮独立上限为 12,000,000 complex RHS、48 小时受控 CPU、24 GPU-hours。按 A–E 独立记账，尚未使用 GPU。每次 full-forward、adjoint、LU、失败和拒绝均保留；原始端点立即落盘。

一次内存中断来自缓存持续强引用旧 dense matrix/LU。六个中断作业的 109,488 RHS 和 3,408 state/LU attempts 均保留。缺失尾部运行时间采用明确标记的保守上界计费。修正后限制 transfer 缓存为 4、逐对象清理并强制一个 worker；640 个在线 context 已恢复，实测约 443–448 MiB。完整恢复审计在 AUDIT 下，不将失败伪装成正常完成。

目前的科学判断是：简单局部 Gaussian reference 尚不能支持本 P1 的可信材料 Bayes 成像；完整 posterior 参考仍是主要障碍。是否存在值得发表的算法增量，要由后续 HMC、校正 proposal、实际风险和完整成本共同决定。没有证据时不宣称强 IEEE TAP 贡献。
