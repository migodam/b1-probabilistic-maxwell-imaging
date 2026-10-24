# P3L-R4 完整研究报告：条件材料输运、风险成像与 Maxwell 自适应计算

## 1. 这轮实际解决了什么

我们完成了独立的训练、物理与泛化实验，而不是把 R3 理论重新整理成计划。结果需要分开看。

材料成像有可重复进步：新方法在 512 和 256 个全新对象上，比历史 R2 CNN 分别降低约 13.10% 和 11.36% 的材料 MSE；比本轮更强、训练数据更多的 CNN 降低约 5.36% 和 4.90%。但普通 Gaussian 产生的全局材料统计与 Spline Flow 一样有效，不能把这一成像收益归因于 Flow 已经学会正确 posterior。

完整 66 维条件密度也有真实的相对进步。Spline 在独立联合样本上的 proper NLL 优于同容量 diagonal Gaussian；固定 Maxwell frame 的三角模型优于随机 frame。有界的 Maxwell 几何 CFM 使材料 energy score 比同谱随机控制改善约 21%。这些是相对分布质量证据，不能抹掉绝对概率失败：固定生成样本的 Maxwell 残差仍比原噪声模型预期大一至两个数量级，重要性权重仍然退化。

物理计算方面，在 8 个全新对象、相同提案和随机数的 8,192 步重放中，收紧数值势区间把 Full 的 71.40 秒降到 24.15 秒，逐步决定和状态路径零差异。它没有改善这个短步长核的混合，也不是连续 Maxwell 的严格认证。普通 GMRES 同样带来收益，A17 在这里没有独立优势。

动态非完整几何没有建立独立算法收益。确实保存了 8 份观测下 111 个新的任务尺度、有限噪声竞争候选；它们来自多种普通和几何方法的共同池。它们不是 111 个 posterior modes，也没有概率权重。四段路径在距离保持与成本上弱于直接方向提案。

## 2. 物理与统计问题始终相同

本轮使用原 P1：24×24 scalar-TM/VIE、256 个材料 ROI cells、8 次照明和 32 个复接收通道。前 64 个独立 Gaussian latent 经相关场矩阵生成材料，后两个 Gaussianized 属性经标准正态 CDF 回到原 uniform 属性。权威 decoder 来自冻结仓库；没有用摘要重新写一个相似模型。

原标准差为 0.00119211798514519。记 Phi=||F(G(x))−y||²/sigma²，posterior 相对于 dx66 为 exp(−Phi−||x||²/2)。复数据沿 source-major 顺序打包，每个 source 先 32 个实部再 32 个虚部；实白化系数为 sqrt(2)/sigma。所有照明共享同一个被动复材料。没有改 noise、gain、采集、材料 prior 或物理后端。

decoder 在两侧精确饱和，多对一映射允许退化材料层；所以 latent density 与材料像素 Lebesgue density 是不同对象。本轮完整密度始终在 66 维原测度上定义，不把 material Gram 行列式误加进 prior。

## 3. 数据与独立性

恢复并复用了历史 32,768 个 joint TRAIN 和 2,048 个 validation 对象。历史昂贵 physics context 只有 512/128 子集，不被记成 32,768 个物理锚点。本轮另外付费形成 1,024/256 的强 context 训练与验证集。

新生成 128 个 development parent，随后在方法、checkpoint、NFE、样本数、失败规则和评分合同冻结后，生成 512 个第一批与 256 个第二批新 parent。两批使用相同方法；第一批结果没有用来重选第二批算法。两批所有风险与密度在线输出、Full PPC 和优化诊断封存后，才开启独立 GT evaluator。

风险面板包含 45 个固定方法实例，密度面板包含 26 个固定模型实例；二者在全部 768 个对象上均无推断失败或 fallback。随机种子不是独立对象。主表先对每个 parent 的三个 seed loss 求平均，再以 parent 做 4,000 次 paired bootstrap；它不是三模型集成图的成绩。预注册主比较使用 98.333% 区间，其他比较属于描述/探索，未宣称 45 方法的全族错误率控制。

## 4. 材料风险：强基线之后仍有收益

|全网格 complex MSE|新512对象|新256对象|
|---|---:|---:|
|独立 Probit physics point|0.01091244|0.01133359|
|冻结 R2 CNN|0.00783202|0.00832320|
|同一 5,842 参数类别重新训练，旧512 context|0.00741971|0.00798060|
|更强 CNN，1024 context、真实物理谱|0.00719170|0.00775824|
|Spline 全局矩 + 局部风险 head，预冻结主要方法|**0.00680623**|**0.00737782**|
|普通 diagonal Gaussian 全局矩 + 同风险 head|0.00682170|0.00739262|

对 R2 的相对改善为 13.10% [11.37%,14.89%] 和 11.36% [8.84%,13.96%]；对强 CNN 为 5.36% [4.02%,6.91%] 和 4.90% [3.14%,6.72%]。Spline 与 Gaussian fusion 仅差约 0.23%/0.20%，区间跨零，不能声称 Spline 在风险上获胜。

主要模型不是“只有 6,418 参数的完整系统”：6,418 是局部风险 head，另有约 0.68M 参数的冻结全局 provider，K=32 材料样本和原 physics context 的取得费用。这个成本需要与单独 R2 CNN 区分。

固定物理谱相对同容量随机空间控制的风险增量约 1.20%/0.83%，小但两批方向一致。它是静态 context 的贡献，既不是动态 weak flow 的贡献，也没有增加观测信息。

## 5. 为什么劣质全局生成均值仍然有用

完整 Spline 模型的直接样本均值，第一批 MSE 约 0.02278，远弱于 CNN。把这个图直接替换重建是坏方法；把它的均值和方差作为四个全局 context 通道，交给带物理基础图的局部 CNN，却有价值。

这与标准条件期望正交性一致：局部 patch 未必包含整个观测的所有任务信息，另一份仅由 y 计算的粗糙全局统计可以补充它。保持相同 head、打乱这四个通道的空间对应关系后收益下降；普通 Gaussian 或直接 field-moment provider 的控制又表明，收益主要来自有用的全局上下文，而不是 posterior 已校准。

按基础图的相区做离线分层，第一批 R2→新方法的全图风险下降中，过渡相贡献 0.00102460，其他单元贡献 0.00000118；第二批为 0.00094420 与 0.00000117。约 99.9% 落在基础图的过渡相。这是分层恒等式与描述证据，不是 Maxwell 独特因果机制或 GT 可达性证明。

## 6. 相阈值条件矩：理论正确，当前算法不胜基线

真实复属性与相事件配对监督，未独立化 E(delta) 与事件概率。先运行 K8/K16 点阈值，再把标签改为精确 h 质量 bins，逐样本求和恰好等于原 decoder；因此去掉了求积偏差。最后只改变为 anchor-induced mass 初始化，保留标签、容量和优化器。

三种版本均未超过强普通 CNN。第一批精确 bin8/bin16 的 MSE 为 0.007742/0.007746；anchor-bin 更差，约 0.010263/0.010019。宿主 plateau 的 mass/logit 梯度存在已验证抑制，但不能由此解释所有纯相及过渡相。停止继续调整 epsilon 或扩大此 head，不把恒等式当成算法成功。

## 7. 完整条件概率：相对改善是实测事实

Spline 更新全部 66 维，显式逆变换及 Jacobian 自动核验；三角模型真正为 q(a|y)q(b|a,y)，不补独立的剩余 prior。所有模型使用原 joint simulation proper NLL 训练。

|独立概率指标，三个 seed 的 parent loss 均值|新512|新256|
|---|---:|---:|
|同容量 diagonal Gaussian NLL，nat/parent|83.07743|83.44856|
|Full-66 Spline NLL|**77.50922**|**77.27755**|
|三角随机 frame NLL|77.95014|77.86049|
|三角 Maxwell frame NLL|**76.59216**|**76.25860**|
|anchor CFM 随机同谱控制 material energy score|0.165310|0.174586|
|有界 Maxwell 几何 CFM energy score|**0.130661**|**0.138822**|

Spline 对同容量 Gaussian 改善 5.568/6.171 nat；固定物理 frame 对随机 frame 改善 1.358/1.602 nat。相同原测度下，joint-test NLL 差估计两个条件模型相对真实 posterior 的期望 forward-KL 差，未知的 conditional entropy 相消；不需要 MCMC teacher。它不提供绝对 KL、全局覆盖或可信区间资格。

CFM 物理权重的 energy score 改善约 20.96%/20.49%。权重在实际网络标准化坐标里有界于 I 和 4I，所有 66 个速度分量保留；相同谱随机旋转控制享有相同容量、数据和 base。这个正面结果不等于 density 已有可用精确 logq。

值得注意：三角 Maxwell frame 虽改善 latent NLL，却没有同样改善 material energy score；第一批随机 frame 0.09964，Maxwell frame 0.10127。不同 proper score 衡量的目标不同，不能挑一个指标把所有结果叙述成一致胜利。

## 8. 为什么仍然不能称完整 posterior 合格

每个 parent 使用预冻结 seed1 的 index0 样本，不能选择最小 residual。原 posterior 的 prior-predictive 检验应给 Phi~Gamma(256,1)，均值 256。实际 Spline 的两批均值为 8,371.95 和 8,662.42，第一批 510/512、第二批 255/256 的 PIT 超过 0.975；几何 CFM 分别为 22,139.70 和 23,557.96，仍严重失配。

对固定前16对象的 K64 Spline 重要性修正，ESS 中位数为 1、最大仅 2.273。MCMC 校正核的正反比正确，拒绝重复状态完整保存，但有限链的 ESS/MCSE 不合格，没有据此输出 FULL posterior 均值或 calibrated mass。

实际训练了三类、成对 physics/supervised 的混合双向 KL 方案。Full Maxwell 力可以降低生成样本 Phi，但 NLL 恶化超过预先登记的复制门，覆盖明显受损；没有继续扩大种子。失败说明当前损失与模型在 response-fit/coverage 间存在冲突，不说明 posterior mean 不如 MAP。

CFM 的数值密度另有独立问题：双精度全迹 DOP853 审核用尽 895 秒上限，两个 tolerance 的某对象 NLL 差 0.136 nat，未过 0.05 稳定门；NFE32 离散值与更精细结果也有 0.342 nat 差。其材料样本 score 仍可评价，但不能把连续 CNF 近似密度当有限步精确 map density 用于 exact-MH。未完成的八个求解标记 NOT_RUN，而不是擅自推断全家族失败。

## 9. 风险图与合法物理世界不是同一个输出

对全部 768 个对象的六种风险图，进行了 4,608 次实际 Full forward。新主要 seed1 clipped 图只有 242/512 和 125/256 满足原 noise-ball；冻结 R2 为 256/512 和 129/256。风险改善没有自动带来更多数据相容图。逐像素 clipping 也不能证明来自原 66 维 decoder。

另对固定第一批前64对象、五个臂进行独立合法材料代表优化，实付 80,640 RHS。所有 320 臂都在预算处停止，未达到普遍驻点；84 个回到原合法 anchor。输出满足数据门的原因包括显式 fallback，不能把它冒充优化全部成功。

原普通 SLSQP 的 callback 只观察大迭代试探点，并非保证线搜索接受点，导致最优已付费合法查询未被选中。原失败保留，零新物理 best-query 修复和新的 Full-J trust-region 诊断在首次 GT 开启前完成。后者证明：data-only TR 可以继续降低全部64对象 Phi，却把 MSE 从约 0.010514 升到 0.016346；保持原 Gaussian prior 的 TR 则几乎不改变风险，约 0.010505。现有 NN 收益不能仅归因于一个严重未优化的 regularized baseline，但也没有全局最优性结论。

材料代表投影保留部分风险收益，Spline raw 约 0.006748 投影后约 0.009154；代表比原图更差。这不是“理想合法 Bayes action 必然更差”的证明，因为投影未充分收敛、表示也受 decoder 限制。

## 10. Maxwell 复用与接受决定：什么真的节省了成本

精确 active-row T-matrix 用正确本构增量、缓存 inverse action 和全照明场。与 Full 数值一致，但有效变化通常涉及 180–256 个 cells，低 latent 维不等于低端口秩。含 setup，Full 7.427 秒、复用 8.374 秒，停止加速主张。

数值势区间的首轮链同样失败过：虽然 RHS 下降，wall 从 229.78 秒升到 355.79 秒。随后只收紧可用逆范数上界并减少不必要动作，未改 prior、kernel、接受随机数和统计目标。前瞻新对象的安静重放得到：

|同一8对象×1024步|完整计费 RHS|全流程秒|
|---|---:|---:|
|optimized Full|65,536|71.397|
|区间 Richardson|4,608|24.153|
|同缓存全局 GMRES，决定驱动|14,691|37.605|
|同缓存 GMRES，固定精度|20,759|46.832|

所有逐步接受/拒绝、潜在区间检查及状态路径零差异。包含每对象 576 列缓存 setup、rebase、残差/低阶线性动作与 fallback 权利；减少 RHS 不是单独可信指标，实际 wall 也下降约 66.17%。结论范围是此 beta=0.002 核的数值重放，不是 posterior mixing 或材料 MCSE 的加速证明。更大 beta 的边界实验显示 fallback 增多，接受和实际移动迅速下降。

## 11. Learned Jacobian 与非完整探索

六次局部 value/JVP surrogate 训练全部未过最大白化输出误差和 Phi 误差 0.05 的门。已有 teacher 成本、普通 field-update、Broyden、frozen-J 和 Full-J 均公平保留。学习动作误差的排序器仍需要候选 Full 场特征，且弱于便宜的 total-field+exact-decoder 更新。它不是完整 LazyDINO 的复现，也不能据此贬低该文献。

四段 Euler 的二阶项混入自身场加速度，不能直接当 Lie bracket。修订为 Heun 后，原任务尺度尝试 288 项中仅 9 项满足2%材料距离匹配，最终全部支持门只有2项；直接方向法有94项，frozen/dynamic 共同阈值描述为59/83。原1% teacher 条件未改写，新的2%门只是离线比较适用域。未匹配平均 GT 误差不能用于算法排名。

共同 archive 按不读取 GT 的 Phi/UID 排序、完整材料 RMS≥0.02 去重，八份观测分别有 15、12、7、16、28、16、15、10 个候选（含基础图），即 111 个新候选。它们满足各自 noise-ball、ΔPhi≤2、prior_delta≤2 和进度门，不代表两个世界的预测距离≤1，不是 exact fiber，也不是 posterior mode 数。

## 12. 成本、结论与下一步

独立资源账本和全部 CLOSED 收据是费用权威。最终资源与失败明细见 [RESOURCE_SUMMARY](RESOURCE_SUMMARY.md)；复用历史训练数据没有把它们的冷生成费写成零。全部原始材料、latent、current/field、预测、模型、chains、失败和配置保留在证据包中。

本轮没有把全部研究押在 MCMC 或某个 Flow 上，三条路线都完成了真实实验和有实质变化的修订。没有启动三维、新 Mie、外部付费计算或未经门禁的48/64网格及分布外扩展。

**科研裁决：继续风险成像与静态物理条件学习；保留限定场景的决定驱动数值求值；停止扩大当前相阈值 head、active-row 加速、动态四段路径和当前 learned-J surrogate。完整 posterior 与强 IEEE TAP 核心创新仍未建立。**

下一轮最重要的不是再换一个网络名称，而是把已有相对概率进步转成同噪声尺度的 response calibration，同时保住覆盖，并检验在更大二维系统上相同决定精度是否仍有端到端成本优势。所有新目标与预算需另行冻结；本轮不再追加计算寻找正例。

## 13. 相关工作与可发表性

完整条件Flow、simulation-based learning、DINO、局部预条件与retrospective接受判定都有成熟先例。本轮不把这些工具拼接后的名称当贡献；逐项原始来源和不支持的主张见[相关先例与主张账本](PRIOR_ART_AND_CLAIM_LEDGER.md)。当前最有价值的是可被独立核验的新实验现象及失败机制，仍需要新的Maxwell专有机制与决定性比较。
