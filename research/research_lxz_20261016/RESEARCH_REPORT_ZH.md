# LXZ：局部–全局Maxwell多解研究，2026-10-16

**本轮真正新增的显著结构材料world：0；原N3附近的合法局部fiber：仍未证实。** 已完成两个真正可运行的实验：四个暴露径向观测的概率采样比较，以及两个公开TRAIN原P0 smooth材料的full-TM局部校准。保留N3强Top1，未将TRAIN已知真材料包装成未知物体重建。

**最明确的正结果是概率修正，而非多解成像突破：** 同4096次proxy likelihood预算，C1全局MH eventMAE .004364，普通local .012039，hybrid .007243。Gain-aware geometry .011720比原geometry .026158好，但未胜过简单global。E3未获分离mode证据，新增fullMie为0。

## 从材料成像问题出发

同一电磁观测可以容许材料不确定性，但图像相似、数据noise相容、原先验合理、posterior mode是不同要求。N3低MAE值得保留，却不是原P0合法anchor的证明。旧T/G任意复像素补偿进入noise球仍破坏固定增量相位。不能从这些材料开始宣称原P0解族。

本轮完整读L/Z/X实际41份文本材料。L正式全文在本包确实存在，但只提供合成代数实现；Z给全局提案设计；X提供标准正确采样法则。统一实现没有强行加Flow或A17。代码分为真实prior、固定target和local/global kernel，corrector单独candidate-only。

## E0及目标隔离

development N3仅2/32过noise球，0/32过原phase必要条件，故无合格P0重建anchor。公开TRAIN smooth250033/weak、250034/ring合法且残差22.3130/22.9696，只能做已知材料校准。JointP2 C proxy和fullMie、原P0和条件smoothP0、未定义P1分别登记，不混比较概率。

## E1：概率机制真实运行

| 方法 | event MAE | 材料CDF grid误差 | recall* | 4096次平均arm秒 |
|---|---:|---:|---:|---:|
| local | 0.012039 | 0.071142 | 1.0000 | 0.4998 |
| geometry | 0.026158 | 0.178926 | 0.9479 | 0.5615 |
| geometry_nuisance | 0.011720 | 0.075228 | 1.0000 | 0.5531 |
| global | 0.004364 | 0.028285 | 1.0000 | 0.5079 |
| hybrid | 0.007243 | 0.042557 | 1.0000 | 0.5231 |
| tempered | 0.016145 | 0.103910 | 0.9896 | 0.0127 |
| prior | 0.011206 | 0.074878 | 1.0000 | 0.4307 |

*固定参考事件质量≥.05、经验估计≥.005；这是本轮false-exclusion诊断阈值，不同于旧P2的.01 discovery。CDF仅保存grid最大差。SMC为64个加权粒子，MCMC为保留拒绝重复的链平均，ESS不可直接互换。


四观测solid/large_cavity/annular_cavity/combined已暴露，四随机seed，预算含64点初始化和8点几何差分。Global用原四维C1normalizedfullMDN的.9密度与.1prior防御混合，加准确MH ratio。Local用Gaussian-CDF prior coordinates；normal/gain-aware metric冻结。TemperedSMC64粒子，所有温度和likelihood调用计费。初始bestprior单链与SMC完整population的差异明确披露。

有correctkernel不等于finitechain精确posterior。所有CDF/mean/SD/quantiles来自同一empirical measure，与旧独立continuous参考分开比较，保存拒绝重复。多事件仍可能在一条连续ridge上，不能据此宣布跨mode成功。原reference有限数值误差保留，没有uniform fullMie/continuum保证。

## E2：真实full-Maxwell局部材料试验

| 校准anchor | 方法 | compatible记录/提案 | 去重compatible | 新结构 | 最大材料MAE位移 |
|---|---|---:|---:|---:|---:|
| 250033_weak | local | 5/5 | 5 | 0 | 0.000354913 |
| 250033_weak | geometry | 4/5 | 4 | 0 | 0.00134309 |
| 250033_weak | prior | 0/5 | 0 | 0 | 0 |
| 250033_weak | continuation | 10/10 | 2 | 0 | 0.000352185 |
| 250034_ring | local | 5/5 | 5 | 0 | 0.000448741 |
| 250034_ring | geometry | 5/5 | 5 | 0 | 0.00115676 |
| 250034_ring | prior | 0/5 | 0 | 0 | 0 |
| 250034_ring | continuation | 10/10 | 6 | 0 | 0.00241842 |


8个原Gaussian参数只有7个材料方向，1个max-normalization gauge。允许小幅noise-compatible连续材料变化，但原多尺度结构门槛新增world0；重复field去重。普通local及geometry的最大compatibleMAE位移约.00035–.00134，continuation另一例最大.00242。这些是对knownanchor的位移，不是N3误差改进或未知物体重建。

第一例weak方向映射到prior参数范数1998.55，高曲率使大步退化。只缩小FD步长10倍，full derivative相对误差从8.097%降到.08097%；原材料Dg整体FD误差仅8.44e−9。说明localweak方向并不自动生成有意义finite material family。停止把一阶nullspace维数当多解能力。

## 失败、完整费用与边界

几何差分index错导致两原臂解释无效，保留原131072次proxy；修复后单独重跑，合计proxy 589824。E2首35state因chart rounding失败并丢失原始endpoint数组，账本仍保留；下一任务缩短至5步并逐endpoint保存。最后5state只做导数核验，总96/96 fullstates、896/3072 complexRHS，GPU0、新Mie0、非同心3D0。失败纯wall未测，成功job加保守失败elapsed上界97.25s，低于3600。未借旧预算或覆盖历史。

## 八个最终问题

1. N3附近非平凡材料不确定性：**未证实**，缺同prior数据合格N3anchor；knownTRAIN只见小局部宽度。
2. A17增量：**NOT_RUN/未建立**，ordinary full-J已经足够便宜，无可信likelihood区间的早筛。
3. 其他重要结构mode：radial有多个有质量事件，**分离mode未建立**；复杂P0无合格竞争见证，非唯一性证明。
4. Physics-context latent：已有C1 density经MH有效，但它只条件于measurement，不证明新增feedback/currentcontext价值。
5. Local–Global：本轮混合劣于global，**没有独立混合优势**。Gain-aware local改进只恢复普通local级别。
6. 瓶颈：complexP0主要是合法inverseanchor、非光滑/狭窄representation及finite-chart曲率；radial是有限mixing与density近似。不是缺更大NN。
7. 下一步：优先合格anchor/有可信reference的新ambiguitybenchmark或active sensing；**暂不训练Flow/Diffusion、不扩ROM、不启动3D**。
8. 强TAP论文：**当前不足**。缺同一复杂物理target上可信重要材料区域及同完整成本的材料成像/风险正收益，不能用proxy误差或图数替代。

候选JSON/NPZ、完整MCMCchains、SMCcheckpoints及失败账本见本目录。后验未知质量与模型外质量保持null/unquantified，不把TopK重归一化成图像概率。

## 可交接的实际proxy概率图

图库另提供四个径向观测的真实global-MH链事件代表，JSON含每事件质量、4seed描述范围、prior质量、条件材料mean/SD/quantile、proxy prediction与完整链样本。展示K仅覆盖至少.95经验质量，剩余事件质量保留，不重新归一化；模型外/漏区域质量不可量化。它们不是新的fullMie候选或已证分离mode。为代表prediction另外支付22次proxy response（不计算likelihood），完整ledger已单列，full物理0。
