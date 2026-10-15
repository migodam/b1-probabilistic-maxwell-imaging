# 理论修订：把补偿机制与材料任务分开

## 1. 正确但不自动有效的原理

局部分析必须先堆叠旧/新白化 Jacobian，再对同一个 nuisance 作 Schur 消元。每个材料世界允许自己的合法 gain；共享仅意味着该世界跨测量历史保持一致。不能强迫不同世界使用相同未知 gain，也不能每频重新抽一个 gain 后乘边缘 likelihood。

实际更新使用 `log L(old,new|material) - log L(old|material)`。解析截断 Gaussian gain 保留 determinant 和截断归一化，27种通道/极端截断组合的归一化posterior与独立joint计算最大差4.44e-12。这个软件核对不构成 Maxwell 成像性能证据。

## 2. 第一轮被否定的算法解释

E0全部七案例中，单独最大化relative-gain compensation-conflict分数都选calibration；完整材料风险/Fisher/EIG都选2.4。冲突分数是一种候选对secant启发式，不等于复合假设风险、minimax regret或完整material Bayes risk。它忽略了任务相关性、候选对代表性和对旧数据的条件约束。后续保留它作为负对照，没有按truth调权挽救。

E1/E2更宽菜单下，普通Fisher和固定3.2已足够强。额外测量改善MAE，与复杂准则优越是两项不同主张；前者得到限定支持，后者未得到支持。当前无需ROM、图网络或神经policy增加复杂度。

## 3. full-Mie非线性见证的严格范围

固定small central r=.08R、large central r=.25R、annular h=.08R，优化被动epsilon与一个有界shared gain。所有12个joint端点经full-Mie验证；对照优化和minimum仍为best-found上界，不是全局下界。

large r=.25在1.6+2.8的full白化距离约3.225，在1.6+3.2约3.176；small r=.08各单个附加频率仍很近。ring在旧频就已明显分开，不能把它的高距离当作解决原歧义的主要成功。

同一epsilon端点下重新分别/共同profile gain，可隔离gain冲突。large3.2的gain冲突平方约.589，ring3.2约7.829；许多其他组很小。不能把所有新增频率收益都归因于shared gain冲突；频谱本身、epsilon补偿限制及已有分支差异同样起作用。

## 4. 从单步见证到有限预算的正确补充

不同频率分别优化的近对可能是不同材料，因此不能仅靠它们宣布三步测量不可能。新增独立6次full审计拟合一个共同epsilon/gain小孔世界，再一次核对全六频、calibration和独立geometry动作。

对这两个固定世界，每个动作的观测是具有已声明相同协方差的Gaussian；策略仅依赖已有历史。链式KL给出 `KL <= (d_old² + 3 max_a d_a²)/2`，随后Pinsker给等设计先验二点分类错误率下界。实际数值与全部字段在 `results/SMALL_CAVITY_HORIZON_WITNESS.json`，三次新增block下的下界约0.440。

这是同一full-Mie数值台架、固定两世界、有限已注册菜单下的minimax诊断；不是完整posterior事件质量、所有材料世界的下界或现实仪器证明。小孔在体积MAE中的权重也很小，结构任务应另外保留。

## 5. 最重要的新限制：支持不完整

625个full-Mie候选的Best-of-K MAE均值仍.062224，Coverage@K(MAE<=.05)只有.40625。频率选得更精细不能突破固定池的oracle下限。

删除central/annular/combined家族后重算动作，确实改变推荐，但残差sentinel常漏报。未知支持质量必须为null；不能把profile residual变softmax后宣称越来越确定。优先修复候选发现，再讨论停止或概率精度。

E2所谓refresh只是把已有2048粒子支持中的样本追加到decision catalogue，没有发现新的后验支持；只有E3独立corrector生成、通过双频full验证后安装的材料才可统计为新增候选。

## 6. 需要继续但当前未获证实的理论

精确的连续复合minimax regret、开放材料空间coverage、一般非同心矢量场景、真实色散/仪器模型和ROM inclusive加速均未建立。本轮不把标准Schur、Bayes risk、链式KL、Pinsker或序贯停止概念重新命名为原创定理。Maxwell-specific贡献候选应来自可复算的物理补偿/失败机制与材料后果，而不能来自这些通用公式本身。

## 7. 完整history参考的数值条件必须独立验证

初始posterior积分稳定不蕴含加测后稳定。新增数据会收缩/移动模态，2048条件样本重加权与固定13/17积分都可能失稳；本轮独立physical-coordinate积分14/16未达原规则，不能拿另一个未收敛估计给原reference背书。实际full字段的候选选择MAE仍有效，但EIG、posterior event mass与expected risk必须保留近似性质。

E3候选补全改善四例材料MAE不等于结构全部正确，也不等于发现完整posterior支持。独立corrector与新增观测同时变化，下一轮需同容量无新测量对照，才能判断信息、初始化与优化的贡献。当前证据优先支持这一机制区分，不支持再造复杂选测准则。
