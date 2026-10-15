# E0：已有多解轨迹的零新增物理诊断

E0已实际执行，新增full-state、RHS、solver构造、重跑corrector及GPU均为0。它解释历史信息，不能计作Research G新发现的材料世界。径向Mie和二维TM分别统计，不拼接成更大的独立样本。

## 已发现多少不同结构

32个径向parent、五方法、1920个proposal中，原静态方法正式archive共有193个材料成员、160个几何代表，分别为每观测6.03125和5.0；multi-start为82个材料、56个几何代表，平均2.5625与1.75。它们是各方法各观测内的计数，不能跨方法/对象相加成全局不同世界。

原静态方法明确用了八个径向保护box，不能直接当G的topology-agnostic连续尺度分层baseline。C1四维径向生成器也不能当A30 full-field模型。A30 Sparse64的35.69%是truth-fit容量，不是在线解或观测相容性。

二维历史四parent×weak/ring共八条件，仅三个weak有合法base。有限编辑多出的一个材料成员不满足新多尺度结构差异；在dS=.04/.08/.16三个预定门槛下，所有archive仍只有三个结构代表，即各合法base。没有A30新增有效结构。

## raw → corrector → validation → archive

|原方法|raw材料代表/观测|修正endpoint代表/观测|full-valid代表/观测|原archive含base|raw新颖但修正回前缀邻域|
|---|---:|---:|---:|---:|---:|
|静态|10.5000|6.53125|5.90625|6.03125|145/384|
|multi-start|11.5625|2.56250|2.56250|2.56250|331/384|
|population|10.7500|2.28125|2.28125|2.28125|341/384|
|MDN|10.5000|2.43750|2.43750|2.43750|305/384|
|adaptive|10.78125|6.21875|5.62500|5.81250|165/384|

邻域折返不是“同一个动力学吸引盆/同一连通分支”的证明。普通径向raw residual只有proxy，因此其full invalid→valid为未知；不能把proxy→full转换当两次full验证。

二维126轨迹、318个非base full trial中，有限编辑32条的full invalid→valid为0。保护格从数组重建后最大偏差精确为0。31/32有限编辑没有形成≥6格封闭数字新洞；唯一形成并保留的新洞属于ring失败条件，最终残差约53.025，仍不相容。这里失败主要是“提案不构成所声明尺度的新洞”或“洞保留但不相容”，不能说corrector把所有洞填回去了。

## archive失败 → 修改 → 相同缓存成本重放

654个合法duplicate端点曾比匹配留存代表残差更小，其中644个优于全部匹配代表。它们没有被证明都是有意义的新结构；first-member的直接问题是拟合质量代表选择。

修改只允许truth-free低残差近邻替换，保持对其余代表的pairwise packing；旧端点进入immutable success ledger，不因替换重获发现奖励。相同192个缓存episode执行369次替换，没有失去原结构邻域覆盖。二维累计残差改善1.92094，但多尺度结构仍3→3；旧6格计数3→4是门槛/packing效应，不计新发现。

该循环已经完成，但只评价相同已生成池的archive选择，无法重建改变parent后未观察到的adaptive提案。E1另以新增物理、相同状态预算比较first与quality两规则。

## 历史是否帮助找到下一结构

按32个parent分组五折，小logistic固定模型与训练fold归一化：主标签是corrected/full-valid endpoint相对全部严格前缀full-valid ledger的物理几何新颖性，不是当前结果或latent差异。A含测量、当前材料、动作和预算；B加前缀历史；C为匹配观测历史置换；D仅历史。

结构B−A log-loss=+0.003786，描述性parent区间[-.031108,.040507]；B−C=+.003490，区间也跨0。材料重复/合法性有预测关联，但不能推出可选择更好的动作。旧adaptive proxy调用为static约3.79倍，原材料/结构发现仍较少。没有logged propensity和反事实端点，不做off-policy价值估计，也不据此训练PPO/GRPO。

详见SEARCH_HISTORY_SIGNAL_DRAFT.md与结果JSON；首次无阻尼Newton数值饱和保留，只有优化器backtracking修复，没有按test表现挑模型。

## 连续表示与边界

固定软件phase-field的一条连续系数路径可从无孔变为有孔，不输入孔数标签。它只说明连续表示允许阈值拓扑变化；有限空间带宽、固定支持与训练分布仍限制可达结构。软件图绝不称Maxwell相容成像成功。完整材料场、固定复contrast尺度与最小6格特征用于G评价，负实contrast也被计入。

原始逐trajectory表、完整资产引用、敏感性、缺失字段与独立回执见E0_DIAGNOSTICS.json及results/E0_*。T当前无已提交handoff；后续若使用历史后端，须保持同target并独立计费验证合法base。
