# 相关先例与本轮主张边界

这不是一次穷尽文献综述。此前实际ScholarQA检索与身份核验保留于AUDIT/literature_collect.json、literature_verify.json；429限制也保留，不能把检索空缺当原创性证明。下面使用原论文或作者版本。

|成熟工具/先例|原始来源|本轮不能据此单独声称的创新|
|---|---|---|
|Flow matching|[Lipman等，2022](https://arxiv.org/abs/2210.02747)|用CFM回归全维速度不是新Maxwell理论|
|Simulation-based conditional flow matching|[Dax等，2023](https://arxiv.org/abs/2305.17161)|用joint模拟标签训练条件posterior不是新概率合同|
|Simulator feedback posterior flow|[Holzschuh与Thuerey，2024](https://arxiv.org/abs/2410.22573)|加入物理/模拟器力不是本轮首次提出|
|电磁条件flow逆散射|[Zhang等，2022，IEEE TMTT](https://doi.org/10.1109/TMTT.2022.3205890)|电磁概率Flow本身已有直接先例|
|Amortized Gaussian/Flow variational inference|[2024原作者版本](https://arxiv.org/abs/2305.20004)|混合方向KL/监督与后验拟合需对照其已有机制|
|预条件normalizing flows|[2021原作者版本](https://arxiv.org/abs/2101.03709)|引入局部GN尺度不自动形成新算法贡献|
|DINO|[JMLR 2025](https://jmlr.org/papers/v26/24-0745.html)|本地小型value/JVP模型失败不等于复现并推翻DINO|
|LazyDINO|[JMLR 2026](https://jmlr.org/papers/v27/25-0858.html)|本轮不是完整LazyDINO复现，不能宣称成本胜过该文献|
|误差界驱动的retrospective接受决定|[2014作者机构全文](https://eprints.gla.ac.uk/96718/)|用概率阈值决定是否加精不是本轮新原则|
|Delayed acceptance效率分析|[原作者版本](https://arxiv.org/abs/1506.08155)、[adaptive DA](https://arxiv.org/abs/1509.00172)|普通近似或训练误差模型不能自动叫certified，setup也不能省略|

原P1/R2/A17-X/F3的具体本地来源与固定commit见SOURCE_INVENTORY.json。有限resolvent、Schur、Woodbury、条件期望正交性、proper-score身份和Lie bracket基础均作为工具，不以首次提出的形式叙述。

## 当前有证据的增量

1. 两批真正新parent上的material-risk改善，并能由强CNN、Gaussian/空间打乱控制分离归因。
2. 原66维测度下的相对proper-NLL与样本energy-score改善，同时保留Full-PPC失败，避免把风险与概率混同。
3. 同缓存、同完整核、同随机数上的Maxwell数值势区间费用收益，含setup和失败边界。
4. 对phase-bin参数化、物理训练coverage冲突、Euler路径污染和动态bracket距离失败给出可复核反例。

这些是实质实验现象。尚未证明它们足以形成新的Maxwell特有理论与算法，尤其缺少response-calibrated完整条件分布、与最强同成本方法的统计效率优势，以及较大二维系统的独立证据。最终强IEEE TAP结论仍为尚未建立，而不是把工程规模当新颖性。
