# R13 最近先例与原创性边界

检索日期2026-11-12（新加坡）。这是有限范围的**先例抽取与边界审查**，不是穷尽检索或最终新颖性结论。核心数学/科学裁决由父线程完成。

## 检索协议与错误

按 ScholarQA 工作流运行一个六query、每query最多4篇/2段的bundle，原始结果为 `scholarqa_collect.json`。六个query依次针对 *Searching with iterated maps*、*Matrix product constraints by projection methods*、*An Analysis of Constraint-Relaxation in PDE-Based Inverse Problems*、NFFT-SOM、CRM、Anderson/DR。**12个请求中的11个返回429**；唯一成功的paper search找到Elser canonical record及其重复项。保留全部 `operation_errors`，没有将失败写成“未发现先例”。

随后从包内明确标识的种子以一次 `scholarqa-cli verify` 批量核验7篇，`scholarqa_verify.json` 显示7/7 resolved，无 unresolved IDs。批量核验只证明书目身份。再逐篇读取原作者arXiv/会议PDF等主来源填补全文/摘要证据；没有改成无来源的搜索摘要综述，没有作大范围引用网络扩展。NFFT会议文献无已验证DOI，保留直接会议PDF而不补造DOI。

预印本年份与期刊年份不同：Elser预印本2016/期刊2017；van Leeuwen–Yang预印本2024/期刊2025；CRM预印本2019/修订2021/期刊2022；Themelis等预印本2020/修订2021/期刊2022。不把Semantic Scholar首次年份强改为期刊年，也不因年不同判作错误论文。

## 逐因素最近先例

|因素|独立核验主来源与阅读等级|直接重合|本轮需要额外建立|
|---|---|---|---|
|最近点对法向漂移、累积不一致|[Elser, *Matrix product constraints by projection methods*](https://arxiv.org/html/1601.01003v1)，§4全文A；[期刊DOI](https://doi.org/10.1007/s10898-016-0466-9)|§4明确给出陷阱附近每次增加投影差以及discrepancy accumulator解释；§3.3已有rank-one SVD投影|不能以“DR保留误差/反射逃离AP陷阱”作原创；应比较有界共享标量的首事件公式及可用后果|
|早期difference map背景|[Elser–Rankenburg–Thibault2007](https://doi.org/10.1073/pnas.0606359104)，书目已批量核验；主文搜索抽取可见，直接PMC CAPTCHA|迭代投影作非凸搜索已有清楚前史|精确normal-ray先例由上行Elser§4全文支撑；不依赖本次无法完整打开的PMC页面|
|完整current/NFFT与DPCS释放|[Wei–Chen2018会议原文](https://usncursi.org/archive/aps-ursi/2018/papers/1927.pdf)，两页全文A|Eq(3)为`J=J_s+F alpha`，完整Fourier校正可改变deterministic current；结论明确讨论修正其噪声误差|解锁不是创新。移植到同一个G19材料块、beta、PR与审计日程可作关键因果控制，但不能冒充原始NFFT方法复现|
|profile/receiver Gram/伴随representer|[van Leeuwen–Yang](https://arxiv.org/html/2403.15292v2)，Theorem1、Corollary1全文A；[期刊DOI](https://doi.org/10.1088/1361-6420/adab89)|消去扩展场后得到材料相关receiver-Gram残差metric，扩展源在接收伴随场张成空间|R13的完平方恒等式不是新的variable-projection理论。Maxwell AP gap与事件联系可作为明确特化，价值须由条件/预测/算法实用性承担|
|Anderson/准Newton反射加速|[Themelis–Stella–Patrinos](https://arxiv.org/html/2005.10230v2)，AssumptionI与§3.3.2全文A；[期刊DOI](https://doi.org/10.1007/s10589-022-00366-y)|明确给出Anderson多割线方向、有限记忆和DRS线搜索|平滑第一项条件不适用于任意两个indicator；可以做明确标注heuristic控制，不能移植不满足条件的收敛保证|
|CRM几何外推|[Dizon–Hogan–Lindstrom](https://arxiv.org/abs/1910.04384)，作者摘要B；[期刊DOI](https://doi.org/10.1007/s11228-021-00626-9)|非凸反射外推与局部加速，数值误差累积敏感性|不能泛称首次几何外推；本worker未审完其全部定理，摘要不证明CRM对当前Maxwell全局有效|
|最新有限DR外推|[Collard–Lindstrom2026原PDF](https://arxiv.org/pdf/2611.02694)，§4.1 Theorem16、§4.2和§5相关部分全文A|两个仿射子空间的DR差分空间、有限步恢复DR极限、LSM/LeSE/DR-GMRES关系；QR减轻Gram平方病态|Pro当时全文不可得，本次已可读。其仿射前提与非凸材料盒事件不同，但否定宽泛“首次DR有限外推”。尚不构成对材料首事件公式的逐式否定|
|Maxwell敏感度|[Darbas–Heleine–Lohrengel2019](https://arxiv.org/abs/1911.02409)，主来源摘要B；[DOI](https://doi.org/10.1080/17415977.2019.1588896)|固定频率三维电参数灵敏度与小幅扰动重建已有先例|不是全局有限对比可辨识性先例证明；标准Maxwell derivative/adjoint不能作为本轮主要创新|

## 与Pro包相比的更新及未解决空缺

**新增可用证据：** Collard–Lindstrom完整27页PDF本次可访问，原协议的“只能读摘要”限制可以更新为本次已审相关定理。Theorem16在两仿射子空间设定中，以差分空间维数s构造LSM_s精确恢复固定点；§4.2明确实际数值秩检测与病态可能妨碍有限步恢复。其结果没有把一般非凸投影变成全局收敛算法。

**仍未解决：** 未完成“有界共享标量rank-one锥的首事件公式”在全部低秩投影/参数化集合/DR文献中的系统查重；没有因当前检索未见而认定原创。Pro35篇清单中的近期CIE两阶段、NIE2016、LSL等本worker未重新逐条全文验证，不能用继承书目的阅读深度冒充自己的阅读深度。它们可作为父线程需要时的明确种子。

**特别接近的内部前史：** R12已经推导共享材料二次斜率投影、2×2 Gram退化与分支切换；Formula Atlas已有精确双边材料导数与Hadamard Gram；R8已有current投影损失与anchor修补，R9已有完整current松弛。因此本轮新意最多是这些成熟构件中一个更窄、可计算并具后果的有限事件规律，而非组合名词本身。

## 可供父线程采用的保守定位

潜在贡献可以写为：在固定实各向同性共享材料盒及指定状态度量下，给出法向射线上的首全局投影分支事件，区分精确跳步与近似proposal，以Maxwell接收伴随结构解释gap，并独立检验其预测和候选生成价值。

目前证据**不支持**“首次DR反射机制”“首次DPCS释放”“新的变量投影基本定理”“全局Maxwell凸化/唯一恢复”“预测box flip即预测GT”“三例开发续算已胜既有DR”或“已验证端到端加速”。父线程应保留正例、负例和条件收益，同时以同空间/同付费审计/同成本的随机或一般外推控制检验新增几何是否有实际价值。
