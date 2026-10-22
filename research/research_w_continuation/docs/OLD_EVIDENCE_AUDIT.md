# 旧 P0 证据审计：已知结果、离散接口与续接边界

审计日期：2026-10-22。审计对象是相邻只读目录 `research-w-p0/`。本次只读取文字、源码、小型 JSON 与账本；没有载入 NPZ 大数组，没有运行物理解、参数搜索、测试或新文献检索，没有改写旧目录。唯一写入是本文。

**旧结果必须继承为：P0-A PASS；P0-B `P0_B_UNRESOLVED_NUMERICAL_QUALIFICATION`；P1-A/P1-B HOLD。** 旧记录没有合格二维 calibration/joint 端点，没有搜索成功率、posterior、递进算法或加速结果。前向资格失败不能改写成物理唯一性或歧义不存在。核心数学、当前续接 gate 与最终科学解释由主线程裁决。

## 1. 先分清两层来源

|来源层|已经做过什么|不能混用的时间语义|
|---|---|---|
|原始 Pro 包 `source_package/`|一维无限平面、无损、非磁板的完整 Maxwell 解析实验，理论与来源审计，P0/P1 协议建议|`README_ZH.txt:29–30` 与 `P0_P1_PROTOCOL_ZH.txt:2` 的“未执行 P0-B/P1”描述的是原始 Pro 报告，不否认后续旧 P0 的二维资格计算。|
|旧 P0 实施目录|独立复现 P0-A；执行两个不规则有损墙体的二维多频全波资格、源接口修订、有限四块诊断；因 Gate1 未通过停止|已执行二维离散求解不等于执行见证搜索。搜索/endpoint 阶段费用均为 0；候选数为 0。|

原始 PDF 为 52 页。已读取其现成文字抽取 `docs/SOURCE_REPORT_EXTRACTED.txt` 的相关段，并与 LaTeX 的实验/判断章节交叉核对；没有重新渲染或重编译 PDF。PDF 阅读器物理页码比正文印刷页码大 1：物理第 38–40 页对应正文第 37–39 页，包含“已执行/未执行”、P0-B 门槛和 P1；物理第 45 页/正文第 44 页保留后续步骤与停止规则。原始 LaTeX `experiments.tex:1–5,22–27` 将已完成的一维回归与建议的 P0-B/P1 区分开。[S08, S09, S10]

## 2. 旧 gate 与精确已知数值

### 2.1 旧 gate，不作本次重新判定

|旧 gate|原记录状态|证据位置与含义|
|---|---|---|
|0，一维解析复现|`PASS`|完整参考面、PEC 级联、幺正性与带宽负结果独立复现。|
|1，二维前向资格|`UNQUALIFIED_SOLVER`|匹配源平均修订后，最细相邻网格的完整数据白化差最小 0.9508595591879654、最大 1.2291015897405526；要求 0.1。|
|2，calibration 机制|`NOT_RUN_UNQUALIFIED_SOLVER`|没有启动候选搜索；没有有限墙体合格见证。|
|3，target-present joint|`NOT_RUN_UNQUALIFIED_SOLVER`|没有未知目标/仪器的联合补偿实验。|
|4，带宽/扰动鲁棒性|`UNRESOLVED_2D`|只有一维指定补偿对的带宽负结果；二维扰动、孔径变化未运行。|
|5，算法/P1 价值|`HOLD_P1_NOT_AUTHORIZED`|P1-A/P1-B 均 HOLD。NN、3D 均 `NOT_RUN_NOT_AUTHORIZED`。|

以上逐字状态来自 `results/GATE_DECISION.json:2–24`；其中 `search_executed=false`、`candidate_count=0`、`budget_exhausted=false`。[S01, S02, S03]

### 2.2 一维证据：正确，但范围有限

两个板 `(n,d)=(2,1/4),(3,1/6)`，共同前表面和后参考面 `D=1`；PEC 位置 `1.5` 与 `1.5−1/12`。独立 P0-A 的中心最大反射模为 `1.974546727233598e−16`，同位置双程响应差模 `1.0000000000000007`、相位差 `1.0471975511965983` rad，位置补偿后联合差 `5.071544317185128e−16`。完整 PEC 边界与级联的最大差为 `1.9860273225978185e−15`。[S11]

|相对全带宽|墙校准白化距离|指定位置补偿联合距离|
|---|---:|---:|
|0.2%|0.3548765484820979|1.216503750515063|
|0.5%|0.88713775173482|3.038246262489969|
|1%|1.773892641321239|6.055071158130078|
|2%|3.544726568359131|11.94185858505982|
|5%|8.808752726311752|27.215436721305647|

每档 9 点，每实/虚正交分量噪声标准差固定 0.01；距离为完整复欧氏范数除以 0.01，不再除以通道数平方根。每档频率重新布点，非嵌套。0.5% 带宽已破坏该组指定联合补偿；这只排除这一个 pair，不能证明一般联合唯一性。位移 `1/12≈0.08333λ0` 还小于二维预先声明的 `0.1λ0` 任务阈值。[S04, S11]

原 Pro 归档的中心反射/联合差分别为 `1.6328623988631375e−16`、`2.449293598294703e−16`。它们和独立 P0-A 都在机器精度量级，但来源必须分别标为“Pro 归档原结果”和“P0-A 独立回归”，不能把原数值写成本次结果。[S05:168, S12]

### 2.3 二维实际资格与离散诊断

固定物理为二维 scalar TM、`exp(−iωt)`、出射 Green VIE，单位线源谱；`D=[−2,2]×[0,3]`，同侧 4 发 12 收，宽带 `{0.8,0.9,1,1.1,1.2}`、窄带 `{0.9975,0.99875,1,1.00125,1.0025}`，共 9 个不同频率。跨频共享 `ε∞` 与 `σ0`，`εr(f)=ε∞+iσ0/f`。查询域 `T=[−1.2,1.2]×[1.6,2.8]`；calibration fixture 将 T 置为空气，实际资格的 target 状态只有同一个固定 probe，不是未知目标推断。[S04:56–66, S06, S18]

修订版 `nx128→192`（即 `128×96→192×144`）完整多频数据白化差如下。[S13]

|seed|宽带 cal|宽带固定目标 total|窄带 cal|窄带固定目标 total|
|---|---:|---:|---:|---:|
|2026102201|1.174762158610265|1.2291015897405526|0.9508595591879654|1.0195011713328501|
|2026102202|1.1275042536096458|1.2163797692448675|1.015356922125296|1.0344399929216679|

这些量超过 0.1 门槛约 9.5–12.3 倍。经验网格差不是连续 PDE 误差上界；此处连预设的一致性要求也没有通过。

|已记录检查|准确结果|只支持什么|
|---|---|---|
|材料方向差分/实伴随|四步长最佳相对差 `3.5098153664790935e−10`；dot 缺陷 `1.9118587546358895e−15`|当前离散接口一致；不能替代前向空间精度。|
|源收互换修订|原缺陷 `3.488290448943158e−7`，平均入射修订后 `7.203885968749396e−16`|当前离散实现的互易数值性质；不能证明连续精度。|
|均匀圆柱解析参考|原 nx32/64/96/128 相对误差 `84.0869%/22.7884%/9.90988%/7.03479%`；平均入射 nx96/128 为 `9.84519%/7.00361%`|外部解析参考下的综合有限网格误差，未满足资格所需精度；误差来源未拆分。|
|内存修订 parity|同材料相对数据差与白化差均为 0；旧 peak RSS `2054520832` bytes，修订 job `572342272` bytes|分块/原位 LU 工程等价，不能改变物理 gate。|
|四块有限诊断|8 项；nx64；2 seed×2 frequency×2 contrast；21 个任务支撑端口；最大直接全波相对差 `8.6675025159728e−16`|同一离散系统的公式/反馈诊断，非整个 T 的签名、未知墙反例或加速。|
|Born 目标诱导近似|上述 8 项误差 `100.5518%–337.4853%`|这些固定强目标中省略内部反馈的代价，不是全目标类统一结论。|
|标准库测试|旧日志 `Ran 20 tests ... OK`；保留一次缺少 pytest 的启动失败|测试通过不等于物理见证成立。|

来源：导数 [S14]；互易 [S15]；圆柱 [S16]；parity/失败 [S20]；四块源码及数值 [S19]；测试日志 [S22]。

两个固定 seed 在最细修订网格的 cal 距离宽/窄带分别为 `37.984019069980654/35.94431409381102`；相同固定目标 total 距离为 `37.913405476322936/35.65521213101644`，目标诱导差为 `3.6946008796890477/4.685430606366778`，两世界目标质心差为 0。它们是保存参考世界，不是搜索端点；既不近数据，也没有目标任务位移。target-induced 差依赖额外公开的 calibration fixture，不能偷偷变成 joint 在线输入。[S21]

## 3. 源码核对：确定观察与待澄清处

本节提出可复核的实现观察，不把它们认作已测量的主导误差原因，也不作新的数学有效性判决。

### 3.1 “源平均修订”只改了入射表示，不能直接称完整 Galerkin

`vie.py:109–118` 显示：

- 传播矩阵 `G` 是源单元积分、观测点位于单元中心的平移不变核；自项是中心到真实正方形单元的积分。
- 接收算子 `S` 是点接收位置到源单元的积分，带 `k²`。
- 原入射是自由空间 Green 在单元中心取值；`source_cell_average=True` 后，入射改为该 Green 的单元面积平均，除以 `h²`。
- `G` 的观测端没有同时改成单元平均，也没有双重测试/源单元积分；`solve` 仍使用单元常数对比和 `I−diag(χ)G`。[S17]

所以代码确实实现了旧报告所说的“匹配入射/接收平均”，并测得互易缺陷下降。可是“互易性达到机器精度”和“变量究竟表示中心值还是单元平均值”是两个问题。续接若更换离散方案，应明确未知量、入射、传播核、材料投影和测量算子的同一约定；不能仅复制这个 flag 就宣称已经采用一致的面积平均 Galerkin。

`vie.py:12–13` 的模块说明仍只写“incident at collocation centers”，与可选平均分支并存，属于需要在新文档中消歧的接口说明。旧源码保持原样。

### 3.2 圆柱 Mie 比较是有效的连续参考，但有限网格材料不是解析圆盘本身

`qualification.py:31–36,53–65` 的事实为：

1. 解析参考是半径 0.3、中心 `(0,1.5)`、`εr=4+0.03i` 的均匀尖锐圆形界面；4 个平面波入射角在 1.3–1.8 rad，频率只有 `f/f0=1`，观测点使用原 rx。它不是墙任务所用线源采集的解析参考。
2. 数值端用 `disk_fraction` 计算正方形网格与圆盘交叠面积，再令每个完整方格的系数为 `1+(εr−1)×fraction`。切割单元成为中间介电率的整格材料；场/电流再按现有单元常数近似求解。
3. averaged 版本的平面入射乘两个 sinc 因子，明确改成方格平均；G 仍沿用上一节的中心观测/源单元积分。解析参考仍为连续锐界面圆柱在点接收器的场。
4. 解析展开固定 `m=−28,...,28`，旧文件没有单列模式截断收敛量；本次也没有新做该检查。[S18]

这是针对同一个目标连续圆盘的网格逼近检查，而非同一有限离散材料的恒等式。7.00361% 可以作为当前整体近似未够准确的外部负证据；它混合了界面/系数投影、单元场近似、积分及参考截断等可能来源，不能单独归因于 LU、入射平均、圆盘面积积分或某一种“solver bug”。修订前后只有约 0.03118 个百分点的相对误差变化，旧记录不足以宣称源平均解决主要误差。

### 3.3 墙体系数平均与网格差也尚未分账

墙体 `material_grid` 使用每单元 `4×4` 个中点采样，跨频从同一连续 seed 生成系数；这不是对不规则界面的精确面积积分。网格加密同时改变切割界面系数的表示和场空间。`geometry.py:27–43` 已注明“coefficient discretization, not a continuous-PDE error certificate”。旧报告 `RESEARCH_REPORT_ZH.md:95–97` 也明确没有量化几何面积平均与常数单元场各自占比。[S04, S18]

圆柱的 `disk_fraction` 与不规则墙的 4×4 系数采样不是同一投影方案。不能把圆柱误差降幅机械外推为两墙的误差率，也不能把相邻墙网格差当作只含场离散误差。

### 3.4 辐射条件、任务覆盖与公平比较的边界

出射 Green 在无限均匀外部表达 Sommerfeld 辐射条件，因此没有 PML/人工盒边界参数；旧协议用此解释替代了原 Pro 建议中的 PML 扩展，明确 PML 为 `NOT_APPLICABLE`。这说明边界建模选择，不能把“没有 PML”记作 PML 验证 PASS，更不能推出内部空间精度合格。[S06:5, S07:15–21]

24 个声明 probe 中，墙网格资格只跑 `probes()[8]`；四块只跑 `[8,9]`（同位置、两对比），不是整个目标类覆盖或构建/验证 probe 的完整分离核验。四块 setup 需 21 个端口 RHS，而一个目标 direct 仅 4 个源 RHS，且另有小核因子与 4 RHS；这些记录没有证明相同总费用下获益。[S18:109–135, S19]

## 4. 费用与失败：三个快照不能互换

旧 prose、旧 `COST_SUMMARY.json` 和当前旧目录账本结束位置不同；计数相同不代表尾项计时相同。本次仅用标准库读取 JSONL 并加总字段，未调用物理模块。[S04:124–128, S23]

|来源快照|job wall / s|process CPU / s|charged / s|说明|
|---|---:|---:|---:|---|
|旧中文报告正文|约 297.6|约 288.0|约 897.6|报告当时快照，文字已声明尾项继续追加。|
|旧 `results/COST_SUMMARY.json`|352.72315450338647|294.114513|952.7231545033865|包括原 Pro 解析复跑，但尚未包括最后的解析复跑发布包装。|
|本次对旧账本只读加总|359.0694514203351|294.617934|959.0694514203351|共 1473 行；最后 cost 为 `p0a-rerun-publication-001`，`2026-10-22T09:28:21Z`。|

charge 额外含 `600 s` 未分项早期工作保守分配，状态是 `CONSERVATIVE_ALLOCATION_NOT_MEASUREMENT`，不是测得的物理 CPU；包装工具按 wall 计费也不能改称 process CPU。最后 6.346296916948631 s 包装尾项使 summary 比账本少一项。[S23]

三个快照的物理计数一致：`1184 complex RHS`、`266 factorization`（含 8 个实际目标小核）、`258 native states`。validation 为 `820 RHS / 178 factors`，repair 为 `364 RHS / 88 factors`；search、endpoint 均为 0。RHS 角色为源前向 1032、接收转置 36、四块端口 84、小核 32。旧总预算是 `8000 RHS / 2000 factors / 7200 charged s` 任一先到即停，并未耗尽。

失败继续保留：系统 Python 缺 NumPy、一次 pytest 启动失败、原离散源收不对称、旧内存峰值超工程目标、两个 seed 的资格失败，以及首次公开回读检查器失败。发布回读成功只是资产可访问记录，不改变科学 gate；本次没有在线验证旧网址。[S20, S22, S24]

## 5. 续接阶段的 plans-only 边界

|名称|旧文件实际定义与状态|本审计允许继承的边界|
|---|---|---|
|P0-C|旧 Pro/P0 包没有定义这个阶段编号或相应授权。|如果续接使用此编号，必须由主线程另行定义目的、费用、停止规则与产物。本文只是静态审计，不是 P0-C 物理 PASS，也不构成运行授权。|
|P1-A，概率传播等价|旧 `P0_TO_P1_DECISION.md:7–9` 为 HOLD。先要 qualified joint 端点，再考虑可独立积分的低维联合 prior reference。|事件/CDF 误差≤0.02、同准确度总费用减少≥20%均只是未来审核目标；没有已建立后验或压缩证书。|
|P1-B，任务递进价值|旧决策 `:11–13` 为 HOLD；未来要做单点/多解×joint/progressive 的 2×2。|同数据、同 prior、同费用；若只有多解增量而无递进增量，删除递进优势。旧 P0 没有这项实验。|
|P2|旧包没有定义 Research W 的 P2 执行协议、预算或授权。|源文件中的“P2”多指既有 B1/Joint P2/材料 fiber 历史证据，不是 Research W 下一阶段已完成。任何新的 P2 在本审计中仅为计划。|
|P3|旧包没有定义 Research W 的 P3 执行协议、预算或授权。|“P3-L”是原 Pro 明确未取得完整原稿的既有方向，不能充当 Research W 的已证理论/执行成果。任何新的 P3 在本审计中仅为计划。|

具体当前 P0-C/P1/P2/P3 的定义与 gate 以主线程续接方案为准；本文只限制旧证据可以支持什么，不替主线程新增科学阶段。

按旧决策，接下来首先是同一两个世界、同一采集/噪声/0.1 条件下的前向资格；四种搜索和 joint corrector 当时既未完成实现，也未运行。只有资格通过才有恢复搜索的依据；若之后 qualified 证据仅窄带 CAL 成立而宽带 JOINT 消除，才讨论 `REDIRECT_TO_ACQUISITION_IDENTIFIABILITY`；如果有 qualified JOINT 歧义，才考虑 P1-A。旧二维结果尚未支持任何一种分支。[S03, S06:13–17]

## 6. 可继承表述与禁止升级

可以继承“指定一维完整 Maxwell 解析见证已独立复现”“二维实际离散全波资格失败”“当前离散伴随/互易/四块诊断通过”。不能升级为“有限复杂墙体已发现歧义”“宽带已消歧”“二维机制 NO-GO”“已有可靠 posterior/mode probabilities”“递进、A17 或 ROM 取得同成本优势”“整个查询域四块签名已被覆盖”。

理论审计中，纤维因子化、附条件微分结论、非交换 resolvent、保持联合先验的 TV 界均按原条件保存；它们不等于新颖性、统一稳定常数或有限墙后验证书。本文没有新增 prior-art 检索、数学证明或最终 novelty 裁决。[S05]

## 7. 精确来源指针

以下链接均指向本地旧目录；`Sxx:行号` 用于上文压缩引用。行号是本次读取所见，不以发布副本替换旧证据。

|编号|来源与关键行|
|---|---|
|S01|[START_HERE.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/START_HERE.md#L3)，3–13。|
|S02|[GATE_DECISION.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/GATE_DECISION.md#L5)，5–14；[机器 gate](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/GATE_DECISION.json#L2)，2–24。|
|S03|[P0_TO_P1_DECISION.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/P0_TO_P1_DECISION.md#L3)，3–19。|
|S04|[RESEARCH_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/RESEARCH_REPORT_ZH.md#L3)，3–25、29–50、56–118、124–148。|
|S05|[THEORY_AUDIT.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/THEORY_AUDIT.md#L7)，7–13、15–131、133–186。|
|S06|[PHYSICS_PROTOCOL.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/PHYSICS_PROTOCOL.md#L3)，3–17。|
|S07|[PHYSICS_VALIDATION.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/PHYSICS_VALIDATION.md#L5)，5–25。|
|S08|[Pro 执行提要](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/source_package/P0_P1_PROTOCOL_ZH.txt#L2)，2、9–11、13–72、74–108；[Pro README](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/source_package/README_ZH.txt#L16)，16–30、45–50、65–73。|
|S09|[Pro experiments.tex](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/source_package/report_parts/experiments.tex#L1)，1–5、22–96、98–136；[judgment.tex](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/source_package/report_parts/judgment.tex#L45)，45–68。|
|S10|[原 PDF](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/source_package/Research_W_Report_ZH.pdf)，物理第 38–40、45 页；[现成 PDF 文字抽取](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/docs/SOURCE_REPORT_EXTRACTED.txt#L1897)，1897–2088、2250 附近；[evidence.tex](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/source_package/report_parts/evidence.tex#L30)，30、92，P3-L 原稿缺口。|
|S11|[独立一维结果](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/p0a/results.json#L25)，25–99。|
|S12|[Pro 一维归档](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/source_package/analytic_results/results.json#L19)，19–26、86–91。|
|S13|[修订墙网格记录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/qualification_source_average/WALL_GRIDS.json#L10)，10–42。|
|S14|[修订导数记录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/qualification_source_average/DERIVATIVES.json#L2)，2–21。|
|S15|[互易修订记录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/qualification/RECIPROCITY_REPAIR.json#L2)，2–13。|
|S16|[原圆柱记录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/qualification/CYLINDER.json#L2)，2–37；[修订圆柱记录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/qualification_source_average/CYLINDER.json#L2)，2–19。|
|S17|[VIE 源码](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/src/research_w/vie.py#L12)，12–13、30–57、99–125、169–207。|
|S18|[qualification 源码](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/src/research_w/qualification.py#L11)，11–36、53–68、82–95、104–135；[geometry 源码](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/src/research_w/geometry.py#L27)，27–47；[冻结配置](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/w_config.json#L372)，372–453。|
|S19|[四块源码](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/src/research_w/four_blocks.py#L17)，17–50；[四块结果](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/four_blocks/RESULTS.json#L1)，1–154。|
|S20|[实现修订说明](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/docs/IMPLEMENTATION_REVISIONS.md#L3)，3–10；[组装 parity](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/qualification/ASSEMBLY_PARITY.json#L2)，2–5；[失败账本](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/FAILURE_LEDGER.jsonl#L1)，1–7。|
|S21|[固定参考 pair](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/REFERENCE_PAIR.json#L2)，2–26。|
|S22|[旧测试日志](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/jobs/tests-002/OUTPUT.log#L36)，36–39。|
|S23|[费用 summary](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/COST_SUMMARY.json#L2)，2–37；[物理账本](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/results/PHYSICS_LEDGER.jsonl#L1473)，共 1473 行，最后 cost；[summary 加总源码](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/tools/summarize.py#L72)，72–86（本文未运行该脚本）。|
|S24|[旧发布回读](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/db7146cae2d74ba231577d42f225dbfc9a389268/research/research_w_p0_20261022/publication/WEB_READBACK.json#L2)，2–5、7–54；仅引用历史记录，不声称本次在线可达。|
