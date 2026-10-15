# Research G E0：已保存轨迹的归因证据

本分析仅重放Research A保存的在线日志和材料数组。没有导入或构造Maxwell solver，没有重新优化，没有读取`__OFFLINE_TRUTH`或locked/test资料。径向Mie与A30二维TM始终分表；所有候选均为candidate-only，不赋posterior概率。

## 已发现的材料与结构不能混为一项

32个径向确认parent、每个5方法、每方法12条proposal，共1920条轨迹。原体积加权材料L2/.03去重得到静态193、multi-start82、population73、MDN78、adaptive186个archive成员；这些是**各方法跨parent计数**，不能再相加宣称748个不同物体世界。对应原r/h几何分离计数为160、56、52、51、154。原物理验证及archive可逐项精确重放。

|方法|raw材料代表/parent|修正endpoint材料代表/parent|full-valid endpoint材料代表/parent|保留archive/parent|raw相对prefix新颖、修正后落回prefix的轨迹|
|---|---:|---:|---:|---:|---:|
|静态|10.5000|6.53125|5.90625|6.03125|145/384|
|multi-start|11.5625|2.56250|2.56250|2.56250|331/384|
|population|10.7500|2.28125|2.28125|2.28125|341/384|
|MDN|10.5000|2.43750|2.43750|2.43750|305/384|
|adaptive|10.78125|6.21875|5.62500|5.81250|165/384|

Raw阶段与endpoint阶段的材料不同不保证最终是新结构。这里collapse只指冻结原材料metric的邻域折返，不是证明进入了同一动力学吸引盆或同一连通解分支。径向raw residual只有proxy，因而普通1920轨迹的full invalid→valid必须为null；JSON另外标记proxy→full转换，绝不将其冒称两次full验证。单独已付费8个no-comp控制见原文件记录。

证据：`research/research_a_20261015/results/CONFIRM_CAMPAIGN.json#parent=confirm_00..31/method=<arm>/attempt=0..11`；新转移表每行保存完整source引用。旧归因8个独立补偿初始化端点距原endpoint最大2.83877e-5，小于原.03；497对结构不同候选中210对gain后白化预测距离≤1。其原文件为`STRATIFICATION_ATTRIBUTION.json`、`PAIRWISE_AMBIGUITY.json`。

## A30二维失败不应描述成“corrector把洞填了”

4parent、8条件，共318个非base full trial state，形成126条轨迹（finite32、multi32、weak32、no-comp30）。三个有合法base的weak条件及五个无合法base条件完整保留。finite archive计数仍为[0,0,2,0,1,0,1,0]，其它multi/weak为[0,0,1,0,1,0,1,0]。所有方法的all-case Coverage@4=1/8；无合法候选条件的MAE是null，整体8条件平均MAE仍为null。

finite32条轨迹：full invalid→valid=0；22条raw与base存在至少6格旧mask差异，22条修正后仍存在这一级别差异。重新从保存center/radius构造受保护mask并检查全部原始/中间/末端数组，保护格最大缺陷精确为0。它们没有被修正器填回。

G-multiscale-support-v1采用复contrast幅值，σ=[0,1,2]格constant-zero平滑、固定occupancy变换，dS=.04/.08/.16并同时要求raw mask XOR≥6格。三档均只有三个archive结构代表（就是三个合格base），finite额外材料不是新结构。原6格定义和新metric均单独保存，没有覆盖历史数字。

新增post-hoc digital-hole诊断使用foreground4/background8、最小6格封闭背景组分：finite31/32没有符合该定义的raw新洞；仅`p250077_ring_a1_snr20/trajectory0`形成一个raw新洞并保留到末端，full residual=53.02497，仍不相容。multi/weak各3个raw新洞在修正中消失，但这些是global/weak起点，不能用它们断言受保护有限编辑填洞。受限阈值和有限尺度下没有洞不意味着连续材料空间没有洞。

## First-member确实丢了拟合质量，但没有据此新发现世界

发现654个合法duplicate endpoint的full residual比至少一个匹配留存代表更低：径向652、二维2；其中644个优于全部匹配代表。绝大多数径向改进极小。二维两例均在`p250069_weak_a1_snr20/finite_history/trajectory2`：原s0 residual22.8953093被保留，s1=20.9757067和s2=20.9743657被判duplicate。

已登记并执行同缓存池archive-only修复：在同原材料metric下，duplicate若full residual更低且替换后与其它代表仍pairwise≥.03，可替换最近代表；immutable success ledger仍保留全部端点，重插无新发现奖励。192个episode执行369次替换，旧结构邻域覆盖丢失0。静态active193→193，multi82→83，population73→72，MDN78→78，adaptive186→186；这些active数量变化是packing顺序变化，不是新物理世界。

二维finite两次替换累计残差下降1.9209436；旧6格几何计数3→4，但G多尺度结构仍3→3，说明这种微弱阈值计数不能称新结构。径向adaptive旧r/h packing154→157，而多尺度182→182，同样不能宣称新发现。此次只重排固定已生成池，不能重建更换parent后未观察到的adaptive提案，更不构成新的search policy优越性试验。

## 产物与限制

- `results/E0_TRAJECTORY_AUDIT.json`：分平台、分阶段、逐parent/方法完整计数、尺度敏感性及证据路径。
- `results/E0_TRAJECTORY_TRANSITIONS.csv`：1920径向+126二维proposal→endpoint转移。
- `results/E0_FULLSTATE_ARCHIVE_TRANSITIONS.csv`：318个二维full state的valid/archive状态。
- `results/E0_FIRST_MEMBER_RESIDUAL_IMPROVEMENTS.csv`：654条质量丢弃证据。
- `results/E0_ARCHIVE_REPAIR_REPLAY.json/.csv`：同缓存、同费用的archive修复比较。
- `code/e0_trajectory_audit.py`、`code/e0_archive_repair.py`：可复算代码。

每个timed执行回执都保留，累计费用见同目录`RESOURCE_ROLLUP.json`。历史solver/constitutive diagnostics未由G重新验证，缺失trust-region-exit等字段保留null/MISSING。grouped-by-parent历史预测由独立分析处理，本脚本不生成RL有效性判断；最终科研裁决由主Codex负责。
