# FBG-R3-E 独立最终科学审阅：主指标、次级增益与归因分开

状态：**全部5748次唯一推理执行及声明后的component控制已闭合评分；不再进行新数值实验。旧原型、新native-affine与两批独立确认分别归档。**

审阅对象：已闭合并离线评分的 E1 开发实验、独立新确认集、E2 条件多根实验，以及 48 个真实 576-cell E3 成熟基线任务。全部结论以有限协议为界；论文主张与最终新颖性判断仍由父 Codex 线程负责。

**已评分 anchored-N1 原型的结论：其任务方向选择、软释放与有限种子，在当前冻结预算内没有显示相对于强基线组合的新增恢复覆盖；N2 展示了不同的条件子问题根，但没有展示全数据根或新增正确恢复。独立确认集保留了这一负面结果。成熟 CIE 的真实 576-cell 成像则取得了明确正面结果，必须作为文献基线成果单独归属，不能包装成 N1/N2 的成功。所评分原型与 N2 的 E3 科学扩展门仍为 `NOT_RUN`；新 native-affine 的独立确认保留一个clean材料/数据次级增益，但主joint覆盖没有增加；post-disclosure normal/soft也能达到该案例，不能归为task评分或A17特异优势。**

## 表示族匹配修正与 native-affine 分支

父线程在旧分支确认完成后识别到一个不能仅由总 DOF 排除的结构混淆：6×6 anchored 原型的核心含16个可观测方向与约2个 dark 方向，而成熟 CIE 可优化20个 dark/APCS方向，且带固定 DPCS affine offset。即使某些总实自由度接近，也没有建立相同 affine current family、相同 dark 容量和相同优化机制。旧结果对所运行原型依然有效，不能据此裁定整个 normal-selection 家族失效。

父线程另建并冻结了 `native_tsom` affine 分支：固定数据决定的 `jD`，15个 APCS core 加1个被选择 APCS，共16个 complex spatial APCS＝192个实 source-current自由度（6源），与 native CIE 的 mf2 容量匹配；使用同一个 beta=.2、source-normalized joint objective、prior及共同corrector，并与 fixed Fourier GN、energy/random、task hard/soft 和 constant-beta/mf2 CIE 对照。该协议的660次受控开发执行现已闭合评分（下节）；当前观察到一个所有finite规则共享的恢复增益，仍没有task/normal selector特异优势。成熟完整 CIE 继续保留。

该变体属于新实现，不能把其结果追溯归入已确认旧分支，也不能将已经曝光的旧确认集重新称作它的 fresh confirmation。若继续，只能在公开 DEV 选择与方法冻结后使用另一个未读取的独立确认集，并单独报告成本、失败和修复。最终不能给出“所有 N1 均无basin增益”或“新 native-affine 已建立普适优势”的结论。N2 的条件根负面证据不受这一 N1 结构问题改变。

## native-affine 受控开发结果

本批660次新增执行＝11方法×12个6×6观测×5起点；另准确复用180个旧成熟CIE/fixed-feedback/Native参考行，不重复推理或收费。新分支和matched CIE为192实current DOF；所记录N1 mode替换容量均192→192，零执行异常、零物理计数超限。所有新方法4s/start、120LU/2400RHS；旧强参考的原始4s协议不全有相同物理caps，因此同时报告其作为较强参考的范围。

| 方法 /12观测 | K1 joint | K5 joint |
|---|---:|---:|
| 普通Fourier、energy、random、normal、task、soft | 各6 | 各6 |
| matched CIE beta=.2/mf2 | 6 | 6 |
| 旧mature CIE dual参考 | 6 | 6 |
| 旧fixed feedback参考 | 2 | 6 |
| finite task | 6 | 7 |
| finite Fourier | 5 | 7 |
| finite random | 7 | 7 |
| finite Krylov | 6 | 7 |

K5四个finite规则增加的恰为同一案例、没有损失：`obs_5ec587d4bc0898c9`，islands、幅度18、jittered。材料NRMSE约.0027789、IoU=.83333、原始数据残差约.00049316，data/material+data/joint均通过。普通Fourier、matched CIE、成熟CIE、固定feedback均未恢复；旧direct-material-seeds、GlobalDual、reciprocal/data也未恢复该案例。它是一个真实有限协议内的新增恢复，不能继续声称新native-affine分支完全无basin增益。

但该增益不是task方向特异：finite random在start0已恢复，task在start2、Fourier/Krylov在start4的在线best-Phi端点恢复；四者达到基本相同材料。K5finite worker-wall sum为137.17–141.68s/12观测，普通选择器66.87–70.71s；matched CIE112.99s，旧完整CIE155.44s。G/QR和实际种子另计；两次预算中断缺失partial辅助计数，账目显式记为NOT_RECORDED，不能当0或声称等总FLOPs。

这是已经曝光旧E1合同上的controlled development，并非新盲验证。第二个专门native-affine确认协议随后另行冻结、生成与完整闭合后才读取新private真值；其结果见下节，旧anchored确认没有被当作它的确认数据。

证据：[e1_native_tsom_scores.json](../results/e1_native_tsom_scores.json)、[共享finite增益案例图](../figures/e1_native_tsom_finite_gain_case.png)。

## 第二批独立确认：保留次级正面结果，不改主joint门

新native-affine方法冻结于2026-11-01 06:45:53.258952 UTC，随后06:45:53.577579开始生成另6个变化对象（3形状族）；2布局×clean/noisy，共24观测。9方法×3公共起点0/1/4＝648次预测于06:49:31.827973全部闭合，评分线程收到明确授权后才首次读取新private目录。全部648次零执行异常、零LU/RHS超限；封存预测与实际clean/noisy Y的残差记录差最大0。评价为`FRESH_BLIND_NATIVE_AFFINE_CONFIRMATION`。

下表每个条件有12观测，单元格为 **data / material+data / joint**；在线K3按原始物理Φ选择，不用GT选起点。

| 方法 | clean K3 /12 | noisy K3 /12 |
|---|---:|---:|
| 成熟CIE dual | 6 /6 /6 | 6 /6 /6 |
| 固定feedback | 5 /5 /4 | 5 /5 /5 |
| direct material seeds | 6 /6 /6 | 6 /6 /6 |
| native-affine普通Fourier | 6 /6 /6 | 6 /6 /6 |
| native-affine task local | **7 /7 /6** | 6 /6 /6 |
| native-affine finite task | 3 /3 /3 | 3 /3 /3 |
| native-affine finite Fourier | 6 /6 /6 | 6 /6 /6 |
| native-affine finite random | 6 /6 /6 | 6 /6 /6 |
| Native物理优化 | 4 /4 /4 | 4 /4 /4 |

主joint覆盖：task与普通Fourier成功集合完全相同；finiteFourier/random也相同，finiteTask在clean/noisy各损失3个既有成功案例。没有主joint增益。feedback clean 的在线K3 joint4，而封存终点oracle可达5，选择差距仍保留。

K1 joint clean/noisy为：CIE5/6，feedback0/0，material-seeds4/4，普通Fourier5/5，task3/4，finiteTask2/2，finiteFourier6/6，finiteRandom6/6，Native0/0。finiteFourier/random的K1增益付出了约两倍普通Fourier的worker-wall；至K3其覆盖没有增加。旧开发集所有finite规则共享的7-vs6增益没有成为这批新对象上的新增K3覆盖。

### 独立确认中的一个真实clean次级增益

`obs_7ad98676df61dbb6` 是curved、幅度16、uniform、clean。task在start0/1达到材料NRMSE约.000310、原始残差约8.9×10⁻⁵；直接material-seed、成熟CIE、固定feedback、普通Fourier和其他预声明finite规则均未恢复该案例。task的clean材料/数据成功集合因此相对强基线组合**新增1例、没有损失**。这是独立确认后可保留的正面观察，不能被“全N1负面”的说法抹除。

该对象GT自身Otsu对冻结support的IoU=.625，task的IoU亦为.625；所以near-perfect材料/物理恢复仍按原定门为joint FAIL。GTself诊断只解释门槛，不改变mask/阈值，不把该case重记为primary PASS。对应noisy观测`obs_b3b3447fe12c27c7`没有被task恢复，另一布局也没有同等增益；不能将其表述成跨对象/噪声稳定优势。

paid日志显示task选的是mixed候选池里的`random_native_3`，不是某个专属task家族向量。两起点都替换12实current DOF、总容量192→192；固定DPCS无自由变量，core15个APCS。task与Fourier都各交给**1个共同corrector种子**，均执行80 joint iterations。

| clean起点 | 方法/所选mode | 日志proposal_gain | 日志paid physical_gain | 全run LU/RHS | common-corrector LU/RHS | wall/s |
|---|---|---:|---:|---:|---:|---:|
| 0 | task /random_native_3 | .0104704 | .00906241 | 44 /906 | 39 /858 | 1.3400 |
| 0 | Fourier /fourier_native_old | 0 | 0 | 23 /444 | 18 /396 | 1.3656 |
| 1 | task /random_native_3 | .00647560 | **−.00703225** | 27 /532 | 22 /484 | 1.3857 |
| 1 | Fourier /fourier_native_old | 0 | 0 | 23 /444 | 18 /396 | 1.3805 |

两task最终材料都精确对应封存`common_corrector_0`，Fourier共同corrector残差仍约.13496。共享配额不等于相同实际工作；不能把选mode、80步joint rollout和后续corrector的整体效果归因于一条局部公式。尤其start1的proxy为正而paid physical_gain为负，仍最终恢复，已经阻止了“该分数单调预测basin/保证下降”的说法；日志中的`Born_no_error_bound`范围必须保留。

[新确认完整评分与逐对象/布局分层](../results/e1_native_confirmation_scores.json)、[分层CSV](../results/e1_native_confirmation_stratified_scores.csv)、[新clean材料增益与预算图](../figures/e1_native_confirmation_task_material_gain.png)。全部9-arm面板保留所有对象、布局、噪声与失败。6个seedless真值NPZ及24观测mapping在闭合评分后才发布于`data/native_confirmation_truth_after_freeze/`。

### 声明后的24次组件诊断：不是另一批确认

父线程随后只在上述clean及noisy两个观测、3起点上运行未改变的ordinary random/energy/normal/soft，共24次；标为`POST_DISCLOSURE_GAIN_CASE_COMPONENT_DIAGNOSTIC`，不加入原648个fresh率，不作统计泛化。

normal在clean start0达到与task指标一致的材料端点；soft在clean start1达到NRMSE .0003115、残差8.91×10⁻⁵；二者仍IoU=.625/joint FAIL。ordinary random/energy在这两个观测均未恢复，四种组件控制都未恢复noisy版本。因此该clean basin**不专属于task评分**，normal/soft也可触达。

normal start0也选`random_native_3`；soft start1也选该mode。本次ordinary-random在三起点恰都选了`fourier_native_old`，只能代表冻结的一个随机选择实现，不能当成覆盖全部随机候选的实验。组件诊断支持一个限定的current-mode/rollout机制观察，不能证明A17特异因果、新颖性或广泛优越性。

证据：[posthoc24评分](../results/e1_posthoc_mode_controls_scores.json)。至此全部物理数值运行停止；没有把N1/N2升级到576-cell优越性实验。

## 1. 评分对象、门槛与真值边界

材料误差为 `RMSE(m−m_GT) / declared amplitude`。开发集幅度为 8/18，确认集为 10/16，历史 E3 为 12；没有按每幅预测重新缩放。材料门为 NRMSE≤.15；图像门同时要求重建自身 256-bin Otsu 前景与冻结 GT mask 的 IoU≥.7；原始全 Maxwell 数据相对残差门为 clean≤.001、noisy≤.003。报告分别保留 data PASS、material+data PASS 和 joint PASS。

评分使用封存材料数组、预测数组和已支付物理求解的残差/候选日志；没有在评分阶段新增 Maxwell 求解。确认集 576 个封存预测与实际 clean/noisy 公开观测数组的残差逐一一致，最大记录差为 0。E3 直接使用已支付的原始全 Maxwell receipt，未另作 forward 验证。

旧 E1 存在真实访问偏差：E2 首次授权评分时，评分程序曾急切加载旧集全部 12 个对象真值，而实际只评分 E2 的 3 个对象；当时旧 E1 尚未闭合。未打印旧 E1 真值、support 或评分，且设计决定事先冻结，但数组读取本身已发生。旧 E1 只能称受控开发/机制集。详见 [TRUTH_ACCESS_DEVIATION.md](TRUTH_ACCESS_DEVIATION.md)，该偏差没有被后来的 lazy-loading 修复追溯抹除。

新确认集使用另一个 private 目录。方法冻结于 2026-11-01 06:00:04.534 UTC，新对象生成随后开始；全部 576 个预测于 06:02:51.626 UTC 闭合；本线程随后收到明确授权，才首次读取新目录并评分。该批独立保持 `FRESH_BLIND_SMALL_GRID_CONFIRMATION`，不继承旧 private 集的访问偏差。公开 truth release 发生在闭合评分后。

开发集是 3 个固定连续几何 ×2 幅度×2 物理网格形成的 **12 个材料合同**，并非 12 次独立形态抽样；布局、初始化、噪声版本都是重复条件。主要描述性 bootstrap 按 3 个形状族聚类，置信度很低，不作人群统计显著性推断。新确认集有 6 个变化后的对象，但仍仅 3 个形状族。

## 2. 原始粗网格 IoU 的实际限制

离线 GT 自身 Otsu 与冻结 support 的比较没有改变门槛、mask 或物理调用数：旧 4×4 curved 两幅度、6×6 islands 两幅度共 4/12 个材料合同，即使材料精确等于 GT，其 Otsu IoU 也只有 2/3。原始 curved-4×4 的近乎完美材料恢复因此可能 joint FAIL。

初始严格 CIE single/dual 均为 13/24 material+data PASS、10/24 joint PASS。其中两个 curved-4×4/幅度8案例 NRMSE≈10⁻⁵，IoU=2/3；另一个 islands-4×4/幅度18/jittered 案例为真实重建边界偏移。原门槛继续保留，但不能把这些都解释成物理/材料恢复失败。

新确认集 GT 自身 Otsu IoU 全部≥.8889，故此粗网格门槛问题不能解释新确认集中的 N1 负面结果。诊断记录位于 [gt_gate_diagnostic_scores.json](../results/gt_gate_diagnostic_scores.json)。

## 3. 成熟 CIE 基线已得到实质性恢复

[CIE_TSOM_BASELINE_AUDIT.md](../CIE_TSOM_BASELINE_AUDIT.md) 与 [CIE_SOURCE_AUDIT.md](CIE_SOURCE_AUDIT.md) 记录了 DPCS/APCS、源归一化目标、原生 Fourier 系数 PR-CG、共享实材料 readout、mask/beta continuation、原始物理 readback 和共同 corrector。

该实现是真实 **2D scalar-TM 改编**，并非 2019 年 vector 3D 论文数值实验的逐字复现；小网格使用显式 Fourier 列，也不能宣称复现原论文大网格 CG-FFT 的速度/内存标度。2011/2016 原全文本轮未独立重新取得，宽泛先验技术覆盖与新颖性结论不由这次窄基线审计支持。

DEV 中原论文类型的系数停止规则只运行 22/23 步，原始数据残差 .32659/.078185；相同 beta/mask 配置改用严格停止后运行 305 步，约 1.19/1.16 秒，残差降至 .000427/.000266。由 22/23 步端点称 CIE 为弱基线是不成立的。严格版本在后续比较中保留。

## 4. 开发集多初值、预算和状态权重控制

初始 504 行包含 216 基线与 288 初始 N1/材料种子结果。初始 N1 只有很小 current backbone，公开数据的最小可表示尾误差约 .554–.901；它被保留为 `IMPLEMENTATION_CAPACITY_FAILURE`，而不是对 normal 理论的反例。公开 receiver rank 修复在读取初始 E1 评分前声明，后续修复集有 264 行；其中两个 TRF 初值边界 roundoff 异常保留，没有重跑擦除。

全多初值协议包含 14 方法×24 观测×5 个共享起点=1680 行，其中 336 行是准确复用的封存 start0，实际新增 1344 次执行。K 顺序为 0→1→2→3→4，在线组合按共同原始物理目标 Φ 选最小值；oracle 只在离线提供覆盖上界。没有 K=8/16 外推。

| 开发集 K5，每方法/24 观测 | data / material+data PASS | joint PASS |
|---|---:|---:|
| 严格 CIE single / dual | 13 / 13 | 10 |
| 固定 scalar feedback | 12 / 12 | 10 |
| GlobalDual | 8 / 8 | 8 |
| 普通材料梯度种子 | 6 / 6 | 6 |
| CSI-state N1 finite task | 7 / 7 | 5 |
| CSI-state N1 finite Fourier | 7 / 7 | 5 |
| CSI-state N1 finite Krylov | 9 / 9 | 7 |
| CSI-state N1 task / soft | 4 / 4；5 / 5 | 4；4 |

有限 task/Fourier/Krylov 新增 data PASS 均未超过 CIE 的案例覆盖。局部 task/soft 相对单独 CIE 多出的一个案例，也已被固定 feedback 找到；相对强基线组合没有新增案例。组合比较只用于检验增量覆盖，不能冒充与单一方法等总成本的竞争。

840 行固定总预算控制采用每起点 .8s、120 LU/2400 complex RHS 的共同上限，K5 名义总预算为每观测 4s。全部 840 行实际计数没有超 LU/RHS 上限；67 次物理上限停止保留（40 LU、27 RHS）。K5 CIE dual 为 12 material+data /9 joint；固定 feedback 12/9；材料种子6/6；N1 task、finite task、finite Krylov 与 Native 均2/2。共同上限不等于相同总 FLOPs，G/QR 工作另计。

600 行 CIE-state 控制使用冻结 beta=.2、相同 observable pool/rho/prior，并包含正确 DP·e 导数。K5：Fourier与 task 均10 material+data、7 joint；soft 为11/8；finite task6/4；finite Krylov8/6。CIE state 权重改善了部分局部 CSI 对照，但没有显示 task 选择优于 Fourier，或新增成熟 CIE 案例。Krylov 的逐运行 data PASS 为21/120，高于 CSI 的17/120，但 K5 案例覆盖反为8 vs9；两种量不能互换。

修复后的 N1 current 自由度为 4×4 网格156→168、6×6网格204→216 个实自由度。同组 N1 selectors 可作等容量比较；成熟 CIE 的固定 DPCS 与阶段性 free APCS 不同，不能称其与 N1 严格等 current DOF。

固定 feedback 的离线候选存在选择差距：全多初值逐运行 joint 为17，候选 oracle 为21；K5在线10、全部候选 oracle12。不能把 oracle 成功当作在线成功。

证据：[e1_multistart_scores.json](../results/e1_multistart_scores.json)、[e1_fixed_total_scores.json](../results/e1_fixed_total_scores.json)、[e1_cie_weight_scores.json](../results/e1_cie_weight_scores.json)。

## 5. 独立新确认集给出的主判断

6 个变化对象、2 布局、clean/noisy 两条件，共24观测；8方法×3个共享起点(0/1/4)=576 次运行。每起点 wall 上限4s、LU120、RHS2400，K3名义上限为每观测12s。全部运行零执行异常、零已记录 LU/RHS 超限。

| K3 在线 joint PASS | clean /12 | noisy /12 |
|---|---:|---:|
| 严格 CIE dual | 4 | 3 |
| 固定 feedback | 6 | 6 |
| Native | 2 | 2 |
| CIE-state N1 Fourier | 4 | 4 |
| CIE-state N1 task | 4 | 4 |
| CIE-state N1 soft | 4 | 4 |
| CIE-state finite Krylov | 2 | 2 |
| CSI-state finite Krylov | 1 | 1 |

该表的 data、material+data 与 joint 计数相同。所有新 N1 成功均包含在 CIE+固定feedback+Native 的联合覆盖内（clean7/noisy6）；N1 对单独 CIE 的一个 clean/一个 noisy 增补案例，也由固定 feedback 覆盖。task/soft 未优于 Fourier。该有限确认集不支持已冻结 anchored-N1 原型的新 basin 或方向选择优势。

N1 并非全部未执行：local CIE Fourier/task/soft 各72次 `EPISODE_LIMIT`；实际 finite seed 插入为 CSI-Krylov120次、CIE-Krylov88次。配额与实际插入数分别记录。clean 的幅度16合同中，只有成熟 CIE 找到1/6个布局观测 joint PASS，N1与固定feedback均0；这保留了成熟 CIE 在困难案例上的真实价值，但幅度与几何变化耦合，不能当作纯幅度因果实验。

证据：[e1_confirmation_scores.json](../results/e1_confirmation_scores.json)。已确认 anchored-N1 的定义与记录保持冻结；后述新 native-affine 实现必须另作协议和证据，不能作为对旧方法的追溯调参。

## 6. E2：有条件根不等于全逆问题根

240 次真实运行覆盖6方法×4匿名观测×2个k×5起点；这些观测只有3个材料合同、2个形状族。120 个 charts 得到33个合法 subset roots，selected-channel 残差≤4.27×10⁻¹³，但它们全数据残差为 .822–1.251，无一个原始或校正后的合法 root 达到全数据门。12个同 chart 根对的材料 RMS 差为 .875–3.143，证明它们并非全部重复，却不证明全材料问题的根枚举。

| /40 次实际运行 | data / material+data PASS | joint PASS |
|---|---:|---:|
| Algebraic | 2 / 2 | 2 |
| block TRF | 4 / 4 | 4 |
| block multistart | 4 / 4 | 2 |
| random | 0 / 0 | 0 |
| full multistart | 4 / 4 | 4 |
| matched-count block multistart | 4 / 4 | 2 |

Algebraic 两次通过均来自空 root 集时的 anchor fallback，不能归因于 root proposal。其实际 corrected seed 分布1/2/3为31/6/3；matched control 精确匹配该分布，LU为1900而 algebraic2015，且在线目标15胜/3平/22负。所有已封存候选 oracle joint 与在线计数一致。数值 quotient/path resolution 不构成严格 all-roots certificate，更不是连续 Maxwell 全问题根完备性证明。

证据：[e2_scores.json](../results/e2_scores.json)、[e2_mechanism_summary.json](../results/e2_mechanism_summary.json)。

## 7. E3：成熟基线的真实 576-cell 正面结果

这是已公开历史目标的 `DISCLOSED_HISTORICAL_MATURE_BASELINE_DIAGNOSTIC`，不是新盲测。12个目标为数字1–9、ring、two_islands、curved_bar；每个最终材料向量576维，四个成熟方法共48次独立执行。保持2D scalar-TM物理、45s队列配置与既有平台。N1/N2没有E3预测，不生成伪空白图。

| 方法 /12 | data PASS | material+data PASS | joint PASS | median NRMSE | median IoU |
|---|---:|---:|---:|---:|---:|
| Native TRF | 0 | 0 | 0 | .39924 | .01209 |
| 固定 feedback β=.2 | 8 | 8 | 8 | .10823 | .75556 |
| CIE single + corrector | 9 | 9 | 8 | .10416 | .76924 |
| CIE dual + corrector | 10 | 10 | 10 | .09413 | .81530 |

single 的 digit5 为 material+data PASS，但 IoU=.67826，继续 joint FAIL。固定 feedback 失败于 digit1/4/6、curved_bar；single 另失败于 digit7、ring、curved_bar；dual 失败于 ring、curved_bar。Native全部失败。所有失败图与逐行记录保留。

digit2 pilot 中，CIE raw-current 端点已经有良好图像（single NRMSE .09567/IoU .80952；dual .09759/.80952），但原始数据残差 .01650/.02530 仍 FAIL。共同全材料 corrector 后降至 .000311/.000296 并 joint PASS。不能把校正后的结果全部归给未经校正的 auxiliary state。

原始receipt中，CIE single/dual共同corrector各12/12以`WALL_SECONDS`停止；Native为7/12 wall停止；fixed feedback有7/12 `LU_CAP_WITH_PAID_CLOSURE`及2/12 wall停止。因此不能把48个worker统称纯wall截断或忽略实际cap停止。CIE current-stage与corrector stopping分别保存，不改写原始记录。

CIE single 仅3/12完整完成预定阶段，9次截断；dual仅1/12完成、11次截断，且只有9/12到达mf10。single12/12到达mf10。因此这些45s结果是预算内性能，不是成熟 CIE 的收敛性能上限。dual在本协议较强，但没有全局恢复保证。

### 576-cell 成本与图像约定

| 方法 | actual worker-wall sum /s | full LU | actual complex RHS |
|---|---:|---:|---:|
| Native | 296.8810 | 1192 | 71520 |
| 固定 β=.2 | 421.3715 | 20391 | 489024 |
| CIE single | 540.1542 | 15926 | 381660 |
| CIE dual | 540.2338 | 11825 | 283272 |

worker-wall 合计 **1798.6405s**；两批2-lane队列 elapsed 合计 **930.279s**，已含setup。二者不能混称。共享 setup 合计约3.47665s、6LU/72source RHS、3次operator upload；580 G/S basis列及SVD按各实际队列记录，共享一次而非逐任务重记。比较不改写成跨 CPU/GPU 或跨网格的固定加速倍数。

原始历史 raster 来自 PIL row order，因此 E3 全部 GT、重建、绝对误差使用 `origin='upper'`，标注 row0 at top。先前 lower-origin 图已修正；仅显示映射改变，材料与评分没有改变。材料色阶统一[.02,20]，绝对误差统一[0,20]。三张主图均已视觉检查：

- [数字1–9](../figures/fullwave_digits_1to9.png)
- [全部12个目标与失败案例](../figures/fullwave_all_shapes.png)
- [数字绝对材料误差](../figures/fullwave_digit_errors.png)

原始记录与评分：[e3_all_scores.json](../results/e3_all_scores.json)、[e3_all_scores.csv](../results/e3_all_scores.csv)。这些成果属于成熟 CIE/固定 feedback 基线与实际工程适配；它们没有把 N1/N2 的 E3 门改成 PASS。

## 8. 数学成立范围与算法结论必须分开

[MATHEMATICAL_AUDIT.md](../MATHEMATICAL_AUDIT.md) 中经条件推导的 pure-normal ±readout parity，具有同号 O(t²) leading readout、两侧材料差 O(t³) 的局部范围。这不是对混合 N1 action、非线性 rollout 与所有有限种子的因果证明。不能用该公式直接宣称“已经证明本次 N1 失败的原因”，也不能从局部公式成立跳到实际 basin 优势。

当前证据允许保留经审计的局部导数、normal/Schur 机制与条件代数根作为数学/机制成果；算法证据是混合的：主指标大多没有新增覆盖，但native-affine存在一个clean材料/数据次级增益；其归因与稳健性受上述对照限制。它们不构成理论整体无效的证明。相反，成熟 CIE 的正面结果属于一个已知方法家族的真实成像表现，不能被重新归属为本次新方法新颖性。

## 9. 成本、失败与最终放行建议

[COST_LEDGER.csv](../results/COST_LEDGER.csv) 当前有5748个唯一付费 inference execution；按 `(execution_source, run_id)` 去重，336个复用start0与E3 pilot没有重复计费。初始失效N1与修复N1的同名运行来自不同执行目录，均真实保留。另列20个DEV/audit条目、2个owner-native DEV/weighted-reference条目、2个公共生成条目、2次E3共享setup与1个未测量启动失败汇总，共5775行。owner native DEV为.815661584s/2LU/30RHS，weighted-reference为1LU/6RHS且wall未知；新native确认生成另记12LU/72RHS及已记录wall .008774708025157452s一次（audit/NATIVE_CONFIRMATION_GENERATION.json）。

已记录账目下限为391744 LU、7616788 complex RHS；这是跨不同矩阵维度的记账和，不是统一算力单位或精确完整总成本。原始 G/QR/SVD、small solves、seed插入/拒绝及各phase成本以不同字段保留，没有把它们当作免费。全源初始化/导入与被覆盖的早期DEV记录存在明示缺口。

必须保留的未知项：第一次 representability 诊断已做完24个forward后，因np.int64 JSON序列化失败；可恢复收费24LU/144RHS，physics wall未测量，整次tool .243s含导入，不能代替physics wall。远程失败launch的物理/时间成本为 `UNMEASURED`，不是0。owner weighted-reference的wall也未知；native确认生成wall已有receipt，已正确计入.008774708025157452s；预算中断丢失的partial G/QR等辅助计数记为NOT_RECORDED。后续成功诊断与成功GPU任务独立收费，没有为修复记录而重跑公式检查。

**最终独立放行建议：发布完整有限证据、成熟CIE真实576-cell正面结果，以及native-affine一个clean材料/数据次级增益；保持主joint无新增覆盖和post-disclosure normal/soft也能恢复的限定。该结果不足以启动N1的576-cell优越性扩展，N1/N2 E3仍为NOT_RUN。下一项科学工作应是跨未见对象与噪声稳定确认材料步增益，而不是从一个对象推导A17特异或通用算法优势。** 旧分支的新确认集支持其冻结协议内这一判断，但只有6个对象/3个形状族，不能升级为所有形状、连续PDE或所有先验上的否定。不能声称宽泛TAP级新算法优势、文献方法收敛上限或全逆问题all-roots。父线程的最终论文与公开解释必须保留上述归属、门槛、访问偏差和成本缺口。
