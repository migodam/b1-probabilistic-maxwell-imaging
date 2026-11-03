# R4 最终静态科学审查

日期：2026-11-03。范围：当前 `RESEARCH_REPORT_ZH.md`、`ABC_VALIDATION_ZH.md`、`THEORY_STATUS.md`、推导第 8 节、既有 confirmation 评分、post-confirmation addenda、两份 failed-endpoint continuation 回执，并在本审查过程中读取已生成的 `DEVELOPMENT_CONTINUATION_SCORES_V1.md`。另核对成本说明与 confirmation 汇总字段。仅静态核对／对已评分行计数；未运行求解、评分、梯度或物理检查，未打开 GT 数组／私有映射，未改算法、冻结合同或原报告。

## 结论及发布前需同步的事项

**未发现要求新增科学实验的 blocking 问题。** 主结论与当前已评分证据一致：A 只有特定旧开发对象的恢复正例，fresh 六对象未超越成熟 CIE；B/C 已测路径没有新恢复；32 网格两臂均恢复、平移平均两开发端点失败、真实相位控制仅支持该实例范围。报告没有将低物理残差、达到数值一阶容差或重复旧正例扩写成普遍成像优越、局部极小、唯一性或 3D 证据。

发布前仍需做两项**事实／文档同步**，不需新增门或实验：

1. **最终成本账尚有过时 pending。** 当前 `COST_LEDGER.md:1,13,17,23` 仍保留 confirmation／Windows 待集成措辞，且候选总账尚待纳入完整 Windows 正例及两份续跑等最终回执。主报告第 143 行将其称为“最终总账”，须在合账完成时同步；若仍有未知计费或 bootstrap 时间，应保留 null，并把状态写成“已记录成本下界，个别字段未知”，不能把已完成的科学实验称仍 pending。父线程已在进行此项合账。

2. **刚落盘的续跑材料负证据应进入主报告。** 第 114 行已经正确报告 PG／FTOL 与残差，但尚未提到 native 在 Phi 继续下降时材料误差进一步恶化。追加简短事实及 `results/DEVELOPMENT_CONTINUATION_SCORES_V1.md` 链接即可：native NRMSE `0.743118035 -> 0.750890266`，onehop `0.761367699 -> 0.760862298`，二者仍 FAIL。不得以新端点替换原 18 个冻结端点或更新原 4/6。

以上属于最终公开材料完整性，不是对核心结论的否定。此审查不代替父线程最终科学接受或发布访问核验。

## 主统计、候选池与负结果

对 `CONFIRMATION_SCORES.json` 的 18 行直接计数：

- onehop/global/native_extended 各 6 个对象，各 4 个 joint PASS；三者成功的四对象完全相同。两弯条均失败。
- 18 行全部 `complete_auxiliary_schedule=true`，全部原 corrector 为 1,000 步上限停止；未见 UNKNOWN／FAILED 被排出分母的情形。
- onehop 相对 global 的材料 NRMSE 是 1 低／5 高；IoU 是 3 高／2 低／1 相等，与主报告相符。六个对象是比较的描述单位，不能把 18 个方法端点当 18 个独立样本。本报告只作逐例描述，没有显著性或总体优越性声称。
- 原 Phi secondary pool 仍为 4/6，42cb 选 onehop，50fc 选 native；两者都恰是该对象三臂中材料误差最高者。主报告保留了这一不利证据，没有按 GT 改选。
- 42cb 的 onehop Phi `0.00219301446` 低于 global `0.00410523163`，NRMSE 却由 global `0.741539073` 升到 `0.761367699`，报告的残差／材料区分准确。

`DEVELOPMENT_POST_CONFIRMATION_ADDENDA_V1.json` 与报告相符：平移平均在旧 R2 正例及 42cb 均 FAIL；42cb NRMSE 降到 `0.641395037`，但相对残差升到 `0.143490953`，不能称成像或数据拟合成功。Windows 重放旧正例 PASS，NRMSE `0.041603262`、IoU `0.951807229`、相对残差 `0.000968350517`；其约 `376.20s` solve 不是性能优势证据。这三行被明确标成已暴露 development／跨平台重放，没有加入 fresh 分母。

主报告已同时保留 24-grid 正例、32-grid 比较优势消失、fresh 无优势、pi 失败、平移平均失败及 B/C 失败。未发现只展示正例而隐去这些已知反例的主结论问题。本次没有重新读取 phase／32 原始数组及完整运行回执，对它们仅核对当前报告间一致性，不能称再次独立数值验证。

## 续跑的科学解释与成本

两份 `failed_endpoint_continuation/*/RECEIPT.json` 与主报告及续跑评分一致：

| 方法 | 实际额外步数 / 额度 | 停止 | 原尺度投影梯度 infinity norm | 相对物理残差 |
|---|---:|---|---:|---:|
| native_extended | 3361 / 4000 | PGTOL | 7.45734496e-9 | .0886718336 |
| onehop | 251 / 4000 | FTOL 相对下降 | 3.86901419e-6 | .0655837432 |

native 确实达到声明的 `1e-8` 数值 bound-projected 一阶容差，但远未到噪声球；所以“该基线失败不能只归因原 1,000 步不足”是有范围的合理结论。它仍不认证局部最小、全局最优或数据不可辨识。onehop 的 FTOL 退出不满足同一 PG 精度；保留“数值停滞而非认证收敛”合适。两者有相同最大额外额度，但实际工作不同，不能将该对照称为精确等算力。

续跑顶层物理成本分别为 native `3385 LU / 81264 complex RHS`，onehop `267 LU / 6432 complex RHS`，共 **3652 LU / 87696 complex RHS**。只能加顶层 `physical_cost` 一次；其 nested corrector 计数、已有端点和此前六次 PG 诊断不再重复加入。该表是续跑额外成本，不能混入原六对象 campaign 的 18 个端点成本来悄悄改变比较口径。续跑独立评分本身记录 0 新物理工作。

## 理论与机制表述

第 8 节的四 parity 分块平均公式正确：水平／垂直邻边在四分块中出现概率 1/2，角邻边 1/4；边界采用非周期 partial/singleton 分块。线性 residual transport 的平均等于给出的 Gbar，但不等于四个平方目标的平均。凸范数界保证共同 q<1 时状态零集不变；该充分条件本身不保证读出更好、成像成功或优化路径稳定。

两个已暴露端点失败，仅否定所测的“固定无调参 parity 平均足以修复迁移”方案，不否定所有空间 transport。主报告第 112 行把否定限制到这一具体解释，适当；不得再扩大为证明所有 partition-origin 影响均不可修复。

第 7 节的满列秩／SPD／零缺陷局部 LS 假设及实际一步需 `C^T M(I-Pi)b != 0` 已在主报告保留，不再将“存在某个可改变读出的 b”偷换成“实际 b 必然改善”。第 8 节也未依靠该局部式给实际 prior／box／非零缺陷／有限内层／双 beta 迭代作收敛保证。

`THEORY_STATUS.md` 的 proved／derived／empirical／NOT_RUN 主要分界与证据一致。新颖性陈述承认 preconditioned CSI、NIE/CIE、T 重整化及加权 LS 先例，没有用有限检索未找到同式实现证明 first。二维 scalar-TM 明确没有外推为未跑的 3D vector 反演。

## 成本口径的已核对部分

confirmation 三臂任务顶层小计为 `18647 LU / 445860 complex RHS`；共享生成一次 `6 LU / 72 RHS`，合计 **18653 LU / 445932 RHS**。JSON 的 `all18_plus_shared_once` 内同时保留 task subtotal 和明确命名的 `campaign_physical_totals_including_shared_once`，报告使用后者，数值正确；不能再把 shared setup 另加一次。

报告列 onehop/global/native 的 worker wall 累计 `1147.09 / 577.95 / 2174.97s`，和确认评分相符。`3900.01s` 是 worker wall 总和，`1179.03s` 是既有远端 orchestrator elapsed，二者未混同；约 `3804.34s` CPU 的 import／bootstrap 范围有说明。小块操作、576／1024／16 维物理工作分列，不按相同 forward 单价混加。金额和缺失连续项目计时维持未知，不影响上述已记录科学结论，但最终合账前不得把候选 ledger 当完备项目账单。

本轮审查建议止于上述两项文档／账目同步，不追加任何实验、参数搜索或成功门。

## 父线程发布前处理记录

2026-11-03：上述两项同步已完成。`COST_LEDGER.md` 现在以 `COST_LEDGER_FULL_V2` 为最终已记录成本清单，包含 Windows 重放、结构修正、驻点诊断及两份续跑；未知计费和未记录启动时间仍为 null。主报告第 5 节已写入 native NRMSE `.743118 -> .750890`、onehop `.760862` 和独立续跑评分链接。原确认 18 个端点与 4/6 分母未变。本节只记录文档处理，不改写上面的独立审查，也没有新增科学运行。
