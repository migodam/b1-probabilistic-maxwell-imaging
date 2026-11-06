# 已完成历史/DEV成本账本审查

主账本为`COMPLETED_COSTS.json`，按源收据逐条保存CPU、wall、状态、LU/RHS类别、矩阵维数和nested视图。这里只整理已经完成或明确中止的历史/DEV工作；**没有运行fresh scorer，没有打开fresh sealed_truth，没有补跑物理**。fresh二维168任务和fresh三维的生成、运行、传输成本全部留待最终终态账本，未使用运行中估计混入本表。

## 已记录的顶层计时窗口

| 阶段 | 范围 | 已记录CPU秒之和 | 已记录各进程wall秒之和 |
|---|---|---:|---:|
| 历史energy分析 | 69历史runs/449保存状态，含2失败及1中止尝试 | 170.648895 | 175.431338 |
| 固定状态更新 | 32完成读出＋失败尝试已付16读出 | 0.967325 | 1.005183 |
| 完整DEV分母候选 | 6原候选＋2保留自项，共8条 | 301.672837 | 310.076151 |
| onehop optimizer-only正式控制 | 2条修正dispatch的LBFGS8完整轨迹 | 74.342174 | 75.993405 |
| onehop错误dispatch尝试 | 2个自身中止进程，实际inner60，不是有效inner8控制 | 119.202706 | 125.547438 |
| 历史离线评分与图 | 原六候选评分、材料/energy图、self-retention及onehop既有评分 | 4.339666 | 4.924779＋未记录部分 |
| 24状态shifted谱诊断及图 | 诊断CPU已记录；图CPU/parse失败CPU未记录 | 0.271964＋未记录部分 | 0.711922＋未记录部分 |
| self-retention导数wiring | 已有有限差分/伴随接线收据 | 0.015851 | 未记录 |
| 文献检索 | 7份检索artifact，含ScholarQA 429失败 | 未记录 | 未记录 |

上述**已记录CPU和为671.461418秒**，已记录各进程wall和为**693.690217秒**。这两个数字都不是完整实际支出的声称：未记录值保持null，文献工具服务端/客户端CPU没有凭空补零。wall为若干调用窗口之和，可能有并发重叠，不是日历总历时。账本重编译本身的极小CPU/wall另存`BUILD_COSTS.jsonl`，未重复并入科学阶段数值。

DEV inverse收据是既有runner的记录窗口；没有外部完整child GetProcessTimes测量，不能将其升级为覆盖全部启动/导入/退出开销的进程寿命。8候选与R9历史参考有同Mac主机记录，可作历史描述性cost对照，但并非同期同负载重跑；fresh同队列成本比较仍需后续独立终态收据。

## LU/RHS必须按算子及维数分开

| 计数类别 | 维数 | 已知LU | 已知复RHS |
|---|---:|---:|---:|
| 原物理forward/adjoint相关 | 576 | 15914 | 380532 |
| 历史小3D离线物理forward | 1536 | 38 | 2280 |
| shifted辅助分母 A_G | 576 | 24 | 288 |
| 候选辅助小分母 | 28 | 135475 | 3251400 |
| 候选辅助小分母 | 19 | 136390 | 3273072 |

原物理复RHS计数包含各原模型收据记录的forward/adjoint列；小分母表是小三角系统RHS。576维shifted辅助LU虽然维数相同，**不是原物理forward LU**。19/28/576/1536维次数不可相加后宣传为同一成本单位或等量计算；真正wall/CPU耗时仍看其原计时窗口。

历史energy中止的vector进程只留下ps观察到的23.24 CPU秒/24 elapsed秒，其可能已产生的额外forward计数无法恢复。因此物理LU/RHS总数也是已知下界。self-retention wiring没有记录物理LU/RHS字段，账本保留null；不因任务名称就伪造实际测量。

## 防止重复收费

每个独立调用只将一个最外层**已记录**CPU/wall窗口计入和。各候选的setup、auxiliary、auxiliary execution、corrector、physical operation timing互相嵌套；全部保留为nested证据，不再加到候选总CPU上。固定状态basis setup和材料readout微秒计数同理。原energy聚合总账与其分收据没有各收费一次；重用R6/R9的旧inverse结果也不作为本轮新inverse支出收费。

失败并未消失：固定状态失败尝试的20物理LU/288 RHS已计；错误inner60的两个onehop进程合计952 LU/22740 RHS、119.202706 CPU秒已计；energy的dtype/shape失败及中止尝试已计；谱诊断的parse错误保留日志，未执行任何LU但没有测得CPU，故CPU为null。ScholarQA验证的429失败保留检索artifact，检索CPU、wall和货币费均未知。

## onehop控制的独立源码/记录核对

`ONEHOP_CONTROL_SOURCE_REVIEW.json`保存了核对结果。两条corrected dispatch的实际`NUMERICAL_SCHEDULE.material_inner_maxiter`均为8；各3000条有inner迭代信息的history记录中，观察到的最大`material_inner_nit`也均为8。

实际callback明确注入原`research_max_struct_r4.src.cie_onehop.onehop_material_readout`，该函数使用原onehop目标/导数、普通CIE种子和L-BFGS-B，没有调用Möbius/MM。schedule中剩余的`aplus_arc_candidates`等字段不被这个callback读取。两条正式控制的629/1769物理LU和32.019755/42.322419 CPU秒已单列，原未修正的两个inner60尝试不能混入有效控制结果。

root已有离线评分显示正式onehop+LBFGS8控制未复现G19对这两个历史DEV实例的材料结果。这只限制“仅换inner优化器已经足够解释改善”的假说；不是对成熟A路线的全面反证，也不是fresh确认结果。本账本没有再次评分或据此修改任何候选/预算。

本次到此为止，继续等待root在168任务全部终态后明确postseal评分授权。
