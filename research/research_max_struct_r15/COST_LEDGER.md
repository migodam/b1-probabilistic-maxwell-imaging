# 费用与真实计算量

本页区分OS整进程CPU、内部body计时、并行wall、物理LU/RHS和辅助算子。原记录及逐任务表在[最终费用汇总](analysis/cost_final_enriched_v2/COMPACT.md)。没有把多个线程wall之和叫并行历时，也没有构造不同计时口径混合的“精确全项目总CPU”。

| 组 | 实际CPU口径/秒 | Maxwell LU | aggregate complex RHS |
|---|---:|---:|---:|
| 主264反演 | OS child 31,372.375 | 203,157 | 3,717,630 |
| 同状态fork56 | OS child 102.375 | 172 | 2,064 |
| 接线、旧端点及所选终点验证 | 各独立进程另列 | 149 | 7,500 |

149 LU包括integration13/156 RHS、旧48端点72/6,000、16材料读出20/288、44终点44/1,056。这里不把source_RHS/adjoint_RHS再加到aggregate。逐次Maxwell因子维度为576；原counter维度累计和不是一个巨型矩阵的维度。

主辅助小分母LU共2,985,359次，三角RHS84,966,012；fork新计22,441次/1,032,624 RHS。它们是6/11/19等小系统，不能折成同量Maxwell LU。background H Cholesky主组96次、fork新24次，其setup/iteration三角RHS、FFT动作和source DOF均在独立命名空间。fork只计cumulative−inherited，新旧precomputation不重复收费。

| 方法，均24对象 | OS CPU秒 | Maxwell LU | complex RHS |
|---|---:|---:|---:|
| G19 | 2,836.656 | 16,948 | 312,804 |
| S19 | 2,889.500 | 18,297 | 329,310 |
| U11 | 2,896.250 | 18,551 | 330,732 |
| U11＋邻域 | 2,887.766 | 17,983 | 326,940 |
| 路由DR融合 | 2,896.516 | 17,599 | 322,098 |
| 路由AP | 2,866.094 | 16,937 | 307,110 |
| 全空间周期DR | 2,859.609 | 17,021 | 308,082 |
| PR reset | 2,890.531 | 18,041 | 327,966 |
| CIE | 2,837.250 | 19,409 | 355,680 |
| Dual-CIE | 2,841.234 | 19,827 | 365,406 |
| 独立DR流程 | 2,670.969 | 22,544 | 431,502 |

主队列实际并行历时3,046.360秒（约50.77分钟），controller CPU11.328秒；31,372.375秒是所有child进程CPU之和（约8.71 CPU小时），不能将其称为wall。

120秒是候选资格的整进程上限，不代表每次都恰好使用120秒。结束序列化有付费尾部，晚完成候选不进入预算内；原状态保留和新缓存也收费。等CPU、共同LU/RHS候选池和固定迭代干预是三个不同实验。

独立DR首个完成corrector审计CPU中位43.25秒；CIE68.13、Dual67.13、stateful臂约73秒。这个数包含第一次审计计算，不是exact corrector start。18个DR对象完成固定projection日程，6个碰proposal CPU；八stateful臂均碰72秒proposal期限。不能拿raw最优候选出现时间充当阶段结束。

## validation、评分与IO

新增导数接线和完整state接线只查新接线，历史Atlas恒等式不重跑。新GT评分14.56 CPU秒包含绘图；失败的首次评分import接线4.55秒单独保留。无GT完整状态分析远端88.44 CPU秒，全部480cycle pairs与36有限响应；新最终图两次绘制合计16.92 CPU秒。各内部子计时不再次叠加。

原始结果8,723,094,459 bytes分16卷无损ZIP，7,364,794,746 bytes。第一次归档因缺远端complex缓存失败，OS CPU531.625秒、wall567.35秒、0卷写出；修正后的V2内部body139.14秒、wall217.66秒，snapshot139.15625秒为嵌套口径，不能相加。审阅导出CPU162.95秒/wall285.82秒另列；缺58个secondary导致116个缺文件槽位，是候选不合格的真实缺失，不是116次运行失败。

最初detached SSH子进程在会话关闭后未存活，未启动任何反演，CPU未知；之后保持前台SSH完成所有任务。一条2.389GB审阅包拷贝在836MB处改为精简副本并保留partial；精简734.8MB副本仍保留全部候选材料池、全部m/Y视图与评分元数据。完整current/field/PR另在全部raw分卷。停止拷贝的CPU无法精确追溯，填null，不当免费。所有大文件仅在外置盘和任务远端。

## 仍未知的费用

货币账单、LLM实际计费未由工具提供，均null。该批算法没有CUDA或NN训练；不把CPU工作报告为4060计算。内存为原运行sample/end RSS，未伪称精确峰值。部分部署/SSH/SCP与发布阶段CPU不能从现有工具可靠拆分，保留[缺口表](analysis/cost_final_enriched_v2/GAPS.csv)。上传、构建及访问回读有另外的publication receipts，不并入120秒算法比较。

## 已完成的发布IO费用

全量raw及metadata28附件已上传到同一草稿release，7,377,892,085 bytes，0失败/重试/覆盖；远端上传body CPU30.953125秒，wall2,424.647秒，进程snapshot31.03125秒为嵌套口径。关键state选集另188,513,757 bytes，remote包装CPU1.625秒、wall1.740秒；外置盘SCP进程树CPU2.599秒/wall188.921秒。原文件未删除。公开发布、快捷Pro包和网页另有实际回执；这些是IO，不归入120秒算法表现。

快捷Pro包构建CPU2.429645秒/wall4.470720秒，成品430,038,635 bytes；6个本地Pro/关键state及metadata附件上传CPU3.399604秒/wall52.662秒。包内已声明快照时间，构建后父线程验收及发布回执以Git最新版为准，不重复构建或重复计数。
