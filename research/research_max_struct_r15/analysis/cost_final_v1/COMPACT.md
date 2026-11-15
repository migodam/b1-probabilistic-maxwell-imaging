# R15最终机械费用汇总

仅读取本地JSON费用/队列元数据，不打开材料/GT数组，不运行physics、优化、评分或远端任务。费用不代表恢复质量或科学结论。

主264任务外部OS整进程CPU重算：31372.375000s；fork56：102.375000s。两组控制器记录CPU分别11.328125s、3.718750s，与子进程分离。保留11方法×24完整分母、全部终态/失败/未知项。

Main实际Maxwell LU / aggregate complex RHS：203157 / 3717630；fork：172 / 2064。source_RHS、adjoint_RHS及operation角色是aggregate RHS子项，不能再相加。

已记录独立validation物理总费：149 LU / 7500 aggregate complex RHS。不叠加其case/arm子计数；分区/常规模型费用分别见原回执。

full_factorization_dimension是原记录的累计维度之和；逐run除以LU所得实际单因子维度另列为576，不将sum当单矩阵维度。denominator_small_LU与其dimension_sum/triangular RHS、FFT调用/变换、background H Cholesky/setup/iteration RHS各保留独立counter命名空间，不折算为Maxwell LU/RHS。峰值bytes/rank仅取max/snapshot，不相加。

Fork仅按cumulative−inherited计新auxiliary/background费用；继承setup不再收费。rank缓存每个实际operator_group的计算次数只记一次，24个复用execution不当成24次SVD。integration恢复历史counter不重复计；无法从原回执拆开的pending replay新aux保留unknown，但其CPU与物理费用已在真实进程总费中。

CPU_BY_CATEGORY_SCOPE按host/category/OS整进程、process snapshot和body分开。PROCESS_COSTS每个独立进程仅选一个首选CPU，内部body只作嵌套参考；不把case/arm/plot子计时重复叠加。Wall是不同且可能并行的区间，不求单个总和；不比较Mac与Windows绝对速度。scoring+plotting无法单独拆分时保留组合类。

失败评分attempt01与首次打包外部531.625s（0 ZIP）都保留。V2原PACKAGE_RECEIPT若尚未copy则pending，父线程报告数字仅列未验证provenance，不进入verified sums；补copy后用--packaging-receipt和新--out只重做metadata汇总。停止的bulk SCP、失败detached launch及未计量部署/上传CPU不填0。currency与LLM费用未知保持null；GPU_used=false只限已记录CPU程序。

| Method | Tasks | OS CPU seconds | Maxwell LU | Aggregate complex RHS |
|---|---:|---:|---:|---:|
| g19 | 24/24 | 2836.656250 | 16948 | 312804 |
| g19_s19 | 24/24 | 2889.500000 | 18297 | 329310 |
| twofold | 24/24 | 2896.250000 | 18551 | 330732 |
| twofold_hop | 24/24 | 2887.765625 | 17983 | 326940 |
| twofold_dr | 24/24 | 2896.515625 | 17599 | 322098 |
| twofold_ap | 24/24 | 2866.093750 | 16937 | 307110 |
| twofold_full_dr | 24/24 | 2859.609375 | 17021 | 308082 |
| twofold_reset | 24/24 | 2890.531250 | 18041 | 327966 |
| cie | 24/24 | 2837.250000 | 19409 | 355680 |
| dual_cie | 24/24 | 2841.234375 | 19827 | 365406 |
| dr | 24/24 | 2670.968750 | 22544 | 431502 |

原始/派生counter：COUNTERS_RECORDED、COUNTERS_BY_CATEGORY、逐11方法METHOD_COUNTERS11；进程/CPU范围：PROCESS_COSTS、CPU_BY_CATEGORY_SCOPE；缺口：GAPS；源核对：RECONCILIATION。全部CSV并有JSON副本。
