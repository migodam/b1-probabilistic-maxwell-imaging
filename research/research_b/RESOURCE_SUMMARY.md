# Research B 实际资源汇总

独立授权为本次Research B的4096个full-Mie material-frequency attempts、E1/E2受控wall7200秒。E0新Maxwell必须为0。历史P2/B1/C1/A余额只读、没有转账或重置；GPU0、一般非同心3D新campaign0、硬件测量0。

| 实际物理用途 | 新attempts |
| --- | --- |
| E0缓存排名/噪声模拟 | 0 |
| 新四频625节点bank | 2500 |
| 独立64材料×4频n14审计 | 256 |
| 前8材料×4频n18与独立递推 | 64 |
| 固定noise锚 | 4 |
| 新32parent×6频truth | 192 |
| 非线性共享补偿端点（缓存复用后） | 22 |
| 同一个小孔世界六频审计 | 6 |
| 四个真实A包未来field/新truth/双频corrector/参考audit | 240 |
| **合计** | **3284 / 4096** |
| **剩余，不再为消耗额度继续运行** | **812** |

物理调用失败0、未完成0；所有失败/审计/truth调用先记BEGIN、统一计费。kernel本身wall=11.187805秒，不等于完整选测计算成本。

主要COMPUTE_LEDGER共16个job，E1/E2受控wall=370.412760/7200秒；包括E0/统计/导出等全部已入账数值wall=377.550891秒，CPU=373.482364秒。一个adaptive报告导出job失败，已完整算入；仅修正摘要键，积分没重跑。

另有单独回执的缓存/adapter/早期合同测试/绘图（未进主ledger）：wall=6.165147秒、CPU=5.520746秒。去重后已知数值总wall=383.716038秒，CPU=379.003110秒。即使把全部这些任务保守归入7200秒上限，仍低于它。没有将重复TEST_RECEIPT_01..03或plot stdout再加一次。

按后端metadata共有92896个boundary/interface线性系统与92896次factorization记账；recursive实现早期遗漏计数通过追加ACCOUNTING_CORRECTION补1584，原END行保留。这是此Mie实现的内部系统计数，不能叫额外硬件测量，也不是通用PDE source-RHS成本。

环境为现有Gaussian/.venv_nn Apple Silicon CPU，BLAS threads1，受控job nice10，串行物理。没有安装基础设施、启动GPU或干扰并行A。记录只覆盖受控数值任务，不假称LLM、写代码、阅读、写报告或Git上传延迟已计为数值CPU。

物理台架一次共享所有菜单预测与模拟setup，再比较策略；逐方法的measurement blocks是真实模拟消费数，portfolio墙钟不能当各独立policy的优化运行速度。calibration/外部geometry各算一个block，没有真实切频/时间/能量报价。

原始追加账本：results/PHYSICS_LEDGER.jsonl、COMPUTE_LEDGER.jsonl；主摘要RESOURCE_SUMMARY.json，去重附录RESOURCE_COMPLETE.json，递推内部系统说明RHS_ACCOUNTING_APPENDIX.json。E3的FINAL_COST.json为当时快照，阶段最终值以上述最新主ledger为准。
