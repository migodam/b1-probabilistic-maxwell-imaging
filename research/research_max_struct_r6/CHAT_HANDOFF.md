# MAX-STRUCT-R6 给Chat的入口

## 首先读取

- [中文完整报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r6/RESEARCH_REPORT_ZH.md)
- [网页和全部图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r6/)
- [数学结构](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r6/THEORY_DERIVATIONS.md)
- [逐例数据与费用](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r6/RESULTS/analysis/ALL_CASES.md)

## 已实际完成

21个材料×CIE/Dual-CIE/Global-Dual/A-material＝84个主任务；另两个预指定对象×Dual/A＝4个真FFT任务。A是普通global current+仅材料块one-hop，继承R5精确重复求值缓存。固定mf3/6/10各500＋mf12共1500辅助步，共同原Maxwell corrector上限3000；576独立实未知量，无GT在线选择。

| 分层联合成功 | CIE | Dual | Global | A |
|---|---:|---:|---:|---:|
| 历史数字9 |8|9|8|8|
| 历史困难6 |1|1|1|3|
| 新确认6 |4|4|4|4|

A额外恢复两个已公开曲条：7b8d NRMSE .04686，b20d .04128；不是新盲测突破。42cb/63ab及历史环仍失败。数字5 A因IoU .6667未过联合合同，NRMSE .1142和数据拟合合格，不能只说“数字5完全失败”。

新共同失败：82ac偏心薄环与40ac开缝壳。真材料可行Phi分别8.70e−6/8.55e−6，算法约.00485/.00381，残差约GT的91–103倍。说明没到达已知可行低目标区，不是本轮已证明的数据不可辨识。薄环四臂达到预算上限，开缝壳四臂梯度停止；不可混淆。

## FFT与NIE

可真正build upon FFT-TSOM：native Fourier action/adjoint＋非周期零填充Green，A读出不变；接收和最终物理LU仍dense。数字1、旧b20d的A CPU分别59.02→49.76秒、70.58→54.44秒，均保留恢复；同FFT的Dual仍更便宜。四个工程终点非逐位相同，不宣称通用纯速度等价。

NIE2016作者全文已取得，原文已有FFT-TSOM实际方法；本轮没有完成独立成熟NIE代码复现，因此主比较是可验证CIE/Dual，不用self-only冒名NIE。

## 成本与判断

主A总CPU1504.32秒，Dual932.81秒（1.61倍）；主157418LU/3773976RHS，FFT6607LU/158376RHS，setup/失败另列。全部新增已知物理164051LU/3932556RHS（其中12LU是16维接线小模型）。主四worker队列运行合计1140.37秒、另修复间隔41.26秒，FFT55.36秒；不要把重叠task wall之和叫历时。

结论：A是有局部互补价值、已实现FFT的候选插件，还不是稳定优于Dual的默认方案；没有3D或实测新成像结论。当前最有价值的后续问题，是P dr+(dP)r如何把材料切向正交补的残差投到有效更新方向，以及四cell传播为何未改变跨腔错误结构。先对成功7b8d/失败开缝壳做状态固定的公式诊断，再新材料确认；不靠扫β或增加网络。

R4 B/C没有重跑，原否证保留。R1–R5报告/数组未改写。新六对象的GT现在随交付公开，之后不得重复标盲测。所有路径、失败与复现入口见同目录README/REPRODUCE/COST_LEDGER/FAILURES及audit。不要把本交接中的待验证假说当结论。
