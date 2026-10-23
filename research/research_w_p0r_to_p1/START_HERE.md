# Research W：P0-R → P1（执行中）

这是独立的新阶段，保留旧 P0/Q/B 及 P3L-R2 历史。当前仅报告已完成证据。

- 真值闭环：两个seed × 两带宽 × cal/task，共8组原无噪声均值的相同接口重放差均为0。新噪声残差约21–23。
- 48×32直接面积投影距离6.57–9.41，不能直接当作最小可达表达误差。
- Oracle-near 多步受约束Full-J TR：两个wide对象达到0.96757、0.98798。数据相容不等于驻点；几何资格另列。
- 测量驱动的FWI仍在运行，不能将Oracle结果当成在线能力。
- 最小特征原文只有0.1标量；已完成精确说明与反例诊断，不事后改写旧先验PASS。
- B2、P0-C、P1-A、P1-B尚未形成结果，状态保持NOT_RUN/条件待定。

[当前实际数组报告](report.html) · [冻结合同](configs/EXECUTION_PROTOCOL.json) · [计算benchmark](BENCHMARK.json) · [新授权](configs/AUTHORIZATION.json)

主程序严格区分 ORACLE_DIAGNOSTIC、PHYSICAL_EXISTENCE_SEARCH、ONLINE_INFERENCE。候选概率默认null。没有NN/GPU/3D或外部付费计算。
