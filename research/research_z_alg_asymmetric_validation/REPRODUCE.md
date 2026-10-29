# 复核与复现实验

本轮CLOSURE已封账，Book拒绝再次启动。不要删closure/ledger，也不要覆盖raw或input重置预算。

## 零新物理审查

在还原完整证据、安装相同NumPy/SciPy环境后，使用独立审阅脚本：

```bash
python audit/results_independent/audit.py --output /tmp/zalg-independent-review
```

具体脚本相对本研究目录；读取保存的A/RHS/current/J和中央账本，重算材料/数据/成本与概率资格。它不调用forward、不构造Book。图脚本只读已有数组。原输入报告和源码保持只读。

## 原执行顺序

E0解析/两cell→序列化缓存恢复→非对称probe→P0 A/B、P1 B/A、P2 A/B→小区间图册certificate→E3换基→E2多启动与往返数据continuation→元数据恢复→旧E4新32/缓存64→新局部coarea积分→70dps独立检查→公平完整trace成本控制→current代数trace重放→738端点伴随rank诊断→封账。所有失败也计费。

若要重新支付物理进行独立复现，须在新的输出目录、独立新授权与账本下保留相同EXPERIMENT_CONTRACT和源版本；不能把当前封账额度转给新项目。env限定OPENBLAS/OMP/MKL/VECLIB线程1，float64/complex128，单物理worker、一个live物理LU。库版本及SOURCE_MANIFEST保存。没有NN/GPU/cloud依赖。

原作者E4默认32/64不能直接重跑后只按10RHS记账；新适配只算32，明确复用64。闭式响应也按source-RHS等价单列。物理数组是无损压缩NPZ，可直接np.load(allow_pickle=False)。完整归档逐字节回读，阅读包与完整过程包分开标注。
