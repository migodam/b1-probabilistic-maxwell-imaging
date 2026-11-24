# 复算四条成熟 Dual-CIE 轨迹

在本目录使用安装了 NumPy、SciPy、threadpoolctl 和 matplotlib 的 Python。原传播算子及输入沿用同一工作树。计时与证据复查可在单独工作树执行。

```sh
python -B diagnose.py
python -B score.py
python -B report.py
```

`diagnose.get_case` 返回原用例设置和公共训练数组。`diagnose.py` 从原 `engine.run_cie_arm(..., 'dual', ...)` 的每次返回读取未修改的当前 J 与总场 E，按引擎墙钟最多 30 秒，四案实际均跑完三段。每案 `updates.jsonl` 的首行是 DPCS 初态，其余每行对应一次实际更新。

每案 `geometry.npz` 的 `E` 是入射场 B，另存原接收 S、观测 Y、Green 频谱和接收掩码。`step_0000_dpcs.npz` 及 `step_NNNN.npz` 的 `E` 是总场，`J` 是电流，`incident_B` 是入射场，`m` 是共享材料。形状为格点×入射源。逐格能量和违背图均为原尺度，功率从 `J.conj()*E` 形成。

光学矩阵以 J*E 的 Hermitian 虚部定义。对角和其余元素分别记 Frobenius 范数和元素 RMS。全矩阵归一使用两个场的 Frobenius 范数之积；成对元素使用对应两个源各自电流和总场列范数乘积的几何平均。原尺度的非对角包含共轭重复项。

`score.py` 在四案均有终点和线上封存后才读三个合成参考与 Fresnel 留出端口。各源训练通道按原索引提取，留出按有效复数读数评分。物理重算采用各案反演网格上的完整离散算子，另外保留辅助电流直接得到的训练残差。
