# N3与本轮目标

N3保留原始A26单图材料风险基线，未重新训练、未用作初值、未应用新physics corrector。本轮真正的在线目标是原conditional P0八个uniform坐标中的full-physics likelihood优化；它不是N3增量训练，也不返回posterior。

| 条件 | 初始复MAE | physics final复MAE | 缓存N3复MAE | physics RMS |
|---|---:|---:|---:|---:|
| ring | .03218134 | .0007283783 | .03870807 | .002242914 |
| weak | .03218134 | .001554472 | .03147809 | .004913129 |

这是同一parent两个已暴露背景条件，physics额外已知family和global geometry，N3条件不同；不能把两例优势写成同信息量的普遍算法胜出。所有GT指标在F1结束后独立计算。N3缓存与这两例的GT数组逐元素相符，但没有其full-forward prediction缓存，所以N3 Maxwell residual保持MISSING，不为补齐它追加求解。

N3实际训练含材料MSE与SSIM，不能直接称posterior mean。Explorer若未来成立，应覆盖合法材料不确定性/有限传输，而非强制降低单张图MSE。本轮没有posterior event mass、Best-of-K、在线选GT或新候选多样性结果。
