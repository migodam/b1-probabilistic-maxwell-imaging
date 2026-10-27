# G1 — 有限joint概率与压缩

状态：**有限缓存数值模型的准确性重放完成；无新PDE速度证据。** 12已暴露parent×2noise，2950原子counting prior，已知target intrinsic/gamma，2×2不表达孔洞。

全部方法输出先封存，再reference/GT评分。M0独立exact、M4完整joint、M5 prior-conditional lift、M6普通cache逻辑的joint TV为0；所有case概率及位置风险门通过。M5在24/24保留118/118簇，未发生压缩。M4在24/24查询2950/2950原子，默认未知likelihood[0,1]直到穷尽才得到有用遗漏质量界。

M1 confidence MAP-wall与M3 task-medoid选择相同，joint TV均值.113268，仅13/24通过joint门；但位置风险门24/24通过。M2是有限bank secant角度单点消融，非成熟FWI基线，joint门0/24。所有方法都独立计2950逻辑查询；读取完整posterior或任务table的工作没有藏到reference名下。

对单墙Q，TV(P,Q)≥1−P(W=g)。最佳组为MAP-wall；这仅否定一般单墙Dirac替代，不否定保留随机材料残差的精确Green换背景。保留目标边缘而破坏wall/target cross-covariance的案例实际存在，额外held-out评分在输出封存后执行。

**原F64区间FAIL必须保留。** 原子posterior越界达2.55e−15、log证据下界差7.11e−15，没有用容差抹除。一次定向高精度归一化修复在168检查点没有观察到越界，范围仅为给定F64 potential的数值归一化，不涵盖FEM/potential形成/连续PDE。额外计算须归方法，未省查询。

[方法表](results/G1_METHOD_TABLE.csv)、[parent表](results/G1_PARENT_TABLE.csv)、[事实汇总](results/G1_FACTUAL_SUMMARY.json)、[signed audit](results/g1_offline_extra/INTERVAL_AUDIT_SUMMARY.json)、[修复](results/g1_normalization_repair/DECISION.json)、[封存记录](results/g1_replay_001/ALL_ONLINE_OUTPUTS_SEALED.json)、[factor历史](FACTOR_HISTORY.json)。全部原子权重、协方差、证据/质量区间、动作trace在Release。

统计以parent为单位，先平均其两个noise；2000次配对bootstrap仅描述开发不确定性，不升级盲测或论文级PASS。物理posterior保真仍NOT_CERTIFIED，不能把G0 .1白化门当成posterior TV .02保证。
