# 商坐标下的独立物理驻点

原Research F Gate2失败事实保持原样：ring .0452734、weak .00200099，对原.001标准均FAIL。本轮是独立新坐标、新停止准则，不回写历史。

从两份原F1已支付端点继续，未从背景重新初始化，未读取N3或GT。使用原六shape uniform坐标与比例t，共7维。每步有界线性最小二乘，真实full-forward回溯1/.5/.25，单条件最多8新states，所有回溯计费。

冻结新标准：七坐标符号正确的projected-gradient无穷范数≤.001（数值active-set1e−8，同时报告exact边界），有界线性GN材料RMS≤.00527555/1000、预测Φ降低≤1e−6，且原P0合法、原noise球相容、材料/响应rank7、后向误差<1e−10。最后实际点均到达**精确比例边界t=0**，不是仅用近边界容差删除梯度。

| 条件 | 新state数 | projected gradient | 可行GN材料RMS | 预测Φ改善 | 实际残差 |
|---|---:|---:|---:|---:|---:|
| ring | 3 | .0005792604 | 5.29818e−8 | 9.30783e−10 | 22.239554758 |
| weak | 1 | .0005356274 | 6.76656e−8 | 4.56623e−10 | 21.512112256 |

两者通过新标准，允许C1比较。总4states/160RHS，包括4份新receiver Jacobian。原缓存初始state/J被读取，不重新求解；其历史40RHS/条件另列amortized来源。

这是likelihood局部有限精度驻点，**不是posterior MAP**。在t=0处，rho梯度分别+.01198/+.11067，方向符合下边界KKT；但商先验密度为0，posterior密度为0。零密度不等于材料不合法，也不等于其任意邻域概率为0。

优化只把原端点微调很小：原Φ已经接近受约束最优，过去“没有完整合格anchor”的障碍在本轮解决；这没有创造新的exact null direction。两个材料及响应的数值rank仍7。

完整梯度、谱、每trial真实Φ与原始材料均在results/F15及endpoints。独立只读重算见evidence/INDEPENDENT_F15_STATIONARITY.json。未报告全局最优或数学KKT精确证书。
