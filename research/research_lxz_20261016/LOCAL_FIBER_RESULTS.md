# 原P0条件局部校准

| 校准anchor | 方法 | compatible记录/提案 | 去重compatible | 新结构 | 最大材料MAE位移 |
|---|---|---:|---:|---:|---:|
| 250033_weak | local | 5/5 | 5 | 0 | 0.000354913 |
| 250033_weak | geometry | 4/5 | 4 | 0 | 0.00134309 |
| 250033_weak | prior | 0/5 | 0 | 0 | 0 |
| 250033_weak | continuation | 10/10 | 2 | 0 | 0.000352185 |
| 250034_ring | local | 5/5 | 5 | 0 | 0.000448741 |
| 250034_ring | geometry | 5/5 | 5 | 0 | 0.00115676 |
| 250034_ring | prior | 0/5 | 0 | 0 | 0 |
| 250034_ring | continuation | 10/10 | 6 | 0 | 0.00241842 |


已知TRAIN条件中存在多个非零但很小的noise-compatible材料位移；依原多尺度结构规则新增world为0。两anchor的8参数material Dg秩均7，1个common Gaussian-weight scale是同一材料的gauge。绝不能按latent维数报告8个物理解自由度。

最弱材料归一化方向的白化response singular scale约12.8487和297.0758；第一例转换为prior坐标的方向范数1998.55，第二例65.52。故单位材料线性位移可能需要巨大参数步，并没有可延拓的大局部fiber保证。continuation受.5 prior-coordinate trust半径后重复到相同端点，已在material-field空间去重。

第一例相对full中心差分误差8.097%，减小步长10倍后为.08097%，与二阶截断一致；精确含peak归一化Dg与原FD整体差仅8.44e−9。这说明是有限步高曲率放大，不是凭小整体Jacobian误差便可保证weak方向有限延拓。最后5个状态只做此诊断，没有重跑搜索。

E2不是N3附近反演结果；N3的合法anchor缺口仍然存在。5步链不足以估计条件posterior质量。Best-of4/SelectedMAE相对N3：NOT_RUN，而不是0改善的独立试验证明。
