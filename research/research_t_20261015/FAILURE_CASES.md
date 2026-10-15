# 保留的失败与限制

- 原N3相容失败30/32，其中固定子集6/8；Research A一步base5/8失败。它们是无效base，不是“无其他解”证明。
- full-space GN平均MAE0.167560，高于N3 0.041756；最低残差base选择0.147656。修复相容没有修复成像。
- v1：p250069_ring/edit1最终27.691464>tau；v2六个、v3一个失败；40个raw都失败。逐端点路径/数组和物理费用未删除。
- v2减少额外空气却降低成功；v3换base不消除补偿后材料误差。所有Best-of-K无改善，不继续调参制造正结果。
- 32truth无严格真空单元。生成空气核是另一种数据相容解释，不能描述为恢复真实孔洞。
- E0有4/32 truth-noise球外实现；这不等于model mismatch，原sigma/tau未改。
- 离线最终CSV首次刷新因重复字段TypeError失败；原始异常和wall=null/UNRECORDED在offline/FAILED_REFRESH_RECEIPT.json。后续修复没有新物理。
- 物理求解数值异常0不等于无失败：raw、失败中间端点、最终不相容材料均已计费。未记录的历史setup/单次后处理wall保持未知，未写成0。
- 没有uniform output error bound、continuum mesh验证、完整覆盖、same-cost A17增量、独立生成泛化或3D证据。
