# 失败归因与不得混淆的结论

|分类|证据|解释与限制|
|---|---|---|
|INVALID_BASE|原N3 2/8、Research A一步GN基础3/8数据相容|5/8失败只能说明当前基础反演不足，不能证明替代世界不存在|
|TOPOLOGY_MISINTERPRETATION|figure05显示增量模；首shell实际7组件0孔洞且无真空|低增量、family名与空气孔洞不是同一概念|
|OPTIMIZATION_FAILURE|A用16/32维DCT/Haar chart、每轨迹3 full endpoints；原N3经过一步GN|容量有限且停止过早是合理待检假设；缺完整KKT，E0不能断言只是优化失败|
|REPRESENTATION_FAILURE|Sparse64改善35.6896%为离线真值目标优化；线上chart与此不同|说明存在更好表示的有限容量见证，不保证观测逆问题能在线识别或恢复这些系数|
|POSSIBLE_AMBIGUITY|部分候选数据相容但材料不同；需独立绝对结构核验|相容端点数不等于模式/连通分支数或概率质量，不得把单像素增量差视为真实成核|
|MODEL_MISMATCH|现有同数据clean/state/constitutive identities一致|目前没有缓存证据支持此归因；真值噪声残差超τ是噪声事件，不是模型不匹配证明|

旧A30补偿表中的`nominal_F/compensated_F`是两材料预测间距离，不是对观测y的重建残差。本目录另保存原表子集，不能与观测残差混并。历史Pairwise 8.84类数值亦必须依原归一化解释。

Sparse64的35.69%：64实Haar系数支撑在读development真值前由N3/TRAIN选择，但随后的系数优化显式最小化真值材料误差。它不是在线data-likelihood反演；没有保存相应full endpoint predictions，因此不做数据相容或多解认证。support选择与连续系数估计、材料误差与数据拟合都应分开。

本次读取不修改旧Gate、旧ledger或缺失记录。A原E2出现5/8无有效候选、额外候选仅1增量掩膜像素变化且无补偿已存在的结论保留。零新增物理审计无法判断全空间是否有合格新结构。
