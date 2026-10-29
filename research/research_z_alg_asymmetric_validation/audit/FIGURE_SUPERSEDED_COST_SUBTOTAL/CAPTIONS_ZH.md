### 01_geometry

非对称固定标签小系统与三个公开参考材料。6实材料、4实观测；不是Full66。

Scope: New six-cell contract; three deterministic references; no population inference

### 02_material_response

三参考上下文的材料位移与独立全前向响应误差。叉号保留有材料端点的失败；无合法端点失败见覆盖图。

Scope: No continuous guarantee between samples; no global fiber coverage

### 03_paired_coverage

同一目标逐项配对：两方法覆盖一致，失败未删除。固定grid均49/49；radial边界终止不代表全域边界认证。

Scope: Sampled target coverage only; no additional branches or global completeness

### 04_cost_windows

两个费用窗口分别保留。优化equal-trace总费用仅约13.6%改善，未过20%门；solver加validation约28.6%不替代总费用。FrozenJ同时展示。

Scope: Timing repeats are not independent parents; original and optimized windows not merged

### 05_branch_search

33起点中31个合格、仅找到1根，2条数据往返路径回到原根。循环途中目标变化，非固定数据纤维，不声称monodromy。

Scope: No bounded-domain uniqueness or branch-completeness certificate

### 06_probability

旧二单元有限噪声后验的正确事件概率约.67516447，错误均匀分支为.5；新六单元只计算局部tile条件核，正确约.50395231，错误权重为.5。两者不可互换。

Scope: New local tile mass is an unnormalized coarea density integral, not a probability or noisy evidence; global missing mass unbounded

### 07_representation

完整正交基只换表示，原prior经推前保持；三系数投影改变材料与支持。728条封存端点含重复，不是独立样本或后验抽样。

Scope: No forward responses computed for truncation; no fiber/gauge or Full66 extension claim

### 08_certificate_scope

保留所有较宽区间失败尝试：仅冻结binary64系数的窄局部graph获证；.05数值积分tile没有整体区间证书。真实Hankel系数包围、连续PDE与全局质量均未认证。

Scope: Exact frozen-coefficient local graph only; no global missing-mass bound; failed bounds are not nonexistence proofs
