# R3交接：从条件响应阶数寻找具体新机制

先读RESEARCH_REPORT_ZH、THEORY_DERIVATIONS、RESEARCH_LOOP、FINAL_INDEPENDENT_REVIEW和cost。用户特别要求公式优先，反对复杂算法拼接和等价换写冒充贡献。

已完成11个同一公开R2弯条的开发控制，全部未恢复。不可写0/11独立盲测。噪声DPCS floor真实存在但clean也失败；loss降低residual却NRMSE从.511恶化到.817；同Phi尺度修正后PGinf9.92e−9仍r=.100734。强模态、GlobalDual、单cell二次根、全部165600swap、邻接以及full-wave-interaction paired rays都没有实质逃逸。

新可用认识：固定材料的硬锁状态floor和固定CN分布；被动rank-one响应圆弧；无prior的rank1内部局部数据最低点已是其圆弧最低点，但box边界能陷住真实Hankel模型；更新算子rank-r决定r维响应包和≤4r−2数据导数多项式阶数。3D isotropicvoxel是rank3，不是1。工具是经典的，首次性未认证。

真实onecell m20→m1.2解析跳跃有效，但成熟full-current CIE readout也1.2。此结果解释baseline的一项能力，不支持beyond-CIE。3D8cells三源的dyad/isotropic响应SVD范围清楚，无3D成像或连续精度主张。

下一轮只在明确公式预测上选择行动：例如从O C^k Q中寻找确实可见的协同反馈，并给一个具体有限路径为何超出现有rank1/2字典的证据。先比较这个方向与同成本普通方向/CIE；不要自动增加block尺寸或枚举更多根。另可研究固定材料噪声floor如何扩展到unknown/data-selected材料的合法界，但当前GT谱不能在线选门。

新算法若出现材料收益，冻结后再运行已准备的fresh generator（DRAFT、目前未运行），使用新entropy/独立GT评分。所有本轮材料都已公开；不能重新称blind。R1/R2不得改写；新增工作另开隔离checkout。发布已获用户授权，沿用B1仓库和直接raw中文入口。
