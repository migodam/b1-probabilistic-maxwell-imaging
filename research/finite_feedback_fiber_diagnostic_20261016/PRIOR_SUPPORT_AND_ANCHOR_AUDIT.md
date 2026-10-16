# 先验支持与anchor来源

E0新物理0；原缓存和此次E1真实测量严格分开。

1. **Physics-only**：同观测BASE_N2，背景adjoint初值+一次R-ROM梯度box投影，无NN/GT初始化。E0只找到材料图，无独立完整F/J缓存。此次付费full-state/J后残差120.4323；先验phase失败。不是已充分收敛的物理重建。独立迭代full-GN/TSOM缓存MISSING。
2. **N3对照**：冻结A26已暴露validation，非fresh验证。E0已有完整状态/J，缓存残差28.0347、N3–GT预测距15.5526；新付费重建用于统一factor/receiver接口，明确重新收费。phase失败，不是原P0 anchor。
3. **Known P0**：原LXZ固定两份TRAIN smooth材料，原generator参数逐点重放；残差22.3130/22.9696，均过原95%界。只用于数学合法路径校准。没有将TRAIN已知材料作为development重建。
4. **PHYSICS_REFINED_N3**：A17/N3_APPLICATION_R1 full-GN/ROM均从N3起点，不能冒充独立physics-only。旧结果与MAE不覆盖本轮anchor资格。

原P0有四个family与固定增量相位/原峰值归一化，不等价于任意合法复pixel box。rho=0是原背景，不是绝对空气。smooth条件chart是原P0的条件子支持，不同于整个先验分布；没有新增P1。

首例truth noise残差24.090925>tau23.785449：这是历史公开噪声尾部，维持原sigma/gain/阈值。不通过挑新noise或放宽阈值获得合格GT。GT直线路径除最后endpoint外无原先验证书，统一OFF_PRIOR_ORACLE_PATH。

已读取“逆散射中的Dyson级数”原会话的finite-resolvent/Neumann/rational-Schur讨论，只将其视为理论上下文；私人会话正文不随仓库/Release发布。原始closure_final与A26/A30/LXZ源码证据见evidence/CONSTITUTIVE_CONTRACT.md和CACHE_AUDIT.md，未改冻结论文/历史账本。
