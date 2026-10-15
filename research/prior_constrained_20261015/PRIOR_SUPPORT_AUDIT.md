# E0 原先验支持审计

材料规则从 A26 实际 config/生成器复原：`chi=chi_b+v*rho`，`v=(1+.2i)/sqrt(1.04)`，原幅度为 .3/1.2，a1=1.2。weak背景为 `.1+.04i`，ring在 `.2<=r<=.35` 为 `2+.1i`、其余weak。rho=0回到背景；绝对chi=0空气孔不属于此增量规则。

原 P0 是四类等权 TRAIN 生成分布（shell/asymmetric_piecewise/smooth/mixed_scale）及原参数support、硬rasterization、峰值归一化。64个独立TRAIN parent对应512条件行。`delegated/prior/TRAIN_REPLAY.csv` 记录512/512逐元精确重放，最大误差0。原生成器合法参数且端点重放一致才是P0_EXACT；PHASE、BOX、低denoising误差均不能替代。

新增256个独立material-only parent，每类64，使用新master seed20261015；weak/ring端点均通过必要相位/范围检查。另有256个既有TRAIN图像插值诊断，非P0_EXACT且不进候选reservoir。保守计512个计算材料数组，低于2048上限；新增paired labels、Maxwell states/RHS、GPU训练均0。完整参数/seed在`delegated/prior/MATERIAL_ONLY_PROVENANCE.json`，support检查256/256通过。

2835审计记录中1671满足PHASE必要条件、1164为PRIOR_UNKNOWN。这里包含别名/重复，不能按记录数声称独立候选。N3原图32/32、T raw40/40及其所有corrected记录、G raw72/72与所有保存correction/endpoints均不满足原必要相位/范围/峰值规则；这些特定端点被必要条件排除。未知候选未能retract不构成不属于P0的证明；保留PRIOR_UNKNOWN/RETRACTION_NOT_FOUND。历史PHASE合格pool没有自动成为P0_EXACT。

结构校准只用TRAIN材料与栅格扰动，不读y、clean response、development truth或MAE。rho/a的excursion阈值.1/.25/.5/.75/.9，尺度0/1/2 cells；centered shape、occupancy、boundary下限分别.08796801446026997/.02549189814814815/.03399884259259259，另需至少3个阈值的跨尺度occupancy差>=.05。DC单列。细节见`delegated/prior/STRUCTURE_CALIBRATION.json`；最终联合门槛由主线程冻结，不可按pair结果调低。

证据限于生成器复原、TRAIN重放、material-only支持与历史必要条件审计，不建立多posterior区域或质量。
