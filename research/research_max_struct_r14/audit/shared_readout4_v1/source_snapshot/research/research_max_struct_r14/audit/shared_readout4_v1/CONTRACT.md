# Shared material inner program：四个曝光对象机制控制

首次运行前冻结：0002/0017/0032/0037 × 00/01/10/11，共16条。本组只排查material列同时改变targetmetric与原内循环程序的混杂，不是新主算法，不作独立24/全48恢复结论。RESEARCH_LOOP最后一节给出目的；原冻结算法和独立24算法均不修改。

四组合的current/PR/sourceweights/stage/mask/beta/box/prior沿用原core.step。仅通过`R13 full_current._local_function`复制step私有globals，分别把`mobius_material_readout`与`_baseline_readout`绑定为共同χ-coordinate L-BFGS-B读出；无全局monkeypatch。M0材料beta=b0，M1材料beta=b1，current指标族的选择完全不改。

每次material读出都从同一个current m出发，原χ坐标、原box/prior，统一maxiter60、ftol1e-12、gtol1e-8、maxls25。没有scalarseed、Arc、MM、frozenmetric或额外候选。M1调用原`onehop_material_value_gradient`完整dP/dm。M0接收core已计算的J@Ls/Q@Ls与onesweights，使用`R=β(Jw−mQw)/(1+βm)`，`dR/dm=−β(Qw+βJw)/(1+βm)^2`，目标`.5*rho*||R||²+prior`与完整实梯度。若最终material目标非有限或增加，保留旧m；所有目标/梯度/失败/回退均计费。

使用原公开uniform初始化，不读GT。`run_shared_material.py`逐行复制当前`run_experiment.py`，唯一变化是StatefulDualMetric的import替换。故120绝对整进程CPU、最多72 proposal、随后原corrector使用到绝对120；不是3D strict48合同。所有paidcandidate/physicalaudit/history/primary/secondary/raw/corrected/PRcheckpoint/费用和超限排除沿用原runner。shared solver标签在manifest task、nativehistory和完整checkpoint frozenflags中记录；原Arc配置字段保留为历史输入但不会被私有替换后的readout执行。

只在现有AppleSilicon CPU上运行两个本地children，每child单BLASthread；不声称独占/绑核，不与Windows绝对CPU混合。原public文件路径复用，不复制/读取truth。source/config在首child前归档；16条全部终止后封存预测并保留任何missing/failed/outbudget，之后仍需主线程另行授权GT评分。
