# 给下一轮GPT Pro：先看反例与完整state

## 核心结论

R14完成576次全程反演。曝光48四pair00/01/10/11为29/35/30/31；独立24为8/12/8/11，成熟CIE/Dual各9、G19 11、DR14。A+相对CIE/Dual各增加3个且无丢失，成对同CPU及共同LU/RHS池仍12比9。新11独有正例r14_new_0021（高contrast patchwork，NRMSE .087/IoU .925），其余全部旧/新对照均差；01却在0012/0015明显胜11。

01→11从机制提出、在独立预测前冻结，但仅保留01同12集合，无新增恢复。不要把它写成已经解决收敛的新算法。看到0021后才加11的portfolio属于后续开发，不是本轮独立确认。

## 公式不能混淆

- M0=P0*P0，P0=b0/(1+b0m)；M1=P1*P1，P1=(I+DaGo)Dp。Go是原disjoint四cell offdiag，不是fullG。
- b0=1；b1按1/.2/.2/.2。00不等于原Dual-CIE；全material导数保留，P0也材料相关。
- 理想T_m需要精确根/无prior/固定活动面/独立current切空间/精确两块BR。native有限PR、冗余Fourier、guard和非根状态没有自动满足这些条件。
- 非根g_metric=ρ Re⟨rΩ,(DP[h])e⟩在根处O(||e||²)，实际56个first-step上范数常与构成力同阶。完整数组保留，不能只拿一个norm相关性当恢复定理。
- 所有复变量须正确realify；源metricΩ与观测whitening W分开。

## 优先读取

1. `RESEARCH_REPORT_ZH.md`、`THEORY_DERIVATIONS.md`、`BASELINE_AUDIT.md`、`THEORY_STATUS.md`。
2. `analysis/independent24_v1/POSTHOC_GAIN_LOSS.csv`、`DESCRIPTIVE_PAIRED_BOOTSTRAP.csv`、`FAMILY_REGIME_AFTER_SEAL.csv`；特别看0012、0015、0021、0009、0020。
3. `analysis/shared_readout4_v1/REPORT_ZH.md`：0017去Arc后仍只有01好。
4. `analysis/self_only16_v1/`、`analysis/stage_control16_v1/`：非对角Green与早期次序。
5. `analysis/recovery_dynamics_v1/`、`audit/material_force_decomposition_v1/`、`audit/local_closed_loop_v1/`及`coupled_mode_trajectory_v1/`。
6. `analysis/pairwise_cost_independent24_v1/`、`COST_LEDGER.md`、`FAILURES.md`与原始候选/状态。

## 数值证据怎样读

关键完整checkpoint为0002@128/512、0017@128/512/2500、0037@128/512。每个复制后四metric×PR保留/清空，保存step0/1/8/32。`StatefulDualMetric`用NPZ+typed JSON恢复：m、J、Q、A、old_grad、old_direction、old_coefficients、DPCS/源metric、stage与β、FFT及onehop cache都在；不要把相同m当相同state。

恢复接口见`src/dualmetric_core.py`与`src/intervene_checkpoints.py`。`state_dict`定义的是完成block step后的数值状态；runner的CPU时钟/paid-candidate pool另存，不能声称跨设备逐CPU秒bitwise重放。所有关键固定步重入已核对。corrector的内部L-BFGS历史没有伪装成已保存outer PR状态。

`raw/primary/corrected.npz`的J/Q是完整物理审计场；solver checkpoint中的J/Q是辅助状态。raw可来自较早的best proposal，final checkpoint不一定对应raw材料。请按candidate/iteration/source字段连接，不要串错两类场。

根实验使用64材料，native64系数只有48独立current维；ideal在独立空间算。材料-only32步收缩不是渐近反例；完整兼容方向在原guard尺度出现128次nit0，也不是渐近稳定证明。原±方向、三个固定扰动尺度、模式数组、正规方程误差分解、counter和失败均保存。

非根力数据为56×576完整向量和J/Q/m/Ω等字段；`FIRST_STEP_MATERIAL_FORCES.csv`给出const/metric/prior对真实Δm的work。它不依赖GT。构成力上坡38/56、full上坡0/56，说明冻结metric投影解释不足。不要将chain rule本身算新定理。

## 尚未回答的精确问题

1. 能否从一个非根完整state预测“该用01还是11”，同时解释0012/0015与0021的相反结果？当前根谱、短期Φ和单个signed-work统计均不够。
2. 有限current正规方程误差传递、非根metric导数力、表示profile偏差与PR记忆，哪项决定早期分流？已有same-state干预可构造反事实，不需重做大扫。
3. 一个有用预测应能反对错误选择：例如别把新0002的IoU过线当成大basin突破，别因一个好Φ方向就宣称材料方向更正确。
4. 五严重失败GT有更低Φ，但多数端点还有微小下降。如何区分病态长谷、可访问集合受限和真正分支障碍？现有证据未证明严格局部极小或材料唯一。
5. 如何保留多个合法候选而不由GT选？原公开01→G19→DR为14/24且额外付费，不能拿离线并集冒充统一算法。

## 费用、物理及文件边界

成功solver实际OS CPU64294.503713秒，失败16次另.514267；setup/validation/scoring/归档/上传另外计。二维576阶与三维648阶LU分开，DR背景Cholesky/RHS独立；cached四cell block RHS不算Maxwell RHS。

三维是216实χ/648 vector current、60合法外部dipole，full analytic material inner，非R7Arc原样移植。3个旧几何4条件仅有限DDA证据，继承Mie/网格误差；新Fresnel/Austria/NN四metric未运行。

原始fullfield ZIP与source/public、offline truth、curated Pro包分开。material-only review NPZ只保留同值m，完整fields仍在对应raw archive；任何恢复都要按README/manifest映射根路径。原seal不因下载mtime变化而覆写，本地复算另建派生seal并审查授权。原ZIP/Chat内容仅作研究来源；当前用户任务及冻结合同才是执行依据。

完整数据入口见DOWNLOADS与PRO_PACKET_CONTENTS。本轮所有新代码/报告在独立R14目录，R4–R13不改写。
