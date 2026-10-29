# 真实图件索引

所有图由保存的数组、封存评分或已付费Jacobian生成。没有理论预期曲线冒充实验。

| 图 | 内容与解释边界 |
|---|---|
| [目标材料总图](figures/R0_TARGET_MATERIAL_OVERVIEW.png) | 三个独立parent、首噪声、真值与三主方法；共同实介电率色标，包含成功与失败 |
| [未知墙实部](figures/R0_WALL_EPS_OVERVIEW.png) | 所有六格均在线未知；空气缺块连边界，不是封闭内孔 |
| [未知墙损耗](figures/R0_WALL_LOSS_OVERVIEW.png) | 与实部独立估计；不能用几何成功替代损耗准确性 |
| [质量与成本](figures/R0_QUALITY_COST.png) | 同parent位置/IoU/实部误差及FE/RHS/墙钟；两CPU缺测明确留空 |
| [全部优化分支](figures/R0_ALL_BRANCHES.png) | 全数据目标函数与离线目标误差；最低目标函数选支，未按真值选支 |
| [实际迭代轨迹](figures/R0_OBJECTIVE_TRAJECTORIES.png) | 9个主运行、各四起点；非收敛保留 |
| [控制臂](figures/R0_MATCHED_CONTROLS.png) | prior-only、wall-first、known-wall及普通joint；未完成时图中留空，按最终图版本读取 |
| [同状态切向与目标Schur后果](research/delegated/sources/CONSEQUENCE_FIGURE/same_state_consequences.png) | 已付dev Jacobian的离线分析；33记录含重复来源，不是33独立物理样本 |

各图同名PDF可用于放大。每个实际重建目录还包含纯在线输出 `reconstructed_target.png`、`reconstructed_wall.png` 和独立评分后的比较图。图中位置单位为参考自由空间波长λ₀。二次噪声的图与表如已生成，以最终`EXPERIMENT_RESULTS.csv`与目录为准；不能用缺失行补零。
# 新增关键图

- `exploration/digits/truth_vs_reconstruction.png`：真实数字8、自由49格材料反演及误差。只有8实际运行，图不代表数字1–9均完成；两个孔洞未保留，两次优化到迭代上限。
- `figures/R0_FULLGN_CONTOUR_ZOOM.png`：相同0.5λ₀视窗的连续边界，兼报半轴、质心相对尺寸和评价格细化。用于区分视觉偏差与IoU栅格误差。
- `research/delegated/sources/CHART_MECHANISM_RUN/FIGURES/public_2x2_mechanism.png`：同一参考chart下墙变化与强目标插入的独立/交互作用。强目标自身变化已足以导致严重切向失配；不是一般TSOM失效证明。

