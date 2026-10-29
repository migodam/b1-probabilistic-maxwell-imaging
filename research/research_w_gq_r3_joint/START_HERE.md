# Research W-GQ-R3-JI

已实现**未知六宏格墙体＋连续未知椭圆目标的二维全波联合成像**。主输出由测量目标函数选择，真值只用于封存后评分。

- 三parent首噪声：Full-GN位置RMSE **0.03234λ₀**；一例同时通过位置/形状/实介电率门（0.01160λ₀、IoU0.84069、实ε误差0.00656）。目标损耗及另两例形状/材料仍有失败。
- 自适应chart恢复固定chart的部分质量，但未优于普通Full-GN，不能宣称加速或新算法成功。
- 新增数字8的真实110变量全波反演：测量白化RMS0.9801，但孔洞未恢复、复contrast相对误差0.7793，两起点未收敛。1–9中只有8实际反演；见`exploration/digits/REPORT_ZH.md`及真实对比图。
- 已完成主比较 **10/18**、控制 **9/9**；全部槽位和缺失保留。
- 当前总预约 **12440 FE / 187732 RHS**；两个CPU/RSS记录缺失，没有补造。

阅读顺序：[中文报告](REPORT_ZH.md) → [41组推导公式](FORMULA_COMPENDIUM_ZH.md) → [单命令运行](RUN_RECONSTRUCTION.md) → [逐例结果](EXPERIMENT_RESULTS.csv) → [分项裁决](GO_NO_GO.json)。

图件见`figures/`及各运行目录的`reconstructed_target.png`、`reconstructed_wall.png`。真值/重建比较图是独立评分产物，在线重建本身不读取真值。原始分支、复观测、迭代与费用全部在`results/`。

本轮范围是受限二维材料模型，不能称任意未知墙、完整复杂形状或连续PDE认证。A17X已有理论与标准SOM/TSOM/Schur均明确归属。真实执行的portable复现状态为 **PORTABLE_EXECUTION_COMPLETED**；全新依赖安装未运行。没有自动发布。
