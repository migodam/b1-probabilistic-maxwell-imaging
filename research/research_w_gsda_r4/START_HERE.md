# Research W-GSDA-R4：真实九数字全波反演，成像尚未通过

本轮已实现并实际运行未知墙与未知自由材料图的二维全波联合反演。数字1–9全部完成B0、B4与普通Full-GN补充，不是字形生成或标签识别。**运行交付成立，准确材料/形状成像未成立；当前不支持GSDA同质量加速。**

| 方法 | 实际数字数 | IoU均值 | 复对比度相对L2均值 | 质心RMSE / λ0 | 训练白化RMS | 留出白化RMS | 完整成像门通过 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B0 Full-FWI+TV | 9 | 0.2870 | 0.9494 | 0.1360 | 0.9606 | 0.5635 | 0/9 |
| B0+Full-GN+TV | 9 | 0.3000 | 0.9396 | 0.1413 | 0.9554 | 0.5632 | 0/9 |
| B4 原生CSI+风险筛选 | 9 | 0.0000 | 0.9574 | 0.1335 | 1.6568 | 1.4104 | 0/9 |

先看 [九数字直接对照图](DIGITS_1_TO_9_CONTRAST.png)、[损耗图](DIGITS_1_TO_9_LOSS.png)、[完整实验报告](REPORT_ZH.md)、[网页端handoff](WEB_CHAT_HANDOFF.md)。全部图源自保存的材料数组，并经实际视觉检查。

实现正面证据：原生固定GS/GD contrast-source确实运行；精确吸收、背景照明及返回接收路径在同离散完整Maxwell下成立；可逆吸收/重开有实际控制运行。负面证据：九数字B4全部KEEP_ACTIVE，没有发生自动吸收；current数据拟合与材料全波预测存在明显差距；单纯GD模态截断失真。

墙体未知范围是公开3×2容器内12个材料参数；目标为7×7、98个自由实/损耗参数。未使用真墙、真实平均墙、字形模板或有限目标标签。不能推广为任意未知壁形/孔洞数。

复现入口见 [RUN_ALGORITHM.md](RUN_ALGORITHM.md)。原研究封账后拒绝新增物理调用；新部署需显式给独立本地时间额度。只重放已有数组和评分无需重新求解。没有自动commit、push或网页发布。

范围与失败：[GO_NO_GO.json](GO_NO_GO.json)、[PHYSICS_VALIDATION.md](PHYSICS_VALIDATION.md)、[DIGITS_1_TO_9_FAILURES.md](DIGITS_1_TO_9_FAILURES.md)、[CONTROL_RESULTS.md](CONTROL_RESULTS.md)、[PHYSICS_COST.json](PHYSICS_COST.json)。本轮结束后停止，不自动开启下一研究方向。

计时口径补充：已有wrapper CPU从BudgetBook/批分区开始，早于此的Python启动导入未覆盖；没有回填估计。完整说明见 [COST_ACCOUNTING_SCOPE.md](COST_ACCOUNTING_SCOPE.md)，不据此声称严格端到端加速。
