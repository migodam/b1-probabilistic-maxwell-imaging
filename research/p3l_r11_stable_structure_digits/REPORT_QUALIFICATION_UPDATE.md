# 报告资格补充：数值拟合通过不等于数字可辨识

这是主线程在打开真实P数字图后作出的更严格视觉资格补充，新增Maxwell RHS=0、LU=0。原始Release与数据不覆盖。

原生成协议的ROI占比MSE≤0.025且IoU≥0.75门保持不变。九个P对象都通过了这项**数值**门；然而P目标3的中间横笔在本轮离线拟合中消失，实际真值更接近右括号。依据任务要求，其语义视觉资格为 **DECODER_REPRESENTABILITY_FAIL_VISUAL_THIS_FIT**。P目标5也有明显上部填充失真，但仍可读作5，记为DEGRADED_BUT_LEGIBLE。

这不证明Gaussian64 decoder整个材料像无法表示3：本轮拟合有1200步上限且未全部收敛，潜在表达能力与构造优化失败尚未分开。该标签只说明**本轮没有构造出清晰P数字3**。

![实际模板与P真值](figures/DECODER_3_5_VISUAL_AUDIT.png)

- 九个原生N数字仍全部保留，均做了四个后端真实反演。
- P九个预定对象也全部保留，没有用更好看的替代真值，没有删除P3。
- 72份预测、材料MSE、IoU、原数值fit门、噪声球门及费用全部不变。
- P3材料误差仍相对它实际的合法decoder真值计算，不能将其解释成对清晰数字3的恢复质量。
- 原始完整tar.zst和reading kit保持不变；补充包只修订报告资格、中文视觉说明与派生状态，不改GT、测量、latent、预测或物理账本。

执行数据固定在科学提交7d4d1fa2930afacdd8f683c3025e5538b4f3be8a。报告修订身份由R11_REPORT_UPDATE_MANIFEST.json另外记录。资源协议仍为FAIL_CAP_EXCEEDED。
