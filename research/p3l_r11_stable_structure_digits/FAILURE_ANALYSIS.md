# 失败、纠错与未执行记录

1. **FAIL_CAP_EXCEEDED。** 跨阶LU上限12,000，实际28,727。局部账本只保护RHS；多个对象与小阶增阶把LU累积漏过全局约束。72个预测已生成时由Sol-B审计发现，Main即时叫停三个分支。原始账本和初始合同保留，未将allorders改成LU256，未扩预算，也未继续物理计时重放。费用不合规是整轮限制，不是可以因数值好看而忽略的小问题。
2. **主线程时间标签错误。** 初始EXPERIMENT_CONTRACT手写20:07时间标签晚于真实文件mtime与首个inverse。审计后保留原字节于audit/INITIAL_CONTRACT_FUTURE_CLOCK_LABEL.json，修正为原文件实际mtime；阈值、方法和数据不变。INVERSE_CONFIG实际先于首个inverse。不能靠修正隐藏原错误。
3. **B/C嵌套计时重复。** CURRENT_SHAPE_PROJECTED已包含在EVALUATE_COMPLETE，原physics_kernel_wall重复计。原回执不覆盖；派生corrected_kernel_wall只求父级事件+最终精确验证。全wall、RHS、LU和数值输出不变。
4. **评分依赖缺失。** 初次评分因未安装skimage在读GT前导入失败；保留SCORING_DEPENDENCY_FAILURE.log，使用明确7×7 uniform-window SSIM公式和scipy滤波，未改二值阈值、主MSE或模型。无新增物理，无训练。
5. **Sol-C瞬时容量失败。** 首次Sol-C agent返回modelcapacity错误，重试同一个Sol最高强度后完成；没有Astra或DeepSeek替代冒名。
6. **固定总场失败。** 注入固定仅移除weighted支持障碍，CA不保持固定Krylov空间。anchor128仍73/192真实超0.01；Green128在C组全部真实超门。并非仅上界保守。
7. **全局耗散界未加速。** b*>0且便宜，但真residual不够小、端口/势门保守和精确J成本共同消除收益；不把理论界PASS改写成计算PASS。
8. **A核与整流程分离。** 稀疏truth上强正面核结果不能代替未知支持反演；N全流程较慢，P只单窗口观察节时。没有全系统普遍加速主张。
9. **图像与先验。** N80J未充分收敛；Ptruthlatent极不典型；unknownproperty初始化恰在truth中点。模糊与丢孔不能单因归结为物理不可辨识。9个设计对象不足大样本泛化。
10. **NOT_RUN。** 专用P1有限支持出生、真正TSOM电流材料联合反演、受控低阶Jacobian、A+B、A+B+C、大样本确认和Full66 posterior。不得以公式或设计替代结果。

11. **非结构性实现开销。** B/C精确总场Jacobian另付输出灵敏度恢复乘法；global谱常数主反演每对象冷建，未跨九对象摊销。主wall比值不能全部归于子空间理论。100% trial回退与真实额外RHS仍是本实现的否定证据。A对R10 rank64未完成统一窗口逆成像计时，不能用不同核窗口冒充。

12. **数值表达门不足。** P3虽然MSE与IoU过数值门，实际真值中横缺失，标DECODER_REPRESENTABILITY_FAIL_VISUAL_THIS_FIT。初始归档保留原数值门和图，当前补充更严格的视觉资格，不改GT与模型。本轮拟合未收敛，不能据此证明整个decoder数学表达不可能。
