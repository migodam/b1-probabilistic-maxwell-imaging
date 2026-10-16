# 是否已经值得学习有限传输？

**目前NO-GO；本文件给出不执行的接口与必要证据，不能作为训练完成说明。**

当前不存在新生成的可复现非平凡finite transport teacher pairs。F0的坐标反例说明：网络可能只是在补偿坏坐标，而非学习昂贵Maxwell机制；F1又指出真实prior边界。此时训练会有显著no-op或记住反演终点的风险。

未来只有多个独立合格anchor出现任务尺度的合法相容移动后，才比较确定性corrector与学习摊销。优先输出材料度量中的可行法向补偿系数/有限传输初猜；输入包括合法当前状态、方向、目标材料距离、局部metric/边界状态和实际physics context。N3的GT点回归架构不作为默认目标。

teacher不能用GT选方向/停机；必须保留非零进度、所有失败和总RHS。网络摊销的是算法corrector，不是唯一“真实材料演化律”。需要多anchor/方向/距离的teacher与独立transfer验证；本轮没有这些。

若用作posterior proposal，应输出可求正反密度的随机核，或建立可逆体积/augmented-state修正。直接投影的many-to-one corrector不能套原proposal密度。候选发现与posterior采样继续分开；未展示/模型外质量不可由residual softmax补造。

完整成本比较应计teacher、训练、context、refresh、失败及最终full validation，达到同一材料风险/覆盖标准后才谈break-even。当前小dense solver已很快，NN是否有任何摊销空间仍未知。
