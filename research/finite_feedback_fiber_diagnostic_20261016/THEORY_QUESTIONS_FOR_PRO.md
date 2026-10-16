# 下一次Pro理论研究应回答什么

先读本次实测数据，再讨论新理论。不要把经典resolvent或常秩结论当作原创贡献。

## 可由本轮约束的问题

1. **先验内有没有精确局部fiber？** 对合法商chart g去掉同材料gauge后，若DF Dg在某邻域满列秩，精确水平集局部为离散点；不能把全pixel复空间的巨大nullity搬过来。本轮两份7维材料tangent数值满秩，只支持有限噪声相容tube。需要声明噪声tube半径、材料度量和均匀余项，不能从一个小singular value推断有限族。
2. **二阶可延续的障碍在哪里？** 等数据曲线需A v=0，并存在a使A a+D²(F∘g)[v,v]=0。后二阶必须包括decoder曲率，且其在Range(A)外的部分不能由normal修正消除。若只近似Av≈0，应给出先验合法步长下的一阶漂移、未补偿曲率和高阶余项共同界。该经典必要条件不是新定理，本轮没有计算其uniform bound或宣称精确fiber存在。
3. **哪种条件能把“GT弱能量高”连接到可恢复性？** 必须先有独立physics-only、原先验内且观测合格的anchor，再度量GT差量是否位于同prior chart的finite可达区域。本次R基线不满足理想前提，N3 off-prior，二者不能验证该假设。需要比oracle能量比例更可观测的敏感分量恢复诊断。
4. **何时学习tangent field，何时学习finite transport？** 当Dg曲率远大于Maxwell反馈曲率，训练tangent vector本身可能无助有限步。先检验material-native retraction/商chart和合法decoder；若几何仍弯曲，才论证有监督或无监督finite transport的收益与物理成本。
5. **A17能增加什么？** 原A17是固定二次材料任务下current representation fidelity，不负责产生物理解。下一步若使用rational Schur加速，应首先比较完整finite响应精度和总成本，而非再次只改善frozen GN gap；当前小型dense solver已稳定廉价，未建立A17独立增量。
6. **分离分支与高阶反馈是什么关系？** 本次resolvent预测精确但没有新增结构，不能把高阶抵消视为跨分支机制。要先有同一冻结prior/观测的合格竞争材料区域，再问全局proposal或posterior kernel。否则先建立ambiguity benchmark或active sensing问题。

## 现阶段停止与保留

保留：finite反馈数值oracle、正确prior/质量度量、decoder-vs-Maxwell曲率分解、原N3 Top-1基线。

暂缓：新NN/Flow/Diffusion、非径向3D、多解成功叙事、通GT路径保证、A17新增成像收益。不得用剩余额度追加正结果。

决定性下一实验应是：从背景独立反演获得同原先验内合格anchor，验证其敏感分量，然后在同material距离/同full成本下比较普通与material-native曲率延续。须新授权；本轮不执行。
