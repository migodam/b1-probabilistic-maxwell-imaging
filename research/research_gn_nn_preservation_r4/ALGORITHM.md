# 可运行算法与接口

输入固定G/E0/S、匿名复观测Y、公开噪声等级和材料box；图像类别、字体、支撑、幅度及GT从不作为网络输入。主模型为576个实材料的完整scalar-TM VIE。

## 共同GN输入

`src/anchors.py`复用R3解析Full-J及8→12→24材料参数化，每阶段仍解576阶物理系统。每个主锚点45秒/500LU/30000复RHS，停止容差1e-12，无数据合格门提前结束。保存真实电流j0、总场Q0、材料m0、公开观测Y和背景Jacobian回投影gb。数值收敛与预算停止分列，full24只表明576材料已释放。

discrepancy及weakreg1e-5为独立后缀基线，不可替换NN的严格GN输入。discrepancy门为clean1e-6/noisy.0011，并在完整24阶段检查。所有失败与重新准备费用保留。

## 六个两级共享权重Cascade

27个输入通道：[m,gb,ReQ的12源,ImQ的12源,noise]，TRAIN公共锚点拟合归一化。编码通道16/32/64，24→12→6→12→24。新模型均从零输出头初始化训练，不使用R3权重。两级只更新材料输入和累计辅助电流；gb、Q及noise保持原锚点，不增加未计费物理刷新。

voxel头输出576维材料修正；current头输出576×12复修正，分别应用完整、接收零空间、G弱、联合暗、软暗投影。所有current头参数相同，但允许的实电流自由度不同，不能与576材料自由度混称。每级通过共享实材料row-LS读出并保持box。

`train.train_model`的alpha、loss_name及seed明确进入checkpoint；选择指标统一为VAL材料L1+0.1边缘L1，不使用TEST。AdamW lr.001、batch64、100轮上限、15轮无改善停止。stage1/stage2材料及辅助电流随预测存档。

## 物理后处理与在线排序

```text
public operators,Y,noise -> strict GN anchor -> NN full material proposal
保留anchor
分别取0.5和1倍完整材料提案，并在box内形成候选
完整VIE计算两份raw材料；不把辅助电流残差当物理残差
每份必要时用有界解析Full-J TRF校正，最多3秒/8次函数评估
总增量最多10秒/100LU/6000复RHS，包含特征/推理/试探/验收
从物理合格的非anchor候选中选距离完整NN提案最近者
无合格提案时回退合格anchor；再无则FAIL
先保存候选与在线排名，再离线读GT
```

clean/noisy物理门为.001/.003。另保存相同候选的旧data+TV排序、全部best@K和数据合格best@K，后两者只可离线计算。质量联合成功和严格预算内联合成功分别报告。

## 场约束消融

先比较材料监督，再在J-full及VAL选定的受约束方案上分别增加.01倍本构缺陷、.01倍内部场变化。两项分母固定于每个公开GN锚点。可选精确完整VIE微调只有达到既定VAL提升门及剩余时间门才启动，并有相同起点/更新次数的纯材料微调对照。没有运行的模块不能写成实验结果。

## 真实误差能量

历史入口`experiments/e0_energy.py`分析已披露R3输出。新TEST入口`experiments/test_energy.py`在封存后重算raw材料和top1的真实Maxwell电流。固定Green强/弱、接收强/弱/核、联合暗/明、顺序Two-fold及数值交集分别保存；不把auxiliary head代替物理电流，不把重叠投影强行加成100%。

## 传统参考与成本

严格GN、继续TRF、弱正则、discrepancy与成熟结构CIE/Fourier-TSOM分别列出。CIE复用冻结原实现的3→6→10 Fourier掩模及单/双beta结构，采用全576实材料和完整原物理校正；显式Fourier列不等于已实现FFT复杂度。其冷55秒口径与NN的GN45秒+增量10秒并列，实际LU/RHS和当前自由度另报，不以相同时间上限冒充同等实际计算量。

本轮无概率后验、全根枚举、全局唯一性、GT必达或transverse全局收敛保证。
