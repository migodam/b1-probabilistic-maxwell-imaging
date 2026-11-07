# 可运行算法与接口

输入固定G/E0/S、匿名复观测Y、公开噪声等级和材料box；图像类别、字体、支撑、幅度及GT从不作为网络输入。主模型为576个实材料的完整scalar-TM VIE。

## 共同GN输入

`src/anchors.py`复用R3解析Full-J及8→12→24材料参数化，每阶段仍解576阶物理系统。每个主锚点45秒/500LU/30000复RHS，停止容差1e-12，无数据合格门提前结束。保存真实电流j0、总场Q0、材料m0、公开观测Y和背景Jacobian回投影gb。数值收敛与预算停止分列，full24只表明576材料已释放。

discrepancy及weakreg1e-5为独立后缀基线，不可替换NN的严格GN输入。旧 discrepancy 在完整24阶段检查 clean1e-6/noisy.0011；新增 `discrepancy_allscales` 在8、12、24每阶段检查同一门，避免粗层先拟合噪声。达到粗层门仍继续 prolong 到下一层，最终仍使用576材料和完整物理。该普通停止策略修正先于 TEST 评分冻结，旧计算及所有失败费用保留。

## 六个两级共享权重Cascade

27个输入通道：[m,gb,ReQ的12源,ImQ的12源,noise]，TRAIN公共锚点拟合归一化。编码通道16/32/64，24→12→6→12→24。新模型均从零输出头初始化训练，不使用R3权重。两级只更新材料输入和累计辅助电流；gb、Q及noise保持原锚点，不增加未计费物理刷新。

voxel头输出576维材料修正；current头输出576×12复修正，分别应用完整、接收零空间、G弱、联合暗、软暗投影。所有current头参数相同，但允许的实电流自由度不同，不能与576材料自由度混称。每级通过共享实材料row-LS读出并保持box。

soft alpha大于零时允许电流仍为完整13824实维，alpha是方向缩放而非硬幅度限制。没有独立明方向范数约束时，不预设它能够保护真实强方向；必须测量实际aux/raw/final更新。

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

## 追加的共享权重反馈补充

[recurrent_feedback.py](src/recurrent_feedback.py) 在电流补偿和材料读出之间加入固定G/S及伴随的物理下降层，不计算完整L逆。[理论及条件](theory/A17_RECURRENT_FEEDBACK.md) 区分固定材料二次下降、读出后的剖面下降与整个NN映射的稳定性。原六臂保持固定锚点Q输入；补充的refreshed-Q、二次feedback1与噪声容许feedback_hinge使用当前辅助Q，网络容量和两次权重共享相同，rho固定1。原J-full复用为无反馈控制，不重复训练同一个模型。三种补充均已从头训练并完成VAL；独立TEST、16步深度诊断以最终运行回执为准。GN初始化与最终验收仍求解完整 Maxwell，不宣称全管线无LU。
