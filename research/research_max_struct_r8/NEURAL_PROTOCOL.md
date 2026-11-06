# Neural warm start：冻结的数据和比较合同

## 研究对象

继承 R7 的实材料 scalar-TM disk-pulse 有限模型、576 个独立未知量、原 smooth/TV prior、[.02,20] contrast box，以及成熟 3000 步 CIE/FFT-TSOM + 最多 3000 步共同原物理 corrector。训练模型输出材料初值，不替代 full-wave 求解、不缩减未知材料空间、不改变最终 Phi。

1536 个训练对象、192 个验证对象使用不同随机种子及独立几何参数。训练家族均衡：椭圆、矩形组合、多圆盘、圆环、嵌套材料、十字、随机 Fourier 二值区域、光滑多峰材料。contrast 2.5/6/12/18，复 RMS 噪声 .001/.01/.03，6/12 外部源及半/全孔径独立抽取。NN 合法使用训练/验证标签。测试已在原权重和选择规则冻结后生成：24 个同分布对象、12 个未训练家族对象（Austria 变体、梳状、星形、多孔）、6 个 36→24 跨网格对象；不将这些不同总体混成一个普遍成功率。

训练验证可用相同家族，材料对象不重合；未训练家族只用于最终泛化。cross-grid 仍是同离散物理族，不称独立连续 PDE 求解器验证。训练先验会影响初始化，NN 收益不是无先验的信息创造。

## 模型与成本

小型两尺度 CNN，16 通道，监督目标是按每个材料峰值 contrast 归一的 MSE。十个公开物理特征：复 Born 对角 LS、复 DPCS 本构 LS、current/total-field 能量、孔径/源数、空间坐标。DPCS 只由已知 receiver S 构造，rtol=.001。特征提取不求 Maxwell LU；其 SVD、Green/receiver 动作与时间单独收费。训练只按独立验证误差选择 epoch，完整120epoch记录保留。Born-only 同结构控制只去掉 DPCS/场能量四通道，不查看测试 GT。

材料初始化后，默认仍使用成熟 DPCS 初始电流和零 APCS。额外联合状态消融先求该材料完整 Maxwell 电流，再最小二乘投到同一个 DPCS+Fourier 仿射可行空间；同一初值的已支付 forward 被复用，projection SVD 独立计费。不能把 NN 材料对应的物理电流直接当满足 DPCS 的辅助电流。

## 比较及选择

主比较保留 CIE、Dual-CIE、Global-Dual、原 A、A+、sourcewise A+、原物理 L-BFGS-B 的均匀初值。NN完整接入CIE、Dual-CIE、Global-Dual、global A+、sourcewise A+和直接物理；另有NN+global A+合法辅助电流投影，以及DPCS解析初值+直接物理。初始化消融含m=.1、m=1、Born/DPCS解析、固定随机材料、NN、Born-only NN。原750项为630=42×15主臂及预定8对象120项额外消融。所有候选均实际运行后续求解；不因NN初始Phi高而删除。

sourcewise/global补充比较按用户要求在测试GT未读时加入，时间和原清单保留。之后由公开开发投影诊断提出固定`Jd+(I−Pd)jphys(m0)`偏置，另冻结126=42×NN Dual/global A+/sourcewise A+，没有覆盖原750项。8项最小pilot及24项其余验证扩展属于开发；新42的两套预测都完成封存后才统一评分。详细公式、正负结果和选根顺序见THEORY_DERIVATIONS及WARM_DEVELOPMENT_ANALYSIS。

最终物理任务按整个对象分配硬件；七个尚未开始的完整对象在Mac释放资源后由i9迁回，未读取GT或按结果选择。原分配与HARDWARE_REALLOCATION变更都保留，最终主750为Mac405/i9345，偏置126为Mac66/i960。远端队列对迁走任务只写RELOCATED_BEFORE_EXECUTION控制回执，不把调度确认计为一次物理求解；真正结果由Mac产生。

最终候选按原 Phi 的已支付 incumbent 选择，未用 GT 选初值/中间步/根。记录原 Phi、数据残差、全部停止原因、CPU/wall、LU/RHS、FFT/材料读出费用。预算对比使用已支付轨迹，训练/特征/推理计费分别报告单次和多次使用摊销。更低 residual 不代表更低材料误差。

测试预测封存后离线报告幅度归一 NRMSE、预测自身 Otsu IoU、box 合法性和噪声差异。光滑材料的阈值 support 只作为声明的形态诊断，不据此宣称唯一物理 support。主要使用连续配对误差与实际成本，成功阈值沿用 R7（NRMSE≤.15/IoU≥.7）以便比较，不能触发整轮停止或筛除对象。

来源：R7 Möbius/MM 已冻结；2019 CS-Net→SOM/TSOM/CSI 已存在直接先例，本轮不将 NN warm start 概念主张为首创。
