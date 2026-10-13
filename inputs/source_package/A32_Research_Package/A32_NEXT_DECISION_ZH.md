# A32 下一阶段研究裁决

## 主线

**有限尺度拓扑分层材料后验 + Maxwell 模型修正的双侧材料似然比签名。**

目标不是更漂亮的多图生成，而是：给定测量、已声明物理先验和计算精度，计算有孔/无孔等材料事件的后验概率，区分观测歧义、先验敏感与数值未决。

## 最重要的已推导连接

在共同实白化 residual、共同 proper prior 与冻结材料坐标下，若

`J1 - J0 = Y Gamma^{-1} N`，

则两后验的精确 density ratio 只依赖

`U(s) = (N s, Y^T J0 s)`。

U 给定后的材料条件分布不变。该结论允许非高斯/混合/固定物理截断先验；Gaussian 情形恢复 A12/A17 的 paired material support。拓扑概率是否改变，取决于拓扑事件与 U 及重加权的条件相关性，而不仅是 forward error。

一般条件重加权不是新数学；真正需要建立的是 Maxwell-specific 结构、跨拓扑非线性余项及可支付的任务级误差控制。

## 第一阶段：参考概率优先

新建3D同心分层介质球 vector-Mie 小问题，两拓扑分支、每分支初始仅两个连续材料变量。先核对完整矢量散射，再用独立积分细化核对 evidence。明确最小孔洞尺度和 proper topology/material prior。

可预注册研究容差：不同积分细化的 P(H) 差不超过0.02、log Bayes factor差不超过0.1；这不是未经证明的严格误差上界。建议先批准最多384个新增材料–频率评估，包含重复、失败与控制。不能建立参考即标 REFERENCE_UNRESOLVED，不扩张 sampler 比较。

## 第二阶段：最小机制检验

对同一参考后验比较 full F、single-anchor calibrated L 与一个冻结近似。测量真实 likelihood defect、双侧签名预测、非线性剩余、事件–签名耦合及实际概率变化。

若签名不能解释有限拓扑变化，或其估计成本接近完整 posterior，就只保留局部解释，不声称部署加速。若 full F 对主要事件概率没有实质作用，就允许采用 calibrated L。

## 第一版算法

无大型网络：分支 SMC/MCMC + 真实 likelihood；material metric 只用于 proposal/preconditioning；弱可观测细节保持随机尾部。标准 delayed acceptance 是必须比较的已知基线，而不是待申报的新颖性。

## 高风险备选

Topology-conditioned nonlinear likelihood-ratio transport，研究保留条件材料分布的低维重加权及条件 evidence 误差。只在可信参考与明确结构后考虑小规模 learned proposal。

## 停止/不启动

不扩展 current-selection/OPM；不将 Sparse-64 oracle 当 posterior 支持保证；不把 CNN candidate softmax 当校准概率；不把 Poynting 恒等式作为独立数据；不重复全 box uniform Neumann certificate 主流程；不在参考后验未建立前训练大型 flow/diffusion。

## 论文成立所需的组合

精确的材料/拓扑事件 + 非平凡 Maxwell 机制定理 + 保留支持的算法 + 三维矢量概率参考 + 同概率精度的全部成本比较。现在是 theory/pilot GO；优越性、复杂三维校准与认证加速尚未建立。