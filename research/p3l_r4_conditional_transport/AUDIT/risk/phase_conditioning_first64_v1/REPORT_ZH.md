# 最终相矩诊断：有限梯度不等于易优化

固定最初64个TRAIN父样本，不拟合、不读取DEV/locked标签，不新增物理。原512/128训练context为48+48，与实际revision2六fits一致。比较uniform exact-bin与eps=1e−4 anchor-mass初始化，保持同一精确箱标签。

对每个阈值箱j及质量logit k，∂M_j/∂l_k=B p_k(1{k>j}−tail_j)。当anchor平台给出近one-hot质量时，多数质量p_k为O(eps)，tail_j接近0或1，故有限但很小的导数可能限制更新。该机制是softmax条件数现象，不是梯度NaN；也不能仅由它证明所有实际风险差来自初始化。

|K|anchor相区|联合moment-Jacobian幅度比|积分材料头导数幅度比|实际箱loss logit梯度幅度比|
|---|---|---|---|---|
|8|host|0.000424228|0.000499957|1.81373e-06|
|8|transition|0.808304|0.711894|0.349139|
|8|pure|0.872311|1.16308|0.107405|
|16|host|0.000583042|0.000699935|2.66551e-06|
|16|transition|1.02917|0.958707|0.438875|
|16|pure|1.22809|1.68493|0.155217|

比值先按每parent非空相区取均值，再在parent上平均，未把像素当独立试验。第一项是完整K×(K+1)Jacobian Frobenius，第二项是1/K积分后材料头导数L2，第三项为同一64TRAIN batch真实加权箱MSE对logits的梯度L2。保留逐parent尺度和空相区，不用GT调eps。

结合父调度器已观察的revision2验证失败，停止phase家族扩展。此诊断不授权新初始化、学习率或额外训练救结果；finite-gradient PASS只证明导数数值有限，不能证明学习良态。
