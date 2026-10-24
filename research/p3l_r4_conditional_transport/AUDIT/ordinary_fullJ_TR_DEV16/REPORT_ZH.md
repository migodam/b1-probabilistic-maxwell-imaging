# 新Full-J信赖域普通控制：DEV16封账事实

预注册同一原full66 probit解码与数据白化，两臂lambda=0的Phi和lambda=1的Phi+.5||x||²。每臂每父对象16新状态/256RHS，所有accepted更新都新付J32，每个Full候选8RHS，初始J复用旧付费context。没有mu、GT、模型/原代表元改动、超参搜索或L1运行。

|目标lambda|RHS|Full状态|接受更新|拒绝Full试探|Phi变化中位|目标变化中位|初始GN decrement中位|末端GN decrement中位|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|0|3808|116|90|26|-1.52137418|-1.52137418|5.42526159|4.92419501|
|1|3832|183|74|109|-3.83987447e-05|-1.54601803e-06|0.00180674825|6.19221959e-05|

合计7640/8192RHS，原Full8部分2392、accepted新J32部分5248，全部计费。受控5.956秒，峰值RSS228786176字节。原cached锚点Full+J冷取得40RHS、96+64上下文优化拥有另披露，不当作免费冷启动。

保存每次原Full current/field/prediction、accepted的完整Jx/梯度、实际与预测下降、信赖比、半径和TR子问题KKT。只读重放验证单位、梯度、接受条件和费用全部一致。GN decrement是gᵀ(JᵀJ+lambdaI)^daggerg平方根，曲率为Gauss-Newton而非真实非线性Hessian；截断rank与被舍弃梯度分量也披露，不追认全局最优或后验资格。

是否“MAP已经驻定”应按逐父对象原/末端精确梯度与decrement量级判断，不用目标微小下降代替数学证书。旧unscaled SLSQP第一步与callback-only返回问题另见ordinary_slsqp_diagnosis；best-query DEV附加候选及原320个L1case都保持不变。此新控制只能作为预GT诊断补充，任何L1强控制需父线程先审阅并另登记，程序不自动运行。
