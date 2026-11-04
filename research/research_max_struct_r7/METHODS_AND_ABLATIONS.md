# 方法、消融与公平边界

## 主比较的实际实现

|方法|Current|Material|最终原物理校正|
|---|---|---|---|
|CIE|继承sourcewise PR；同β|普通共享实材料readout|原L-BFGS-B最多3000accepted|
|Dual-CIE|继承sourcewise PR；currentβ=1|materialβ=1→.2|同上|
|Global-Dual|global PR；currentβ=1|普通共享实readout|同上|
|A-material|同Global-Dual|真实四cell one-hop，完整dP梯度，原L-BFGS readout|同上|
|A+|同Global-Dual|Möbius状态二次驻点候选+原prior-MM接受；最多2sweeps|同上|
|physical L-BFGS|不使用辅助current|原Maxwell解析adjoint直接优化576材料|最多6000accepted总额度|

共同材料空间576独立实数，不使用GT的mask/family/fonts/support。初始m=.1；重复初值实验m=1。box、prior和原Phi相同。前五均mf3/6/10各500+mf12为1500；所有源一起保留。FFT实现完整action/adjoint与零填充Green，接收和最终物理LU仍dense。最后mf12有576native系数，但DPCS affine约束使独立APCS current秩少于576；不能把它称无约束Full-J。实际receiver_rank和state/source计数见各次stage/receipt。

主R7 FFT循环只扩展R6支持的arm与readout callback，没有改PR、β、DPCS定义、Fouriertransport或原物理corrector。成熟sourcewise控制的数值路径受FFT roundoff影响，不能承诺与旧dense逐位相同；这不是新算法收益。A vsGlobal-Dual对材料机制的因果对照更干净；A vsDual还混有current的global/sourcewise区别。

所有候选在线只用原Phi：初始host、真实stage端点及原corrector的已支付incumbent。辅助更低loss不自动当原物理成功。停止原因、约束、数据残差与材料指标分开。

## A+实际开发链

四个披露对象7b8d/b20d、开缝壳40ac、正常近双岛00f1，各5版本完整运行，共20run：原A、精确quintic候选MM、projected8、继承gtol的quintic-MM、最后arc-candidate A+。所有版本保存，没有在新30上重选参数。开发端没有GT和新30结果输入。

Möbius残差坐标仿射是精确代数；最终arc-candidate只对**状态项**取精确驻点，再评价完整prior-MM。它不是正则化条件全局最小，不是逐cell除法，更不是全反演线性化。旧点在候选集中保证理论MM非增；实际浮点loss检查及原seed/old回退保留。容差来自原material_inner_gtol，不是新增科学成功门。quadratic求根确实替代quintic，根数与时间都记录；近似目标更易下降不代表更准确材料。

开发20次的固定电流loss上升记录、公开Phi/CPU/LU/RHS见delegated/aplus/FULL_TRAJECTORY_COMPARISON.md，最终配置见configs/aplus_final.json。projected8丢失b20d的负结果不隐藏。最终A+的新30结果和旧九数字材料评分均在全部主预测封存后计算。

## 有理搜索

同一个公开失败锚点比较：原anchor继续corrector、rank2有限cell对、rank8两四cell块、同动作族从零点局部精确优化。一次结构修正到rank32的两幅度跨区动作，保留gradient/translation/弱物理Hessian支持及随机对照。实际rank不等于幅度参数数量。

Exact finite evaluator保留完整背景多次散射；候选来自材料状态、数据残差、源/收可见性和合法box，无GT。原prior一项不改，端点原Maxwell验证和同corrector再比较。有限grid+local不是全域全局求解证明。最低原Phi不是最低材料误差。局部原点/随机是动作控制，未执行独立全盒multi-start/完整GN/HVP搜索比赛，故没有证据声称超过一切全局化方法。

## 不采用冒名基线

NIE原文与最接近CIE/FFT-TSOM已有方法完成文献审查，但没有本轮独立成熟复现，不把self-only/A叫NIE。LUCIE仅primary abstract已确认，不冒充已运行baseline。三维ordinary/dual-like明确是vector-adapted joint controls，与成熟二维CIE分开；实测只运行Dual/global/A，未运行A+的实测版本。NN/MLP/LNO、LSL均未训练/实现，不能列性能零或理论失败。

## 成本比较

equal iteration/CPU/LU/RHS是四个不同问题。主报告同时给完整endpoint费用和逐对象成对共同CPU/LU/RHS的已支付checkpoint结果。CPU包含辅助工作；LU/RHS不包含大量材料求根，不能代替实际CPU。所有shared-precompute、拒绝候选、失败、数据生成、开发和独立验证另收费；已冻结R1–R6引用不是新增运行。串行rational预算不作独立matchedCPU。CPU时钟受本地调度/频率和路径长度影响，操作数与配对分布同时给出。
