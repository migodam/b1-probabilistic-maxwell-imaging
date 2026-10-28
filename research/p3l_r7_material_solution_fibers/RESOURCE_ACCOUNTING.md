# R7实际费用

总新增37,456 complex RHS /50,000上限；未用12,544。forward10,600、derivative26,856；state/LU各1,325；19个物理job、物理动作失败0。拒绝试探、恢复及修正重放全部计费。

|用途|RHS|
|---|---:|
|首次实际probe|200|
|最终DEV4|7200|
|保留的DEV修正前/重放|2952|
|八份新truth|64|
|八对象anchor、探索、路径|26336|
|普通DEV20态profile|576|
|同材料与交叉copula物理控制|128|
|合计|37456|

物理动作wall14.7676秒；包含持久化/setup的物理job wall72.2100秒、CPU67.4749秒；物理进程峰值189,972,480 bytes。E1另7.18秒与490MB峰值。编辑、文献、绘图、打包和交付不伪称都被物理job钟覆盖。没有训练、GPU或云服务。

本轮fresh anchor是direct-Probit64，目标是几何；不能与R6不同优化流程的1.73秒换算成成像系统加速比。R6的6,860,984 RHS不重计；九次合法T/J复用的cold-equivalent288 derivative RHS是历史cache权利，不是新支付。

中央追加账本在完整Release；snapshot及独立逐行replay在仓库。C的24 RHS setup按四档整体计一次，缺失端点不得构成加速主张。
