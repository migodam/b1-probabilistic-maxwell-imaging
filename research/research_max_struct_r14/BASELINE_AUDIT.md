# R14 实现、强基线与因果解释审计

## 研究实际改变了什么

本轮没有更换 Maxwell forward。二维仍用已冻结的 disk-pulse scalar-TM VIE、576个独立实χ、原非周期FFT Green动作、原生Fourier APCS和噪声选取的receiver DPCS。材料box、平滑/TV prior、初始χ=.1、源权重、原物理验收与min-paid-Φ选择均继承。M1只使用原互不重叠2×2四cell块，保留对角自项分母与非对角一次传播；不是full-G、two-hop或任意稠密反馈。

本轮四组合都使用global PR；同一material列内切换current family只改变current残差变换及其正确伴随。同一current行内切换material family同时改变了实际目标和原有有限readout程序，因此主四组合本身不足以排除“普通读出与Arc程序差别”。额外四对象共享χ坐标L-BFGS60控制专门处理这个混杂，没有用冻结metric近似替代完整导数。

| 本轮标签 | Current | Material | 含义 |
|---|---|---|---|
| 00 | P0(b0) | P0(b0) | stage内一致scalar目标 |
| 01 | P0(b0) | P1(b1) | 原global A+-Arc数学更新 |
| 10 | P1(b1) | P0(b0) | 反向不一致 |
| 11 | P1(b1) | P1(b1) | stage内一致one-hop目标 |

二维b0=1；b1=1/.2/.2/.2。原Dual-CIE的scalar material后期β=.2，因而不等于00。把两个scalar函数都叫“局部metric”不能消除它们的β差异。原A+的P1依赖材料，material full derivative保留；其根处理想线性化中dP乘以零缺陷消失，不允许据此离根漏导数。

## 原生更新和完整状态

`StatefulDualMetric.step`继承R8实PR+、精确二次current线搜索、R6 NativeFourierAPCS以及R7 Möbius/Arc两sweep读出；非局部current接线有独立伴随核对。01在0017的32/501/1001步与原程序m/J/A逐项一致。checkpoint保存原PR历史梯度和方向、Fourier系数、DPCS、源metric、阶段、旧stage系数、J/Q/m及FFT/块缓存；不把相同m称作相同完整状态。

跨stage重入验证覆盖实际stage边界，三对象七checkpoint的干预从同一个完整state复制；retain/reset PR分别记录。checkpoint使用数值NPZ+类型JSON，不使用pickle/object arrays。无GT在线输入。完整状态不是仅用于图片的材料终点。

## 强基线与停止规则

旧48上的CIE、Dual-CIE、G19、DR、A+-MM/A+-Arc保留原已封存主结果，明确标历史参考。新24另实际重跑成熟CIE、Dual-CIE、global-Dual、G19、DR及direct，并与四pair和单次stage切换候选一起冻结；没有用普通GN替代CIE，也没有把00冒充成熟CIE。

复用的C2数值runner保留原文件；统一预算版本只有一处corrector deadline改变，精确差异见`audit/BASELINE_BUDGET_DIFF.txt`。Direct从相同χ=.1直接做原物理L-BFGS，允许完整120进程CPU；不人为限制成72+48分段。它的raw是初值、corrected是物理搜索，不能混入“辅助proposal之后corrector”的平均效果。

原计划与主四组合实际预算的偏差已在主评分和新24预测前记录：proposal至绝对72 CPU秒，之后可使用至总120的剩余CPU。旧C2的corrector最多48秒；因此新四臂与旧结果不能声称停止完全相同。主比较另报strict48已支付候选池保守/宽松截取、新四臂同CPU前缀及componentwise同LU/RHS池。它们不是重跑后少付费，所有实际费用仍记原值。

所有方法的setup、物理审计、固定分解、失败、最终写盘和超时操作均计费；超时操作不能把候选提前放入预算内。候选时间是原进程CPU，外部OS完整child时间另记；未知的逐候选外部时间不伪造。

## 其他控制

- self-only 01/11将原块非对角Goff置零，保留原Gii、β、其他代码；同组有原global-Dual控制，区别self分母、β与空间传播。
- 01→11只在原mf12入口切一次；反向11→01在曝光16作次序对照。它们不是按GT选metric，也没有新的残差阈值。
- 共享LBFGS16的四臂都使用同χ坐标、起点、box、prior、60步上限与容差。其本地算力使全部完成3000步；不能把与旧Windows Arc的差别全部归因readout。
- 三维四pair共用原full analytic L-BFGS12材料程序。另有实际sourcewise CIE/Dual参考。三维不是R7实标量Arc的未经推导移植。

## 本轮没有声称的事情

没有重训NN，没有全局metric导数冻结，没有改变数据或先验来增加恢复率，没有用GT选根、填补缺失或决定corrector是否执行。没有把离线成功并集称为统一算法。没有新Fresnel/Austria反演、连续PDE认证或完整tensor材料反演；已有外部数据和NN证据只作明确标注的历史背景。
