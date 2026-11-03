# R5基准、物理合同和实际运行边界

## 实际完成的二维研究

全部为exp(−iωt)、λ0=1、εb=μb=1、二维scalar-TM disk-pulse Hankel VIE。材料m=εr−1是576个独立实数，box[.02,20]；源和接收位于未知域外。公开算子/观测来自R2及R4，未从GT重建基。

本构j=X(Einc+Gj)，F=Sj；物理对比电流为−iωj（归一化背景）。所有成像使用同Phi=.5||F−Y||²/||Y||²+R(m)，smooth/TV各1e−4、宽度19.5、TV epsilon=.01。该1/||Y||为原冻结目标归一化，不冒充单位噪声白化；噪声discrepancy另按公开σ和复噪声维数评价。复CN噪声的σ为复幅度，实虚分量方差σ²/2；公开99%归一化噪声半径为25.165422196633，不搬用旧clean .001门。

12个源、48个接收，576个复观测。source metric Ω_s=diag(1/sum|Einc|²)用于辅助current/state，和物理数据归一化、噪声σ不同。共享材料跨所有源，不逐source拟合材料。常数m=.1初始化；材料没有字体、支撑或类别自由度缩减。

CIE mf3/6/10分别36/144/400个Fourier/APCS独立复系数每源，12源对应864/3456/9600个实状态系数；DPCS固定于数据，不计作独立待优化材料/状态参数。实际receiver rank17记录于HISTORY；basis来源是背景Green/receiver和预定Fourier模式。材料始终576实DOF。此实现显式形成IDFT列，费用记录FFT=0，不能按FFT后端假计时。

四材料：旧弯条b20d、已披露确认失败42cb/63ab、环100b。主重放各global/one-hop；self/反相仅旧历史阶段态。新7次干预：旧current-only、旧material-only、三对象material-only扩展、旧两次自然阶段切换。全部为公开开发机制诊断，不是7个独立新盲测。

最终候选按相同原Phi选；封存后才读GT评分RMSE/12、Otsu IoU及历史joint合同。历史门只用于端点描述，不决定整轮研究是否结束。保存实际停止原因和投影梯度，不把run completed当scientific pass。

## 九数字：完整历史引用，不冒充R5运行

FBG-R3成熟E3使用公开576自由材料九数字及三个其他对象，已发布结果中dual-CIE九数字1–9全部恢复，总体10/12；ring/curved_bar失败。45秒是历史任务预算，不是CIE能力上限。

[完整九数字single/dual端点表](RESULTS/HISTORICAL_DIGITS_REFERENCE.md)和对应JSON逐例列NRMSE、IoU、residual、终止及实际LU/RHS/wall。不只展示成功图，也不重新评分。R5global/one-hop/material-only/Full-J在九数字上的新运行状态全部NOT_RUN。

不启动九数字大批比较的具体理由：本轮简化未修复两个关键反例；已有九数字基线位于全恢复结果，重复它们不能区分当前失效机制。下一次扩大应先冻结一个能在独立对象上预测作用的材料方向，再设置等原物理预算、真实wall和充分基线续跑。不能把现有方法名称排成空白benchmark表就称已完成。

## 三维：真实矢量小模型，不是三维成像

复用R1 TensorDDA：G=V Goff+(-1/3+ik³V/(6π))I，complex128、全6×6resolvent；exp(−iωt)，固定μ1。两个体素中心(±.02,0,0)，体积.04³，ε=3I；两源沿+y传播，分别x/z极化，严格横向。外部六近场接收位于±.5坐标轴，输出总/散射电场与内部E/j。只一次LU解2RHS。全6mode谱另包括未实际入射激励的y与反对称mode，不能把它们称额外测量。

沿体素连线的纵向与横向Green实部符号相反，使正相残差metric的coherent权重分别大于/小于1。[实际结果](RESULTS/vector_pair/README.md)。这是离散矢量Maxwell结构反例；没有tensor反演、rank-one材料恢复或新增照明比较。

R1独立Mie球ε3,r=.12已有h.04→.03、far vector误差.10468→.03655的有限资格，本轮不重放。该资格不赋予两体素连续精度证书，更不覆盖H(curl)、任意高对比形状或完整三维CIE成像。没有新增三维材料优化，因而没有把未实现兼容的3D mature CIE记为算法失败。

## Fresnel和eispy2d

官方来源、原论文和字段合同见[BENCHMARKS_PUBLIC](BENCHMARKS_PUBLIC.md)。二维TM纯介质Fresnel较接近当前物理类别，但必须取得原文件/header、有效通道mask、真实频率/几何、入射参考、一次性标定及长度单位转换；当前下载和标定未闭环。三维PP/TP饱和零必须mask，不能当物理零；无逐pixel真值时不能制造NRMSE。

eispy2d是公开合成框架，采用与本项目相反的时谐约定；未证实其附带官方Fresnel实测集。本轮没有安装/替换solver，也未用新合成图冒充实测结果。

实际实测反演：NOT_RUN。其原因是物理/数据合同未完成，且本轮没有足够可迁移算法证据值得为此增加适配；不是任意数值gate。下一轮需要先解决明确数据/标定问题，再研究方法收益。

