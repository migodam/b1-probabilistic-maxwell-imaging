# 首个材料步的非根梯度分解

固定7个已有checkpoint × 四metric pair × PR保留/重置，共56/56枝完成；失败0。只读step0的m0和首个完成step1的J1/Q1、实际beta/Os、缓存onehop块。没有GT读取、forward、Green/FFT传播、current更新或optimizer。

材料梯度在实际readout开始位置m0、current更新后的固定J1/Q1处计算。g_const冻结P，只对e=J1−mQ1求导；g_metric定义为原完整gradient−g_const−g_prior；P0/P1都计算。先验直接复用原Regularizer的默认构造，与原intervention的MaxwellModel.reg一致。

g均表示梯度分量，下降力是−g。CSV列work_on_actual_dm=dot(g,m1−m0)，负号表示与该负梯度分量同向。记录三项和完整梯度的L2/Linf范数、实际材料步与−g的夹角、各梯度分量间夹角，以及已有paid physical-step gradient angle；零向量角保持null。

样本是7个已曝光完整状态的56个依赖分枝，不是56个独立恢复任务；不使用truth/recovery标签，不拟合成功，不扫尺度、不重写完整导数、不重复机器精度恒等式验证。新增block动作/源权乘积/先验计算全部作为validation收费；历史checkpoint计数不重收。

第一步分解用于描述ideal根线性化之外的非根metric项。是否支持进一步数学解释由主线程判断，不能由本图推出全局恢复、渐近稳定或新链式法则。

本次正文CPU 3.760280s；receipt前整个进程CPU 3.789770s；新LU/RHS/optimizer均0。精确完成后OS子进程CPU另见LAUNCH_RECEIPT.json。

机械概览（无恢复标签）：metric/const L2 ratio范围 0.1002569–2.0673428；g_const、g_metric、g_prior、g_full对实际Δm的正work枝数分别为 {'const': 38, 'metric': 6, 'prior': 44, 'full': 0}。这只是该56枝的有符号方向计数，不是优化目标单调性或恢复预测。

完整OS子进程CPU（含仅重绘既有56行的刻度整理）合计 4.734225s；0失败。源码初版与原图保留SOURCE_AT_ACTIONS.py和FIRST_LAYOUT；最终小图仅稀疏化symlog刻度。

计数口径：material_inner_evaluations=28是P1原完整gradient函数计数，不是optimizer；P0另有28个diagnostic_scalar_material_gradient_evaluations。material_prior_evaluations=28为P1逻辑请求，实际先验evaluate共56次，切勿相加为84。onehop block complex RHS是块作用列数，不是物理LU RHS。详见COUNTER_AND_SOURCE_SCOPE.json。
