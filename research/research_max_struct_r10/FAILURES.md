# 失败、反例及未覆盖范围

本文件将实现故障与科学反例分开。所有数值失败、部分状态和费用保留；修复入口错误不等于修改科学算法以追求成功。

## 已发生的工程故障

|事件|原因及处理|费用与保留位置|
|---|---|---|
|历史energy第一轮dtype错误、确认形状错误|诊断读取类型/形状不兼容，修正分析器后重做必要分析，未改原反演|delegated/energy_audit，全部attempt计费|
|历史vector energy序列化问题及自身中止|数组/长诊断序列化，保留已产生文件及ps观察计时；未知物理计数不补零|energy总账，至少23.24 CPU秒的中止窗口|
|固定状态40ac缺独立operator文件|原PUBLIC_TASK依赖未一并读取；保留前16次已付readout，按真实公开依赖修复|RESULTS/denominator_fixed_state_failed_missing_operator_01|
|谱诊断parse错误|日志/字段解析，修复后一次完成24状态；失败CPU未测|delegated/denominator_spectrum/FAILED_ATTEMPTS.json|
|onehop优化器控制错误传参|R8 arm_config覆盖顶层8步请求为60；中止自身两进程，通过正确的aplus_options通道修复；正式实际schedule为8|RESULTS/onehop_solver_control/FAILED_ATTEMPTS.json，119.202706 CPU秒/952 LU/22740复RHS|
|首次三维CPU探针回执组装失败|已有forward后字典字段重复，保留1 LU/60 RHS；完整失败CPU未知|vector CPU探针attempt记录|
|三维Windows包漏已有R9依赖|缺vector_extension.py，首个v04生成在denseG组装后、0LU/0RHS时退出；只补既有源文件，归档失败后执行原合同|remote_v04_transport_attempt1_missing_dependency及GetProcessTimes记录|
|初次CUDA等待外部视频任务|CUDA context创建前等待并停止自身等待进程，转同物理CPU；用户通知GPU释放后重新使用尚未启动的完整cohort|等待778.84秒、主进程4.03125 CPU秒及查询2.0 CPU秒单列；不是算法失败|
|ScholarQA/出版商访问429|保留失败检索，用已有全文和可访问的学术索引/原始资料补充，不高频重试|audit/literature，工具费用未知|

失败计时窗口、正式运行与嵌套明细见[已完成DEV成本审查](delegated/cost_audit/REVIEW.md)。未记录费用是unknown，不是0。上述中止只涉及本任务自身进程，没有停止其他用户任务。

## 已建立的科学反例和限制

- 固定状态辅助目标下降不保证原Maxwell Phi下降；完整corrector的Phi下降也不保证材料更准。
- 原A+在b20d有实质收益，在40ac仍失败；FFT oracle表示余量很小不能消除这个差别。
- 仅将旧one-hop换为8步L-BFGS没有复制新矩阵分母的两个DEV结果；这一控制不等于否定所有成熟A实现。
- Green19、等秩Green28和Two-fold28在两开发对象相近，尚无Two-fold接收信息独有收益的证据。
- 保留原Gii的控制仍成功，不支持将此前失败简单归为自项写错。
- 正介电分式归一化不保证目标全局凸；一维反例见理论第5节。没有全局basin或唯一GT恢复定理。
- 公开G压缩耗散最小浮点特征值近舍入尺度，不能据其符号发严格浮点证书。精确解析盘脉冲模型的可逆性定理与实际矩阵条件数、浮点误差、原Maxwell成像质量分开。
- β=.2时简单全box实部coercivity充分条件失败，不代表实际状态奇异；未据此改变冻结β或停止确认。
- 2D新集为同有限离散生成/反演；3D跨18→12虽不同网格，仍共享CM+RR物理近似。二者都不是连续Maxwell精度证书。
- 新低秩矩阵分母未做三维移植、复杂介电材料或实测验证。三维本轮比较的是继承原A+与CIE，不将两条方法的结果混合。
- 本轮未新增NIE、NN训练或Fresnel实测。原R9已有经典/实测/原尺寸三维负证据继续有效。

## 最终确认中的负证据

- 168条二维主反演均完成，但Green19平均误差/CPU没有优于CIE；高对比环获救与桥损害同时存在。开缝壳A类低梯度坏拟合不能统归预算不足。
- 三例PGP→GP修正仍只保住环，桥材料更差；较佳低秩Green逼近不足以保证较佳成像。三例701.64 CPU秒全部计入。
- 新二维Two-fold实际秩30，历史DEV秩28。早期文档误沿用28；修正标签和图，保留原标签图及重绘费用，没有改实际basis、阈值或反演结果。
- 条件确认0001首个worker因输出路径保护拒绝，在任何模型/forward前退出；0.125 CPU秒、0 LU/RHS保留。修复允许路径后运行相同配置的首次真正求解，不是第二个科学初值或按GT重跑。
- 三维原18条及续跑3条都到各自80步预算；续跑证明部分预算不足，但未使A+与CIE形成明显可重复差距，不能称全部已收敛。
- 两个三维跨网格clean端口差超过声明3%噪声；matched包覆/空心体则存在GT与估计同在实际噪声球的材料差异。模型误差和有限噪声歧义分别报告。
- 最初3D FFT/GMRES后端确实过慢。公开m0成本比较后转兼容CUDA retained-LU，原慢运行、停止及完整v08旧终点不删除；4.41倍是特定任务工程收益，不是A+算法增量。

科学数值队列均已闭合。原生辅助系数停止、物理数值停止、预算上限和工程启动错误分别记录，运行结束不等于科学主张成立。
