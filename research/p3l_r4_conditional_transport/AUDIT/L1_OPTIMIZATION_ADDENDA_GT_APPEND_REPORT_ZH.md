# 第一64个L1：封存优化补充的隔离GT附录

父释放在全部risk/density两cohort、原physical/PPC和两项choices封存之后。此不同PID先核对chronology/全部封存身份，写固定比较计划，才第一次读取本进程heldout真值；未重新选择、缩放、迭代、拟合或查询Full。统计单位为相同64个父对象，两个lambda臂不是128个独立样本。这里仅固定64子集补充比较，主45risk/26density评分由主队列拥有。

|固定动作顺序|full MSE均值|full MAE均值|ROI MSE均值|
|---|---:|---:|---:|
|original_anchor|0.0105139335|0.0417291503|0.0236563504|
|original_ordinary|0.0105139335|0.0417291503|0.0236563504|
|fixed_spline_raw|0.00674821238|0.0288963313|0.0151834779|
|fixed_spline_clipped|0.00673323729|0.0283991488|0.0151497839|
|fixed_r2_raw|0.00775691932|0.0320888566|0.0174530685|
|fixed_r2_clipped|0.00773437422|0.0314370089|0.017402342|
|original_rep_spline_raw|0.00915413526|0.0376018487|0.0205968043|
|original_rep_spline_clipped|0.00916059213|0.0376315999|0.0206113323|
|original_rep_r2_raw|0.00945203517|0.0384642947|0.0212670791|
|original_rep_r2_clipped|0.00945645334|0.0384741116|0.02127702|
|TR_lambda0|0.0163462236|0.0514595662|0.0367790032|
|TR_lambda1|0.0105051517|0.0417330467|0.0236365914|
|best_query_ordinary|0.0105098703|0.0416940922|0.0236472082|
|best_query_spline_raw|0.00915413526|0.0376018487|0.0205968043|
|best_query_spline_clipped|0.00916059213|0.0376315999|0.0206113323|
|best_query_r2_raw|0.00945203517|0.0384642947|0.0212670791|
|best_query_r2_clipped|0.00945645334|0.0384741116|0.02127702|

RESULTS.json给48个固定对照的按父对象配对bootstrap；95%区间仅探索性描述，不是主方法confirmatory或多重检验控制声明。全部64保留，无根据GT过滤对象或选择最优lambda。原callback-only代表元、固定Spline_s1/R2 raw/clipped和原锚点均原样作为参照。

成本拥有：原320查询池80640RHS，新Full-J TR 30752RHS；best-query新增Full0，原source成本不消失。公共64context实际拥有240304RHS，cached初始Full+J冷参考40RHS/parent，均不重计或当作免费冷启动。原raw诊断/PPC费用另保留原receipt。本评分新增Maxwell/LU/RHS/训练为0。

callback major-iteration事件不等于已通过线搜索的下降步；原预算停滞存在选择语义与未缩放方向问题。Full-J TR降低Phi不等于降低GT风险，lambda1小GN decrement也不是全局最优、MMSE或后验均值证书。best-query角色是已经验证的查询候选，不是优化收敛。两项均为方法冻结后/GT前登记的探索性强控制，不改写primary结果；最终科学解释由主线程保留。
