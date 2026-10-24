# 几何有限距离门禁：保留全部尝试的离线复核

按父线程明确指定的 abs(actual RMS/target−1)≤.02 做距离适用门禁；这不是把原算法/教师的1%条件改成2%，原始 radius_matched 与 teacher_eligible 全部原样保存。统计仅原前8DEV；四Euler已知自身场缺陷依旧不作为交换子理论反证。

|研究/类|计划尝试|完成|2%距离合格|合格率|合格且noiseball|合格且DeltaPhi≤2|合格且进展≥.5r|所有支持门禁合格|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|dynamic_geometry_first/Euler|288|283|4|1.389%|1|0|1|0|
|dynamic_geometry_first/direct|288|288|288|100.000%|207|95|287|95|
|dynamic_geometry_first/dynamic|144|144|144|100.000%|137|83|144|83|
|dynamic_geometry_first/frozen|144|144|144|100.000%|113|59|144|59|
|dynamic_geometry_heun_r1/Heun|288|287|9|3.125%|7|3|6|2|
|dynamic_geometry_heun_r1/direct|288|288|288|100.000%|207|94|287|94|

未完成/部分/跳过都保留在原计划分母。分维数、秩、半径的完整分数见 GATE_FRACTIONS_BY_RANK_RADIUS.csv；重复配置不是新增独立父对象。noiseball直接用原几何白化半径，合同Tau另列；DeltaPhi≤2、进展≥.5r、prior_delta≤2、残差≤1e-10来自原Heun教师阈值。对首轮没有通用教师资格规则的类仅作同阈值描述，不追认原先已通过教师门禁。采样路径和所有付费Heun predictor的noiseball标志分开，均无连续证书。

先保存不读真值的全部逐尝试门禁及固定直接控制配对资格，再读取原隔离离线GT评分。ELIGIBLE_ONLY_PAIRED_GT_SCORES.csv只含双方均满足2%距离、或更严格全部支持条件的同父对象配对。选择不使用GT，缺失/不合格分母仍在计划表；这是已暴露开发数据的事后条件描述，不能推广成盲验证。

原离线表中的Heun均值.0229/.0231与直接法.0200不构成同材料距离比较；严禁用那些不匹配平均GT误差排名算法。未满足距离门禁本身就是该算法在所请求半径的失败证据，不能通过离线缩放/挑近点救回。原始GT评分和失败结果都保留，没有新物理或rescale。
