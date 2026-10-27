# G1：有限prior概率实现与未运行项目

正式G1未运行。G0还没有新物理资格，且没有G1资源授权。没有计算新的2950原子posterior、没有方法成功率、没有PDE加速数字。

## 已实现并测试的有限概率语义

`src/bgga/intervals.py`不接受reference、truth或future data，且没有文件读取函数。未知原子的likelihood界为[0,1]，保留prior所有支撑；已查询近似误差不会自动归零。forward/potential经验误差标为EMPIRICAL_INTERVAL，不会因归一化被提升为数学PDE证书。

归一化用事件和补集分式上下界。log域保存微小likelihood和prior尾部；极小值不能因为exp下溢就被当作无质量。查询规则是预先声明的`prior × likelihood区间宽度 / cost`启发式，不使用full reference risk，也不声称最优风险策略。

同一观测的模型版本更改只替换旧factor；从base prior和当前factor重建以恢复旧零支撑，拒绝同一observation ID以不同factor重复加入。完整joint pushforward可保留条件核；使用prior条件核的canonical lift默认NOT_CERTIFIED，不能自己把joint TV写成0。

任务聚类明确区分`max_probe_l2`和预加权Euclidean任务。前者取每个probe白化范数的最大值，不把多probe拼接L2冒充max。计算簇完整diameter，拒绝仅靠single-link链式合并的推论。浮点directed rounding与统一Maxwell稳定常数都未认证。

## 原始数据已定位，未误用

2950原子定义、全波means、24个finite-prior含噪case均存在；它们与连续A0另24case、旧W另24target probe不是同一东西。见[缓存与来源审计](research/delegated/source_audit/SUMMARY.md)。已有全表可供逻辑开发，但重放身份必须为REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE。

M0–M6正式方法矩阵均NOT_RUN；目前只交付概率核，不冒称confidence/angle/task/multi-Green/quotient方法均已完成。未来所有reference由独立评分器持有；在线query只能访问prior、data和自有cache。任何读取reference的依赖必须进入部署费用。

## 是否存在明确节省机制

完整118背景×2,257ports的eager准备为266,326 RHS，已经远超G1建议20,000。全表direct是11,800 RHS且有旧9288 CPU秒记录；旧数值只作量级预测。必须先找到对同一任务有效、普通ROM无法免费获得的按需节省，再申请G2。保守[0,1]区间若要求接近全表查询，结论应为CERTIFICATE_COST_NO_GO，不能删除低prior原子凑成功。
