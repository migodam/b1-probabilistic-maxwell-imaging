# R9 fresh32可能材料集合：统计方案说明

状态：代码与草案已交root冻结，尚未执行评分、尚未读取任何GT。本文仅解释已经写入统计草案的方法，不改变冻结源码或历史负结果。

## 对象与主要问题

唯一独立评价单位是32个fresh parent。三个方法为FIXED_Q8_RAY、DYNAMIC_Q8_RAY、RANDOM66_FULL。每法有16个非anchor方向预算：四个半径.02/.05/.10/.15、每半径两个antithetic pair；顺序为第一pair依半径+/−，再第二pair依半径+/−。K0/1/4/8/16均指方向预算前缀；失败、重复或未达门的方向消耗原slot，不用下一条成功方向替补。

主资格采用原绝对95%噪声球Phi≤282.8737921622988、原decoder合法、实际材料RMS距目标半径误差≤2%、signed progress≥.5 radius。旧DeltaPhi≤2、DeltaPrior≤10只作secondary诊断，不筛主cloud。每方向在无GT条件下选择合格候选中的最低真实dataPhi。实际材料RMS使用原576格norm/24。

共享anchor仅在噪声球及原decoder合法时进入主qualified cloud；不合格anchor不悄悄混入主集合。anchor baseline始终单独报告；主cloud为空时，部署动作可显式fallback到anchor，全32个parent均进入动作风险。yield排除anchor，避免anchor让qualification成功率恒为1。

## 材料覆盖的严格定义与宽度

固定epsilon=.075为主，其余.03/.05/.10/.15为secondary。对每parent/方法/K，只看已封存的实际合格unique材料，定义C_epsilon为这些材料的完整RMS球的并集。GT命中等价于离线nearest RMS≤epsilon。nearest只评价离线coverage capacity，不能选择任何在线成员、权重、方向、停止步或模型。

不把逐像素min/max框当成材料集合：这样的框可能拼出从未产生的材料。报告pair diameter、mean pair distance、RMS dispersion、离anchor最远距离、epsilon tube外径diameter+2epsilon、球重叠连通分量和greedy packing count，以避免“覆盖率变高只是集合无界变宽”的空洞结论；不估计高维Lebesgue体积，不声称完整连续材料支撑。

unique规则是在冻结在线顺序上greedy保留首个representative，材料RMS差≤1e-6视作近重复。它是有限集合的去重约定，不是精确物理quotient或posterior branch mass。原始noise-qualified visited states另作variableK诊断，不替代固定预算主cloud。

## 四个预注册主对比

1. FIXED_Q8_RAY−DYNAMIC_Q8_RAY：epsilon=.075、K16的paired parent覆盖率差。
2. FIXED_Q8_RAY−RANDOM66_FULL：同一覆盖率差。
3. FIXED_Q8_RAY−DYNAMIC_Q8_RAY：离线best-of-K材料MSE差；只衡量集合capacity。
4. FIXED_Q8_RAY−DYNAMIC_Q8_RAY：包含anchor、setup分配、全部失败/拒绝尝试和I/O的完整方法wall差。

四项使用相同的20000次parent-paired bootstrap索引和98.75%双侧percentile区间（Bonferroni4）。不把方向、seed、源、半径当成独立replicate。n=32的小样本bootstrap仍是近似不确定性估计。yield fixed−dynamic的−.10绝对非劣margin是预声明secondary95%诊断，不能包装为严格小样本非劣证明。其余epsilon/K/radius/任务曲线均标pointwise secondary。

每个覆盖率同时给95%Wilson、95%Clopper–Pearson及单侧95%exact lower bound。“大概率含GT”若被使用，只能指预注册fresh-parent总体的empirical hit probability：给定epsilon和delta（草案delta=.10），只有exact单侧下限≥1−delta时才支持该总体频率说法；它不表示当前某一测量的posterior confidence。

## 空集合、失败与动作

所有32个parent保留在coverage分母：空qualified cloud记hit=false、availability=false；nearest error未定义/无限，明确记录，不以complete-case删除parent、不发明有限惩罚值。若主nearest-MSE比较中出现这类缺失，输出NOT_ESTIMABLE，保留其余覆盖/availability结果。

uniform与exp(−Phi)归一化的likelihood-weighted动作在同一去重qualified cloud上计算；空集合fallback到anchor并标记。它们没有prior/proposal/Jacobian权重，只是可部署approximate ensemble动作，不是posterior mean或posterior samples。每parent完整complex material MSE单独保留，再求总体平均和parent bootstrap区间，不用GT选择哪种动作更好。

材料ROI均值、≥.75像素计数、至少64/256像素达阈值的事件、以及sigmoid((Re material−.75)/.02) soft occupancy分别记录任务覆盖和范围宽度。pixelwise soft范围只作marginal诊断，不构成主材料tube或空间联合概率。

## 费用与时间

所有forward/J/失败/拒绝调用按runner实际counter保留，reduced RHS256不等价为full RHS576。32-parent共享原operator setup的264 full RHS历史费用显式按1/32分配到单方法部署视图，但不重复计为本轮新物理支出；历史整个setup研究2184 RHS另保留来源。setup wall若没有封存可靠记录，完整cold wall仍为未测，不能用组件和冒充。

study视图将每parent anchor/runtime/frame与共同overhead分摊到三个方法一次；standalone部署视图给每法完整anchor/runtime/frame，共同reporting overhead按1/3分配。不同部署视图不可相加当research total。主比较使用同parent/proposal窗口内旋转方法顺序，loader验证window/order metadata。

到首个returned qualified alternative的时间包含先前失败方向和共同acquisition/setup；没有合格替代时在实际消耗预算处right-censor，不能只对成功样本平均时间。未记录的within-trial首次资格时刻不会被捏造。

## 读GT的唯一入口

`src/possibility_statistics.py`不import Maxwell/online engine。`load_sealed_predictions`只读公开PREDICTIONS_SEALED、COMPLETE、ANCHOR、SELECTED、VISITED及原ROI索引；检查全部32 parent、3arm×16slot、public parent/jobs文件在seal内、原资格、实际RMS、无GT选择和同窗口轮换顺序。

私有truth磁盘读取只在`score_sealed_predictions`中。它要求root已批准的冻结统计config、全部32预测seal以及seal之后的单独root评分批准；批准必须绑定statistics id和protocol id。评分开始先写TRUTH_ACCESS_BEGIN，拒绝静默覆写已有offline_statistics；GT接触后不得纠正模型、改选成员或改变统计方法。若评分接口失败，保留失败现场，不能以修复名义重新选择在线模型。

当前代码只做过源码语法解析，没有调用prediction loader、score函数、任何统计实验或GT读取。原R9失败与数值posterior资格限制均原样保留；新目标允许近似ensemble，真实posterior准确率不是这次实验的硬门。
