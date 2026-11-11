# 已封存basin候选池：online与离线oracle

本项为root明确授权的封存后GT诊断，0新forward、0优化。oracle只描述已付费、已保存、bodyCPU≤240的有限池，不是预测、early-stop或新方法。材料/归一化/IoU完全复用原score_basin.py，online仍取最小原Phi。

|对象|方法|savedK / eligibleK|online NRMSE|offline best@K NRMSE|oracle候选CPU|oracle数据残差|
|---|---|---:|---:|---:|---:|---:|
|r10_twod_0001|direct|6/6|0.422880|0.422880|0.188|0.136844|
|r10_twod_0001|multistart|3959/3956|0.748880|0.374538|99.094|0.148631|
|r10_twod_0001|alternating_projection|174/174|0.422880|0.422880|163.344|0.136844|
|r10_twod_0001|douglas_rachford|1779/1404|0.165853|0.165853|239.969|0.00960351|
|r10_twod_0004|direct|6/6|0.515137|0.515137|0.062|0.589798|
|r10_twod_0004|multistart|3791/3791|0.456326|0.452016|209.734|0.0611405|
|r10_twod_0004|alternating_projection|173/173|0.515137|0.514830|1.703|0.617976|
|r10_twod_0004|douglas_rachford|1806/1447|0.113775|0.113775|239.984|0.0096873|
|r10_twod_0022|direct|4203/4200|0.436727|0.336028|1.297|0.168249|
|r10_twod_0022|multistart|4214/4214|0.436730|0.336028|1.266|0.168249|
|r10_twod_0022|alternating_projection|1776/1438|0.351438|0.334774|1.875|0.204332|
|r10_twod_0022|douglas_rachford|925/925|0.106175|0.106175|208.656|0.00974641|

|对象|控制|控制pool oracle|DR online|控制pool是否含比DR online材料误差更小的已保存点|
|---|---|---:|---:|---|
|r10_twod_0001|direct|0.422880|0.165853|False|
|r10_twod_0001|multistart|0.374538|0.165853|False|
|r10_twod_0001|alternating_projection|0.422880|0.165853|False|
|r10_twod_0004|direct|0.515137|0.113775|False|
|r10_twod_0004|multistart|0.452016|0.113775|False|
|r10_twod_0004|alternating_projection|0.514830|0.113775|False|
|r10_twod_0022|direct|0.336028|0.106175|False|
|r10_twod_0022|multistart|0.336028|0.106175|False|
|r10_twod_0022|alternating_projection|0.334774|0.106175|False|

所有12臂都核查候选行号与history一一对应，online材料与已封存240前缀逐元素一致；任何缺失/错位会列于coverage而不补跑。K含重复已付callback，不是独立起点数。raw未审计inner/shadow没有恢复，不能据有限池断言方法从未到过某个区域。

30/60/120/240 bodyCPU前缀和firstK视图分别保存在SAVED_POOL_CPU_PREFIXES及SAVED_POOL_BEST_AT_K。两者均按公开Phi产生online结果，oracle只在离线另列，不与在线胜率混合。

结果应区分：控制的Phi-selected估计较差；控制已保存pool材料oracle是否也较差；DR自己的pool是否含被Phi选择遗漏的更好材料。即使有限池oracle仍差，也不能推出全部未保存迭代或其他restart分布都失败。完整运行超过240的候选均排除，完整费用仍保留。
