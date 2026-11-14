# Historical same16 global-PR Dual-CIE reference

这里只复制原mechanism16评分表和原queue/receipt已经保存的raw/corrected/primary行；没有新评分、参数、候选选择或共同成本截断。global_dual是global-PR Dual-CIE，非physical GlobalDual。A+为原R13 quintic-MM/no-guard。

|method|raw均值NRMSE|corrected/primary均值NRMSE|平均whole-child CPU秒|
|---|---:|---:|---:|
|global_dual|0.1390717147|0.1262009596|103.36914062|
|a_plus|0.1261098080|0.1096624315|120.38769531|

全16保留，包括0045。每个CSV行有原JSON索引、CSV行号、candidate ID、queue job索引和receipt路径。均值是同对象整日程端点描述；实际消费不同，不能宣称等CPU胜负。

## r13_fresh_0045

|method|view|candidate|NRMSE|Phi|source|
|---|---|---:|---:|---:|---|
|a_plus|raw|38|0.0804502436|0.0003899249116|SCORES[1337]|
|a_plus|corrected|899|0.0596554328|0.0001970222005|SCORES[1338]|
|a_plus|primary|899|0.0596554328|0.0001970222005|SCORES[1339]|
|global_dual|raw|26|0.2566261222|0.07128563264|SCORES[1345]|
|global_dual|corrected|872|0.3265423307|0.02315542389|SCORES[1346]|
|global_dual|primary|872|0.3265423307|0.02315542389|SCORES[1347]|

0045的global_dual校正降低Phi同时增加材料误差；原A+从不同raw候选出发后材料误差下降。不同raw路径不能写成同一个corrector起点的比较。

## 同16 direct参考（原缓存primary）

|method|对象数|primary NRMSE均值|平均whole-child CPU秒|
|---|---:|---:|---:|
|direct|16|0.4096614708|55.91992188|
|nn_direct|16|0.2722266400|101.38574219|

DIRECT16_REFERENCE.csv取core48 direct中同一every-third16，匹配原NNdirect16。没有把NN16与uniform48均值相减。NN另带R8训练先验与既有预处理费用；这里仅列原选定终点，不重选共同成本状态。
