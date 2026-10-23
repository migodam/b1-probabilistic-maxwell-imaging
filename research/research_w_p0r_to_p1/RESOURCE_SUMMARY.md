# 最终资源与封账

唯一计费依据为追加式PHYSICS_LEDGER.jsonl，逐job RECEIPT只用于核对，不重复相加。失败、拒绝、额外精度、adjoint与准备工作均保留；所有阶段上限和scope没有转移或重置。

|阶段|complex RHS / 上限|full FE LU / 上限|受控job CPU秒 / 上限|受控job wall秒|
|---|---:|---:|---:|---:|
|R|12,000 / 12,000|1,195 / 3,000|2251.516 / 14,400|2380.219|
|B2|16,140 / 24,000|1,695 / 6,000|8386.091 / 28,800|9021.166|
|C|2,700 / 8,000|360 / 2,000|1992.751 / 10,800|2097.625|
|P1A|13,880 / 24,000|3,470 / 6,000|15283.426 / 28,800|15639.063|
|P1B|15,640 / 24,000|3,491 / 6,000|17767.982 / 28,800|18231.500|

合计 **60,360 / 92,000 complex RHS**，**10,211 次full FE材料–频率LU**；受控CPU **45681.764秒（12.689小时）**，job wall合计 **47369.572秒**。剩余 **31,640 RHS**封存，不是新授权。GPU=0，新三维=0。

每阶段同时核对CPU和受控job wall，均未超其冻结上限。各阶段依次串行，wall不是人为把并行核时当单进程时间。wrapper CPU包含导入、输出和持久化，是保守的完整受控执行费用；solver内层时钟嵌套，不能再加一遍。

最大完成wrapper峰值RSS **4,325,605,376 bytes**；独立采样峰值 **4,309,909,504 bytes**，与冻结6GB过程上限分列。九份早期job receipt没有独立采样RSS字段，原rusage峰值仍保存；未把缺失值填成0。

## 求解与复用的费用语义

full FE LU是材料–频率因子分解，不是独立候选世界数。Green控制另有80次reduced dense-core LU（1228/16/32维），CPU已计入P1B；没有用warm Full-LU=0隐藏小矩阵求解。其1760 RHS/21次full FE LU、setup、核验、残差和instrumentation全部属于P1B总额。详见P1_B_GREEN_COMPUTATION_CONTROL.md。

同一bank的历史精确节点复用另列来源：P1B实际650节点中10个坐标复用已支付P1A端点，200 RHS为历史setup，未重复记为新增；所有新noise parent、连续点控制760 RHS和独立验证都实际付费。P1A/P1B共同teacher、reference、bank不能当免费部署资产。

## 失败与后处理

- `R/r3-graded-optimizer-controls`：CPU 69.725851秒，wall 73.493004秒；失败费用和原目录保留。
- `P1B/p1b-independent-offline-risk`：CPU 0.062040秒，wall 0.792757秒；失败费用和原目录保留。
- `R/r0-truth-closure`：启动前失败，0 RHS/0 LU；0.043587秒为已保守计入CPU总额的launcher elapsed归因，并非实测子进程CPU。wall未独立测量，不写成0。原receipt及ledger均保留。

纯数组独立复算、图库、结构审计和打包在AUDIT/tables中分别保留profile，不新增Maxwell RHS/LU，不冒充从零部署免费。人工代码审阅、文字编写及工具显示没有完整CPU采样，未声明为0。发布与网络用时在发布receipt另列，不挪入物理余额。

本地保留所有现存原始状态。Release包含完整响应、材料、粒子、配置、源码和第一joint见证cal/joint×两世界×五频的明确20份完整复场子集；其余full fields保留本地逐文件库存，不声称全部上传。

封账凭据：[CLOSED.json](CLOSED.json)；逐阶段机器表：[RESOURCE_LEDGER.json](RESOURCE_LEDGER.json)；追加式原账本：[PHYSICS_LEDGER.jsonl](PHYSICS_LEDGER.jsonl)。任何新物理查询均需独立授权。
