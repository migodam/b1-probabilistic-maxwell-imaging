# E0 点似然与缓存pair审计

白化保持512实通道、每实/虚通道sqrt(2)/同一已知sigma；`ell=||prediction-y||²/sigma²`，tau=23.785449004057032。没有候选幅度重估sigma，没有新forward。

2835审计记录中2100缺full prediction，full residual/ell留空；A29 pool的L是背景线性近似，不能填充full Maxwell似然。其余读已付费exact缓存。相同y/sigma下30068个distinct-material比较，230个满足|Delta ell|<=log(10)，其中两端同时PHASE合格为0。没有P0_EXACT竞争见证；近点似然不等于posterior mass。

主缓存screen是根`CACHED_PAIR_SCREEN.csv`：64个独立TRAIN parent，2背景×2幅度×2016 pairs=8064比较。每pair冻结lower-ID首端原始sigma，并共同用于两端；不分别按预测归一化。1632个满足冻结结构规则；0个结构合格pair达到优先目标d_pred<=2。全体最小25.40616862516396（weak,a=.3,parent250037/250051，结构不合格）；结构合格最小44.8267177317206（ring,a=1.2,parent250001/250005）。四格最小见`cached_pairs/SUMMARY.json`。这些不是8064个独立几何实验。

`delegated/cached_audit/CACHED_PAIR_SCREEN.csv`仅为替代诊断：幅度1.2、4032比较、每背景一个固定development median sigma，最小27.105430725344018。它不是主表，不混合两者分母、sigma或最小值。

历史truth有4/32在95% noise ball外，属于预期随机事件；不称材料物理非法、不调tau。旧32条件及250065–250080均为exposed development。G有18/72 fresh validation进入noise ball，但相对其同观测旧base的Delta ell为13.3707至64.4537，均超过log(10)/log(100)；旧base也不能冒充新的P0-compatible reference。

有限缓存未找到pair不证明不存在竞争材料，不称GLOBAL_UNIQUE；posterior_mass=null。本阶段新增states/RHS/factorizations均0。
