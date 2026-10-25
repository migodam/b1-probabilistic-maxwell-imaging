# 成本封账

最终共享 StageA 累计：32705状态预算、261632 forward RHS、30152 derivative/adjoint RHS，总291784 RHS，CPU315.360s，controlled wall374.924s。unfinished=0，2次中断/failure保留。63状态余量封存，不继续消耗。

状态预留一次在实际RHS前中断，故forward RHS比states×8少8。所有父线程曲率、接缝probe、qualification都包含在共享总数。逐job完整数值见`COST_BY_JOB.json`。

V1/V2费用全部保留。没有只给最后成功chart节点计费，也没有把8-source block称一个RHS。没有证明相同可信精度下优于最强baseline20%，所以成本扩大门HOLD。OS峰值RSS未测；缓存/数组设计低于1.5GB，不能写成实测内存验收。
