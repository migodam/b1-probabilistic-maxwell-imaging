# 本轮执行偏差与修正

## 进程隔离偏差

已执行A版本e5bd9d3将dataset生成与反演放在同一OS进程，未满足计划中的进程隔离。反演函数只接收观测、prior与物理接口；调用链未使用truth初始化、方向、排序或停止，但这不能替代OS隔离。独立offline评分是另一个进程。该偏差削弱隔离保证，不追认PASS。

封账后已新增独立generate_material_fiber_observations.py入口：先单独生成观测并保留同一A账本，退出后再启动反演进程。新入口未执行真实物理，不据此改写本轮结果；当前CAMPAIGN_CLOSED继续禁止新调用。
