# R4 资源快照与费用语义

|阶段|complex RHS|state attempts|LU/factorization|受控作业wall秒|物理动作wall秒|MPS秒|
|---|---:|---:|---:|---:|---:|---:|
|0|0|0|0|7.27|0.00|0.00|
|I|6,134,384|207,299|207,299|20778.86|2224.35|8376.42|
|II|574,742|52,424|51,400|2174.50|430.74|0.00|
|III|771,392|55,020|55,020|1400.82|502.71|0.00|
|IV|3,239,234|109,073|109,073|3469.24|939.39|0.00|

当前闭账总RHS **10,719,752**；受控作业wall之和7.731h；其中物理动作实测1.138h；MPS 2.327h。并行作业wall求和，不用日历间隔代替；训练wall和MPS时间并列，不相加当额外耗时。最后打包/绘图/审计仍各有analysis收据，最终以RESOURCE_LEDGER.jsonl及封账快照为准。

硬上限58M RHS /94h受控wall/96h GPU；stage分别独立约束，不借历史余额。所有forward、receiver/source adjoint、LU、reduced小系统RHS、失败、拒绝、Full fallback、cache setup、数据与模型训练计入对应作业。离线读取与打包不产生新RHS，但实际wall另记。

历史32,768+2,048 joint pairs与旧R2模型是复用证据，不声称从零部署免费；新增contexts、768新observations/anchors和所有物理控制计入本轮。多seed共享同一teacher的拥有费只记一次。小系统RHS与full状态分列，不能仅看RHS节省宣称加速。

旧GPU行有序列化字段省略，但after.stages.GPU_seconds已真实累计；仅修复将来字段，不重复补记。clarification见AUDIT/GPU_LEDGER_SERIALIZATION_CLARIFICATION.md。


## 最终科研封账

sequence 1522869，closed jobs 337；10,719,752 complex RHS、423,816 attempts、422,792 LU。作业wall之和7.731043h，物理动作1.138107h，MPS2.326782h。无active job。打包与发布的实际时间另记DELIVERY_LEDGER，不伪装为新增物理或免费计算。


## 交付费用单独计量

完整120分卷封包与上传已完成远端名称/大小/状态核对，共175,313,076,233字节。两次封包（含首次中断）的实际wall合计2015.315秒、CPU合计445.831秒；并行上传累计协调器加流wall为18164.499秒，含失败与安全收尾；此前两次串行上传协调器wall合计655.353秒，另以两倍作为协调器加单流上界。网络等待不混入算法速度表。

[独立交付资源汇总](DELIVERY_RESOURCE_SUMMARY.json)对未逐项计时的父线程/最多三名代理整理时间，另用封账以来区间的4倍给出宽松保守上界，含重叠、不是实际消耗；并预留1800秒最终发布检查。即使用该过度保守计数，Stage IV与全程仍低于原32/94小时上限。原科学sequence1522869与10,719,752 RHS保持不变，不重开旧预算。后续出版读回回执可以晚于此快照。
