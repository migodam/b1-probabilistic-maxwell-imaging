# Research W：先修复可行性，再检验墙–目标补偿与概率递进

**WORK_IN_PROGRESS。已有经过独立网格核验的 Oracle joint 近等价对；尚无在线多解发现，也尚无 posterior 成果。**

两份真实材料世界的目标中心相差0.1个波长、质量相同，完整联合数据白化距离0.798716；两世界joint端点的粗细网格差0.005675、0.015201，低于0.1。两墙calibration距离为3.560512，所以它不是calibration近等价对：墙体变化补偿了目标位移。原24个固定目标probe已全部执行。

原minimum_feature=0.1缺少数学定义，这一缺口不被回写为PASS。已确认的固定支撑材料满足明确的材料范围、被动性、面积、查询域及跨频合同；受控后续实验的范围与原宽材料先验的无条件合法性分开。

原真实世界在相同正演合同下复现均值数据的距离为0；48×32 box/area受限Oracle拟合达到0.968和0.988，但几何全合法性与KKT收敛未证明。测量驱动在线反演仍未达到d≤1。P0-R已经按12000RHS封账。

当前在运行P0-C。两分支受控概率先验已经在任何新观测之前冻结，之后用新独立观测研究P1-A；P1-B仍以可靠概率参考为前提。整个研究没有把Oracle存在性当成在线能力，没有将候选分数当posterior质量。

- [可浏览结果与真实材料图](report.html)
- [Oracle joint检查点](B2_CHECKPOINT.md)
- [R阶段物理与在线裁决](P0_R_GATE_DECISION.json)
- [B2条件裁决](P0_B2_GATE_DECISION.json)
- [新受控prior](configs/W_CONTROLLED_PRIOR.json)
- [独立新授权](configs/AUTHORIZATION.json)

后续结果继续写入同一目录与同一追加式账本；不覆盖旧P0/P0-Q/P0-B，也不重置本轮预算。
