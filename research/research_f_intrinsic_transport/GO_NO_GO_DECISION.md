# 科研裁决

| 项目 | 判定 | 数值依据/边界 |
|---|---|---|
| F0材料几何定义与数值审计 | GO（限定范围） | gauge去除rank7、三差分步+解析D2、坐标重参数化；缺完整Maxwell二阶向量已隔离 |
| 完整三项向量恢复 | MISSING | 旧缓存只有Maxwell norm，不能补造向量 |
| 独立physics-only数据相容材料 | 已取得2条件 | residual22.2396/21.5121，tau23.7854，P0合法；同parent已知家族/全局几何 |
| Gate2完整anchor资格 | NO-GO | 冻结最优性.001；实测.0452734/.00200099 |
| F2配对finite transport | NOT_RUN | 遵守Gate2依赖；不扩大F1分配、不放宽标准 |
| Gate3任务尺度位移 | NOT_RUN | 无传输候选可评价 |
| Gate4/learned Explorer训练 | NO-GO | 非平凡teacher pairs=0，训练本轮不授权 |
| 原P0 exact-local-fiber主张 | NO-GO（当前证据） | 四个测试点均数值满列秩；非全局唯一性定理 |
| 材料原生坐标/边界约束研究 | CONDITIONAL GO | 已证明raw曲率大部分可由坐标项解释；仍缺配对有限传输实验 |
| 强IEEE TAP算法成果 | 尚不支持 | 数学分离有诊断价值，未建立同成本新材料探索收益 |

本轮封账。下一步具体优先级：先在原P0上实现权重比例商坐标/可行切锥与更可靠的有界物理驻点停止；需要另行明确后续执行授权后才继续求解。保留原uniform生成测度，不偷偷定义均匀七维prior。达到anchor资格后，再做事前固定同材料距离、同物理预算的retraction对照。当前不选择训练NN，也不转而放宽P1来救本次结果。
