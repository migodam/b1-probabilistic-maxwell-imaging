# Operator pack v1：已付算子教师与材料动作误差对照，网络未训练

状态 `READY_UNTRAINED`。原几何两轮10,351份Full-J数组全部可读；按同一parent、同一坐标、完全相同q+x去除288份重复，得到10,063个context-anchor→paid-J-state配对。没有按Phi、GT或几何成败筛选，没有读取material GT。II analysis-only准备111.446144秒、CPU81.519901秒、RSS129,024,000 bytes、GPU0；新complex RHS、新Maxwell states、物理LU均0。

## 资产与缺失

- 新Jq和Jcoordinate：10,351份PRESENT，native66/field64/probit66坐标链一致。
- 新current、total field、材料和prediction：在原III paid endpoints保存，配对全部可读。
- 原anchor T0/E0/J0/D0：PRESENT；由T0 E0 D0重组成base J的相对误差通过1e−10检查。
- 新D1：未逐点持久化；按冻结P1 decoder及q/x解析重算，probit两属性乘normal pdf，0物理查询。
- **新receiver T1：全部MISSING。** 不从J1伪逆恢复，不调用576列inverse或32个adjoint补齐。因此原T1排序的独立三项分解仍为MISSING_T1。

可计算的替代是精确ordered telescope：decoder=`J(T0,E0,D1)−J0`；field=`J(T0,E1,D1)−J(T0,E0,D1)`；propagation=`J1−J(T0,E1,D1)`。最后一项来自已付Full-J teacher差分，等于final E1,D1下传播变化的作用，**不是独立恢复T1，也不进入部署features**。保存三项norm与三对cross inner products，不把norm之和当可加因果贡献。

## 材料动作与固定对照

对所有方法使用同一个已付candidate Full residual r，在当前坐标计算

`H=I+J'J; delta=−H^-1(J'r+x1); action=D_M1 delta`。

66为完整Gaussianized prior I66，Dg含du/deta；64为固定原u的Gaussian64条件问题，prior I64。共50,315个real≤66 GN参数系统solve，单列为实参数计算；它们不是complex Maxwell/current-system RHS或物理LU。标签比较确切的局部线性材料动作差，不能推出有限材料成像GT风险或真实优化步收益。

五个matched方法：frozenJ；固定T0/E0且更新真实D1的decoder-exact；固定T0且更新真实E1/D1的field-updated；使用已付Full secant的一次good-Broyden；FullJ reference。前四方法的action error均相对FullJ。Broyden使用候选8-RHS Full响应，仅为一次anchor secant，不声称Jacobian integrability或生产adaptive MCMC正确性。

完整66的operator-test（8个EXPOSED开发parent，先每parent平均，再对parent等权）材料action error为：

| 方法 | mean local material action error |
|---|---:|
| frozenJ | .074738 |
| decoder-exact | .010242 |
| field-updated | .007242 |
| one-secant Broyden | .022863 |

这是固定资产分布上的算子fidelity描述，不是网络结果、posterior收益或锁定成像测试。field-updated平均误差较低，但各parent的pair中位数与尾部排序不同，不能只拿一个均值宣称统一最优。

## 划分与pack

原始parent固定：0–15 train、16–23 validation、24–31 operator-test；全部已暴露开发对象，不称locked。Heun仅前8 parent多出teacher，未来训练按parent均匀抽样，再在parent内抽pair，避免其样本量带来额外独立权重。

| 坐标 | train pairs / parents | validation | operator-test |
|---|---:|---:|---:|
| probit66主任务 | 3448 / 16 | 791 / 8 | 790 / 8 |
| field64独立条件对照 | 3452 / 16 | 792 / 8 | 790 / 8 |

`results/II/operator_pack_v1`约276.8MB，包含全部matched delta/material-action数组、4个action-error标签、J误差、Broyden factors、分解norm/cross、固定source J/endpoint引用。完整float64教师留在原来源，训练feature/action拷贝float32，科学误差统计float64；每parent/坐标另有一个确定性float64 dense-J匹配例子。source清单保留所有重复及MISSING_T1，不销毁已付记录。

## 可训练接口，尚无优化

`src/p3r4/operator_surrogate.py` 定义1176个features和343,620参数MLP：256/128/64、LayerNorm、SiLU、4输出，低于0.8M cap。features仅用anchor已付量、candidate Full residual/field和解析decoder；不读J1、teacher误差或传播residual。candidate冷部署仍需8 RHS，anchor J获取32 RHS仍计费。

目标为 `log1p(action_error/.02)`，输出是经验error/refreshranker分数，不能当certificate或oracle fallback。它不是DINO复现，也不独立声称新颖性。`tools/train_operator.py --prepare` 已写冻结模型协议，**未导入torch构造网络、未运行optimizer或训练**。未来由父线程显式调用 `--train`；仅train统计标准化，checkpoint由validation parent均值选，operator-test只在独立显式evaluate入口打开。三seed沿用2026102501/2/3。4项算子代数核查PASS。

## 费用归属

原III获取总费用459,712 complex RHS、16,060 states/LU保留。Full-J source states的forward8+J32共414,040 RHS，其余非J端点45,672 RHS；重复teacher的真实原获取费不因dedup而减去。共同StageI 96+64 context获取费另在原R4I ledger，不能记0。本次仅做缓存数组与real参数线性代数，没有新物理量。

准备入口 `tools/prepare_r4_operator_pack.py`；`FROZEN_PROTOCOL.json`、`MANIFEST.json`、`SOURCE_MANIFEST.jsonl`及 `AUDIT/operator_transport/BASELINE_RESULTS.json`。最终理论、novelty、部署触发器与下一步训练判断归父线程。
