# 四平移分块平均传播：两例开发预测已封存

状态：`POST_CONFIRMATION_DEVELOPMENT`。两例各启动一次、顺序 CPU1 执行，结束后停止。worker 未读取 GT、私有映射、确认分数或确认预测，未进行材料评分。父线程负责评分与最终科学判断。

新算子把四种 2x2 分块平移得到的**传播算子**平均后再形成一个残差及损失：轴向邻居真实复 Green 权重 1/2，对角邻居权重 1/4，其余非对角元为零，边界不周期 wrap。该结构移除任意分块原点，未增加观测或先验，也不是平均四个目标函数。

正向与伴随动作分别使用 4324 非零元 CSR 及其共轭转置；current/material/adjoint 使用同一个 off_action。物理 G 保持不变。保留标量 CIE seed、source metric、global PR、Fourier/DPCS、三个 500 次阶段、材料 inner 60 次上限、host 0.1、box/prior、原物理 corrector 1000 次与只按 original Phi 选择。

实际邻居绝对值最大行和为 0.03906209771759682，其中轴向 0.027136304510478354、对角 0.011925793207118465。两例的全 box 可逆充分界均为 0.16809971545079827（β=0.2）和 0.03784328908973763（β=1），与实际原始 disjoint block 的解析值差均为 0。

唯一新增稀疏动作接线方向：固定 J=.1E、Q=E，β=0.2、h=1e-6，相对误差 1.6483373247944926e-8，PASS；0 physics。第二例未重复。

| 输入 | Auxiliary original Phi | Final original Phi | 相对残差 | Wall / CPU (s) |
|---|---:|---:|---:|---:|
| obs_b20d5a55677d4eca | 0.0172020699636 | 0.00581543424693 | 0.107660170671 | 42.49006 / 40.67534 |
| obs_42cb233e4f224bb3 | 0.0851915456345 | 0.010320544208 | 0.143490952543 | 49.20717 / 47.85336 |

| 输入 | Auxiliary wall (s) | Scalar divisions | Nonzero neighbor RHS | Physical LU / RHS | Corrector iterations / evaluations |
|---|---:|---:|---:|---:|---:|
| obs_b20d5a55677d4eca | 32.26509 | 1,116,847,872 | 7,835,641,472 | 1043 / 24948 | 1000 / 1036 |
| obs_42cb233e4f224bb3 | 39.85668 | 1,427,764,608 | 9,990,135,008 | 1019 / 24372 | 1000 / 1012 |

每个 sparse nonzero neighbor RHS 是一个非零矩阵元与一个 RHS 列的乘加贡献；不能当作小块 RHS 或完整物理 RHS 相加。没有小块 LU/求逆。两例各完成 1500 次 PR；两个 corrector 都达到 1000 次上限，未证明收敛。完整梯度接线成本、物理成本、时间、历史、三阶段 m/J/Q/coefficients、候选值与 prediction.npz 分别保存在各输入子目录。源与配置均在运行前保存快照。

两例程序总 wall 91.69723 s、CPU 88.52871 s；物理总计 2062 LU / 49320 complex RHS。未自增预算、调参、追加网格/seed、重跑、远程操作、提交或推送。

元数据说明：frozen config 保留了 inactive 的历史 block_alignment 标签，描述原始分块分量；新模块不读取它。实际算子由新模块以及 block_Green/residual/weights 字段明确给出。为保留审计记录，已冻结配置与快照未回改。

这些是确认结果公开后的结构开发实验，不能计为新的确认结果。物理残差、可逆充分界与梯度接线通过均不构成材料恢复或 basin 证书。
