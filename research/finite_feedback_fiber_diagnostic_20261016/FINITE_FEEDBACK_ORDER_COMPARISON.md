# 一阶、有限阶与完整反馈的真实Maxwell比较

本轮新二维证据，非旧L/Z/X合成核验。完整场均共享同一chi，无任意receiver-dark current候选。

L0=I-X0G，R0=L0^-1，E0=incident+G j0。A30的a(chi)=chi，DeltaX=diag(chi1-chi0)。

q1=R0 DeltaX E0；q(n+1)=R0 DeltaX G qn；Delta j=(I-R0 DeltaX G)^-1 q1=L1^-1 DeltaX E0。

每个q项支付8复RHS；rational复用endpoint的LU，再支付8。没有显式R0逆、576列隐藏求解。1/2/3/6阶作用于真实finite材料增量，不是对非线性decoder坐标的Taylor阶数。GT终点额外记录20阶。

| GT全程 | 一阶相对误差 | 二阶 | 三阶 | 六阶 | 二十阶 | rational |
|---|---:|---:|---:|---:|---:|---:|
| physics_R_base | 0.201301 | 0.0496301 | 0.012062 | 0.000174219 | 4.49912e-13 | 3.56e-15 |
| N3 | 0.0825559 | 0.00339382 | 0.00020325 | 2.66101e-08 | 2.94079e-14 | 2.93e-14 |


N3的高阶项与一阶部分抵消，幅度减少约5.18%，但最终响应15.55非零；physics基线的高阶项并未帮助抵消。当前实验没有出现抵消至等散射的GT方向。

## 收敛与稳定性

采用可计算充分上界||R0 DeltaX G||2≤||DeltaX G||2/sigma_min(L0)。GT全程两上界1.129/1.032，检查不通过，只能INCONCLUSIVE；不是实际norm或spectral radius>1的测量。20阶的source-response数值收敛不等于所有RHS一致收敛证明。Rational无需这个Neumann条件，实际endpoint条件数、backward残差均记录。未遇到近奇异失稳，不能声称已测试高非正规极端情况。

GT路径0/.25/.5/.75/1及额外midpoints通过独立full-JVP检验平均Jacobian，Simpson4/8差分别支持数值一致；无等数据、先验合法或路径存在性推论。

A17原三维polarizability是非线性的；其有限DeltaX应为a(chi1)-a(chi0)。这里a=chi恰好线性，不把二维结论冒充原三维测试。A17固定材料的resolvent已含全部多重散射；本轮研究的是材料改变后的finite反馈。
