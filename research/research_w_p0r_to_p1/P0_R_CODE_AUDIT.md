# P0-R 代码审计与边界编码修复记录

状态：**WORK_IN_PROGRESS / 保存结果说明**。本文件只读现有结果、审计与源码文字；新增 Maxwell、LU、RHS 均为0。科学判决与后续实验由主线程维护，本页不修改门或账本。

## 1. 旧搜索与新优化器不能混写

完整旧审计在[AUDIT/old_code/OLD_CODE_AUDIT_ZH.md](AUDIT/old_code/OLD_CODE_AUDIT_ZH.md)。旧错误版本曾把被拒提案当作corrector锚点；其费用保留，科学比较排除。修正版恢复accepted anchor，但仍只有一次提案和一次corrector，且存在不连续支持投影。旧20–50残差不能被简写成“完整FWI已收敛仍失败”。

新通用优化器接口与接受规则见[AUDIT/optimizer/implementation.md](AUDIT/optimizer/implementation.md)。它用实际/预测下降比、盒约束子问题、accepted full-J refresh与独立投影梯度终止。objective floor、call cap、迭代上限、失败与梯度收敛分别记录；这些实现规则不保证非凸问题全局收敛。

## 2. 固定材料导数检查

[results/R2/RESULTS.json](results/R2/RESULTS.json)包含6状态×3步长=18行。对边界使用合法单侧方向，内部状态使用中心差分，二者分开解释。

|case|scheme|最佳FD相对误差|已保存最大伴随缺陷|
|---|---|---|---|
|air|one-sided|2.41549264e-06|1.30750182e-15|
|area_boundary|central|1.38133085e-09|3.99997816e-14|
|box_boundary|one-sided|1.02576371e-05|3.24578537e-14|
|cavity|central|2.12241788e-09|7.9862021e-15|
|high|central|3.08053489e-09|1.04188196e-15|
|weak|central|1.06997569e-09|1.85972772e-15|

中心差分各状态最佳值中的最大值为3.080534886e-9。它支持当前固定离散系数导数与实伴随；不验证top-K投影可微性、支持生灭、形状导数、拓扑覆盖或连续几何合法性。旧审计没有找到支持修改Jacobian符号/共轭的证据。

## 3. boundary roundtrip bug：合法σ在编码后越界一ULP

[AUDIT/boundary_encoding/ROUNDTRIP.json](AUDIT/boundary_encoding/ROUNDTRIP.json)记录：在f=0.8和0.99875，旧编码会得到`f*Im(epsilon)=0.15000000000000002`，严格后端因此拒绝合法上界σ=.15。新编码的对应值为.15，最大物理浮点变化2.7755575615628914e-17；`physical_bounds_enlarged=false`，`backward_solver_tolerance_changed=false`。

[src/wr/physics.py](src/wr/physics.py)的`PixelModel.coefficient`先验证实际σ盒约束，再用实数除法构造虚部；仅对roundtrip向外的一ULP用nextafter向内编码。这里修的是边界表示一致性，没有加宽物理允许区间，也没有放宽线性解残差。

首个[LBFGSB结果](results/R3/online_2026102201_wide_LBFGSB/RESULTS.json)迭代数0、d=45.14926584，停止原因为passivity异常。相关ITERATIONS明确保存1次evaluation_failed；三个旧TR轨迹也保留同类异常（seed2201wide 11次、seed2202wide 13次、seed2202narrow 2次）。不能把这些受异常影响的路径当作修复后的算法公平对照。

独立[LBFGSB_REPAIRED结果](results/R3/online_2026102201_wide_LBFGSB_REPAIRED/RESULTS.json)达到d=17.60304259，17次接受、20次完成评价、0次评价失败，最终受max_calls停止，投影梯度仍为1。它说明这次编码异常已不再中断该路径；没有证明在线反演已解决。

## 4. 资源与没有补做的检查

[R封账](R_SEALED.json)保留12000 RHS、1195 LU及2251.515528 CPU秒，disk/smooth_block controls未完成。本文件没有运行物理程序、重新收费、删除失败或修改历史/新物理源码。控制结果与局部非凸现象见[P0_R_FWI_CONVERGENCE.md](P0_R_FWI_CONVERGENCE.md)。
