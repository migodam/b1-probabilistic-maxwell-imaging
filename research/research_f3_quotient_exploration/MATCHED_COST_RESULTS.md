# 增量费用、冷启动等价费用与匹配范围

总计55新states、712complex RHS、55LU、17.640707751秒受控wall；其中440forward RHS、256adjoint RHS、16方向Hessian RHS。55个endpoint均立即落盘，0丢失。GPU/Mie/新3D/训练均0。

| 阶段 | 新states | 新RHS |
|---|---:|---:|
| F1.5 | 4 | 160 |
| 原pilot余下C1/C2（共享缓存后） | 42 | 480 |
| 追加既有端点路径检查 | 9 | 72 |
| 合计 | 55 | 712 |

原pilot首attempt在8states/224RHS后因1ulp规范权重越界被validator停止；这些成本包含于46/640，不删除。一次数值修复后只读恢复8个state和已付费Jacobian，未重新运行F1.5或已完成A/B端点。所有job的wall、失败与cache事件都保留在同一本ledger。初次closure、recovery closure、最终stage closure均保留，无重置额度。

每条件C1实际费用：A两个端点16RHS；B第一个复用A、第二个新state+刷新J，增量40RHS；C两个端点16RHS；D两个端点16+Hessian8+刷新J32=56RHS。共同anchor与Jacobian准备费用另列，不给D免费准备优势。

C2 native pool为2个C端点复用+6个新端点，实际增量48RHS；MH8个新proposal=64RHS。比较D56与native前7端点时，用**独立冷启动等价成本7×8=56RHS**，包含两个已在C臂付费的端点。不能把“native增量缓存费用”与“D从零费用”直接比较。双方同anchor/方向信息/noise，prefix无GT选择。native前7 weak有2个任务尺度competitive正密度端点，D有1个；ring均0。

匹配的是完整operator/RHS预算，不是实际消耗相同state数，也不是已完成cold/warm wall吞吐benchmark。整个controlled wall含decoder、SLSQP和持久化；没有可信的每方法独立wall数据，不宣称时间加速。

剩余9states/1336RHS没有转用于搜索正例。旧ResearchF24/960、finite-feedback45/1880、LXZ及A17预算只读；原背景G/S/incident及两个旧F1state/J的来源和历史费用披露，未伪装成从零部署免费。
