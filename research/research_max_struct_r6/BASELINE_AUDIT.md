# R6 基线与可比较范围

| 名称 | 电流更新 | 材料读出 | 来源 |
|---|---|---|---|
| CIE | 原sourcewise PR；current/material同β | 原实材料CIE读出及公开prior | FBG-R3成熟实现，仅加入观察器 |
| Dual-CIE | 原sourcewise PR；current β=1，material β按阶段1→.2 | 同上 | 同一成熟实现 |
| Global-Dual | global PR；同Dual β | 同一scalar CIE材料读出 | R2原实现 |
| A-material | 同Global-Dual current | R4一次邻域传播，只用于材料块；R5精确重复值缓存 | R5固定改进 |

每臂mf3/6/10各500步，再自然续接mf12共1500步；没有45秒截断。每个阶段的原始Fourier系数零填充，PR历史在原stage边界重置。随后同一个原Maxwell L-BFGS-B corrector最多3000接受迭代，原ftol、gtol、box与prior均保持。报告实际停止原因；达到上限不称收敛。最终选原Phi最小的已付费候选，不按GT或IoU选早停。

CIE单β是每个阶段两个块同β，阶段日程为1→.2→.2→.2，不是全程恒定β或本轮参数寻优。Dual-CIE是已存在的分块β机制，不属于本轮发明。主四臂不是四个独立数据样本。

原code的DPCS截断为相对阈值.001，所有方法相同；Fourier/APCS容量扩满仍保留DPCS affine约束。最终材料始终576个独立实未知量，不知道数字/shape类别或support。所有算法保持同一原物理、数据、初始化和公开正则。新数据为同离散生成，不能算独立网格或实测验证。

## 实现审查

`src/mature_cie_observed.py`是冻结成熟CIE的可读局部副本；差异与源文件存于audit。`benchmark_solver.py::refine_observed`沿用原R4corrector，只观察已经付费的正演/梯度。一次最小观察器接线检查证明该观察不改变相应有限诊断输出，记录在`audit/OBSERVATION_EQUIVALENCE.json`；它不构成新的Maxwell正确性实验。

`audit/INDEPENDENT_METHOD_REVIEW.md`逐项检查同prior、source metric、beta、native Fourier transport与FFT伴随。无阻塞主比较的数学/公平性错误；已知Global/A没有同粒度逐步辅助operation_counts，所以该细项不填零、不伪造完全统一的算术时间线。完整臂级counter和物理LU/RHS仍可比较。

## FFT及NIE

FFT副实验只在预先指定的历史b20d和b2449上各跑Dual-CIE/A，共4任务；与主任务同数值日程。替换真实U/G action与adjoint，接收和corrector仍dense。它是同观察的工程复跑，不增加独立确认对象。浮点路径差异须保留；若物理解不同，不用GT选择更好后端再汇报。

NIE原论文有实际FFT-TSOM算法；本轮已核查全文但没有另行完成成熟独立NIE复现。故不填一个名为NIE的自制弱基线，也不将缺失源码解释为NIE无效。主CIE/Dual-CIE比较直接对应A的继承框架，参见`NIE_BASELINE_AUDIT.md`。
