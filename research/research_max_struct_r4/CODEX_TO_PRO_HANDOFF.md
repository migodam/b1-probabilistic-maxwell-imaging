# R4到下一轮：公式优先的交接

本轮完成A/B/C最小验证、六对象冻结确认、分块/相位/网格/平台检查及一次结构修正。未得到成熟的普遍优越算法，不把这个状态当整项Maxwell研究终止。下一轮从已确立的结构与反例继续，不能重跑整套恒等式或调β找positive。

## 已有可复用入口

- `src/cie_onehop.py`：冻结四格一次交换；current与材料梯度一致，O(Ns)附加动作，无小LU。
- `src/cie_partition_average.py`：平移平均修正，已失败，不默认替换onehop。
- `src/port_profile.py`：固定eta的精确envelope实现；已测无收益。
- `src/cooperative_geometry.py`：复用完整Maxwell物理Hessian/feedback，不是新forward。
- `src/confirmation_solver.py`：只接public数组，无文件读取；原18配置保持冻结。
- `experiments/replay_public_confirmation.py`：跨机器单任务公开数据重放；输出必须是新目录，不能重新成为blind。
- `results/failed_endpoint_stationarity`：六个失败端点的实际原物理投影梯度。
- `results/failed_endpoint_continuation`：A/native各一次固定追加预算，停止原因、材料退化和全部费用保留。

## 目前最值得推进的问题

精确局部读出关系是 `D_M−C†=(CᵀMC)⁻¹CᵀM(I−Π)`。问题不再是“左变换会不会改变算法”——会改变已经明确。真正的问题是：传播过来的本构不可解释defect，何时对应原物理目标的有益材料方向，何时只推动错误的box饱和结构？

可先在已保存成功/失败阶段状态中，显式分离材料切向和正交defect，记录 `CᵀM(I−Π)b` 与实际readout改变量、原物理梯度的夹角/作用。源metric必须先吸收，实材料须实化，不可把复Hermitian投影自动当实材料投影。prior、活动box、非零defect与不精确内层要分别记账，不能套用根处无prior公式给整个算法下结论。用少量已有状态检验这一解释，比继续任选分块、扩大邻域或增训练数据更有针对性。

若提出新方向，只给一个清楚的机制修正和决定性反例；先冻结，再生成少量真正新数据。现有六对象已暴露，允许作为开发分析但不得重新命名confirmation。保留成熟双βCIE、充分普通模态/迭代以及同原物理corrector。任何GT选方向只能标OFFLINE_DIAGNOSTIC，不能作为在线改进。

## 科学边界

R4是有限2D scalar-TM，并非3D vector或连续PDE。M的H(curl)问题和T的可用异常维数/误差证书仍需真正对应的资格，不可借R4图像通过。R1张量同场族和Z-ALG-E1图册属于信息结构，不能把本轮远超noise的失败端点当它们的多解见证。两失败材料的真实合法模型确实能满足noise球。

EBA预条件CSI（2000）、NIE/CIE、CS-EB、DCTMC和WRI/CC-CSI都有先例；把局部传播称为新发明不成立。可论证的论文价值首先是特定off-shell几何、读出结构和已复现的适用/失败边界；没有成熟优越算法证据时不写普适突破论文。

## 执行习惯

读取实际`AGENTS.md`和operator catalog，复用公共G/E/S和已有CPU/CUDA原语。这个小块SciPy工作负载的Windows replay比Apple Silicon慢，不能因为有4060就声称应全搬GPU；真实GPU适配须另有whole-task证据。不要新建服务/队列框架或改变驱动。

子agent显式配置：Sol medium处理读取/整理/绘图；Sol high处理数值/评分/审查；核心推导与最终科学判断留主线程。费用写实际scope，失误也收费。发布继续使用现有GitHub/site，保留历史证据；不验证SHA。
