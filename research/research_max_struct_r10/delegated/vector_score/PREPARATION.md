# 独立三维后评分：准备完成，未执行

root编写`experiments/score_vector_comparison.py`。只读producer源码、公开合同和少量已保存公开receipt；尚未打开任何fresh sealed_truth，也没有运行评分。

必须有root明确解封授权，并先逐项核验`vector_all18.json`全部18终态及公开seal，之后才读GT。预先冻结的χ>.5、χ相对L2、幅度归一化RMSE、ε相对RMS（字段paper_Err）和区域误差保持不变。全四材料、两跨网格条件均展示；不做四对象总体统计显著性声明。

cross材料GT取同几何的grid12有效contrast；current误差也对grid12的coarse-model reference单列。两份grid12/grid18干净端口差从已有truth字段计算，注明是模型/离散差异，不加入测量σ、也不作为在线prior。

保存端点与40/80、各原生阶段都报告；最终选择沿用已封存公开Phi，离线材料分数不选点。current只使用5184维复矢量场，material1728维不投影进去。只重组一次公开receiver S及一次economy SVD，复用已保存G top32；后者只作captured/remainder，不能称完整strong/weak/kernel。辅助/物理cross先检查完全相同m，不匹配则缺失，不补forward。

不新增任何forward/adjoint。图覆盖六条件的固定中心三正交切片、χ>.5三维体素图及原物理corrector收敛。每臂CPU/window明确是嵌套数值区间，完整进程成本仍由vector_compile_costs汇总，不能把它当完整冷启动成本。
