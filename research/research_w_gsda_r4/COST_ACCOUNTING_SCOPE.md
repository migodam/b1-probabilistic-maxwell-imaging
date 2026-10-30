# 成本口径与已发现的计时缺口

FE/RHS以物理动作之前的追加式reservation为权威，失败尝试也计入。GD/GDH、正向、加权伴随、接收伴随与小dense-Schur分别记录。kernel CPU不再叠加到外层receipt。批总账只补分区未覆盖的CPU，不与所有分区再次相加。

**本轮已有receipt字段wrapper_process_cpu_seconds的精确范围是BudgetBook构造之后，或batch分区开始之后。** 它包含该范围内的网格、装配、LU、求解、日志、I/O、输出和检查，但不包括较早Python启动/依赖导入。科学作业的完整外部CPU不能从已经退出的旧进程恢复，不补猜测值。早期部分只读评分和交互编辑亦未计时。因此全研究所列已记录CPU是已覆盖范围总和，不是机器全部端到端CPU的严格总量。

该缺口在成对控制执行期间发现并保留。后续B0控制启动器如取得RUSAGE_CHILDREN user+system差，则单列whole-child记录，不与已有receipt重复加入。CSI控制外层有真实start/finish墙钟；其启动CPU仍标未测。六次错误解释器的preimport失败为0PDE，但墙钟实测、CPU未测，原日志与命令均保留。

小dense计数覆盖显式记录的模态Schur LU和condition-SVD。其它GN/Gram/投影诊断中的小矩阵求解计入外层CPU，但没有统一的逐调用factorization计数；不能把该字段冒充全部LAPACK调用。FE/RHS全预约账与这一小dense计数范围分别保存。

每方法冷部署必须包含公共参考setup；warm分摊以实际完成case数报告。B4公共setup可同样给普通B5，不属于吸收特有节省。Full-GN从B0终点开始，比较成本必须B0+GN累加。

由于成像质量未通过、旧CPU口径存在启动计时缺口且部分方法未收敛，本轮不提出严格端到端CPU加速结论，也不由较少FE/RHS宣称同质量优势。原始计时值不改写；PHYSICS_COST.json、每job receipt和外部计时文件共同供复核。
