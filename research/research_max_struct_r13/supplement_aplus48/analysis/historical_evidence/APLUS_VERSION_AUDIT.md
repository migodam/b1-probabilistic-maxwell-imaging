# A+名称与实际版本核对

这是来源与执行分支核对，不是版本优劣的因果实验。没有新增逆解。

|版本/范围|冻结配置和接线|实际执行证据|
|---|---|---|
|R7最终主30|`research/research_max_struct_r7/configs/aplus_final.json:49-51`：`aplus_mode=mobius`, `aplus_arc_candidates=true`, `aplus_inherited_gtol=true`|30/30 primary A+ execution receipts都有arc roots与gtol checks，无quintic roots。总arc roots 207295488、gtol checks 269944。|
|R8现存A+回执集合|`research/research_max_struct_r8/configs/solver.json:50-52`同两旗标为true|扫描现存RESULTS内395个含a_plus路径的AUXILIARY_EXECUTION.json，395/395 arc+guard，0 quintic。包含开发、重复、扩展等，395绝非独立材料分母，也不是额外科学消费合计。|
|R13冻结oldC2/机制16|`research/research_max_struct_r13/src/run_experiment.py:198-218`：`solver.arm_config`后直接接`mobius_material_readout`；未合并R7/R8最终json两旗标|16/16机制A+回执有quintic roots，无arc/guard keys。总quintic roots 221212800。|
|R13补充32，与上16组成same48|同一冻结R13 runner；未改算法|32/32新增执行回执有quintic roots，无arc/guard keys。总quintic roots 417300480。|

源码链：

- `research/research_max_struct_r7/src/benchmark_solver.py:13-17`：arm_config只把a_plus映射到R6 a_material并设置arm；不会读取R7最终JSON。
- `research/research_max_struct_r6/src/benchmark_solver.py:34-82`：从common复制，设置sourcewise/global、双beta、原inner tolerances和onehop transform/cache；没有设置两A+旗标。
- `research/research_max_struct_r7/src/aplus_readout.py:50-99`：arc_candidates=true使用状态项二次驻点；false使用包含prior曲率的五次驻点数值求根。
- 同文件`:121-129`：guard默认false；true才检查继承projected material gtol。
- 同文件`:144-166`：候选选择从config取arc默认false，返回`mobius_arc_candidates`或`mobius_mm`。
- `research/research_max_struct_r13/src/run_experiment.py:207-218`：A+的initialize_physical_current为false（只有method.startswith('g19')为true）。

每一份执行回执的路径、关键计数、精确JSON行号见`APLUS_RECEIPT_COUNTERS.csv`；分组汇总见`APLUS_VERSION_RECEIPT_SUMMARY.json`。计数缺失在CSV中为空，不伪装成记录到的数值零；汇总的0表示该类型正计数回执数为0。默认分支的判断同时由源码和非零quintic计数支持，不能仅由缺失key推断。

例：

- `research/research_max_struct_r7/RESULTS/primary/r7_obj_22__a_plus/AUXILIARY_EXECUTION.json`：arc roots 6912000、gtol checks 9000。
- `research/research_max_struct_r13/supplement_aplus48/RESULTS/aplus48_extension32/r13_fresh_0020__a_plus/AUXILIARY_EXECUTION.json`：quintic roots 12337920、sourcewise_cg=false、initialize_physical_current=false。

因此当前48结果的准确称呼是“R13冻结A+ quintic-MM/no inherited-gtol-guard变体”。它在same48 primary平均NRMSE为0.1521001711、common actual CPU为0.1574596474。这个结果不能改称“R7最终arc/guard A+在R13 48对象上已经确认”，也不能仅靠这份核对断言quintic比arc好或坏。
