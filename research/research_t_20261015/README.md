# Research T 复现入口

[独立中文研究报告](RESEARCH_REPORT_ZH.md) · [图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_t_20261015/) · [完整证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-t-20261015)。

代码在仓库src/research_t/，Python使用现有NumPy/SciPy，CPU BLAS单线程。没有新增包安装、GPU训练或重写Maxwell solver。原A30 commit717d3d2、其frozen A29 commit5d5e31f复用；source、原始输入子集和所有NPZ在Release恢复相同路径。

读结果不需要求解：先看METRICS.json、PER_EDIT_METRICS.csv、offline/和T_PHYSICS_HANDOFF/CANDIDATE_INDEX.json。将证据zip解压至仓库根目录后可运行保存数组检查：

```sh
PYTHONPATH=src python research/research_t_20261015/T_PHYSICS_HANDOFF/validate_handoff.py research/research_t_20261015/T_PHYSICS_HANDOFF/candidates/p250065_ring_a1_snr20__minchange_base
```

它不调用solver、不检查SHA、不建立新Gate；只校验真实保存材料与其预测、白化残差和合同。需要重画可用tools/t_offline_results.py；真值只在该离线脚本中。tools/t_analyze.py重算已付费指标，tools/t_export_handoff.py导出已付费端点；不会新建Maxwell状态。

实际已执行顺序（历史记录，不在封口run上重跑）：

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run qualification
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run mechanism
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run mechanism2
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run mechanism3
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.verify_actual
```

各付费阶段代码冻结commit为c752aa6、439773b、1e1b613、a42845f；原协议与收到的Pro包保留。默认Targets通过A30_ROOT指向现有A30或证据中的T_PHYSICS_HANDOFF/vendor_backend/a30_root。任何新物理都需要新的真实用户授权与独立enforcing book，不可删除RUN_CLOSED/重命名job重置预算。

原始cache、truth和offline摘要只作历史/离线证据。G的在线接口加载sanitized背景和默认base，不加载teacher J、truth或oracle step。具体用法见T_PHYSICS_HANDOFF/README.md；它不是另一套Maxwell solver。

证据解压请使用`unzip -n`保留当前checkout源码；旧执行代码快照与当前G接口修复说明见[DELIVERY_NOTES](DELIVERY_NOTES.md)。同Release的interface-overlay.zip可直接应用于只下载主证据的使用者。
