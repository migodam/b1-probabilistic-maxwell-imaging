# 完整资源说明

新增物理0，历史384 BEGIN/384 END和high-water384保持。CPU-only、GPU/NN0、旧locked0、SHA0；未执行额外smoke或完整旧test suite。

主body wall 69.517666s，CPU 67.308709s；补充probes/CLI {"wall_seconds": 5.209958041086793, "cpu_seconds": 5.008944, "polynomial_response_evaluations": 26689720, "new_physics": 0}；报告绘图 {"wall_seconds": 2.0083502089837566, "cpu_seconds": 2.003024, "figures": 8, "new_physics": 0, "scope": "measured report/plot body, excludes interpreter startup and initial imports; native service costs unmeasured"}。

polynomial response evaluations与adaptive次数全部记在results/RESOURCE_SUMMARY.json，不算Maxwell调用。阶段JOB_LEDGER相加一次，不重复相加嵌套总wall。模型服务、人工交互、startup/import及旧历史测量按各自scope保留，不能把未测成本称0或称完整科研wall。

本次所有已测process-body部分合计wall约78.743186s，包含两次绘图和修正容差的限定重诊断；启动/import/交互/模型服务scope未完整测量。full-forward、adjoint、Maxwell RHS、物理acquisition/endpoint均0。
