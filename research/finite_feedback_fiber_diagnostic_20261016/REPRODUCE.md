# 复现与数组

实际物理执行版本：`d74f747`，源码快照 `evidence/source_executed/`；Python见SOURCE_MANIFEST.json。使用已有Gaussian/.venv_nn，不安装新的solver。

历史输入：A30 `inputs/a29_release/inputs/a25/{validation,train,background}`；A26 `results/a26/evaluation/261009/PREDICTIONS.npz`；A17/N3_APPLICATION_R1的同case旧Jacobian；原LXZ worktree两份TRAIN校准NPZ。

E0（只读缓存，无物理）：`PYTHONPATH=src python -m finite_feedback.e0`。

实际物理入口：`PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m finite_feedback.diagnose`。账本已封闭，当前入口自动拒绝重复执行，不得删除RUN_CLOSED/ledger重新计数。新物理复现须独立授权与原始历史保留。

仅分析/画图：`PYTHONPATH=src python tools/analyze_finite_feedback.py`、`python tools/report_finite_feedback.py`、`python tools/plot_finite_feedback.py`。只读已支付数组，无新的F/J求解。

结果数组在原B1公开Release的本轮证据包，含原weak/ring完整已支付background G/S/incident/L_b数组，复制不支付新物理。新授权复现时优先读取本轮inputs/background；仍复用固定版本A30 solver源码。`results/endpoints/*.npz` 45份，逐端点即时保存chi/current/field/prediction/LU/pivots/backward residual。`*_JACOBIAN.npz` 全Jacobian、transfer、质量实化；`*_DECOMPOSITION.npz` 谱与GT系数；`*_FEEDBACK.npz` 1/2/3/6/20/rational实际响应；`*_ORACLE_PATH.npz` 9点Jv与复合积分；`*_CHART.npz` 初始坐标/Dg/弱方向。26份legal路径记录对应13个endpoint/anchor，包含共同的frozen/predictor重复记录，不能叫26个解。

review重放脚本只用保存prediction、Jacobian和原decoder，重现12个updated路径端点，不重新解物理。当前源码修复了方向定号、curvature中心与未来wall deadline；实测数据仍归属原执行源码，未暗中重跑。

没有新probabilities，所有候选仅作diagnostic/candidate-only。已暴露development与TRAIN不称独立验证，原始私人Dyson聊天不公开。
