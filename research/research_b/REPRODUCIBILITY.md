# 复现与费用

现有Apple Silicon CPU环境 `Gaussian/.venv_nn/bin/python`，NumPy/SciPy；没有安装GPU、CUDA、daemon或新框架。`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`。程序运行优先级nice10。

完整原始字段、forward cache、逐步rankings和实验数组随本阶段Release证据ZIP提供；仓库保留报告、代码、必要配置、summary和费用。不能仅凭缺少本地缓存就运行新solver补齐E0。

- `code/verify_theory.py`：用户要求的11项零物理理论/概率核对。
- `code/run_e0.py` 与 `code/e0_candidate_controls.py`：旧缓存上的真实排名和候选更新；硬拒绝新物理。
- `code/prepare_e1.py`：唯一付费forward入口，先预留attempt，失败和数值audit计费。公开结果已完成，不为发布重跑。
- `code/active_experiment.py --stage E1/E2 --job NEW_ID`：读取已保存模拟观测，执行冻结策略，产生实际材料评分，不调用Maxwell。
- `code/coverage_stress.py`：删除后重新选测与独立instrument mismatch，物理0。
- `code/compensation_experiment.py` 和 `code/horizon_witness.py`：非线性合法补偿与同一材料世界的全菜单见证，付费回调统一入账。
- `code/e3_real_handoff.py`：四个真实A包只读接入与新观测/候选refinement；A不被修改。
- `code/analyze_results.py`：从已完成raw结果重建parent级表和固定bootstrap；不要将两份噪声当独立物体。
- `code/history_reference_audit.py`：限定固定子集joint history积分误差，不重跑物理。
- `code/plot_results.py`：标准科学图，读取已保存结果。

新独立复现实验应有自己的授权目录与账本；不得重置或覆盖已发布PHYSICS_LEDGER。研究包原协议保留，本轮授权仅在 `configs/AUTHORIZATION_OVERLAY.json` 生效。既有P2 reference所需raw资产路径与revision见SOURCE_MANIFEST；下载旧P2发布资产用于只读复算即可。

所有权重与概率语义以相应记录为准。A和固定625候选没有posterior weights；reference策略使用连续proxy的有限数值近似，不是exact full posterior。old observation只能条件化一次，repeat必须使用新独立noise ID且同一持久gain。

## 文件布局和无新物理复算

先clone本B1仓库，下载本阶段Release的证据ZIP，在仓库根解压（ZIP自带`research/research_b/`路径）。B的E1/E2数组、所有full响应和E3实际交接输出都在该目录。`code/analyze_results.py`及`build_delivery.py`直接读取保存结果；这不产生新物理。大JSON与全部parent表在Release，仓库另有final摘要切片。

E0/从头物理准备还需要Joint P2 Release的历史原始资产。原运行布局为仓库的同级`b1-joint-p2/research/joint_p2/results/`；具体资产名由SOURCE_MANIFEST与CACHE_MAP列出。E3从头模拟另读同级`b1-research-a/research/research_a_20261015/`的原A观测/离线truth；B保存了四份native snapshots及实际simulator字段，原A阶段也已在本仓库保留。可以复建此文件布局，不用改物理模型或旧资产。

不要把“无法找到外部缓存”作为新付费求解的理由。`prepare_e1.py`/E3的paid入口仅供独立新授权实验；现有结果已完成，不重置旧账本重跑。重新生成所有图表、bootstrap或额外reference积分也是新的计算，应在新独立复算目录记录。

本阶段科学合同27项在`checks/FINAL_CONTRACT_TESTS.json`，源11项在`checks/THEORY_CHECKS.json`。前者检查likelihood、精确材料loss、预算及A/B数据语义，均不调用solver；它们不能替代实验结果。
