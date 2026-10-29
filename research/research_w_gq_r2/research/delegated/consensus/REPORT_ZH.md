# E3 旧 bank 条件目标共识：零新 PDE 回放

本轮结论：此旧 bank 的全体条件目标分布明显分歧，条件共识不能支持省略背景后验相对质量。24 个案例均为 `NOT_CERTIFIED`；未触发冷 PDE。标签为 `REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE`。

在线只开放有限 prior 元数据、TX0/1 主观测 y 和单原子主端口 API。28 个完整曲线由 prior 元数据预先选定，各 support 分支先取一墙，再按 prior 加权最远点填满；另付全部 118 墙的中心 anchor。三种方法共享同一 790 原子记录。条件共识只用已查询行的条件目标分布均匀混合，不计算背景后验权重；它是有意标为未认证的候选。

24 个案例先全部封存，之后才读取 2950 条旧主端口响应生成 teacher 并读取旧 truth_index 评分。新增 FE/RHS 均为 0；在线共享查询共 18,960，teacher 建立 2,950 原子，24 个参考各 2,950 原子等价评估，共 70,800。历史 bank 建立成本仍为 2,950 FE、11,800 RHS、9,288.02304 CPU 秒。此次主任务 wrapper 为 4.147437 CPU 秒、4.898 墙钟秒、峰值 RSS 83,066,880 字节。未重复 M4/M5，未作新物理加速主张。

封存文件的独立回读校对另耗 0.225483 CPU 秒，零新增 oracle 查询；本单元已记账总计 4.372920 CPU 秒。两份 wrapper 收据分别在 `results/jobs/e3-consensus-logic/RECEIPT.json` 与 `results/jobs/e3-consensus-artifact-review/RECEIPT.json`。

| 方法 | 平均目标 TV | 最大 TV | TV≤0.02 案例 | GT RMSE | 相对参考均值的平方 regret |
|---|---:|---:|---:|---:|---:|
| ORDINARY_RESPONSE_INTERPOLATION | 0.074915 | 0.970135 | 18/24 | 0.613416 | 0.021077 |
| ANCHOR_WEIGHTED_CONTRAST | 0.357518 | 0.999984 | 13/24 | 0.755230 | 0.216763 |
| SAMPLED_CONDITIONAL_CONSENSUS | 0.783602 | 0.998349 | 0/24 | 0.734923 | 0.233261 |

全有限参考的 GT RMSE 为 0.595430。先在每个 parent 内平均两次噪声重复，再等权平均 12 个 parent；RMSE 为该 MSE 的平方根。parent 与 family 分解见 `PARENT_RESULTS.csv`、`FAMILY_RESULTS.csv`。empty/nonempty 与 support_count 分组仅在封存后由 truth 产生，未用于查询。

采样条件 TV 直径范围为 0.999998043–1.000000000，全表为 0.999999989–1.000000000；全体 24 个案例均远大于 0.02。直径出现 1+约 2e-16 属浮点舍入，未声称严格舍入认证。未查询背景的实际后验质量从 1.52272e-29 到 0.999998238，因此不能用统一删除未查询质量来救场。

完整响应与目标增量的白化实 RMS 配对距离中位数分别为 2.405921 和 0.902169；背景差在全族增量空间内的投影范数中位数为 19.343964。该全表有限增量矩阵在相对奇异值阈值 1e-12 下的经验秩为 64（实观测维数 64），不能把增量变化较小直接当作残差投影可忽略。四类完整配对距离保存在 OFFLINE_STATIC_DISTANCES.npz 与各 OFFLINE_DISTANCES_*.npz；这些全 teacher 诊断不产生在线查询节省。

默认未知 likelihood 始终使用 [0,1]，全体条件 TV 半径上界为 1；采样直径不覆盖未查询墙。有限 saved-Gaussian 后验区间与所有已付原子均经过封存后回读校对。完整残差投影恒等式最大绝对差为 2.56e-13，两背景差分恒等式为 4.3e-13；一般 Dobrushin 混合证明、逐行尺度反例与未查询原子反例详见 PROOF_AND_SCOPE_ZH.md、ALGEBRA_CHECKS.json。它们是代数与一般概率工具，不是新增 Maxwell 世界或物理证书。

保留边界：此结果只否定本冻结表示和旧离散 prior 上的全体条件共识入口。它不否定所有物理分组或局部条件共识路线；后者仍需在完整表之外可获得且覆盖未见背景的物理界。没有该界，不启动冷 PDE。
