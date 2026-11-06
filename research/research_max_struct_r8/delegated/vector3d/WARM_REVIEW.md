# R8 warm-start 只读审阅

审阅对象：`NEURAL_PROTOCOL.md`、`run_case.py`、`score_warm_starts.py`、R8/R6 solver 与原生 Fourier 接线、公开任务 manifest、硬件分配和 amendment 02。未打开 `sealed_truth/warm_test` 的任何文件，未运行 seal 或 score，未修改运行中的算法。

## 必须处理的实质问题

### P1：公开 seal 未绑定评分时再次读取的关键公开证据

位置：`experiments/score_warm_starts.py:23–33,54–61,66–82`。

当前 seal 仅固定总 manifest、任务回执、预测、预算和 initial 快照。评分通过校验后，仍会重新读取未纳入 seal 的观测 Y/sigma、`test_manifest.json`、`warm_primary_manifest.json`、初值文件、`generation_test.json` 和两个 initializer 的成本回执。上述文件若在 seal 后发生改变，当前校验不会发现，但噪声判据、initial 分数、primary 归属、cohort 或成本结果会改变。因此，“全部预测已封存”尚不等于“整套评分证据已封存”。这不是已发现 GT 泄漏；这是离线证据绑定的缺口。

最小修复只涉及 scorer / seal，不需要改变或重跑求解器：

1. 对评分实际读取的公开文件去重后加入 seal 的 bytes/mtime 清单，并在打开任何测试 GT 前验证：全部 `public_files`、公开 test manifest、primary manifest、generation/feature 成本回执、两个 initializer 回执，以及新增加的本机计时回执。
2. 初值分数直接读取已经封存的 `initial.npz['m']`，代替再次读取外部 `initialization_file`；这是实际求解初始状态，且可少固定一组依赖。
3. seal 同时绑定最后的 solver 配置、协议、amendment 02、最终硬件分配及纯 metric 定义的源码版本。使用现有 bytes/mtime 机制即可，不需要 SHA 检查。
4. seal 前验证 amendment 后的固定 task-ID 集合、唯一性、42 对象和750任务，并核对每份回执/PUBLIC_TASK 的 task ID、case、arm、initialization 与 manifest 一致。当前文件集合已通过下文核验；此项是让最终 gate 保留这一合同，不能只打印任意长度。

### P2：部署前缀混合硬件，但行级 hardware 只标记求解器硬件

位置：`experiments/score_warm_starts.py:44–48,68–74,88–95`。

feature/inference 时间来自本地预先生成的统一回执；目前这些回执没有硬件字段。对于 XINAN 求解任务，scorer 把上述前缀直接加到 Windows solver CPU，并将整行及预算轨迹标为 `hardware='XINAN'`。这个值可以作为明确声明的异构流水线消费账，但不能直接称为“XINAN 本机 NN 部署 CPU”或用于纯同机绝对 CPU 优势。

主线程已提出正确修复方向：保留原 Mac 前缀消费账，另在 XINAN 用相同冻结权重和公开输入独立测量本机前缀，不替换已冻结初值、不读 GT；新计时标出硬件和是否包含模型加载 / SVD。评分中区分 actual hybrid cost 与 same-hardware deployment estimate，并按硬件分层比较绝对 CPU、在同对象内配对。新增计时回执也应纳入最终公开 seal。不需要重新运行当前逆问题。

## 已核验，无需改运行算法

- 最终合同已经从666更新为750：`warm_all_manifest.json` 有750个唯一 task ID、42对象；primary 为630，初始化/特征消融为120。34对象各15臂，8对象各30臂。
- NN current/readout 的完整2×2存在且每格42对象：`cie_dual/nn`（sourcewise+scalar）、`global_dual/nn`（global+scalar）、`a_plus_sourcewise/nn`、`a_plus/nn`（global）。`benchmark_solver.py:13–17` 正确将 sourcewise A+ 转到相同 A+ readout，并显式保留 `sourcewise_cg=True`；未被默认值覆盖。
- 最终执行分区为旧已封存34 + Mac remaining221 + Mac factorial30 + XINAN465 =750，task ID 无交叠、无遗漏。15个 Mac 对象与27个 XINAN 对象不交叠；同对象的全部臂和附加消融保持同一求解硬件。旧 `HARDWARE_ASSIGNMENT.json` 的411/221是 amendment 前快照，最终交付需按已知计划更新；没有把它当作最终分配错误。
- 辅助电流初始化合法：`cie_fft.py:190–197` 使用已支付初始物理电流，计算 `argmin_A ||U A − (j_phys−J_D)||`，随后使用 `J_D+U A`；不是把物理电流直接冒充 DPCS 可行电流。后续 mask 按原 Fourier coefficient 零填充传输。投影 SVD/列数有独立计数。
- solver 不读取测试 GT。公开 Y 输入有私有路径/键拦截；NN/解析初值由公开观测生成，测试幅度从公开 condition 移除。初值不按物理 Phi 先淘汰。
- 辅助日程为500+500+500+1500=3000 PR 更新上限；普通辅助臂追加3000物理 corrector，direct 使用6000物理迭代上限。两者不等同实际 LU 或 CPU，实际成本和停止原因均有独立回执。paid-state observer 读取已支付状态和 incumbent，没有为快照再做 forward。
- material metric 分母确为 `||m−truth|| / sqrt(576) / amplitude`，并同时输出相对 L2；它没有伪装成 `||error||/||truth||`。Otsu 阈值只依赖预测。FAILED 行保留 `material_image_success=False` 并仍进入输出；不能在后续统计时删去失败分母。
- 必须按24同分布、12未训练家族、6跨网格分层报告；scorer 保留 cohort。跨网格没有把 noise-only 门槛冒充完整模型误差预算，且不构造 coarse-GT 的可行 Phi 证人。
- 远端最初8个 import 失败属于额外执行尝试，应保留在 remote audit 并独立收费，不增加750科学任务分母。重试成功不能抹去该开销；主线程已明确保留。

结论：未发现需要停止或修改正在运行的逆问题算法的错误。需要在测试 GT 打开前补齐公开 seal 的最小证据绑定，并修正/分开硬件成本口径；其余核验支持继续当前750任务。

## 审阅后主线程反馈

主线程反馈已将公开 Y、初值、manifest、cost、weights 纳入 seal 的 evidence 并在读 GT 前检查；将补回执 task 身份核对及从 initial.npz 读取实际初值。XINAN 同机部署前缀测量由主线程实施。本工作线程未再次运行 seal/score，未修改这些正在运行的算法。以上问题用于最终整合检查，不能把审阅时快照误当最终修复状态。
