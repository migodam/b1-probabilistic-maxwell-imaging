# 严格 GN 锚点健康审计

状态：**METADATA_AND_PUBLIC_M0_AUDITED**。快照 `2026-11-07T07:17:33.395741+00:00`。仅读取 online manifest 指定的 TRAIN/VAL canonical JSON、这些 ID 的公开 `m0`，以及父线程明确指定的两份中断记录；未读取 truth、TEST 数组或 GT 分数，未运行 GN、forward、Jacobian 或任何公式/误差地板测试。

共 1164/1164 条观测，582 个独立母体；clean/noisy 是相关观测，不能当独立材料。全部严格 primary：ftol/xtol/gtol=1e-12、discrepancy_stop=False、weakreg=0。

## 实际终止状态

COMPLETED_NUMERICAL **618**；BUDGET_STOP **546**（WALL_SECONDS 476，MAX_LU 70）。数值终止不是全局最优、唯一性或材料恢复证书。

|split/noise|观测|数值终止|预算停|最终数据残差 median [min,max]|wall median [min,max] 秒|
|---|---:|---:|---:|---|---|
|train|0|504|504|0|4.15201e-07 [1.25024e-08, 1.29007e-06]|7.08504 [0.786323, 40.8377]|
|train|0.001|504|24|480|0.000963598 [0.000944937, 0.000978858]|45.1263 [16.9735, 46.6874]|
|val|0|78|78|0|4.332e-07 [4.14017e-08, 1.07942e-06]|3.84251 [0.966027, 9.64656]|
|val|0.001|78|12|66|0.000962023 [0.00094873, 0.000977114]|45.014 [16.0947, 45.3619]|

## 8/12/24 阶段与真实工作量

|阶段|事件状态|nfev recorded / missing|njev recorded / missing|forward 操作 median [min,max]|J 构造 median [min,max]|阶段秒 median [min,max]|
|---|---|---:|---:|---|---|---|
|8|{"FINISHED": 1160, "PHASE_BUDGET_STOP": 4}|1160 / 4|1160 / 4|37 [12, 100]|35 [12, 99]|0.980766 [0.327629, 7.81939]|
|12|{"FINISHED": 1020, "PHASE_BUDGET_STOP": 144}|1020 / 144|1020 / 144|106 [9, 200]|93.5 [9, 185]|4.97264 [0.325852, 19.914]|
|24|{"FINISHED": 618, "GLOBAL_BUDGET_STOP": 546}|618 / 546|618 / 546|38.5 [1, 358]|34.5 [1, 307]|17.25 [0.0124225, 42.1641]|

nfev/njev 的完整分布和逐案例字段在 JSON 中；缺失保持 NOT_RECORDED，未以总 forward/J 次数冒充 SciPy 返回值。所有阶段的接受步数均 **NOT_RECORDED**。`accepted_endpoint_available` 只表示有接受端点，并非接受步计数；这里也不把 `njev-1` 自动改名为实际接受步。

全部 1164 条标记释放 full24，最终阶段均为24。宣称 full24 但 full24 **0次真实forward评估：0**。full24 只有 **1次forward且1次J：31**；这些案例没有第二次 full24 目标评估，不能拿‘full24已释放’当作做过 full24 迭代更新的证据。逐ID清单在 JSON。

## 实耗与端点诊断

|量|median [min,max]|总和（仅canonical）|
|---|---|---:|
|LU|191.5 [23, 500]|248540.0|
|RHS|10572 [1332, 28272]|13608480.0|
|wall_seconds|25.4524 [0.786323, 46.6874]|30178.937207917217|
|feature_wall_seconds|0.00057875 [0.000103042, 0.0136437]|0.99054978962522|
|process_task_wall_seconds|25.5578 [0.882292, 46.8485]|30308.411843090085|
|wall_overshoot_seconds|0 [0, 1.68797]|105.56683608563617|

投影梯度/梯度映射：**NOT_RECORDED，1164/1164**。full24 SciPy `optimality` 有 618 条，分布 4.93227e-13 [2.91742e-14, 4.17044e-10]；它不是本轮另行计算的投影梯度映射。没有为补字段多算 J。

box 活动从公开 m0 计算，阈值1e-10：lower活动cells=37204，upper=0，有任意活动约束的观测=581/1164；每观测活动比例 median [min,max]=0 [0, 0.814236]。这些是端点位置诊断，不能替代KKT或恢复证书。

## 对“只是预算未收敛 / 噪声过拟合”的证据边界

所有582条clean观测均数值终止；noisy有546/582预算停。因此不能把整个TRAIN/VAL锚点统一称作“预算耗尽”，而 noisy 中确有强烈的预算受限混杂。是否NN对这些组提升、提升是否仅来自预算截断，需结合已授权的分组结果/同预算传统继续求解对照，当前纯公开metadata没有材料质量证据。

noisy最终残差 <.001：582/582；<=.0011：582/582。这是拟合含噪Y的差异量，并不足以证明材料过拟合；没有读取GT、真值电流误差或另一停止规则的材料质量。数值终止也可能是局部平坦/尺度相关，不能从gtol单独断言原因。

## 先前尝试分账（不双加到canonical）

1. 协议纠正前17个已完成discrepancy尝试：solve=148.919414s，process=150.875218s，LU=2011，RHS=109332。它们已保存为可复用 `_discrepancy` 后缀结果，是另一组既有尝试；这里未再次加到1164条严格canonical合计，也不应在baseline重用时再算一次已发生的尝试。

2. 同次协议纠正终止时最多4个inflight：计数/时间 UNKNOWN_NOT_ZERO。初次detached启动失败记0完整输出、setup未知；无法证明它是另一个独立inflight集合，不额外叠加“最多4”。

3. 本地控制器恢复记录的最多2个本地inflight：obs_id及成本 UNKNOWN_NOT_ZERO，另列；i9队列当时持续运行且未重启，不能把其canonical完成任务当作重新执行成本。

真实总尝试成本因中断和setup未知而 **NOT_IDENTIFIABLE**；本报告不给未知成本补零、不用预算上限替代实耗。来源为 `audit/ANCHOR_PROTOCOL_CORRECTION.json`、`audit/data_recovery_20261107T060129Z.json`。

完整的每阶段计数差、状态分布、公开m0活动比例和逐案例表见 `audit/GN_ANCHOR_HEALTH.json`。

## 另备离线梯度诊断（尚未启动）

`experiments/anchor_gradient_diagnostic.py` 已按父线程给定公式准备，状态 **PREPARED / NUMERICAL_DIAGNOSTIC_NOT_RUN**。父线程协调CPU后可用现成 `Gaussian/.venv_nn/bin/python` 执行；默认仅TRAIN/VAL。每观测仅新做1次complex128 LU及12个adjoint RHS，不解forward、不建Full-J、不优化、不做公式或FD核验。输出单独时间戳目录的 `PER_OBSERVATION.csv` 与 `RECEIPT.json`，标为 **OFFLINE_RECOMPUTED_NOT_SOLVER_REPORTED**，记录真实梯度范数、unit-step box gradient mapping、box Frank-Wolfe gap、1e-8活动比例和实际派发LU/RHS/wall。TEST模式要求METHOD_FREEZE且同样不读任何GT。它不改原始状态或anchor回执；本metadata报告中的原回执字段缺口仍保持NOT_RECORDED。
