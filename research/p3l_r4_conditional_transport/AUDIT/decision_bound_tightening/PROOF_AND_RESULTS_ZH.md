# 可选范数收紧：闭账star数组审计，尚非chain收益

状态 `COMPLETE_DEVELOPMENT_AUDIT_NOT_CHAIN_RESULT`。读取封账 `decision_refresh_star_v1` 的8个已付R0、1024个原保存候选材料与残差历史、原付费Full prediction。没有生成新物理候选，没有读取或修改当前运行链的数组/转移日志，没有GT、inverse、LU、线性solve或新RHS。II analysis-only job wall2.631714秒、CPU1.757453秒、RSS105,725,952 bytes、GPU0。

## 直觉与代数

用Frobenius范数上界逆算子的spectral norm会随维数膨胀：I_576的Frobenius为24，而1/∞组合界为1。当前R0接近单位主对角但有非零耦合，已付数组的Frobenius约24.16–24.70，1/∞组合上界约4.42–6.76。

对任意复矩阵B，精确算术下有

`||B||2 <= ||B||F`，`||B||2 <= sqrt(||B||1 ||B||∞)`。

第二式可由 `||B||2² = rho(B*B) <= ||B*B||1 <= ||B*||1 ||B||1 = ||B||∞ ||B||1` 得到。故 `b=min(bF,b1inf)` 仍为上界。

设 `L1=L0−diag(Delta chi)G`，`eta0=||I−BL0||F`，`d=||diag(Delta chi)G||F`。则

`||I−BL1||2 <= eta0+||B||2 d <= eta0+b d = eta`。

若eta<1，Neumann不等式给 `||L1^-1||2 <= b/(1−eta)`。B不必是精确逆；有限方阵的可逆性来自BL1与I的距离小于1。eta≥1只表示该充分条件未建立，不代表L1奇异。对eta0<1、d固定，`b/(1−eta0−bd)`随b递增，因此收紧b不会缩小该充分条件窗口，也不会放大相同残差的误差半径。

同一个已存近似current的残差rho满足 `j_true−j_hat=L1^-1 rho`，故complex输出误差Frobenius不超过 `||S||F * b/(1−eta) * ||rho||F`。白化乘sqrt2/sigma后得到epsilon；若已存白化残差范数为r，Phi位于 `.5[max(0,r−epsilon)²,(r+epsilon)²]`。原prior-reversible full66 pCN使用同一个logU、已知精确base Phi；这里没有改变proposal或接受规则。真正chain还须正确处理当前势的区间、cache基点累计材料变化和history，不能直接从star推广。

## 浮点资格

以上是精确算术不等式。实现用普通float64绝对值、归约和sqrt，没有directed rounding或validated interval库。标记统一为 `CONDITIONAL_EXACT_ARITHMETIC_FLOAT64_NOT_DIRECTED_INTERVAL_CERTIFIED`。0个观察到的区间漏包或错误判断不是全域certificate。

## 保存候选上的结果

8个R0均存在，没有MISSING。b上界收紧因子3.6330–5.4726，中位数4.0089；相应允许的matrix perturbation Frobenius充分阈值 `(1−eta0)/b` 同因子增大。这不是latent半径或真实链cache寿命的实测倍数。

| 同一保存状态 | Frobenius旧界 | min(F,1∞)新界 |
|---|---:|---:|
| eta<1的候选窗口 | 1022/1024 | 1024/1024 |
| Richardson保存历史能决定的候选 | 1022 | 1022 |
| Richardson未决定 | 2 | 2 |
| GMRES已存终点未决定 | 2 | 0 |
| 观察到的错误resolved决定 | 0 | 0 |
| 观察到的Phi区间漏包 | 0 | 0 |

Richardson第0近似状态Phi区间宽度中位数由66.1925变为11.9906；对共同有效候选逐一计算的宽度改善倍数中位数6.0095，最大106.77。39个候选可以在更早的已存refinement状态做出判断。没有把减少的R0矩阵作用列数写成新RHS节省。

原star并非全局0 fallback：R4DEV_00004有2个旧F窗口失败候选。其Richardson历史在iteration0就停下；本次严格只回放已存历史，没有补造后续迭代。因此新界虽让这2个窗口有效，仍不能凭其单个初始近似状态作决定。GMRES已有高精度终点，新界能在这2个已付终点决定；原GMRES全部preconditioner/迭代与Full审计费用仍保留，不能从这次后处理删除。

## 费用和接口

新增1/∞ setup只需一次entrywise abs及row/column sum，共8个数组累计0.012899秒；新增matrix product、factorization、inverse columns、complex RHS均0。为核对eta0，本审计另做原本就需要的8次B L0矩阵乘，累计0.125889秒；不属于新增物理solve。

可选模块 `src/p3r4/decision_refresh_tight.py` 提供 `inverse_norm_upper_bounds`、`neumann_constants`、`phi_interval_from_saved_residual` 和显式opt-in `DecisionRefreshTight`。旧 `decision_refresh.py` 默认未改，mtime/size在本job前后不变。新类只替换bound_constants，iteration和decision mechanics继承原实现；caller须在rebase/B变化后构造新实例。不能把这个开发模块静默接入正在运行链。

两项diagonal闭式toy测试PASS，使用标量公式直接给出参考current/逆范数，未调用inverse、LU、SVD或solve。含旧eta=1.2失败、新eta=.4有效的例子及四个已知Phi闭式包络核查。

原R0获取费用仍是封账star中8×576=4608 RHS；本次没有将其变成免费。下一步是否真正跑独立revision chain由父线程结合当前first chain结果决定；本报告不宣称链fallback率、MCSE、wall或后验正确性已改善。

复核：`RESULTS.json`、`REPLAY_ROWS.jsonl`、各parent SUMMARY、II job `decision_bound_tightening_r1/CLOSED.json`，执行脚本 `tools/audit_r4_decision_bound_tightening.py`。
