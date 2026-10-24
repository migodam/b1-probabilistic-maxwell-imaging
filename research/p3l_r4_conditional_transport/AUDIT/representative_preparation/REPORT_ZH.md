# Output3代表解：实现与模拟核查完成，尚未启动物理求解

状态 `INTERFACE_READY_NO_ESTIMATOR_REQUEST`。新模块 `src/p3r4/representative.py`、入口 `tools/run_r4_representatives.py` 已实现；只运行 `--prepare` 和4项decoder/mock测试，没有真实Representative、Full、adjoint或LU调用，没有新RHS，没有GT或网络优化。

## 目标与梯度

给定冻结估计器提供的mu_raw或mu_clipped，求一个原decoder支持内的材料代表：

`min_x (1/576)||chi(x)−mu||_complex²`，`Phi(x)<=282.873792`。

x为完整66维probit坐标，chi仍由原P1 C²饱和decoder产生。材料梯度为 `(2/576) Re[Dchi(x)* (chi−mu)]`，两属性的链式导数乘normal pdf。材料目标和导数只调用decoder；Phi与约束梯度调用原Full nonlinear Maxwell Target，不用Born代理或线性noise-ball。

它是Output3局部candidate representative，不是posterior mean/sample，也不新增prior。mu在原decoder像外是允许的；求解不会偷偷把raw重clip成另一目标。

## 有限预算与可行性

只从现有Full-qualified online anchor起步，初始Phi及backward residual须有限且满足原阈值。优先读取已保存有限probit context；没有有限context时不向内clip整个prior。每parent/arm最多16新states、256complex RHS，调用Full前同时检查forward+gradient费用。

默认每新点forward8+source-adjoint8 RHS；receiver-J32为显式可选路径，此时每新点40 RHS。已付anchor receiver T/J可用纯矩阵乘复用，原获取费用仍保留。所有真实预约由原Book/Target收费，不按一次多RHS函数调用算一次。每个Full trial的材料、prediction、current、field会立即由Target保存；梯度失败或trial被拒绝不抹去费用。

SLSQP自由计算材料objective/Jacobian，只为真正约束查询支付Full。最终只选择原anchor、SLSQP接受且经Full核对Phi严格可行的迭代，或新的、经Full核对接受的可行性restoration点。未接受/infeasible trial绝不直接替换最后可行anchor。

在同一16/256 cap内为最多2个restoration点预留费用；若某个已接受SLSQP点材料目标更好但数据不可行，沿原anchor到该点的预声明比例.999、.5核查端点。不能证明整条路径可行，也不超cap二分救结果。都失败或无可行改进则明确 `fallback_to_original_anchor=true`。source-anchor若还需新gradient8，则可用trial数相应减少，不能借名义16states突破256 RHS。

## 对照与KKT

独立ordinary extra likelihood control也从原anchor起步、同样16/256 cap，只最小化原Phi，不把mu放入proposal/selection。它可每parent执行一次并在raw/clipped两目标中共享，费用有唯一owner，mu距离仅事后计算。

最终输出材料objective、原Phi/white residual、anchor差值、优化status/iterations、accepted次数、实际state/RHS/LU、fallback、restoration及KKT数值诊断。KKT含constraint slack、active flag、非负lambda、stationarity∞和complementarity；这是局部一阶诊断，不是全局min证明。最终梯度使用已缓存付费量，不为报告额外突破预算。

## 准备和运行接口

`--prepare` 无manifest时仅保存 `configs/REPRESENTATIVE_PROTOCOL.json` 及本接口状态，不会启动物理计算。

有实际冻结估计器后，manifest须提供request_id、estimator_id、predictions_npz、显式mu_layout、parent allowlist、outputs、frozen_estimator=true、no_GT=true。只允许读取 `parent` 和 `prediction_raw`/`prediction_clipped`，支持显式complex_N576或realimag_N2_24_24布局。`--prepare --manifest ...` 只冻结mu快照、原anchor及资格，仍为analysis-only。缺失/不合格anchor写skip，不补算。

实际付费需要父线程之后单独明确调用 `--run --prepared-request ... --job ...`；本次没有调用。独立offline GT评价尚未实现/执行，不混入在线solver。

## 核查

4项PASS：原真实P1 decoder material gradient（含两probit属性）与有限差分；mock非线性约束中的可行final/拒绝trial保护；mock硬预算耗尽回退；不合格anchor在任何查询前拒绝。mock tau=1只用于通用控制逻辑，不是新物理target。原物理阈值固定282.873792。

已生成预备analysis-only job凭据 `II/jobs/representative_interface_prepare_v1/CLOSED.json`；没有真实代表解可供offline风险评判。下一步实际估计器/预算/运行由父线程安排。
