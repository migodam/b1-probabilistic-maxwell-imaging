# 原输入静态审计（未运行物理实验）

审计日期：2026-10-29。根目录简写 K=`inputs/Research_Z_ALG_Research_Kit`，P=`research/research_z_alg_asymmetric_validation`。以下路径均相对 P。所有数值为交付包已有记录，未在本次审计重跑。未做 SHA256 检查。

## 核心判断和适用范围

原包实际新增的物理证据限于真实 Hankel 系数、圆盘 pulse-basis 中心配点离散 scalar-TM 的两 cell 交换阳性对照及另一个源合同中的 dark-current 反例。E4 为同一对称两 cell 模型下的数值概率一致性阳性对照。非对称 4–9 cell 图册、非对称远处分支、严格全覆盖和普遍加速均尚未建立；不能把包内 CONDITIONAL GO 转述为本轮非对称验证已通过。

原包的 `symbolic/certified` 输出措辞强于其浮点认证能力，README 第五节主动澄清为条件解析结构与浮点见证：真实 Hankel 系数没有区间包围；非零浮点 determinant 不是严格 rank 证书。报告源 `03_maxwell.tex:138–177`、`06_experiments.tex:28–30` 支持此限定。

## 逐项证据表

|ID|原 claim ledger 精确陈述|可定位支持|包内结果与资格|本轮状态|
|---|---|---|---|---|
|C1|Useful fixed-label C2 material symmetry exists in the frozen two-cell TM/VIE.|K/CLAIM_LEDGER.json C1；K/report_source/sections/03_maxwell.tex:1–19,116–166；K/zalgo_maxwell/two_cell_witness.py:13–44；对应 JSON|材料 (1,3) 与 (3,1)，相对预测差 1.2215904044358529e-16；材料距离 2.8284271247461903。条件恒等式 + 已有浮点见证；非连续 PDE 认证。|静态核验，未重跑|
|C2|Two distinct regular branches obtained via symmetric p,q coordinates.|K/CLAIM_LEDGER.json C2；K/report_source/sections/03_maxwell.tex:138–166；two_cell_witness.py:45–55；对应 JSON|p≈4,q≈3，根 3、0.9999999999999992；quotient det=-0.11846655750648813；state sigma_min=0.8087481239013333。完整两根结论依赖商系统可逆、根合法、resolvent 非零等条件，无区间 Hankel 证书。|静态核验，未重跑|
|C3|Receiver-dark current may lack any common scalar material across sources.|K/CLAIM_LEDGER.json C3；K/report_source/sections/03_maxwell.tex:172–177；two_cell_witness.py:63–80；对应 JSON|另设源 (0,-0.7),(0.1,-0.9)；预测差 3.5598230403955905e-17，但 shared-material residual=0.006351035627535262；两 cell normalized minors≈0.00449824,0.00448347。LS恢复材料虚部≈0.1062,-0.06494，不在原实域。|静态核验；不能与 C1 同合同混用|
|C4|Original posterior and correct quotient/fiber integration agree in finite-noise E4.|K/CLAIM_LEDGER.json C4；K/zalgo_probability/e4_two_cell_posterior.py:65–75,118–190；对应 JSON；06_experiments.tex:1–30|[0.5,4]^2，rho=x1/27.5625，8个实分量 sigma=0.001；64阶 Z=0.15957386627129155；两方法 evidence 相对差5.218e-16；32/64差最大证据2.188e-9。数值参考，无严格积分误差界。|已有结果，不是本轮独立复现|
|C5|Uniform branch weights can preserve evidence yet give wrong material tasks.|同 C4，脚本161–169；JSON results.64|正确均值≈(2.386598,1.721835)，事件0.67516447，x1²期望6.254154；错误均匀分支均值≈(2.054217,2.054217)，事件0.5，二阶矩4.914982，证据不变。|已有明确负对照|
|C6|Full66 locally full rank forbids a smaller smooth lossless local forward quotient.|K/CLAIM_LEDGER.json C6；02_geometry.tex；01_framework.tex:19–26；K/SOURCE_MANIFEST.json|标准 rank 链式法则的条件推论；历史数值 rank 非数学满秩证明。冻结 R7 报告被引用但历史大实验未重跑。|来源有限，非新 Full66 结论|
|C7|Full576 real material prior has regular local fibers of dimension at least64 with512 real observations.|K/CLAIM_LEDGER.json C7；02_geometry.tex:174–180|新支持/先验、合法内点、常秩及 decoder Jacobian 条件下的维数定理；未执行新576实验。|条件理论，非物理结果|
|C8|General nonsymmetric Maxwell branch coverage and speedup.|K/CLAIM_LEDGER.json C8；06_experiments.tex:61–83,128|明确 NOT ESTABLISHED；本轮应独立检验。|NOT_ESTABLISHED|
|C9|Exact posterior disintegration, coarea and error/omission bounds.|K/CLAIM_LEDGER.json C9；04_probability.tex|标准理论应用与推导，非原创数学主张；材料 measure 与 recovery Jacobian 不可省略。|理论审计由主线程判断|
|C10|Periodic retracing target path yields identity under unique regular lifting.|K/CLAIM_LEDGER.json C10；06_experiments.tex:106起|定理具正则唯一提升前提；Newton 跳支/奇点不能据此称 regular monodromy 或完整覆盖。|理论审计由主线程判断|
|E0结果|cubic roots / branch weights / QD counterexample / complex monodromy|K/zalg_algebra_e0.py, .json；K/zalgo_probability/analytic_probability_checks.py, analytic_checks.json|QD独立拼接均值2.081368981；cubic权重1/4,1/2,1/4；两复闭路每条800步，residual≈9.687e-14。|原脚本不等于协议要求完整 E0-1至6全部自动检查|

## 来源范围与未交付项

SOURCE_MANIFEST 共11个冻结 R6/R7 文本 URL，均指向 GitHub raw 与明确 commit：R6 `d1bdc8aceae4cd6a843f873b61b019d5eee3cd1d` 四项；R7 `975ebe7f17c2ba837a8cac0a5f8eb4ccf7d4e19b` 七项。这些历史文本本身不随本包对应 `zalgo_sources/` 文件夹交付，因此此审计仅确认 manifest 声明，未访问网络/未核验原报告的当前可得性。

INPUT_AUDIT 还记录旧 Research Z zip 的 library-file 身份及读取范围。原 zip 不在本交付包内；PairTSOM/F/F3 原实验包未独立完整检查，这是原作者自报限制。不能把源 manifest 当作原代码和数据已全部复现。原报告 PDF 与完整 LaTeX 源均已提供；本审计读原脚本、结果 JSON、claim ledger、manifest、协议和相关源段落，不声称逐页视觉复核 PDF。

## 用户任务与 supplied protocol 差异

|事项|原协议|本轮用户/主线程约束|执行含义|
|---|---|---|---|
|优先顺序|E0→两cell→非对称，涵盖 E0–E4 与跨表示、monodromy 等广泛候选|本轮聚焦非对称4–9 cell；原对称例仅回归阳性对照|不得用旧对称成功充当本轮方法贡献；低价值扩展不能挤占核心预算|
|成本|每 run 10k actual complex RHS，2CPU h，8GB；总30k，原文在“进入扩展实验后”措辞留下初始控制记账歧义|总30k、每run10k为硬上限；所有物理评估先纳入统一账本|两cell解析公式和回归验证也先审批预算登记，不能以原包毫秒耗时当免费|
|先验/噪声/几何|新完整域须另设 P_Z；允许另建新采集 contract ID|预冻结后不得为挽救结论修改 prior/noise/geometry|不同初始子问题必须事前声明；后验改合同不能救原合同 PASS|
|检查文件身份|协议要求 SHA256 / preregistration hash|用户既有偏好及主线程本轮指令不加 SHA256 检查|保留原 metadata 字段为输入事实，不运行额外 hash gate|
|复现脚本输出|README命令会原地写结果JSON|inputs只读，本worker只写 audit/input_evidence|任何后续执行先复制到输出目录并用显式 output；勿直接照 README 运行|
|科学结论|最终STOP / CONDITIONAL GO / GO，几何/概率/费用分开|父线程保留科学方向、核心理论和最终判断|本审计不独立宣布最终GO|

## 运行前费用估算

### two_cell_witness.py

静态逐调用计数：8次 `forward`（基点、交换、2变量中心差分共4、dark基点、LS投影后验证）×2源 =16 complex RHS；解析 sensitivity `solve(A,diag(field))` 两列=2 RHS；`inv(A)` 等价求解两列 identity=2 RHS。合计20 actual complex linear RHS；预期10次复杂2×2分解（8 forward+1 sensitivity+1 inv），另1次实2×2 quotient求解。两项独立 disk Green quadratures 不含 Maxwell RHS，但有未计次数的 adaptive Hankel kernel evaluations；嵌套cross quadrature与CPU/峰值内存未instrument。原结果JSON没有保存完整 J/J_swap/全部 prediction 数组，故不满足新协议“全部端点/current/predictions可审计”保存要求；重跑应增加输出而非修改物理。

### e4_two_cell_posterior.py 默认 --n 32

层级 n,2n；box每层2n²（全box及三角事件），quotient每层n²，因此 likelihood=15n²=15,360。响应检查额外 direct5+quotient5=10；共15,370解析物理响应，其中 direct=10n²+5=10,245、quotient=5n²+5=5,125。另数值 forward 5次×2源=10 actual complex RHS / 5分解。没有 sensitivity。

**必须保留两个费用口径：** actual linear-solver RHS=10；解析 forward RHS-equivalent=2×15,370=30,740，合计 conservative physical source-RHS-equivalent=30,750。默认原版若按此保守物理评价预算计，超过整个30k且单run10k。不能直接运行默认E4再称只花10RHS。若中央预算仅按actual solver RHS计，也必须独立报告15370物理解析响应及其CPU费用，不能把其视为免费。

一般n的保守计费=30n²+30，故单run10k最多n=18（9,750）；n=19需10,860，超过。较低积分阶数若未达到原冻结参考精度，应标 reference insufficient，不得放宽门槛。原JSON记录 CPU≈0.01107s/wall≈0.01108s仅是他平台历史结果，不能当本机先验保证或非对称加速证据。当前未执行任何物理正演/大积分。

### 其他脚本

`zalg_algebra_e0.py`为1600次cubic多项式求根，无Maxwell RHS；`analytic_probability_checks.py`为toy一维及嵌套积分，无Maxwell RHS。它们不覆盖协议全部E0项目。figure脚本依README只读现有结果绘图；本worker未执行。所有非物理工作仍有CPU/内存，未视为实验免费。

## n16/32适配与已有收敛信息补充

原JSON只报告32/64，不能从其结果严格推出16/32会通过。box的32/64证据相对差2.783e-15、事件差1.087e-9；quotient证据差2.188e-9、均值最大2.124e-9、事件差3.758e-10、非线性任务差1.213e-8。64阶两方法几乎机器精度一致，这支持32已相当准确，但没有16阶结果，因此16阶是否足够必须标 UNKNOWN_UNTIL_RUN。原脚本默认gate evidence/means/event=1e-7、nonlinear=1e-6；原协议下一轮容差较宽但必须事前冻结，不能看16阶结果后改变。

默认n改为16时，两层16/32总 likelihood=3840，analytic responses=3850，plus numeric RHS10，保守合计7710；再加原two-cell20得到7730。若后续只补**单层40**，三组积分共4800物理响应，2源=9600 equivalent（若重新数值QA再加30）；但直接调用原脚本 --n40会自动跑40/80，代价48030，超预算。适配应明确单层模式并复用已登记检查结果，所有额外费用仍记账。

原E4脚本13、194、195行导入/计算SHA256。父线程已决定使用只读原件的适配副本，身份字段改为BLAKE2b，不能原地运行原SHA路径；本审计没有计算任何hash。
