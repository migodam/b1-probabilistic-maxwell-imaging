# P1 概率实现独立只读审计

审计日期：2026-10-23。状态：`IMPLEMENTATION_REVIEW_COMPLETE / SCIENTIFIC_DECISION_NOT_ASSIGNED`。本报告由独立审计工作线程完成；父线程负责修复与科学裁决。这里列出的修复均由父线程实施，审计线程只写本目录。

本轮发现两个 P2 实现问题、一个 P3 数值展示问题和一个 P3 费用溯源缺口。父线程已修复；修复前反例与修复后源码/纯数组读回复核均已保留。目前没有由本审计反例支持的未关闭 P2。这里的“关闭”仅指列明的局部实现问题，未执行实际推断，也不表示 P1 科学或数值门槛通过。

## 范围与证据

工作树：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1`。以下代码路径均相对阶段根目录 `research/research_w_p0r_to_p1/`。源文件行号为本次读回的一基行号；“修复前”行号及字面片段单独保留，不能据当前文件行号反推旧版本。

只读检查了 `src/wr/posterior_math.py`、`prob_reference.py`、`progressive_inference.py`、`progressive_offline.py`、`prob_handoff.py`，经父线程明确扩展后检查 `sequential_particles.py`、`prob_bank.py` 的规则函数及 `probability_family.py` 的公开坐标/材料接口；读取配置仅限 `configs/P1_PROBABILITY_PROTOCOL.json` 和 `configs/W_CONTROLLED_PRIOR.json`。未读取结果规则数组、新观测、OFFLINE_TRUTH、GENERATION_PLAN 或任何 GT。源码中出现的真值路径只作为访问控制代码阅读，没有打开对应文件。

纯数组检查通过 AST 只提取函数定义，未导入上述源码模块。未调用 `main`、`smc`、`Family.point`、材料解码器或建立 Family/Book/Forward；PDE、RHS、LU、新观测读取、GT 读取均为 0。初次检查使用 numpy/scipy.special 的合成数组；修复后检查使用 numpy 及检查脚本本地的稳定 log-sum-exp，不依赖项目模块导入。

证据文件：

- `ARRAY_CHECKS.json`：修复前检查结果，保留原始两个反例；不会用修复后的结果覆盖。
- `run_array_checks.py`：当时的 AST 隔离检查脚本。它读取运行时源码，因此不能把当前重跑当作历史版本复现；请勿覆盖历史 JSON。
- `REPAIR_READBACK_CHECKS.json`、`run_repair_checks.py`：父线程修复后的局部纯数组检查，以及公开积分规则/坐标接口检查。
- `REVIEW.json`：发现、源码片段、读回位置、关闭范围及待父线程决定事项。

## 发现及修复读回

### F01 / P2：零权重粒子槽位遮蔽完整漏分支

修复前位置：`src/wr/progressive_inference.py:257–276`，关键旧行 263：

```python
missing=[bool(eligible[branch] and not np.any(b==branch)) for branch in (0,1)]
```

判断只检查分支编号是否出现在粒子槽位，未检查该槽位是否拥有正概率权重。合成反例为两个槽位分别位于分支 0/1，权重 `[1,0]`，精确有限参考分支质量 `[.5,.5]`。旧函数报告实际质量 `[1,0]`、两分支均 eligible，但 `branch_missing=[False,False]`、`branch_MISS=False`。这是漏分支诊断错误；分支质量本身计算正确。不能据该反例断言已经发生真实实验漏报。离线集成质量评分还另有检查，不能将每次运行诊断问题扩大为全部离线评分失效。

当前读回：`progressive_inference.py:263–285` 区分 `branch_particle_slot_present` 与 `branch_positive_weight_present`，eligible 分支没有正权重才 missing；分支内 ESS 先除以最大权重，避免极小权重平方下溢。相同反例现报告 `branch_MISS=True`、missing `[False,True]`、槽存在 `[True,True]`、正权重存在 `[True,False]`、分支 ESS `[1,0]`。把第二分支权重改为 `1e-300` 时 ESS 仍为 `[1,1]`，没有把微小正权重误判为精确零。

关闭范围：`CLOSED_LOCAL_SOURCE_AND_SYNTHETIC_READBACK`。正权重存在不等于有效采样充分，故必须保留已有 ESS、分支质量误差与祖先数，不能单凭 MISS=False 宣称精度合格。

### F02 / P2：离线完成标记采用负向排除，允许缺少完成证明

修复前位置：`src/wr/progressive_offline.py:291–298`。旧检查仅排除显式 `False` 和少数状态：

```python
record.get('inference_frozen') is False
str(record.get('status','')).upper() in ('RUNNING','INCOMPLETE','FAILED','NOT_RUN','PREPARED')
if 'observations' in record and record['observations']!=count: ...
```

因此 `{}` 和 `{'status':'INCOMPLETE_NO_CONVERGENCE_CLAIM','observations':32}` 都通过旧 `_completion_record`。这说明检查函数本身不能建立“明确完成并冻结”；不是实际真值泄漏的证据。旧主流程仍会要求两个标记文件存在、32 个完整观察索引和 parent 一致，审计没有执行该流程或读取任何真值。

当前读回：`progressive_offline.py:291–300` 正向要求 inference 标记 `FULL_NODE_INFERENCE_COMPLETE` 或 point validation 标记 `FULL_POINT_VALIDATION_COMPLETE`，且 `inference_frozen is True`、`truth_access is False`、`observations==count`。inference 对应写出端 `progressive_inference.py:473–474` 写入匹配的 status/冻结/无真值字段。`progressive_offline.py:335–341` 在源码中的真值打开语句之前检查两个标记、32 个索引和 parent 一致。

修复后合成检查：空对象、不完整状态、31 对象、truth_access=True 均拒绝；明确完成的 32 对象 inference 标记和 point-validation 标记均接受。没有调用主流程。point-validation 的实际写出器未在本轮允许清单内，故仅验证其消费契约，不额外宣称该写出器端到端已验证。

关闭范围：`CLOSED_LOCAL_GUARD_AND_INFERENCE_WRITER_READBACK`，实际运行状态与离线访问日志由父线程持有。

### F03 / P3：posterior odds 人为分母下限产生伪有限数值

修复前位置：`src/wr/prob_reference.py:112`：

```python
posterior_odds=ref['branch1']/max(1-ref['branch1'],1e-300)
```

若分支概率在 float 中舍入为 1，旧表达式报告约 `1e300`；它并不是从未舍入分支积分求出的可靠 odds。该问题影响 odds 展示，不改变已经计算的后验质量。`prob_handoff.py:72–78` 原本就把零补集情形单列处理。

当前读回：`prob_reference.py:21–27` 从两个分支的 log-sum-exp 直接计算 log odds，只有可表示时才转为普通 odds，超范围为 null 并注明 `FLOAT_RANGE_EXCEEDED_USE_LOG_ODDS`；输出使用处为 118 行。合成二节点、相同先验、log likelihood `[0,1000]` 得到 branch1=1、log odds=1000、odds=null，状态明确，没有伪造 `1e300`。极小 odds 仍可能在普通 float 下溢为 0，但保留的 log odds 给出可恢复信息；不把该显示值当作精确零事件概率。

关闭范围：`CLOSED_LOCAL_SOURCE_AND_SYNTHETIC_READBACK`。

### F04 / P3：离线费用提取只保留顶层，遗漏分臂费用

修复前位置：`src/wr/progressive_offline.py:301–303`。旧字典推导只扫描顶层键，因 `SMC_runs`、`methods` 本身不含 cost/cpu 等关键词，遗漏其内部每个方法/seed 的 CPU 和 wall 时间。原 inference RESULTS 保存了这些数据，所以这是离线费用来源文件的完整性问题，不是费用被真正删除，也没有据此证成已经发布了错误的等成本优势。

当前读回：`progressive_offline.py:303–315` 递归遍历不含费用标签的字典/列表路径；含费用标签的对象整体保留。合成输入保留 `SMC_runs.direct[0].cpu_seconds=3`、`methods.point.wall_seconds=2`、共享 proposal library cost。当前 `COST_PROVENANCE` 写出处 386–389 明确这些路径是原始来源、不能自动相加，根账本是实际计费总量。

关闭范围：`CLOSED_LOCAL_SOURCE_AND_SYNTHETIC_READBACK`。共享建设费如何分配到方法仍是父线程的比较口径决定，见下文。

## 未发现上述类型错误的数学链路

### 1. Uniform GL 测度与截断 CDF

`progressive_inference.py:55–63` 使用 `u=.01*z`、`c=3+.06*z`，每连续轴的 uniform 先验积分权重均为 `w/2`，再乘分支 `.5`。`posterior_math.py:43–53` 与 `progressive_inference.py:72–75` 对 gain 幅度/相位也各用 `w/2`。这是在明确的独立 uniform(a,phi) 先验下积分，不能额外乘复平面极坐标 Jacobian a，否则改变配置先验。

纯数组质量：650 个材料节点权重和 `0.9999999999999999`，两分支质量各 `.5`；gamma 权重和 1、均值幅度 1、均值相位约 `5.93e-21`。先验位置均值 `-.8499999999999998`，先验平方风险 `.0025333333333333323`，与解析值 `.05²+.01²/3=.002533333333333333` 一致。

父线程扩展授权后补查 `prob_bank.py:55–67`：截断区间上界为 `min(.01,cut-center_b)`；质量乘 `(hi-lo)/.02`，保留原始 prior density，未把 CDF 分子局部归一化为 1。3/5/9 阶各自对常数 integrand 的检查，在五个预注册 cut `[-.909,-.9,-.85,-.8,-.791]` 得到 `[.025,.25,.5,.75,.975]`，最大浮点偏差小于 `3e-16`。区间外 cut 与完整规则也检查为 0/1。本项只验证规则函数生成的公开坐标与权重，没有加载已生成响应数组。

`prob_reference.py` 和 `progressive_inference.py:202–226` 使用独立截断积分分子除以对应完整规则 evidence；有限规则的 atomic CDF 则在 `progressive_inference.py:229–253` 单列。有限 reference 的 qualified 条件 245 行同时包含连续 reference、nuisance、CDF coherence、CDF 误差及 branch 误差，而非只检查 CDF。

### 2. 一个共享 gain，所有 source 只累加一次

`posterior_math.py:70–84` 先累加 F、C、Y，再在 87–138 行形成一次 likelihood 并积分一次 gamma。交叉项 `cos(phi)*Re(C)-sin(phi)*Im(C)` 对应 `Re(exp(i*phi) C)`，符号与 `C=sum(f*conj(y))` 匹配。每实分量 sigma=.01，绝对密度常数为 `-2M log(sigma sqrt(2pi))`；没有按通道数再次归一化。

`progressive_inference.py:78–108` 保留 source 轴但每个材料只有一条共享 gamma 轴；133–143 行 joint prior 只乘一次 material_weight*gamma_weight。纯合成 complex F/y 检查，source panel 对显式残差的最大 log likelihood 差 `1.3642420526593924e-12`；与 posterior_math 的 marginal 差 `1.8189894035458565e-12`；整块与分块统计合并差 0。错误地每个 source 各自积分 gamma 会在该反例改变 log target `5.2227863268108194`，当前实现没有走该错误路径。

### 3. SMC 温度、source 增量、MH 比与重复状态

`sequential_particles.py:53–65` direct 四桥依次将四 source 同时加至 `.25,.5,.75,1`；progressive 四桥使用下三角 source 前缀。二者最终 beta 全 1，增量均由相邻 beta 的差定义。211–220 行每桥只加入该增量，证据增量取当前权重下的似然增量均值。系数 0 显式跳过，避免 `0*(-inf)`。

85–113 行 proposal library 提前生成全部公共随机数；拒绝/重采样不会改变后续随机数消耗。137–145 行的局部索引环绕是对有限索引图的对称提案，不是把物理 likelihood 宣称为周期函数。229–240 行 global proposal 采自 prior，MH 比中 proposal prior 因子与 target prior 抵消，故只保留 likelihood 比；local 对称提案则显式加 `logprior(new)-logprior(old)`。混合权重固定且与状态无关。零 prior 提案拒绝。

使用包含 singleton 维度及非均匀 prior 的 24 节点合成网格，枚举 direct/progressive 共 8 个桥 target 的完整转移矩阵：局部提案对称误差 0；最大 detailed-balance 缺陷 `5.204170427930421e-18`；转移矩阵行和误差 `4.440892098500626e-16`；四桥 likelihood 增量和相对完整 likelihood 差 0。这是枚举核代数检查，实际 `smc` 调用数为 0。

205–208、241–246 行每次 mutation 都记录状态；只替换接受粒子，拒绝和自提案保留重复状态。222–226 行非重采样桥保留 identity 索引，重采样时同步更新祖先；trajectory 保存相应权重。254–265 行没有把不同温度轨迹混成最终 posterior，明确 log normalizer 的 log 不是无偏 log-evidence 估计。未从上述路径发现少加 prior ratio、重复观测增量或删除拒绝样本的问题。

### 4. Point、有限概率、连续参考、handoff 的对象区别

`progressive_inference.py:167–199` direct point 是对全部材料节点做 shared-gain profile 后选最小 Phi；progressive point 是每个 source prefix 的两次材料/gain 交替更新，未证明全局最优。该 point MAP 使用连续 uniform 密度的 profile 目标，不乘 GL 积分权重；不能把它误解为有限原子分布的 MAP。SMC 使用正确的 GL 原子权重；continuous 5/9 reference 与有限枚举另行比较。二者描述对象不同，当前 objective 文本 197 行有明确标注。

`prob_handoff.py:44–69` 对 joint material/gamma posterior 计算依赖的预测，而不是 `E(gamma)*E(F)`。合成两材料/两 gain 的 joint 权重反例中，与显式双重求和差 0；错误因子化的预测差范数 `.1039230484541327`。代表材料是既有付费节点中靠近有限 conditional mean 的点，不能替代 posterior mean 为一个物理世界；源码 143–166、184、236 行区分连续 branch 质量与有限条件预测。

### 5. 离线 prior baseline 与 paired 评分

`progressive_offline.py:136–154` 先验动作位置固定 `-.85`，分支概率 `[.5,.5]`，不从观测估计，不冒充 posterior。qualified continuous Bayes mean 另列；不合格参考保留原值但不计入该 qualified reference 的主评分。184–204 行先在每个 parent 内平均两次 noise；缺一次则完整分数为 null，partial mean 另记。207–221 行 CI 按 parent 配对重采样（2000 次、固定 seed），正 gain 表示 loss_a-loss_b 有利于 b；没有把两个 noise replica 当 32 个独立 parent。本轮没有调用带 truth 参数的评分函数，也未运行 bootstrap 或读 offline 评分结果。

## 费用公平性与待父线程决定事项

1. **Direct/progressive SMC 的计算核对照已有明确相同条件。** `progressive_inference.py:429–448` 二者复用同一 likelihood panel、同一 prior、同一 gamma 规则和同 seed 的固定 proposal library；均四桥、每桥四 sweep、2048 粒子。不同 ESS 路径导致实际重采样/状态接受不同，是已记录行为；没有为 progressive 添加信息。需要区分“相同目标及核预算”与“端到端成本相等”。
2. **共享建设费不能从每臂 CPU 相加直接得到。** reference_cost 418–419 含连续积分、有限 target 构建/精确枚举；point_cost 428 是聚合，methods 内还有分臂成本；SMC_runs 内各自记录 CPU，proposal library 成本共享保存。paid full bank 与 validation 另有费用来源。父线程需指定部署是否需要 reference/完整 panel、共享 bank/setup 如何分配、cache rights、是否比较边际重复使用成本。递归费用修复保留数据，但不自动给出这些决定。
3. **协议的 direct adaptive control 与当前实现覆盖范围。** 配置 `direct_adaptive_control` 要求 full3 及独立 truncated order3，与 5/9 同容差/cache rights 比较。`prob_reference.main` 目前计算完整 3/5/9 并保留比较诊断，而不是一个明确的按观测停止、结算 fallback/setup 的自适应部署策略。源码 gate 仍声明 compression 未 cost-qualified。因此本轮只能确认 direct3 诊断及对照材料存在于代码路径，不能把它们自动计为已完成 adaptive 成本优势。父线程需决定先保持诊断名称，还是另立明确停止/费用规则；审计不改配置。
4. **四臂损失比较的解释。** 两个 point 臂返回单点，SMC 臂以 posterior position mean 决策，且 direct point 是全节点 profile、progressive point 是有限次交替更新。因此四臂差异同时包含决策规则与搜索/近似误差；其中 direct/progressive SMC 对照更直接隔离 schedule。是否支持“多假设优于点”“progressive 有额外收益”等表述由父线程结合实际有限参考误差、paired loss 和费用裁决。本审计没有读新观测或结果，不给此类结论。

## 公开 family 接口与几何边界

配置明示受控新 prior、固定 wall support、两固定 wall 向量、无 cavity、branch 中心间距 `.1`，每 branch target offset uniform `[-.01,.01]`、contrast uniform `[2.94,3.06]`；并声明 `strict_original_prior_membership_claim=false`、original minimum_feature 未定义/未解决、outside family mass 未量化。概率质量来自新明确先验，不从一个 witness 的存在或 likelihood 拟合反推。

`probability_family.py:54–74,228–247` 分支只接受整数 0/1，坐标遵守冻结 bounds，仅允许同一边界外一 float64 ULP 映回本边界；公开材料解码接口不把 gain 乘入介质。纯坐标检查覆盖两个端点、一 ULP 接受、两 ULP 拒绝、bool 分支拒绝、越界 contrast 拒绝。未建立 Family，未读取 wall 数组，也未调用 material/point。材料接口内部调用的实际几何解码器及网格合法性不在本轮执行核验范围；不能以本次坐标检查替代完整几何合法性证明，尤其不能改变原 minimum_feature 语义。

本审计交付的是局部实现审计及父线程修复的有限复核。没有读取 P1 新观测或 GT，没有执行推断，没有签发概率合格/最终科学 GO。
