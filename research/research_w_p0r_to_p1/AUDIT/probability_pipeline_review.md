# 概率流水线独立代码审查（零 PDE）

日期：2026-10-23。范围：`src/wr/prob_data.py`、`prob_bank.py`、`prob_reference.py`、`configs/P1_PROBABILITY_PROTOCOL.json`；另核对已有 Family/概率代数/账本接口和授权数额。审查期间主线程修复了两项已即时报告的问题，下文区分已复核修复与仍开放项。

本审查没有修改主源码，没有调用 `Family.freeze/load`、`Book`、`Bank`、`Family.point` 或任何 `main`，没有打开真实观测或 truth 数组，没有生成真实观测。小检查仅通过 AST 提取纯函数，并输入固定合成数组。唯一新交付文件是本报告。新增 PDE、RHS、LU 均为零；本报告不授予 P1 科学资格。

## 已即时报告并复核的修复

### 1. 截断积分可以得到不合法 CDF，却原来仍获资格

初读 `prob_reference.main` 只检查 order5/order9 的 CDF 差，不检查范围或单调性。正 quadrature 并不自动使独立截断节点的分子与完整规则分母构成合法 CDF。例如完整规则上 logL=0，而 prior mass 为0.25的截断规则上 logL=log10，则分子/完整Z=2.5；两个阶数都取2.5时，阶数差为0。即使同一函数在不同节点上取值，这种有限节点混叠也可能发生，不能仅靠相邻阶数差排除。

主线程已加入 `cdf_coherence` 和 reference/compressed 各自的 coherence 标记，并将完整 reference 与 full3/compression 的 nuisance 资格分开。AST 提取复核：`[0.2,2.5,0.8]` 和 `[0.2,0.4,0.3]` 均拒绝；不再由近似方法单独的 nuisance 失败拖垮独立 reference。修复已观察到，未运行真实 reference。

### 2. logsumexp 之前先相加会丢掉 prior 权重

原 `posterior` 先计算 `log(weights)+ll`。当 ll 带巨大共同偏移时，prior 的对数体积因子在相加时已经被舍入抹掉，后面的 logsumexp 无法恢复。实际纯函数反例：full9的162个合法 prior 节点、常数 `ll=−1e20`，原实现返回权重和162、branch1=81；`ll=1769` 时权重和偏离1约1.1e−14。

主线程已改为先做 `ll−anchor`，再加 `log(weights)`，并明确归一化后验权重、检查似然。修改后对四点非均匀 prior `[.1,.2,.3,.4]` 与同样巨大偏移复核，权重和 `1.0000000000000002`、branch1=0.7，符合该 prior。此反例是通用代数稳定性问题；没有证据表明真实 Maxwell 观测已产生 `−1e20` 量级的 log likelihood。

## 仍开放的具体问题

### [P2] 重启后的“新增向量”与缓存命中费用字段不准确

`prob_data.main` 末尾无条件写 `new_full_vectors=count`。在全部16个 parent 的 `COMPLETE.npz` 已存在时，循环不会调用一次 `Family.point`，但完成记录仍声明新增16个 full vectors。这是明确的增量计数错误；建议在实际新调用/成功保存处分开累加 `new_full_vectors` 与 `reused_full_vectors`，保留总 parent 数。

`Bank.point` 缓存命中直接返回历史 row；该 row 的 `incremental_rhs=20` 是原创建时的费用，命中时没有新的预约，却仍返回同一字段。当前 `Bank.main` 只使用其中的 path，未据此再次计费，因此**没有发现账本重复收费或重置**；但该接口会误导之后的成本比较。应区分 `original_evaluation_rhs=20` 和本次 `new_charged_rhs=0/cache_hit=true`。

### [P2] 缓存只绑定坐标，没有验证公开墙向量与获取配置

当前 key 是 `(int(b),u.hex(),c.hex())`，在**同一个不变的公开 family** 内，这个节点身份是充分的；γ 已在外部积分，不应进入 raw F 的 key。相同节点可以合法复用，stage不同也不意味着必须重付已有 F。

但 `INDEX.jsonl` 和缓存命中没有将当前 q0/q1、公开 prior/acquisition/numerics 与原响应绑定或核对；已经完成的 rule NPZ 还会因存在而直接跳过整个生成过程。因此，若公开输入被错误替换、复制了另一运行的缓存，或响应实现修复后沿用旧目录，当前代码会无声接收旧 F。Family 的正式 freeze 拒绝覆盖降低了正常路径下的风险，但缓存自身没有 fail-closed 身份核验。

建议保存并精确比较小型公开 family/config snapshot 与 q0/q1 数组（无须 SHA），缓存命中再核对 `FULL_RESPONSE.npz` 的 branch/u/contrast、完整 `(5,4,12)` shape/finite 与端点的公开参数元数据。整条 rule 的完成标记也需要同一公开 family 身份，不能仅由文件名和存在性决定。这里核对的是公开 prior 节点元数据，不是访问 GT。

此外，`Bank.point` 在校验前就 `int(b)`，会把不合法 branch=1.7 或 True 转成合法1，绕过 Family 的严格 branch 检查。应先调用 `family.coordinates(b,u,c)`，用返回的规范坐标构造 key；这也会把允许的单 ULP 边界回移映射到同一个缓存身份，避免重复付费。冻结规则自己生成的 branch 都是0/1，所以未发现当前这批 rule 已触发该错误。

### [P2] 数据恢复只核对计划参数，不核对实际恢复数组身份

`prob_data` 对已有 plan 只比较 parameters，没有比较保存的 noise/seed/shape；对 parent `COMPLETE.npz` 则直接读取 F/y，不核对该文件保存的 branch/u/c/a/phi 是否匹配当前 plan，也没有 shape、finite 或公开 family 身份检查。随后 `OFFLINE_TRUTH.parameters` 使用**重新构造的 plan 标签**，而 F/y 可能来自另一个完整文件。文件被误放/替换时，会把错误标签与旧观测配对，而不会立即拒绝。

正常同一冻结协议、同一 RNG 版本、文件完好时，重新构造的 plan 是确定的，当前逻辑不会主动更换 parent 或按结果挑选案例；本项不是已经观察到真实错配。建议从已有 plan 读取并核验完整噪声与标签，恢复 parent 时核对已存标签、形状和公开 family 标识。损坏/部分写入的 NPZ 当前通常会报错停下，不会无账本地补做，这是正确的保守行为。

### [P3] 某些诊断证据没有完整留存

compressed CDF 的 `marginal` 返回 convergence history，但当前调用只保留 ok、丢弃 history。于是 `compressed_nuisance_qualified=false` 时，无法从结果看到究竟哪个 cut、哪个阶数未收敛。应将每个 compressed-CDF history 写入 nuisance diagnostics。reference 的全部完整/截断规则 history 已有保留，未发现同类缺口。

另外，`marginal` 的 `[11,21,41,65]`、`2e−4` 与事件/CDF的 `.02` 当前和冻结协议一致；这不是本轮错误，但函数未读取对应配置字段，未来不可只改配置就声称实际阈值已经改变。

## 核验通过的代数与边界

### 正规 prior 因子与截断体积

`node_rule` 的质量因子

`0.5 × (w_u/2) × ((hi−lo)/0.02) × (w_c/2)`

正确包含 branch质量0.5、原始 u 的均匀密度1/0.02和 contrast 密度1/0.12。截断 u 区间的体积保留，未重新归一化为1。

对 order3/5/9，完整 prior 质量分别约 `0.9999999999999999/0.9999999999999998/0.9999999999999998`。在常数 likelihood 下，冻结五个 cut 的 prior/CDF 分别为：

|cut|−0.909|−0.9|−0.85|−0.8|−0.791|
|---|---:|---:|---:|---:|---:|
|正确累计质量|0.025|0.25|0.5|0.75|0.975|

三种阶数都与这些值相符到浮点误差。`cdf_*` 分子减去对应 `full_5/full_9` 的 logZ，使用了完整分母，没有用截断规则自己的 logZ 将条件概率误归一化成1。中央空隙 cut=−0.85 对应 branch0 总质量；后续可额外记录其与 `1−branch1` 的一致性作为交叉诊断。

### 共享 γ、噪声与收敛

`prob_data` 对每个独立 parent 只抽一个 a/phi，构造一个共享 γ 乘完整 `(5,4,12)` F；两个 noise replicas 共用该 parent/nuisance，并各自加入独立 Re/Im 标准差0.01的噪声。它们是32次观测、16个独立 parent，不能在离线覆盖率置信区间中当作32个完全独立 parent。代码已有 parent_ids 和独立 parent 数，当前没有发现生成器按频率或按源另抽 γ。

`prob_reference.marginal` 一次把完整 F 传给既有共享 γ 积分器，未 profile、未独立源 nuisance、未增加噪声。11→21→41→65 的每个节点最大 logZ 变化阈值与冻结协议一致，达到上限仍不合格则保留未合格标记。

只用合成 `f=10`、`y=1.003 exp(.004i) f` 的完整形状做纯函数核验，阶数变化为约 `172.2370、3.68309、3.13727`，65阶仍正确返回 ok=false。这是有意的高集中度负对照，不是物理实验失败，也不能说明真实数据必须使用某个更高阶。相邻阶数一致仍是经验 quadrature 检查，不能宣称严格积分误差界。

### 节点 tuple 与插值排列

完整节点排列是 branch→u→contrast；`interpolate` 按分支筛选后 reshape 成 `(nu,nc,5,4,12)`，与 `node_rule` 的循环顺序一致。以分支分离的双二次复数响应做合成检查，order3响应到order9节点插值的最大误差约 `3.55e−15`，没有发现 branch混合或 u/c 转置。

响应插值的外推、surrogate本身误差和 compression 成本资格仍需真实验证；本检查只验证排列与多项式代数。

### GT 隔离

`prob_reference` 只打开 `OBSERVATIONS.npz` 与公开 rule bank，不读取 `OFFLINE_TRUTH.npz`、`OFFLINE_GENERATION_PLAN.npz` 或 parent `COMPLETE.npz`。`prob_bank` 不读取观测或离线标签。公开文件 y/parent_ids/noise_ids 与离线 labels/F 分开保存。没有发现这三份源码中的在线真值读取路径；这种隔离依赖进程与代码合同，并非文件权限沙箱。

## 数组、预算与重启

冻结规则总计972行 F，其中按现有精确 tuple key 去重为646个公开节点。若全部新算：bank 为12920 complex RHS，16个 parent 另需320 RHS，总计13240 RHS；材料频率 LU 为3310。它们低于当前 P1A 的24000 RHS/6000 LU 上限；CPU能否在28800秒内完成，不能由此静态核验保证。

所有规则的 raw F 合计约3,732,480字节；最大 nuisance 面板为 `162×65²=684450` 节点。每个 observation 顺序处理，未发现 reference 代数阶段必然突破内存的数组规模；这不是完整进程 RSS 上界。物理阶段仍须依赖既有 Forward 的单频流式保存、110k DOF/3.5GB factor gate 和外部6GB进程限制。

账本使用同一 `PHYSICS_LEDGER.jsonl` 追加预约；成功缓存命中不预约，新失败尝试保留独立目录并继续使用原账本。没有发现这三个模块在重启时清空账本或重置 RHS/LU 总额。CPU累计是否完整仍依赖外部 launcher 的过程费用收据；Book 自身的进程计时起点重建不能被解释成总 CPU 归零。上述“新增向量/缓存 incremental 字段”问题影响报告口径，不改变这个账本事实。

## 科学范围

本流水线只对应新冻结的受控两墙、固定解析 support、连续 target offset/contrast 与单一共享 nuisance 的 prior。原宽材料类的 `minimum_feature=0.1` 语义仍未解决；新 family 的概率质量不能转写为原未知墙材料 prior 的质量，也不能把 outside-family 质量写成0。当前代码将该项保留为 null/unquantified，这一边界正确。

本次交付是代码与合成代数审查。真实 reference 资格、compression公平成本、P1B结果与离线校准评估仍由主线程执行和判定。
