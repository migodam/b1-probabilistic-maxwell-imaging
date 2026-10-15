# Joint P2 缓存只读审计（Research B）

2026-10-15；本审计不运行 solver、reference 积分、训练或 acquisition campaign。新增 Maxwell attempts = 0。所有写入仅在本 delegated 目录。未检查 SHA。

## 可直接采用的来源

- 完整数组来源：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-joint-p2/research/joint_p2/results`。
- 源代码 revision：`7ed4faef9bc3dd214d4de689c2203cc383f7495f`，发布 tag `joint-p2-20261014`。
- 本工作区：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-research-b`。当前分支的正式文档可读，但大型 NPZ 在相邻 `b1-joint-p2` 中，不能因为当前 checkout 缺数组就重算旧 campaign。
- 完整 NPZ 键、shape、dtype、文件大小、20 条补偿终点的实际缓存索引在同目录 `CACHE_MAP.json`。
- 可运行的只读例子：`../../.venv_nn/bin/python research/research_b/research/delegated/cache_audit/read_cache_example.py`。只读旧数组和 package 内的已验证 full fields，不构造物理 backend。

## 当前 reference 身份和陷阱

`FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz` 的 target 为：

`joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1`

它对应 **kR=1.6 的原 32 个实通道、连续 coefficient prior 和 proper bounded real gain**。它是该连续插值响应后验的 event/moment/CDF 数值摘要，不是 full-vector-Mie 完整后验粒子，也不是 exact 全材料空间 reference。

实际读取确认：128 observations；64 独立 parents，各两份 noise。125 行 order13，行 12、39、48 升级 order17，另有21阶 cross-check。event mass 归一化最大误差 `6.439293542825908e-15`，保存 CDF 单调性缺陷为0。`parent` 与 observations 一致。

**不能混用的旧字段：** `tools/p2_accept_reference.py:16–19` 只更新这三行的 `CDF`、`conditional_CDF`、`mean`、`SD`、`event_mass` 和 `reference_orders`；accepted NPZ 中的 `partition_log_evidence`、`partition_event_mass` 仍逐元素等于原 `FRESH_CONTINUOUS_REFERENCE_PARTITION.npz`。这两个键是旧 partition 诊断，不能当作三条 upgraded reference 的新 log-evidence 或新 event mass。若需要新的 future-observation likelihood/reference，按原 prior measure 和同一 joint likelihood 单独计算并保存 Research B 输出。

`C1_FINAL_DATA.npz` 中的 `reference_samples`、`reference_log_density`、`reference_event_mass` 仍为原 order13。后者与 accepted event mass 最大差 `0.004960582111149026`。可以将旧样本当作有限 proposal/resampling 支持，但不应把它们改称 refined accepted posterior draws。

`ACCEPTED_REFERENCE_AUDIT.json` 的经验性范围：全部128例 event/mean/SD 13/17核对，三条升级并21阶核对；CDF13/17核对为固定前8加三条升级；direct inverse quantiles仅前8。保存 partition 与原规则最大 log-evidence差 `.0099873311`、event差 `.0055055695`，前8 partition CDF差 `.0009737469`、conditional CDF差 `.0041437708`。这些是数值近似差异，无全支持 uniform certificate。

## 真实 full-Mie bank 与连续 surrogate

`COMMON_RESPONSE_BANK.npz` 包含两频 `1.6`、`2.4`。每频物理坐标 `(central_radius, annular_halfwidth, epsilon_real, epsilon_imag)` 各轴五个 Chebyshev-Lobatto 点，共625个 full-Mie n14 节点。每频原 field 和 coeff 的 shape 均为 `(5,5,5,5,32)`；另有3点 coarse coeff，不能将其当 integration order。

`src/b1_joint/response.py:26 Response` / `predict` 使用纯 NumPy 插值，不导入物理 solver。其 `evaluations` 是 proxy material-frequency evaluation 计数，不是 Maxwell attempt。bank 节点不是离散先验；不能对625节点直接任意等权宣称原连续 posterior。

112个独立 prior probes 的 `INDEPENDENT_FORWARD_ARRAYS.npz` 保存两频真实 full/proxy `(112,64)`。旧数值误差最大 relative field `0.0002985689911131159`，joint whitened discrepancy `0.014870657754625755`。这是有限 probes，不能保证新 action utility/ranking 的均匀精度。原 full-Mie数值资格见 `PHYSICS_QUALIFICATION.json`：七种stress材料两频的 n10/n14/n18 与独立 recursion，n14/n18 保存差为0、recursion最大相对差约4.44e-16；不外推新频率的收敛。

## 最便捷的 candidate-only full-field 输入

七个 `results/packages/*.json` 合计37个代表材料。每个 `hypotheses[j]` 都含：

- `representative_latent`、解码几何、complex epsilon；
- `forward_prediction.unscaled_field`：**已支付且独立 full-Mie n14 验证的双频64维响应**；
- `posterior_gain_mean/SD`：该材料在旧1.6观测条件下的截断gain后验摘要；
- 实际 event 定义和材料/观测距离，但在线 candidate-only policy 不必读 posterior mass。

例：solid、large-cavity 各四个代表，events `[9,0,6,1]`；small-cavity 三个；annular/combined/gain-negative-control 各七个；prior-boundary五个。这些是 event-medoid 材料而非 event 点概率。代表外的模型支持和未发现材料不能删除或赋0。

`packages/*_samples.npz` 提供 whole12-event weighted samples。solid/large-cavity 使用 **full-Mie residual-corrected event masses × proxy event-conditional samples**，其他五例为 proxy reference。原粒子文件还存在 uniform weights 的 quadrature resamples；不要混淆这两层。`solid_FULL_CORRECTION.npz` / `large_cavity_FULL_CORRECTION.npz` 有192/384个材料的 full field，但**仅1.6、shape(N,32)**，没有这些材料的付费2.4批量响应。

`forward_cache/index.json` 2594个参数键，包括 k1.6:1887、k2.4:707；matching2582、recursive12；n14:2570、n10/n18各12。707组参数记录具有双频缓存（包含不同 audit implementation/order，不能全当707个新独立 validated 材料）。部分非cache audit 文件不在index；全局 attempts 仍由旧唯一ledger定义，不能以index大小替代。

## 补偿终点的可用性

20条 `COMPENSATION.json.records` 已按 exact latent decode、epsilon、nmax14、matching 和 k 匹配缓存，逐条文件名见 `CACHE_MAP.json.compensation_endpoints`。

| 材料/优化目标 | 1.6 field | 2.4 field |
|---|---|---|
| small cavity，双频 ε+gain | forward_01302.npz | forward_01303.npz |
| large cavity，单频 ε+gain | forward_01305.npz | 缺失 |
| large cavity，双频 ε+gain | forward_01308.npz | forward_01309.npz |
| annular，双频 ε+gain | forward_01310.npz | forward_01312.npz |
| combined，双频 ε+gain | forward_01317.npz | forward_01318.npz |
| prior boundary，双频 ε+gain | forward_01320.npz | forward_01321.npz |

所以不能把大空腔“单频优化终点”在2.4的 proxy 预测冒称已有 full-Mie endpoint；双频重新优化终点是另一个材料/补偿条件。字段读取为 `np.load(results/'forward_cache'/filename, allow_pickle=False)['data']`。此读取无 solve，直接调用 `MaterialForward.predict` 即使期望cachehit也不属于本审计的安全只读入口。

## 材料、通道与 gain 合同

`src/b1_joint/material.py`：latent `(zc,zs,epsilon_real,epsilon_imag)`，proper uniform product bounds `[-3,1]² × [2.2,2.8] × [.04,.16]`。`geometry` (`:23`) 将 z<-1 解码为 `(.45,.18)*sqrt(1-1/sqrt(-z))`；inactive z∈[-1,1]分别有 prior mass.5。ring中心固定.65R；R=1，rmax=.2925518253、hmax=.1170207301。`decode` (`:63`) 输出同心多层、vacuum epsilon1，其余各层共享 epsilon；passive loss为正。不支持非同心通用3D拓扑；也不能把 coefficient-uniform 改成 radius-uniform。

每频8角度20°至160°、每角两个far-field振幅 `(S1,S2)`，乘 `i/k`（已去已知range phase），先 angle-major flatten16complex，再Re16和Im16。对应 `src/b1/physics.py:82 angular_data`。双频排列是 `[Re16(1.6),Im16(1.6),Re16(2.4),Im16(2.4)]`，不是全双频Re32后全Im32。complex通道j对应实索引j及16+j。角度θ_i的一种偏振p对应j=2i+p。

固定每实通道 SD 分别 `.05044714750236199` 和 `.08064376735208054`；complex总variance为2σ²。噪声锚定固定solid epsilon2.5+.1i，不随候选增大。

gain 是每个世界的**单一实幅度**，proper truncated N(1,.15²)，[.7,1.3]，跨同一测量历史共享。各材料世界允许各自不同合法gain。`conditional_gain` (`likelihood.py:18`) 返回截断前的条件高斯 `(m,v)`，须按 [ .7,1.3 ] 截断后采样/求均值，package的 `endpoint` (`package.py:39`)示例正确。

`log_likelihood` (`likelihood.py:24`) 保留 determinant 和截断normalizer。加入新block后应先 concatenate old/new F、y、sigma，再 marginalize **同一个gain一次**。若从旧 posterior weights 更新，可使用 `logL_joint-logL_old`，不能简单另加用gain先验独立边缘化的新频logL。先验外和数值underflow需显式处理，不能 residual softmax顶替。

## E0 action availability

| Action | 缓存状态 | 必要解释 |
|---|---|---|
| no measurement | AVAILABLE | 无新增likelihood或观测 |
| second frequency2.4 | AVAILABLE | 原C条件仅用1.6；fresh/dev已有独立噪声双频观测及full truth response。若初始条件D则不能再算新增证据 |
| 1.6 repeat | AVAILABLE_FORWARD / NEW_NOISE_SIMULATION | 可复用fullF生成新独立noise；同一episode保留原gain，不得把旧noise原封不动再次作独立观测 |
| 合法角度/偏振子集 | AVAILABLE_RETROSPECTIVE | 完整1.6已用通道不能再添加；先从history删去该subset并重算conditional reference；实/虚索引成对 |
| 参考calibration | STANDARD_RESPONSE_ONLY | `forward_cache/forward_00000.npz`、`forward_00001.npz` 为已知solid2.5+.1i的精确付费双频noise-anchor响应；bank[0,0,2,2]也处于相同物理中点（ei浮点差约1e-17）。没有单独已采集calibration noise观测。可另行声明固定标准、共享gain和有成本的新独立noise模拟协议，不能当成历史真实数据 |
| 新频率1.2/2.0/2.8/3.2 | NOT_AVAILABLE | 不能E0静默调用solver；必须进入独立E1授权账本 |
| 新入射角/通用极化/非同心3D | NOT_AVAILABLE | 当前径向与两个unit polarization far-field通道不等于新任意入射的缓存 |

`FRESH_FINAL_OBSERVATIONS` 的128条noise episodes对每条单独抽gain；同一parent的两条旧noise不是独立物体。对给定episode的未来measurement须保持该条隐藏gain；`truth`、`gain`、futurey和 `F_parent` 属offline simulator/scoring，不得传到onlineactionpolicy。

## 读代码定位与 Research B 注意事项

- `src/b1_joint/reference.py:22 cell` 产生continuous quadrature节点和原prior weights；`:52 bank` 拼12事件、`:58 solve`正确定义C/D。`bank`实际会求两频proxy，但不求Maxwell。使用其weights而非节点均匀权。
- `src/b1_joint/reference.py:108 particles` 是 randomized-stratified quadrature resampling，inactive fibers再按uniform[-1,1]抽样；不能称独立continuous exactposterior draws。
- `src/b1_joint/run.py:137 split_observations` 定义隐藏gain和noise；`:26 physics_bank` 保存dev固定stress数据；不要导入/运行该CLI，因它构造可付费solver。
- `tools/p2_fresh_reference_refinement.py` / `tools/p2_accept_reference.py`给accepted生成逻辑；本审计只读。
- `src/b1_joint/package.py:39 endpoint` 是会调用solver的函数，**不可用于本轮零solver读取**。读JSON现成fields即可。

没有阻塞读取E0的资产缺失。限制是：连续参考主要为proxy；全Mie残差修正仅两例且仅原频；新频率无缓存；A真实材料包尚未由本审计确认。所有后续真实action排名、风险收益、closedloop和新物理仍由主线程执行并单独记账；本审计没有声称完成它们。
