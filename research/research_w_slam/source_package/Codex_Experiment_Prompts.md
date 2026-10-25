# Research W-SLAM 分阶段 Codex 实验任务

版本：2026-10-24 UTC。配套报告：`Research_W_SLAM_Report_ZH.md`。

本文件可独立移交给执行者，不依赖此前聊天。以下任务是未来实验协议，不是已取得的 Maxwell 成果。当前研究仅完成文献与仓库审计、条件数学推导和两个小型 NumPy 代数核验；不得把附带脚本通过写成物理机制通过。

使用方法：先给执行者本节“统一总约束”和对应的一个阶段 Prompt。各阶段在完成输出和 gate 判定后结束，后续阶段只在所列前置条件与预算合同具备时执行。未定义的任务、硬件和 prior 字段必须在查看测试结果前确定；不能为了当前结果更好而事后调整。

## 统一总约束

### 研究目标与证据边界

目标是研究未知墙体与目标的联合 full-wave 反演：多次观测何时打破墙–目标补偿，何时可以暂缓计算，何时必须重开历史墙体，以及在相同概率/任务精度下是否能减少 Maxwell 计算成本。当前优先主线是动态补偿与有限可辨识性；任务证书为小型应用验证。一般在线系统、生成模型训练和三维扩展均不是默认任务。

必须读取并锁定以下历史版本；它们是证据，不是可随意改写的工作目录：

- https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/958a478739d9553269316ded8125fda9df225c0c/research/research_w_p0r_to_p1
- https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/483d260c6b5b111db9fb814b096685fccb362448/research/research_w_continuation
- 配套的 `Research_W_SLAM_Report_ZH.md`、`code/observability_toy.py`、`code/freezing_toy.py` 和各自 `results/*.json`。

用户旧 IS-SLAM 已经包含主角度–广义保留谱、shared-compensation 核、Schur 信息单调、保留率非单调、motion prior 与 rank-stratified 跟踪。新工作不能把这些基础关系重新标为原创。Research W 已有受控联合补偿见证，但 measurement-only 一般未知支撑相容锚点尚未取得；本轮必须独立过这一门。

### 统一执行规则

1. 在隔离研究分支/工作树上工作，建议新目录 `research/research_w_slam/`。不覆盖历史证据，不修改旧结果以追认新合同。保存 Git commit、环境、依赖、设备和随机种子。
2. 首先建立版本化预注册文件。明确采集、几何、材料、噪声、任务、prior、优化与推断表示、成本、指标、接受规则和停止规则。未定字段先记 `null`；依赖该字段的物理运行不得开始，其他可完成的审计照常推进。
3. 第一物理平台为二维 scalar-TM 全波多次散射模型与已合格的贴体高阶 FEM。它不是三维 full-vector Maxwell 认证。新材料、拓扑、频率、姿态和端口仍需自己的数值资格。
4. `G_GRID_0125_V1` 是新离散合同：公开墙体候选 ROI 内边长 `0.125 lambda0` 的宏单元并；禁止2×2对角接触模式；固定一致的4连通材料/空气规则；另行冻结边界接触、材料跳变和墙–目标间距条款。它不等于旧 `minimum_feature=0.1`，也不证明连续几何所有方向厚度至少0.1。旧真值仅作成员检测，未通过或不可判定时保持旧标签，并独立生成新的合法真值集合。
5. ROI 是公共允许区域，不是 GT 材料支撑。主优化/采样进程不得读取 GT 支撑、孔洞数、真实材料、真实目标位置、GT 附近初值或真值误差。真值只由隔离评价进程读取。机制实验可以使用已知局部基点，但必须单独标记，不能替代 measurement-only 主测试。
6. 材料 prior 必须 proper、被动、频率一致，并有明确概率测度。几何变维、occupancy、材料条件分布、坐标变换 Jacobian 与 nuisance 相关性必须显式。算法 LM 阻尼不是物理 prior。
7. 冻结噪声模型与实化/白化约定。首轮新合同默认每实/虚分量标准差 `0.01`，必须与下节源归一化、gain 和数据单位共同冻结；它不意味着未定义源尺度时实验难度已确定。不能因不相容而事后增大噪声。白化均值距离不除以数据维数。
8. 不把两个无噪声均值的距离 `d <= 1` 当作高维含噪数据拟合门。主接受规则用预注册 posterior predictive、独立留出预测或经 prior-predictive 仿真校准的 discrepancy statistic。拟合后自由度不能机械设为原数据维数。
9. 每份真实观测只注册一个 `factor_id`。同数据优化、重线性化、提高精度和重采样只替换同一因子的计算版本，不增加 likelihood 指数。SMC tempering 仅使用增量幂；MCMC mutation 必须保持当前目标分布。已有 posterior 作为递推 prior 后，不得再次乘一遍旧数据。 若旧近似likelihood已用于粒子/分支权重，刷新后要乘新旧likelihood比并正规化，或用以修正后全目标为不变分布的校正核；Gaussian系统移除旧因子贡献再加入新贡献。硬删为零的支持不能只靠比值恢复。
10. 明确区分同数据多算、独立噪声重复、新频率、新收发位置/孔径。新采集必须严格追加并保留旧通道；旧非嵌套窄/宽带比较不构成追加信息单调性的证据。
11. 保留 joint covariance、历史原始 factors 和模式质量。SOFT-FREEZE 是计算调度，不把材料变成精确已知常数；重开旧状态需要可恢复原始非线性因子或充分的精确边缘信息。
12. 分支概率必须由合法 prior、proposal 和 likelihood 计算；禁止把 residual softmax、优化命中次数或任意温度评分叫 evidence。Best-of-K 仅为补充诊断，不能替代在线 Bayes action。
13. 局部 Gaussian 证书只保证其声明的 surrogate。真实 posterior 风险还需 forward 数值误差、非线性余项、模式覆盖和尾部预算。有限 rank 或少数随机探针不自动提供全空间严格稳定下界。
14. 对所有适用方法开放相同 warm start、缓存和普通预条件器。强基线至少包含合格 FWI、batch joint、普通 incremental smoother 和可用的小型 direct posterior。不可只比较 hard freeze 与复杂新算法。
15. 成本统一记录全波 RHS、FE/线性算子分解、CPU、wall-clock、峰值内存及 forward/adjoint 分类；包括失败起点、拒绝 proposals、分支发现、设计计算、历史因子刷新、独立审计和 fallback。不能只报迭代次数。
16. 主要加速比较以相同 CDF/事件/模式质量与任务后悔精度为前提，报告未达门比例，不删除失败世界。reference 可有单独离线预算，但新方法在线依赖它时必须计入部署成本。
17. 所有 gate 使用 `PASS`、`FAIL` 或 `INCONCLUSIVE`，给出具体证据与剩余资格。没有运行的项目标 `NOT_RUN`。数学真、文献新、物理可行和算法更便宜分别判定。
18. 证据标签固定为：S 标准结果，U 用户旧结果，D 本报告条件证明，C 待证猜想，N 实际数值事实，H 未运行算法假设。数值未成功不能用更宽叙述掩盖，局部优化未找到也不能称为全局不存在。

### 可直接实例化的新实验默认值

以下值同步主报告 E0c，是独立新实验候选合同，不是旧补偿见证参数或已通过场景。执行前可以统一修改并重新冻结，之后不可根据算法表现调整。

| 项目 | 默认值 |
|---|---|
| 背景与尺度 | 空气，`mu_r=1`，`f0=1 GHz`，`lambda0=c0/f0` |
| 公共墙 ROI | `[-1.5,1.5] × [0.5,1.25] lambda0`，24×6宏格，全部实际支撑未知 |
| 小 reference 区域 | 公共 ROI 分成2×2或3×3矩形区域，均为宏格并；枚举候选后筛几何合法 |
| 占据 prior | 条件于几何合法性的 Bernoulli(0.65) 区域 prior；小类显式归一化，宽类核查合法采样 |
| 材料 prior | `epsilon_inf ~ Uniform[1.5,6]`，`ell=sigma/(epsilon0 omega0) ~ Uniform[0,0.25]`；先共享两参数，再单独增加区域差异 |
| 起步目标模型 | 半径 `0.12 lambda0` 圆盘，已知实对比度 `chi=epsilon_r-1=3`；未知中心在 `[-1,1] × [1.6,2.6] lambda0` 均匀分布 |
| 未知目标材料扩展 | 半径 `Uniform[0.08,0.18] lambda0`，实对比度 `Uniform[2,5]`；损耗和拓扑另作独立扩展 |
| Tx / Rx | 基础 Tx 为 `(-1.5,0),(1.5,0) lambda0`；16个 Rx 在 `x∈[-2,2] lambda0,y=0` 均匀分布；核验无源接收重合 |
| 嵌套新增动作 | 基础 `{f0}`，可追加 `0.85f0,1.15f0`；扩孔径增加 Tx `(±2.25,0) lambda0`，保留全部旧通道 |
| noise / gain | 每实/虚分量标准差0.01；共享幅度 `Uniform[0.99,1.01]`、相位 `Uniform[-0.01,0.01]`；模拟真值从同 prior 生成 |
| pose | 起步已知；未知 pose 另定导航 prior 和锚点；各向同性点源没有角度敏感性时，不虚设可辨识姿态角 |
| 任务 | 目标中心二维平方损失；另冻结一个半平面占据事件；`tau_floor=(0.005 lambda0)^2` |

**E0 必须消除的剩余阻塞：**外域/PML 的固定数值、电流源归一化、数据单位，以及完整几何边界条款尚未指定；半平面事件边界、rank/discrepancy 接受规则和随机种子仍需预注册。在这些字段冻结并完成相关资格前，不生成正式测试数据。目标域与公共墙 ROI 分离是公开 prior，不是 GT 墙支撑。

### 统一预算与准确度门

| Campaign | 建议首轮上限 | 使用规则 |
|---|---|---|
| E0 + E1 合计 | 4,000 全波 RHS；1,200 FE 分解；2 CPU小时 | Prompt 0 和 Prompt 1 共用同一账本，不各自重置额度 |
| A0 小类与宽类锚点 | 12,000 RHS；3,000 分解；4 CPU小时 | 包括所有起点、支持搜索与失败 |
| E2/E3 开发 pilot | 12,000 RHS；3,000 分解；4 CPU小时 | 不自动启动锁定测试 |
| 锁定测试、E4–E6 | 尚未设定 | pilot 测得单位成本后事前统一锁定；未锁定时仅做实现和协议准备 |

以上是建议预算，不是本轮已执行费用，也不是历史92,000 RHS余额的延续。可在运行前依据硬件统一调整并记录理由；不得看见某算法落后后单独加预算。任一 RHS、分解或 CPU 上限先达到即停止。

已建议的初步门：数值见证白化响应差 `d_num <= 0.1`，且明显小于所主张的恢复间隔；该门不直接证明 posterior 足够准确。概率误差采用目标一维 CDF 最大差0.02、指定事件概率绝对差0.02、重要模式质量绝对差0.02。任务后悔门为 `tau = 0.05 R_ref + tau_floor`，首轮默认 `tau_floor=(0.005 lambda0)^2`；如需改变任务精度，必须在运行前统一修改。

开发 panel 为12个独立父世界，锁定测试为24个独立父世界，每世界两个独立噪声重复；这是可行性规模，不能据此宣布高精度95%校准。SBC 必须从声明 prior 抽样。锁定测试加速门为：达到相同精度门时，累计 CPU 与全波 RHS 的中位数均至少降低20%，并给出按父世界配对的误差区间和失败比例。20%是研究资源决策门，不是统计定律。

### 建议输出结构

```text
research/research_w_slam/
  contracts/
    experiment_contract.json
    geometry_contract.json
    material_prior.json
    noise_model.json
    acquisition_contract.json
    task_contract.json
    budget_contract.json
  manifests/
    public_roi_manifest.json
    development_world_manifest.json
    locked_test_world_manifest.json
    provenance.json
  data/observations/
  ledger/
    factors.jsonl
    factor_versions.jsonl
    computation_costs.csv
    acquisition_costs.csv
  results/e0/
  results/e1/
  results/a0/
  results/e2_e3/
  results/e4_e5/
  results/e6/
  checkpoints/
  reports/
  tests/
```

真值文件必须与优化器可读数据分离；上述目录仅为建议，任何实现都需在 manifest 中给出可复现实路径。每阶段最终输出一个机器可读 `gate.json` 和一份中文报告，附原始 arrays/样本/轨迹、命令、环境和累计成本，不仅输出图片。

---

## Prompt 0 合同审计 几何判定与前向资格

你负责 Research W-SLAM 的 E0。请应用本文件全部统一约束，并仅完成本阶段。

**目标：**把历史 `minimum_feature=0.1` 的歧义替换为事前确定、可判定的新实验合同；保证后续观测、likelihood 和成本可追溯。此阶段不做模式搜索、不训练生成器。

**输入：**两个固定 Git 提交、现有二维 scalar-TM FEM 实现、配套主报告、协议示例 JSON。先核对可取得的原始代码与运行 manifest，不能仅引用摘要代替代码检查。

**实施：**

1. 建立隔离工作目录和 provenance；读取旧合同、网格、材料、频率、端口、PML、噪声和 numerical audit 的实现。
2. 实现 `G_GRID_0125_V1` 布尔判定器，冻结公共 ROI、宏格尺度、对角接触规则、4连通拓扑、外边界接触、材料跳变、墙–目标接触/间距规则。边界条款未确定时不得声称合同完整。
3. 判定器输出合法性、拒绝原因、材料连通分量、孔洞数和面积；为合法/非法/边界例建立有解释的核对案例。优化器只看到公共 ROI，不得到真值拓扑标签。
4. 对旧真值仅做成员资格检测。不合格或无法判定时保留旧条件见证标签，另生成、锁定新的合法开发/测试真值，生成过程先于算法表现评价。不得挑容易成像的世界。
5. 若采用连续 rolling-ball 备选族，单独命名 `G_ROLLING_R_V1`，保留 `r0=0.05 lambda0`、双侧内切球、`C1,1` 和全局分离条件；没有区间曲率/非局部分离资格时不得用少量采样点报告严格 PASS。首轮不要求实现这一备选族。
6. 建立 passivity、frequency-consistency、prior 密度、参数单位与实化/白化检查。噪声合同完整记录每实/虚分量标准差或相关 covariance；未知参数填入并冻结，不能拟合失败后重定噪声。
7. 实现唯一 factor registry 和版本 registry。核验同一 observation 重算只增加计算版本；新独立噪声重复才有新 observation/factor ID。记录 prior、transition、likelihood 的身份，避免重复计入。
8. 小型前向资格分别检查网格/阶次、域/PML、求解残差和第二实现抽检。初始见证门为 `d_num <= 0.1`，并报告与预期信号差的比；posterior 资格另保留 log-likelihood/任务精度检查，不能用同一阈值替代。
9. 实现 RHS、分解、CPU/wall-clock、内存和审计费用计数。E0 与后续 E1 共用4,000 RHS/1,200分解/2 CPU小时 campaign cap。

**产出：**`contracts/*.json`、真值成员资格表、独立新真值 manifest、判定器与核对案例、`ledger/factors.jsonl`、`ledger/factor_versions.jsonl`、成本账本、`results/e0/gate.json`、`reports/E0_Contract_And_Forward_Audit.md`。

**通过门：**几何结果唯一可复核；未泄漏真值支撑；频率材料与噪声定义完整；相关前向误差受控；factor 重复检查通过。若任一门失败，输出已完成部分和阻塞条款，后续不解释物理信息谱或模式概率。达到 cap 停止，不通过放宽几何或数值阈值追认。

---

## Prompt 1 动态补偿理论与有限全波恢复

你负责 E1。前置条件是 E0 对本次输入给出 PASS。请应用统一约束，并仅使用 E0+E1 campaign 剩余预算。

**目标：**验证局部动态补偿条件在真实全波响应中何时成立、何时失败；同时保留严格不恢复对照。这个阶段允许已知局部基点，但不得据此宣称 measurement-only 反演成功。

**代数 gate：**先运行附带 `observability_toy.py` 与 `freezing_toy.py`，保存运行环境及新输出。增加独立目标、共享静态目标、有限运动 prior 三种结构，以及共享/独立 gain、pose gauge 的对照。对每种结构独立计算完整 joint Hessian、直接 Schur 和推导式。

记白化线性模型为 `delta_y = A w + B x`，`C=B^dagger A`、`R=A-BC`、`G=B*B`。在 `G>0`、运动协方差 `Q>0` 时核对

\[
K=R^*R+C^*D^*(Q+DG^{-1}D^*)^{-1}DC.
\]

非单射 `B` 必须先消去 `ker B`：令 `N` 为核基、`V` 为其正交补、`P=I-P_{Q^{-1/2}DN}`，使用主报告 D2，绝不能直接把 `G^{-1}` 换成 `G^dagger`。重复附带4维目标反例；错误公式给出虚假正信息，而完整 Schur 应为零。时间核条件应核对 `Rw=0` 且 `DCw in Range(DN)`。

**物理实施：**

1. 在公共合法几何中建立一堵共享墙、一个弱至中等对比度目标、两套事前固定采集。墙使用2–8个公开区域材料自由度、目标使用2个坐标及少量 nuisance。该区域划分仅服务机制实验。
2. 分别测试不同发射位置/孔径和严格追加频点；共享静态目标、两时刻自由目标和有限运动 prior 分开；已知 pose 与已正确固定 gauge 的未知 pose 分开。
3. 计算伴随/自动 Jacobian，并以合法中心差分或边界单边差分核对。rank 阈值必须结合网格细化与导数误差事前确定，不凭单张谱图手动选择。
4. 输出每次与联合 `K0`、`Keff`、retention 谱、angle sum、绝对信息、最小有效奇异值、材料度量、运动 prior 强度、`R`、`DC`、D2投影和恢复方向。区分 data-only 与 prior-assisted 信息。
5. 在恢复方向上做有限正负扰动的全波响应扫描，并重新优化全部允许 nuisance；不能只比较固定补偿端点。检查线性预测与实际 residual 的差及 resolvent/trust-region 限制。
6. 同时设计不恢复控制，例如独立目标分别补偿、共享仪器自由度吸收信号，或目标隐藏自由度消除动力学冲突。固定目标/墙初值错误、近共振和过大步长导致证书失效也要保留。
7. 新观测条件 innovation 与旧 Schur 更新应在同一固定线性化下核验。新增目标状态先按真实转移积分/消元；若改变模型或线性化点，不能直接引用旧信息单调定理。

**产出：**`results/e1/linear_algebra_checks.json`、Jacobian/Schur/metric arrays、numerical rank audit、finite-step curves、重新补偿优化轨迹、完整成本、`results/e1/gate.json`、`reports/E1_Dynamic_Compensation.md`。

**通过门：**至少一个合法非退化场景的联合增强超过数值和 rank 容差，并在一个明确有限扰动区间得到全波核验；同时有严格不恢复对照。理论实现与物理机制分别判 PASS。预算不足或局部搜索失败标 INCONCLUSIVE，不能据此证明全局无解。若一个预注册 campaign 内未出现独立物理恢复或内容完全被旧理论覆盖，暂停 A 主线。

---

## Prompt 2 真正未知支撑的 Measurement-only 锚点

你负责 A0。前置条件是 E0 前向与合同通过；E1 是否通过单独报告，不允许用已知基点成功替代本阶段。

**目标：**不读取真实墙形、孔洞数或 near-GT 初值，从观测和公共 prior 找到数据/留出预测相容的联合墙–目标状态。这是一般在线地图实验的强制前置门。

**实施：**

1. 设置隔离的优化与评价进程。优化可读观测、公开 ROI/采集信息、proper prior 和材料合法域；真值 manifest 只供评价进程。输出访问边界与复现命令。
2. 从4个公共候选宏区的全部16种占据组合开始，统一应用几何判定器并记录合法数量。连续未知量包括事前冻结的共享或区域介电参数、损耗、目标坐标、共享 gain。先用少量连续材料参数维持可核查 reference。
3. 扩展至3×3区域的至多512种 occupancy，再应用同一合法性筛选。孔洞数由占据模式决定，不预给优化器；小类成功仍不等于一般墙体成功。
4. 宽模型使用全公共 ROI 的多分辨率材料/level-set 或可逆 occupancy 编辑。初始化只来自均匀材料、公开 prior 样本、普通 FWI/TR 或测量驱动 proposal。第一次 backprojection 的 Top-K 不得永久截断 posterior 支撑。
5. 强解法至少包含合格梯度、线搜索/信赖域、多分辨率、合法多起点、joint target/gain 优化，以及支持 birth/death 或交换。使用固定已有频率的 continuation 属于计算，不登记为新采集。
6. 对小类建立独立求积/枚举参考，并用温度桥/SMC 或独立 MCMC 交叉核查。给出每个离散支持的 proper prior、连续密度与 evidence 积分，不用优化最小 residual 直接当模式质量。
7. 接受锚点须同时满足冻结的拟合/后验预测规则和独立留出通道预测规则。先以 prior-predictive 仿真校准统计量；不拿含噪 `||y-F|| <= 1` 当通用门，不用 GT 误差决定停止。
8. 记录全部初值、失败、局部极小值、支持改变、数值资格和费用；将表示不足、非凸优化、prior 排真、support 搜索失败及仪器混淆分开诊断。

**预算：**A0 共12,000 RHS、3,000分解、4 CPU小时，任一先达即停止。不得从 E1 或历史预算自动借用未声明余额。

**产出：**小类合法 support 列表和 reference、宽类 blind search 轨迹、holdout 预测、factor/cost ledger、数据访问审计、`results/a0/gate.json`、`reports/A0_Measurement_Only_Anchor.md`。

**通过门：**主优化未读 GT 且在冻结合同内达到拟合与留出预测标准。小类与宽类分别判定。若宽类失败，停止一般在线系统 C 的成功主张；可继续受控局部理论，但禁止以旧两墙 prior、GT 支撑或 near-GT 初值替换失败后继续报系统 PASS。

---

## Prompt 3 任务安全回访与强 Smoother 基线

你负责 E2/E3 开发 pilot。前置条件是本次 A0 主测试和小型 posterior reference 通过。本阶段只做开发 pilot，锁定测试预算未冻结时不得自动启动。

**目标：**在相同概率表示、数据、prior 和任务精度下，检验任务证书能否比普通 incremental smoothing 更便宜；同时验证旧地图修正是否确实参与目标风险改善。

**实施：**

1. 先固定同一个相关 Gaussian 或声明的 mixture 表示，固定 full-wave 后端和数据序列，仅比较调度。基线为 hard layer freeze、普通可迭代 layer stripping、batch joint FWI、普通 incremental smoothing、周期回访、不含任务指标的 soft-freezing、任务证书策略和低维 direct posterior reference。
2. 至少覆盖三种序列机制：早期墙–目标近补偿；新孔径激活旧弱方向；运动与错误墙相位发生冲突。真实时间戳、状态转移与动力学边固定；另用一个预注册 factor 处理顺序随机排列作数值敏感性对照，不能通过洗牌改变 moving-target 的真实时间顺序。
3. 对候选冻结块 `F`、其余联合变量 `U`，计算局部 `S_F`、`r_F`、任务传递 `T_F`。区分二次目标损失 `0.5 r_F* S_F^-1 r_F` 与平方任务后悔 `||M^1/2 T_F S_F^-1 r_F||^2`；不能把二者互换。
4. 实现对偶残差区间 D8。若只有经验最小 Ritz 值，没有可验证 coercivity 下界，结果标“经验审计”，不得写严格证书。统计 cache、dual solve、历史因子刷新、全局 fallback 的全部费用。
5. 明确局部 Gaussian regret 与真实 posterior regret。将 forward、非线性、模式覆盖和尾部误差加入独立预算；均值不需更新不等于可以丢弃 covariance。每次 relinearization 只替换原 factor 的版本。
6. 保留三条对照轨迹：后期数据进入但禁改旧墙；允许改旧墙且当前目标使用同等优化预算；完整联合 reference。输出后期 innovation、历史墙变化、目标 task-risk 和旧因子刷新数，不把单例相关自动称因果定理。
7. 概率指标至少含 occupancy Brier/log score、材料条件分布、区域 covariance、墙–目标 cross-covariance、模式质量、目标 CDF、事件概率、credible interval 覆盖与宽度、在线 Bayes action。
8. 开发 panel 使用12个独立父世界、每世界两个独立噪声；若做 SBC 必须从声明 prior 抽样。reference 数值误差须低于比较容差。目标 CDF、指定事件和重要模式质量门均为绝对误差0.02；任务门为 `0.05 R_ref + tau_floor`，首轮 `tau_floor=(0.005 lambda0)^2`，必须事前冻结。
9. 画风险–RHS、风险–CPU、概率误差–CPU曲线。记录未达门比例，不能只比较成功样本后隐藏失败。所有方法可使用同等缓存、warm start 和预条件器。

**预算：**开发 pilot 共12,000 RHS、3,000分解、4 CPU小时。之后仅基于单位成本和开发数据拟定锁定测试预算，不自动执行24世界 panel。

**产出：**调度策略实现、证书区间与实际误差配对、完整 posterior/任务轨迹、校准与成本表、所有失败和 fallback、`results/e2_e3/gate.json`、`reports/E2_E3_Revisitation_And_Calibration.md`、锁定测试预注册草案。

**通过/停止门：**有合法误差控制并显示同精度成本节省迹象，才进入锁定测试。最终方法收益须在同精度下 CPU 与全波 RHS 中位数均至少降低20%，附配对不确定区间。若普通 smoother 同样好或证书比完整更新更贵，结论为无独立调度收益。最多进行两次事前定义、对全部方法公平的合理工程修正；仍无优势则停止 B 加速主张。

---

## Prompt 4 宽合法材料类的模式发现与传播压缩

你负责 E4/E5。前置条件是相应物理合同、A0、小 posterior reference 通过，且本阶段预算已经依据 pilot 单位成本锁定。未锁定预算时只完成实现、预注册和小型纯代数检查，不启动物理模式 campaign。

**目标：**保留真正墙–目标联合模式，评估错误早期分支能否被后续观测修正，以及传播压缩是否能在不牺牲任务概率的前提下降低成本。

**实施：**

1. 忠实回放旧 conditional joint witness 的原合同与独立数值审计，另报新几何成员资格。它只作机制起点，不作为全部 prior 支撑。
2. 在新合同下使用独立合法墙体和未知支撑搜索。分开固定旧对、允许连续材料/目标重新补偿、允许合法支撑改变三种实验。固定旧对距离增长不证明宽类唯一。
3. 先允许错误分支附近初值作为明确标记的机制测试；主结论使用 measurement-only 起点。记录分支变形、墙体承担目标位移、目标 credible set 和在线 Bayes action，禁止 GT 选分支。
4. 比较 single MAP、独立 marginals、相关 Gaussian/Laplace、mixture、Rao–Blackwellized particles、direct reference 及传播推前表示。所有方法使用同 prior、支持和数据，不能给新方法免费真值附近模式。
5. 实现 Birth/Survival/Reweight/Merge/Reopen/Retire。有密度 proposals 使用合法 SMC/MH/importance 校正；变维时明确占据条件测度和 Jacobian。局部模式质量不能只按峰值 residual 或发现次数计算。
6. 小 reference 中核对模式发现与 posterior mass；宽类只能报告已发现质量、覆盖诊断和风险限制，不能声称穷尽。保留小质量高任务损失模式的压力测试。
7. 若实现精确场消元，使用 `F(W,X)=d_W+R_W(K_W+Delta L_X)^-1 f_W`，保存 `p(K_W,f_W,R_W,d_W,X)` 的联合分布。不得把它分解成独立墙传播与目标 marginals；任务直接询问材料/孔洞时，保留对应标签条件分布。
8. 冻结压缩所覆盖的目标域、频率、pose、激励和未来查询集合。当前端口等价不意味着未来全部端口等价；超出合同必须解压/重开。
9. merge/retire 分别报告有界损失预算或相应矩条件、mode mass 与未来 likelihood-ratio 限制。当前低权重不提供永久删除保证。预测 TV 界只支持声明的平均决策结论，不能自动宣称每个具体新数据后的 posterior 都近似。
10. 与普通 Galerkin/ROM、direct quadrature/direct posterior 比较全成本。A17 是可替换后端，只有在同 current space、同物理精度控制后仍有净收益才另作贡献。

**产出：**`results/e4_e5/` 下的模式样本、归一化质量与证据检查、合法性/全波审计、压缩查询合同、merge/retire 误差记录、历史恢复轨迹、成本、`gate.json` 和中文报告。

**停止门：**模式权重不可核对时只报告候选集；重要模式覆盖无法评估时不报告校准 posterior；传播压缩与普通方法同价或更贵时不声称净收益；新查询超合同时立即撤销旧安全标签并恢复原始分支。

---

## Prompt 5 严格嵌套的主动新增测量

你负责 E6。前置条件是明确允许新增采集的合同、可核查当前 belief、小型 reference 和已锁定本阶段预算。不得把同数据重复计算混入采集收益。

**目标：**比较测量设计能否在计入真实采集费用和设计计算费用后降低任务风险，并检验旧补偿模式能否通过新测量下的重新优化存活。

**实施：**

1. 在 acquisition contract 中分别定义同数据重算、同几何独立噪声重复、追加频率、追加收发位置/孔径。采集成本、时延和允许动作集合未定时，不执行主动实验。
2. 旧通道永久保留在数据集合中；每个真正新增 observation 有新 factor_id。设计器只能读当前 belief 与公开动作信息，不能看真值和未来测试数据。
3. 在低维 reference 上实现嵌套 task EVSI：外层对新观测预测平均，内层对共享墙、目标、pose、gain 和模式变形重新推断。允许材料和目标重新补偿，不能只比较固定旧端点。
4. 比较随机、固定均匀宽带、总 EIG、最大 singular value、angle sum、generalized spectrum、Gaussian task proxy 和 exact task EVSI。共享 nuisance 的结构和 prior 对所有设计一致。
5. 报告每次动作的 acquisition cost、design-computation cost、推断成本、expected task-risk 改善与实际风险分布。比较同总采集预算下的风险；不能把免费多频方案称作算法加速。
6. 分离信息增加和数值逼近改善。重复同一激励的独立噪声可以降噪，但精确 `F_a(theta1)=F_a(theta2)` 时不能破坏这对物理等价。
7. 估计内外层 Monte Carlo 误差。若前两名动作差距低于设计估计误差，按预注册平局规则处理，不用 GT 选更好的动作。
8. 旧窄/宽带比较如果不是严格包含关系，仅作背景。对本轮新测量分别报告固定对距离、连续重补偿最小值和宽支持搜索结果；局部搜索不能证明全局唯一。

**产出：**动作集合与费用合同、每步候选预测及不确定区间、采集账本、嵌套 factor registry、重补偿轨迹、同费用 risk 曲线、`results/e6/gate.json` 和 `reports/E6_Active_Sensing.md`。

**停止门：**采集合同不完整、reference 模式覆盖不足或估计误差不能区分动作时，报告 INCONCLUSIVE；如果简单固定/随机设计同费用达到相同任务风险，不把主动策略重新命名为成功。

---

## Prompt 6 证据裁决 论文路线与三维启动门

你负责最终审计，不新增未预算的 PDE campaign，不为取得正结果继续增加模块。

**目标：**根据原始证据决定推进 Theory First A、Algorithm First B，或停止有关主张；集成系统 C 不作为默认结论。

**实施：**

1. 汇总所有版本化合同、gate、原始数组、posterior 样本、分支质量、失败世界、数值审计和成本账本。核对报告中的每个数字是否能回溯到 commit、配置和运行输出。
2. 审计 factor ledger：同 observation 是否重复指数；tempering 是否使用增量幂；使用旧 posterior 后是否又乘旧 likelihood；墙重开是否误作物理变化。发现重复使用数据时撤销相应概率结果并重跑预算内可复核部分。
3. 对照旧 IS-SLAM 结果和主报告文献矩阵，逐项标 S/U/D/C/N/H。共享 Schur、主角度、普通 smoothing、一般任务 OED 和 Gaussian covariance 恒等式不单独作为新贡献。
4. A 路线只保留通过真实全波恢复/不恢复见证、目标隐藏自由度处理和有限非线性限制支持的具体命题。若只是标准代数加数值例子，降低或撤销原创定位。
5. B 路线必须在相同 posterior/任务精度下超过合理 incremental smoother/direct inference；计算端的 CPU 与全波 RHS 均满足预注册收益门，并报告失败比例、配对区间、参考建设费用和部署费用。
6. 明确区分 prior-assisted 恢复与纯数据可辨识；measurement-only 初始化与已知基点机制；固定两点 Bayes 下界与宽材料 posterior mass；局部非线性有效域与全局唯一性；二维 scalar-TM 与三维 full-vector。
7. 检查两次事前定义的工程修正之后是否仍无任务证书成本优势。若没有，推荐 ordinary joint smoother/direct posterior 并保留负结果，不再用复杂后端掩盖。
8. C 的前置条件包括宽支持 A0、可核对模式覆盖、任务校准和回访成本门；任何一个未过，都不报告一般未知墙在线成像成功。
9. 仅当 E1 物理机制、A0 与至少一个 E2/E3 成本–概率门通过，才另拟 E7 小型3D vector Maxwell 协议。它需要重新认证 curl-conforming 离散、散度/边界、极化、源接收和网格误差。本阶段只输出协议，不自动执行三维。

**产出：**`reports/FINAL_EVIDENCE_AND_ORIGINALITY_AUDIT.md`、机器可读总 gate、claim–evidence 矩阵、A/B/C/STOP 裁决、完整未完成项和受限结论、可复现执行清单，以及仅在符合启动门时编写的 E7 协议。

**最终停止规则：**无法证明新增数学内容、无法取得合法有限全波机制、无法完成未知支撑 measurement-only 锚点、无法核对模式概率，或同精度下不比强基线便宜，均应按对应范围结束，不将“实现了更多模块”当成科研成功。
