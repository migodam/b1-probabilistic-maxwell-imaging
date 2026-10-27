# Research W-BGGA：固定来源与本地缓存审计

审计范围：只读两个指定提交的报告、配置及代码，并核对原始文件位置、文件大小和 NPZ 成员/NPY 头。没有导入旧研究模块、运行新物理、重算概率、修改旧文件、做哈希校验或进行新推导/全局新颖性裁决。完整文件表见 `SOURCE_FILES.json`，路径与数组 schema 见 `CACHE_MAP.json`。

下文 `W:` 指提交 `958a478739d9553269316ded8125fda9df225c0c` 下 `research/research_w_p0r_to_p1/`；`S:` 指提交 `52351aee6ddd41b09851fe5e52d1a8861891ef13` 下 `research/research_w_slam/`。行号均取 `git show` 内容，未用后续分支替换证据。

## 1. 固定来源与实际位置

- 任务指定的旧 W 10 份 md/config/JSON、3 个 wr 程序，以及 W-SLAM 报告、E1/A0/E2/E3/E4、动态说明、6 个 wslam 程序、2 个工具与最终 bank 费用均已读取。另补读 E0 gate、物理协议、prior、bank 构造和旧授权/封账等，共登记 48 个固定来源。
- 旧 W 完整原始根：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-w-p0r-to-p1/research/research_w_p0r_to_p1`。
- W-SLAM 原始根：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/Research_W/research-w-slam`。发布目录里的同名路径保留源码/报告/receipt；核心大数组应从这个原始根读取。
- 对登记来源做普通逐字节比较：47/48 个原始本地文件等于固定提交；唯一差异是 `S:COST_SUMMARY.json`。固定版 31818 事件、E0E1 CPU 754.130853 秒；原始本地版 31831 事件、E0E1 CPU 769.385534 秒，其余 stage 的 used 相同。指定科学报告/源码和 `JOB_COST_FINAL.json` 相同。此尾项差异不能默默覆盖固定版。比较没有计算哈希。

## 2. 三种“24”必须分开

| 数据 | 原始根内路径 | 实际存储与身份 |
|---|---|---|
| 2950 原子定义 | `results/e2e3/finite_bank/online/FINITE_GRID.npz` | support `(2950,4)`、pose `(2950,2)`、prior_weight `(2950,)`；14 support、118 wall groups、25 pose |
| 全波均值缓存 | `results/e2e3/finite_bank/online/FINITE_BANK.npz` | complex128 `data=(2950,4,16)`；2950 个 `paid_rows/0000.npz` 至 `2949.npz` 连续存在；仅 0000/1475/2949 三个代表全场另存 |
| 有限实验 24 case | `results/e2e3/finite_bank/online/OBSERVATIONS.npz` | `data=(12,2,4,16)`，12 parent×2 噪声；base/outer 各 `(12,2,2,16)`；静态同一目标逐步增加通道 |
| 连续 A0 另 24 case | `inputs/online/A0_OBSERVATIONS.npz` | `base=(24,2,16)`；配套 `inputs/offline/A0_HOLDOUT.npz`、`CONTINUOUS_TRUTHS.json`、`CONTINUOUS_CLEAN_MEANS.npz`；不是有限 prior 样本 |
| 旧 W 的 24 probe 任务 | `results/B2/validation/joint_seed2_narrow_fullj/probe*_PAIR.npz` | 24 对固定目标插入响应，属于旧 W 根；不是上述 24 个含噪案例 |

有限模型的 gamma=1、目标半径/contrast/loss=.12/3/.02 已知；2950 是离散 prior 定义，不能叫连续 A0 prior 的收敛求积。证据：`S:src/wslam/finite_bayes.py:31-68`、`PRIOR_AND_LIKELIHOOD.md:19-29`、`tools/finite_bank.py:94-124,176-201`。连续 A0 则未知有界共享 gain，477 个原型只是提案库：`S:A0_ONLINE_RECONSTRUCTION.md:7-11`、`tools/a0_prepare.py:12-54`。

完整参考权重也在原始根 `results/e2e3/inference/POSTERIORS.npz`：reference `(12,2,3,2950)`、method `(12,2,3,11,2950)`；E4 权重为 `(384,2950)`。这些是已暴露输出，未来算法进程不能把参考答案当自己的在线状态；旧 cohort 不能改称新盲测。证据：`S:E2_BELIEF_MAP.md:13,26,93`、`SCIENTIFIC_DECISION.md:42`。

## 3. 白化、排序与静态/动态身份

- 实现的 pack：先按 C 顺序展平复数据（有频率时为 frequency/source/receiver），再全 Re、全 Im，统一除 `.01`；packj 同样把 Re/Im 上下堆叠。不除通道数平方根。所选实向量的 Sigma 为 `1e-4 I`。证据：`W:src/wr/physics.py:15-17`；`S:src/wslam/physics.py:16-33`。
- 旧 `W:configs/ORIGINAL_P0.json:112` 写的是逐 source/frequency 的 `[Re(receiver),Im(receiver)]`，与当前 pack 的字面排序不同。iid Sigma 下只是固定置换，不改变范数；新合同仍须明确采用哪一个排列，并让 y/J/Sigma 一致。旧 W 4TX×12RX×5 频率，不能与 W-SLAM 4TX×16RX 单频表混用（`W:src/wr/argmax_validation.py:78-81`；`S:RESEARCH_W_SLAM_REPORT_ZH.md:21-25`）。
- 有限三因子依次为 TX0/RX8、TX1/RX8、TX0/1 的其余 15RX；目标静止。base 全块与 t1/t2/t3 是同一数据的两种表达，不能同时乘 likelihood。每个 repeat 单独分析。证据：`S:tools/finite_bank.py:108-124`、`E2_BELIEF_MAP.md:13`。
- E1 的 pose0=(.2,2)、pose1=(.3,2.05)、位移 prior 与 D1/D2/D3 属于另一个动态/代数实验；不可作本轮静态目标的新科学证据。证据：`S:DYNAMIC_COMPENSATION_THEORY.md:7-21`、`E1_MAXWELL_OBSERVABILITY.md:34-47`。factor 更换版本应替换旧 likelihood；registry 与 log-ratio 已有实现，见 `S:src/wslam/factors.py:20-40,52-69`。

## 4. 缓存能复用到哪里

- `S:public_mesh` 显式插入 `circle_vertices(pose,radius,n)`；即使 wall-only 求值也经同一网格入口。Evaluator backend key 是 `(pose_x,pose_y,radius,f)`，变化即 clear/remesh；last-world key 还包括精确 support/材料 bytes、contrast、loss、include_target。证据：`S:src/wslam/forward.py:9-37,54-85`。
- 同一 live FEM backend 的 LU 只在 quadrature epsilon 逐值完全相等时复用；新材料先清旧 factor。即便复用 LU，实际源/转置动作仍计 RHS。证据：`S:src/wslam/fem_backend.py:292-295,459-504,541-557`。本地保存响应/全场/基向量，不等于有可跨 target mesh 使用的已序列化稀疏 LU。
- `results/e0/four_blocks/BLOCKS.npz` 实际存在：A `(16,2335)`、B `(2335,2)`、C/mass `(2335,2335)`。它是 coarse `.22/128/grid1`、固定 pose/radius 的同离散恒等式；准备 2335 port RHS 对比直接 2 source RHS。没有跨目标重用或加速资格。证据：`S:tools/e0_four_blocks.py:17-40,53-61`、`E1_MAXWELL_OBSERVABILITY.md:7-11`。该工具顶层会执行物理，不能为了读取定义而 import。
- 旧 W 的 `P1B/green_reuse_control/{BASELINE,ROM_SETUP}.npz` 与 678 节点 P1A bank、650 节点 `P1B/full_material_rule/FULL_RULE.npz` 均在原始根。Green 对照只含单频、一个 wall branch、固定中心(-.8,2.15)、半径.12、contrast 2.94–3.06；代码 `SIGMA=.02` 是目标损耗，观测 quadrature sigma 仍 .01。证据：`W:src/wr/green_reuse_control.py:31-36,180-198`。普通 ROM 同享付费 R0V；没有实测完整 standalone cold 部署或 likelihood interval。见 `W:P1_B_GREEN_COMPUTATION_CONTROL.md:17-29`。

## 5. 旧结论与费用不能升级

- 旧 joint 见证保持 oracle/fixed-support 身份：d_joint=.7987161019689071、d_cal=3.5605119112114565、max24 task=4.81819338103861；它不提供 calibration omega_T(1) 的见证。probe6 自己的粗细误差 .00705354/.01317174，不能借 probe8 精度。证据：`W:P0_B2_JOINT_COMPENSATION.md:3-29`、`P0_C_NONLINEAR_MODULUS.md:3-9`、`P0_B2_FULLWAVE_VALIDATION.md:7-19`。
- 既有状态保留：E0 `QUALIFIED_DISCRETE_FULLWAVE`，E1 `INCONCLUSIVE`，A0 `A0_ONLINE_NOT_REACHED`，E2 仅 finite-prior PASS，完整 E3/E4 仍有限定/缺口。A0 主72拟合及补优化均触评价上限。证据：`S:SCIENTIFIC_DECISION.md:7-19`、`A0_ONLINE_RECONSTRUCTION.md:25,45-56`。
- 2950 表真实 setup：11800 complex RHS、2950 FE、9288.023040 CPU 秒、9496.344884 合计 job wall 秒，包含失败启动；不能把数组阶段 0 新 RHS 写成免费总成本。证据：`S:results/e2e3/finite_bank/online/JOB_COST_FINAL.json:3-8,33-49,104-109`。

## 6. 预算授权边界

- 旧 W `configs/AUTHORIZATION.json:7-33` 的 R/B2/C/P1A/P1B 是 20261023 旧 campaign 分阶段授权；`W:RESOURCE_SUMMARY.md:13,35` 明写剩余 31640 RHS 封存且任何新物理需独立授权，`CLOSED.json` 标记 `CAMPAIGN_CLOSED`。`ORIGINAL_P0.json` 是历史 P0 配置，不能拿它的旧 not_authorized 字段替代后来的 AUTHORIZATION，也不能把后来旧授权搬给 BGGA。
- 旧 W-SLAM caps 仅为 E0E1=4000/1200/7200、A0/E2E3/E4 各=12000/3000/14400、E5=8000/2000/10800（RHS/FE/CPU秒），见 `S:src/wslam/accounting.py:23-37`、`EXPERIMENT_CONTRACT.json:315-354`。ledger 拒绝不兼容旧账本，按 stage 预约且失败不退款（`:44-62,104-138`）。
- 原始本地账本当前 A0 仅余 84 RHS；E2E3 余72 RHS/18 FE。E4 的12000未用 RHS 仍属旧 E4 scope，不能默认借给 BGGA。固定版与本地版 caps/used/remaining 全列 `CACHE_MAP.json`，不调用会写 ledger 的 BudgetBook.snapshot。
- 附件自身 `source_review/PRO_REPORT_EXTRACTED.txt:1454-1458` 把 BGGA 数量称为建议设计，明确不是既有授权。**本只读审计没有在旧来源中发现新 BGGA 物理资源授权。** 若人类在当前主会话另有直接授权，由主线程单独登记；本审计不替主线程判断当前授权。

## 7. 未完成的核验与所需主线程判断

文件存在、大小及 NPZ header 已核对，未逐一检查所有 payload 数值或重演历史数值结果；因此不把 schema 存在等同于全数据正确。12 点 bank 独立网格审核的最大 .0207587 与20个不同后验端点审核是有限检查，不覆盖2950状态的统一 PDE 误差。证据：`S:E2_BELIEF_MAP.md:15`、`RESEARCH_W_SLAM_REPORT_ZH.md:29-33`。

没有对核心公式或全局新颖性作裁决。主线程仍需决定本轮静态合同、参考/在线进程隔离、新预算是否已被人类明确授权，以及新算子是否真正拥有独立于 target mesh 的状态空间。
