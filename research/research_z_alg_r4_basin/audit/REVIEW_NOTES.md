# 独立审查工作记录

更新：末节已追加 V3 真实结果复审；首版记录保留。

日期：2026-10-30。首版范围：只读数学与实现；无新数值运行；无理论/code/results 修改。拥有文件仅为本记录、`FINAL_SCIENTIFIC_REVIEW.md` 和根目录 `PAPER_NOVELTY_REVIEW.md`。

## 已读资产

- `MATHEMATICAL_THEORY.md`
- `theory/RESOLVENT_GEOMETRY.md`
- `theory/BRANCH_NAVIGATION.md`
- `FORMULA_AND_CODE_REUSE.md`
- `BASIN_GEOMETRY.md`
- `src/basin_solver.py` 与 `src/rational_rays.py` 核心路径
- `src/physics.py` 指向的冻结公共后端中 forward/J/Gram；共同 `Regularizer` 与 material basis
- `configs/DEVELOPMENT_CONTRACT.json`、`audit/INITIAL_RESEARCH_DECISION.md`、`audit/literature_bundle.json`

代码正由其他代理并行修改，本审查是可见快照，不是不可变版本证明。没有为此做 hash 或重复源码比对。

## 已向主代理即时通报

- P1 机制范围：校正前 `candidate_rms=.25` 合并并非只合并相同种子，可能剔除分离低谷。
- P1 命名范围：默认 ray corrector 是 L-BFGS-B+adjoint，full-J 主要用于 GN ray 方向；不能写每低谷 full-J correct。
- P2 一致性：构造候选末尾 clip 必须限于浮点级盒边缘修正；若材料改变显著，raw ray response 与实际材料不再相同，最终 forward 单独决定物理评分。
- 非致命精确措辞：同槽策略不是每候选严格相同实际秒数；`corrected` 不是局部最优证书；best-only 是 detected-top1。

## 隐藏条件核对摘要

两侧 resolvent 因式分解、精确余项、共同幅度分母、Schur 导数、Krylov 二阶残差和 Frobenius 椭球证书均未发现乘法顺序错误。范数证书仍依赖真正上界；数值估计不能替代。局部收缩依赖精确根及邻域正则性；非零残差与正则项改动必须另算。

两个一材料反例处于盒边界；额外 tied-material 内点例子改变参数域，不可扩展为全独立材料空间结论。光滑两变量 least-squares 多低谷反例正确但非物理。它们可作解释，不能作本次 Maxwell 性能证据。

## 文献 ledger

| Evidence | 来源 | 可支持 | 不可支持 |
|---|---|---|---|
| `Tier C` | 现有 ScholarQA bundle 的 Habashy–Abubakar 2004 DOI 10.2528/PIER03100702 | 记录身份与相关题名 | 特定机制是否既有、最接近先例、全球新颖性 |
| `ERROR_429` | 现有 bundle 多次查询 | 检索未完整执行的事实 | 没搜到等于不存在 |
| `LOCAL_PRIMARY_READ` | R2/R3/R10 与本轮 reuse map | 项目内直接公式/概念前身 | 全球文献优先权 |
| `PARENT_VERIFIED_PRIMARY_ANCHOR` | CSI 1997、recursive linearization 2017、arXiv 2403.15292v2 | 主代理可据其已核对片段继续最终 claim ledger | 本审查未重取全文，不能冒充自行取得 Tier A |

本轮独立审查没有新检索、没有新增外部引文结论，没有把现有 ScholarQA 失败补写成完成。已阅读 `scholarqa-research/SKILL.md` 并依照其证据分层复用 bundle。

## 真结果到达后的有限追加审查

核对 multi-valley 与 detected-top1 的差异候选是否可追踪；检查较高 raw score 的低谷是否真的贡献更低 post-correction 原目标，以及是否满足残差/RMSE/IoU 联合预声明门。计入所有失败校正和谱 setup。比较同 corrector 普通多起点以区分结构化种子与一般重启，缺少该对照时明确 attribution gap。

保留历史 GT exposed 状态、冻结时间、样本数、停止原因、预算 overshoot 与未完成方法；禁止从 leading-32 能量推断全576良态满秩。真实结果不足以解除 HOLD 时，保留失败结论。

## V3 post-results 追加记录

依主代理限定，在五分钟内定向读取：`results/fresh_v3/SCORES.json`、`configs/FRESH_CONTRACT.json`、`audit/FRESH_PREDICTIONS_FROZEN.json`、`ALGORITHM.md`。从既有 JSON 摘取逐行指标；未启动物理求解或任何额外检索/校验。记录声明：18/18 预测于 03:15:35 UTC 关闭，03:19:07 UTC 评分，首次真值读取在预测关闭后；独立审查只评这些声明与结果，没有审计隔离进程。

历史首审问题处置：V3 pre-dedup=.0000001；最多四个 amplitude-spread 检出低谷并记录截断；每 seed 有 wall/LU/RHS 配额；ordinary seed1 使用已付费 anchor；L-BFGS 为有界仿射坐标与精确链式梯度。当前算法说明已准确称 L-BFGS 物理 corrector，未继续冒称 full-J 校正。原先 .25 pre-dedup 问题已不作为 V3 现存缺陷。

裁决的直接依据：六方法 joint 0/3；bridge 在 7/8/2 的数据残差 .01542465/.11146390/.04576007，但 top1 RMSE 10.890845/10.593344/11.327416，oracle best-K 仍为保留 anchor 5.050260/6.272147/6.337567。best-only 三例残差分别 .76223514/.78080918/.91138354。ordinary multistart 三例 .05463381/.13702526/.04153385，在 digit2 优于 bridge，亦无 GT 恢复门成功。

主张状态更新：`GT_IMPROVEMENT_NO`、`JOINT_RECOVERY_FAIL_ALL_METHODS_0_OF_3`、`FINITE_DATA_ACCESSIBILITY_LIMITED_POSITIVE`、`MULTISTART_SUPERIORITY_NOT_ESTABLISHED`、`NEW_THEORY_NOT_SUPPORTED`。全球新颖性仍 HOLD。低残差候选不是精确同数据根、物理分支或概率解集。

必须保留：bridge digit8/2 为 RHS 60000 的 BUDGET_STOP，只校正三个候选；ordinary multistart 三例 2500 LU 的 BUDGET_STOP。没有 overshoot 与完成全部候选是不同命题。首版 leading32/full576 条件边界完全不受本次性能结果改变。
