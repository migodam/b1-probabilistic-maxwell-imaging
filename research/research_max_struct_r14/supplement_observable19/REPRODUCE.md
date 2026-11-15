# G19 Observable19 补充复现

本补充不是R15，不改旧R14科学结果。冻结主实现为提交 `fbcc230ef88b7a21b79ca753d11397dba85ccb78`，随后只补离线分析、错误记录及公开网页。最终代码和下载见本目录 `evidence.html`，原R14公共输入来源见[原R14下载页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r14/evidence.html)。复现不需要重新生成材料或训练NN。

## 环境和来源

实际计算：Windows i9，12个已核实逻辑核0/2/4/6/8/10/12/14/16/18/20/22，8P+4E；每对象六臂同核顺序，一BLAS线程、NORMAL优先级。没有新增CUDA实现、修改电源或驱动。主队列在2026-11-15 UTC完成144/144，全部原始数组、完整current/field、分母基、readout快照、候选池、日志和费用保留在本补充release，未经材料-only缩减。

仓库完整checkout提供R8/R10/R11依赖，补充自己的runtime提供冻结R13接口。原R14的runtime/equal_budget_C2_runner.py、frozen_C2_runner.py、src/queue_independent.py及scoring_history应保持原样。程序会在实际导入时记录依赖路径，不制造历史目录。准备公开输入：把原R14 source/public archive恢复到父目录，使 `R14/data_public/independent24/r14_new_0000.npz` 至 `0023.npz` 可读。本补充config通过相对父路径读取它们，不复制GT进公开观测。

以下R14_ROOT、S_ROOT和REPO_ROOT表示你实际checkout的绝对路径，不是声称已存在的通用路径。原始远端实际分别为 `D:/AI/MAX_STRUCT_R14_DUALMETRIC_20261115`、其 `supplement_observable19` 子目录、`D:/AI/MAX_STRUCT_R11_DEV_20261107`。

## 实际主命令

```text
D:/python/python.exe R14_ROOT/src/queue_independent.py --root S_ROOT --legacy-repo REPO_ROOT --manifest S_ROOT/configs/observable19_exposed24_v1.json --cores 0,2,4,6,8,10,12,14,16,18,20,22
```

单任务入口如下，task_id必须来自144行冻结config，不自己重新选样：

```text
python S_ROOT/src/run_independent.py --root S_ROOT --legacy-repo REPO_ROOT --manifest S_ROOT/configs/observable19_exposed24_v1.json --task-id r14_new_0021__g19_observation_galerkin --cpu 0
```

同一输出路径不覆盖；如重跑，复制源码/config/公共父输入到新的研究输出根，保留旧回执。CPU以整进程时钟设置120秒上限，超时操作仍收费、超时完成候选不能被选。外部OS测得的少量尾部越界和内部候选CPU分别保留，不能把120说成每个进程精确在120.000退出。

## 封存、评分与解释性诊断

程序运行不读取真值。所有144任务终态后才生成seal，父负责人复查六臂各24、费用和预测清单，再另存authorized=true的批准副本；seal没有SHA操作。重新恢复到不同路径时应新建本地seal，不改远端原封存。

```text
python S_ROOT/src/analyze_results.py --root S_ROOT --make-seal S_ROOT/predictions/PREDICTION_SEAL.json
python S_ROOT/src/analyze_results.py --root S_ROOT --seal S_ROOT/predictions/PREDICTION_SEAL_APPROVED.json --out S_ROOT/analysis/observable19_v1
python S_ROOT/src/fixed_state_analysis.py --root S_ROOT --legacy-repo REPO_ROOT --seal S_ROOT/predictions/PREDICTION_SEAL_APPROVED.json --out S_ROOT/audit/fixed_state_v2
```

第一条只生成未授权清单，不自动批准。后两条的新输出目录不能预先存在。本轮第一版fixed-state诊断的新增一致性输出误取了物理state里的m_after，造成21个G-PGP诊断记录不完整；修正为原保存common快照后完整重放，v1和费用仍保留。它未影响144主反演、候选或选根。

评分复用原R14/R13函数。若只分析已评分结果，读 `analysis/observable19_v1/ALL_VIEWS.csv`、SUMMARY、RECOVERY_GAIN_LOSS、NORMALIZED_*、CORRECTOR_EFFECTS即可，无需打开GT或重复正演。重做离线评分所需真值是父R14独立生成/评分目录里的material_coarse，不是在线输入；原生成出处、seed和合同见原R14 `generate_independent24.py` 与REPRODUCE。不要以GT重选候选或修改阈值。

fixed_state_v2只用冻结七对象×三次readout的21个实际保存共同状态，做126次原LBFGS8条件读出与全物理验收。无GT、没有新基SVD；不加入主候选池，不以诊断结果调整六臂。本轮不扫rank、不换anchor、不新增算法循环。

## 完整原档与网页

本次release中的raw两分卷保存5005个原文件，共885133907原字节；原RESULTS、queues和config逐字节恢复。release同时提供PACKAGE_PLAN/PACKAGE_RECEIPT，源代码及离线状态/分析另作附件。恢复沿用父R14 `src/package_evidence.py restore`，分卷下载到外盘再恢复到新目录；普通ZIP成员/字节检查，不做SHA。跨机器的CPU截断轨迹可能因吞吐和尾部时序不同而略有差异，应复核实际花费与候选资格，不能保证位级相同终点。

```text
python R14_ROOT/src/package_evidence.py restore --plan RAW_PACKAGE_PLAN.json --archive-dir EXTERNAL_DOWNLOADS --destination NEW_EXTERNAL_RESTORE
python S_ROOT/src/build_site.py --root S_ROOT
```

建站仅呈现已存在结果，无新科学计算；所有下载地址必须在公开release匿名回读后写入DOWNLOADS.json。

静态Git页面不重复承载约50MB的PAID_HISTORY.csv及其95MB JSON；原大表和全部verbose cost/selector元数据完整放在release分析分卷，数据没有丢弃。其余逐对象CSV和13张实际图直接网页可读。
