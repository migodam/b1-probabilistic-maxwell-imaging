# R5 隔离副本评分入口修复

2026-11-03，Sol medium。修复未来复现入口；已交付预测、评分、seal和receipt均未修改。此次GT读取0、物理LU/RHS 0，没有重跑科学实验或改变算法/指标/阈值。

## 缺口与修复

原`offline_score_stage_i.py`依赖预先存在的`audit/STAGE_I_SEAL.json`和`audit/OFFLINE_SCORING_RECEIPT.json`，但原复现说明要求隔离副本不携带这些历史输出，缺少可执行的准备步骤。

新增`experiments/prepare_stage_i_scoring.py`，只用Python标准库读取固定配置和完成回执，不导入物理或读取GT。按冻结run列表枚举replay、diagnostics、两个legacy固定状态组和第一批factorial输出，并保留Stage I配置、protocol/theory与当前src/experiments源码的bytes/mtime。它要求这些已有任务处于原完成状态，仅是文件准备的前置完整性检查，不增加科学gate。

prepare在全新隔离副本创建两个原scorer所需文件，初始`first_GT_read_utc=null`、`truth_reads=[]`、物理费用0，并使用独占新建保护拒绝已有任一seal/receipt。默认ROOT来自脚本；`--root`仅用于指定另一个隔离R5目录，不改变科学设置。

`offline_score_stage_i.py`只新增首GT读取时间的初始化：第一次真实`np.load`打开GT文件成功后，若receipt时间仍为空，立即填入该次打开时间并持久化，再读取m/mask数组。准备时间不再充当GT读取时间，已有历史时间不被重写。其余评分逻辑不变。

## 原执行来源保留

修改前将原执行版完整复制到`audit/scoring_source_original/offline_score_stage_i.py`，并写`PRESERVATION_RECEIPT.json`。这是事后保留原源码，不是历史执行时生成的新快照，不代表再次运行评分。原科学结果与原费用记录继续有效。

## 正确接线

只在全新隔离副本，先完成8条replay、174个诊断和最初2个factorial，再执行：

```bash
python research/research_max_struct_r5/experiments/prepare_stage_i_scoring.py
python research/research_max_struct_r5/experiments/offline_score_stage_i.py
```

extension和handoff仍由各自scorer独立封存。不能在原交付目录删除seal/receipt来强行重跑。

当前`src/trajectory_diagnostics.py`已包含ordinary/self transport解析零与raw_roundoff保留、prior曲率分离修正；当前源码新诊断直接产生该语义。`ZERO_PHYSICS_CORRECTIONS_RECEIPT.json`仅是历史后处理证据，不是复现必须另调的CLI。诊断runner本身仍有物理费用。

## 验证与边界

- 对新增prepare和修订scorer做内存语法编译，PASS；未调用scorer。
- 在自动清理的独立临时fixture中，仅创建最小配置、6个完成回执、protocol/theory及占位源码文件。prepare成功封存12个fixture文件，生成所需两个JSON。
- 确认初始GT时间为null、truth_reads为空、seal/receipt时间一致、每条文件记录有path/bytes/mtime_ns、GT读取标记为false、物理费用0。
- 第二次prepare明确拒绝已有seal/receipt，PASS。
- fixture不包含GT、真实数组或物理算子，没有运行scorer、求解或新优化；不能据此称已完成全新环境的端到端物理复现。

未作hash计算，未推送或发布；修复应作为后续独立提交保留历史研究提交。
