# 复现本补充

## 1. 来源与固定条件

工作目录是原仓库隔离worktree，分支 `codex/r13-aplus48-20261112`。原R13 C2执行bundle、公开观测、原16 A+、基线48和NN材料均复用，不重新生成、不重训。物理、材料单位和噪声约定继承上层 `DATA_CONTRACT.md`；运行120秒整进程CPU、72proposal/48共同corrector、BLAS1、原八P核、原min-paid-Phi。

真实执行入口是原 `src/finite_queue.py` 与冻结C2 runner。当前R13仓库runner还包含C2之后的参考方法扩展，**不能用它冒充旧bundle**。本补充另保存：

- `src/frozen_C2/run_experiment.py`：实际原C2执行文件；
- `src/stack_runtime/run_experiment.py`：实际112固定版本/NN接口runner；
- `audit/STACK_NUMERICAL_SOURCE_DIFF.patch`：二者只差指定任务合并五个R7最终选项；
- `src/diagnostic_instrumentation.py`：仅诊断复制/记录副本的构造器；
- `configs/`：原缺失32、8对象40诊断、112版本/NN任务及比较规则。

执行源包 `aplus_stack112_execution_sources.zip` 有20个实际source/runtime/config文件。它仍需要原R13/R11冻结依赖与原公开输入，不是独立的新物理后端。原bundle还原方法、依赖和原始下载入口见上层 `REPRODUCE.md`、`ARTIFACTS.md`。本补充的新增结果/分析包由本目录 `DOWNLOADS.json` 与 `ARTIFACTS.md` 给出。

## 2. 本次真正执行的远端队列

三个任务专用目录分别为 `D:/AI/MAX_STRUCT_R13_APLUS48_20261112`、`D:/AI/MAX_STRUCT_R13_APLUS48_DYNAMICS_20261112`、`D:/AI/MAX_STRUCT_R13_STACK_20261114`，共同legacy依赖为 `D:/AI/MAX_STRUCT_R11_DEV_20261107`。

实际有限队列调用形式如下；运行清单分别换成extension32、dynamics8、stack112，完整已执行记录在各终态queue receipt。这里展示stack112真实命令：

```text
D:/python/python.exe D:/AI/MAX_STRUCT_R13_STACK_20261114/src/finite_queue.py --root D:/AI/MAX_STRUCT_R13_STACK_20261114 --legacy-repo D:/AI/MAX_STRUCT_R11_DEV_20261107 --manifest D:/AI/MAX_STRUCT_R13_STACK_20261114/configs/aplus_stack112.json --cores 0,2,4,6,8,10,12,14
```

复跑应使用**新目录**：还原原C2 runtime、放入相应实际runner和配置、从原公开输入拷贝所需NPZ，只将manifest中的绝对执行路径映射到新公开目录，不更改数值/种子/候选规则。若重新跑出的计时停止点不同，应当记录为新重放，不能覆盖本文结果。现有队列、旧结果与seals不得改写。

## 3. 封存、验收和离线分析

先使用上层实际 `src/audit_predictions.py` 检查公开初始化/候选/成本，再由 `src/seal_results.py` 记录终态结果。父研究者审阅后授权，分析程序才读真值；预测选择索引先于GT指标序列化。封存文件中的绝对本机路径和mtime是原执行来源，不是跨机器的加密校验。

移机分析时，先按下载包还原原结果与公开输入，完成同样公开验收，再创建独立的本地路径映射/新seal；保留原seal作为来源。不要只复制原authorized=true来绕过缺失结果，也不要修改原结果时间戳以伪造匹配。此流程不检查SHA。

以下命令中的 `R`、`S` 表示用户还原后的绝对R13及补充目录，不是已经存在的系统变量。实际入口已随代码提供：

```text
python S/src/compare48.py --root R --supplement-root S --authorization S/audit/COMPARISON_AUTHORIZATION.json --out S/analysis/NEW_comparison48
python S/src/compare_stack.py --root R --supplement-root S --stack-root S/stack_run --authorization S/audit/STACK_AUTHORIZATION.json --out S/analysis/NEW_stack_comparison --include-direct-reference
```

两条入口参数已与实际源码核对，详细说明见 `--help` 和 `src/README_COMPARE_STACK.md`。使用原NN材料，不重新选择初值。scorer只读已选候选，不运行Maxwell正演或训练。

物理诊断入口是 `src/analyze_physical_prefixes.py` 与 `src/analyze_native_states.py`，参数与状态语义见各自说明。它们新增的同离散forward/adjoint已单独收费，不叫独立验证。`src/plot_dynamics_summary.py` 只重画已授权CSV；修改图像不重新反演。

## 4. 文件对应与公开范围

三个完整结果ZIP保留raw、corrected、候选池、history、执行回执与queue；动态包保留完整有限状态和原生PR历史。分析包另含gradients、dense scalar表、评分、索引和成本。原NN权重/公开观测复用上层发布附件，训练成本不计零。

本机为通用seal建立了指向原公开输入的 `D:` 硬链接别名，只解决Windows绝对manifest的本机定位；别名树不进入发布包。私有SSH脚本、密钥路径、传输日志和任何新daemon均不属于公共复现代码。

本补充没有新增CUDA移植或E/HT吞吐优化；耗时比实验更久的基础设施改造没有开展。已有CUDA算子目录仍可用于其兼容任务，不能不经适配替换当前冻结算法。
