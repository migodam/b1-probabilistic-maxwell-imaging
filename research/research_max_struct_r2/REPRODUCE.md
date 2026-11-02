# 复现

依赖既有B1仓库内的R1 vector backend、FBG-R3 mature CIE/corrector、R3/R4 scalar physics；没有新安装求解框架。请使用本次发布对应的完整仓库，保留研究目录相对位置。Python、numpy、scipy、matplotlib、threadpoolctl足够；实际使用本机既有Python3.13环境。

原始结果已经保存，不必重跑科学实验来阅读结论。若要重算，在**另一个复制目录/checkout**执行，保留本次原始目录；runner会跳过已有预测，复跑前只在新副本移走 `results/blind/` 与 `predictions/obs_*`，仍保留公开数据和原配置。

从 `research/research_max_struct_r2/` 运行：

```sh
python experiments/run_blind.py
python experiments/diagnose_endpoints.py
```

第一条依次执行3对象×4臂；单例可加 `--case obs_1c40f624b9334a6d --arm native`。使用同一已保存Y即可精确复现执行数据；运行耗时依机器变化。程序包含一次新profile导数接线检查，已有回执存在则跳过。不要以再次生成不同几何替换本次观测。生成器默认使用新熵，是生成**新的**盲测的入口；本次seed/材料/观测均已在封存评分后公开保存。

三维各阶段是独立进程：

```sh
python experiments/active_family.py plan
python experiments/active_family.py acquire
python experiments/active_family.py fit
python experiments/active_family.py score
python experiments/active_family.py coverage
python experiments/plot_active.py
```

`plan`不读真值；`acquire`是模拟采集端，允许生成真值；`fit`只读取public新增观测；`score`在预测写出后访问GT并做一次细网格复算；`coverage`是明确的模型遗漏开发控制。它们会覆盖副本同名结果，故不得在本次封存原目录运行。`score`的细化与反演同类DDA，不能称独立连续求解器。

`src/cie_global.py`复用成熟Fourier表示与readout；`src/profile_current.py`是标准变量投影控制。`experiments/run_blind.py`中的三阶段、正确器预算和候选规则是实际执行配置，不需要猜测历史入口。所有零/未知费用的含义见COST_LEDGER；不计算源码/配置SHA。

下载论文全文只作本地查阅，GitHub不再分发第三方完整PDF/txt；DOI、原文定位、检索/核验收据及自己的推导随包提供。Release研究ZIP需要放回同一B1仓库的research目录才能使用历史依赖，不称完全独立单文件软件。
