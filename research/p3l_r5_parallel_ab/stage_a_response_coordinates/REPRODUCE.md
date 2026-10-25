# 复现与禁止误重放

已经关闭物理阶段；复现物理脚本会消耗同一累计硬上限，**不要自动重跑**。源 prior/forward 与全部成本由共享 adapter 统一拥有，历史账本和来源只读。

可零新增物理重做：

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/r5_a_mpl "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" code/math_units.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/r5_a_mpl "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" code/test_v2.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/r5_a_mpl "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" code/finish_analysis.py
```

以上命令从本StageA文件夹运行。`host_direct_readback.py`只读取已有V3 endpoints、不构造physics Target，零新增物理。它的pure-layer underflow边界必须保留。

真实物理源码：`code/run_maxwell.py`(V1)、`code/run_v2.py`、`code/run_v3_strata.py`、`code/remaining_spaces.py`、`code/final_response_rank2.py`、`code/final_baselines.py`。固定种子、边界与版本见3份FROZEN_PROTOCOL。endpoint马上持久化，仅compact χ/prediction/q/residual；没有巨型LU输出。精确χ字节相同才复用F；不把response-only cache冒充paid LU。

主输出：18行CSV、final cost/QA JSON、3幅图、完整报告、各版本raw JSON/NPZ与失败日志。未做SHA-256检查。
