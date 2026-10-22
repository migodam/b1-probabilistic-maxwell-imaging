# 复现：缓存结论与真实全波分开

## 1. 环境与原始数组

在本目录创建独立Python环境并安装`requirements.txt`。原运行是Apple Silicon CPU、BLAS单线程，未使用GPU。精确版本和硬件见`docs/ENVIRONMENT.json`；Python 3.13。不安装CUDA、调度器或新后台服务。

大型复场数组位于同仓库Release **research-w-continuation-20261023**。下载`research-w-continuation-arrays.zip`，在本目录解压，恢复`results/.../*.npz`。`ARRAY_MANIFEST.json`列出路径、字节数及用途；小型JSON/CSV、代码、运行日志直接在Git中。没有数组时绘图应明确失败，不能用合成数据替代。

完整复观测保存了实际运行的所有频率。较大的网格场、材料积分点和复总场保存代表性的中心频率；其他频率的完整场未全部持久化，不能将Release说成每次solve的全部中间张量档案。其观测、参数、求解残差和费用仍在原始结果与账本中。没有为补齐展示数组而重复求解。

```sh
python -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
export PYTHONPATH="$PWD/src:$PWD/vendor"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
```

## 2. 零新增物理的结果复核

下列命令只读已有复数组/标量，不执行Maxwell求解。会重新生成派生图表或判决，因此请先复制下载目录到独立复现目录，保留原发布数据。

```sh
python -m unittest discover -s tests -v
python tools/qualify_q.py
python tools/plot_q.py
python tools/plot_b.py
python tools/final_report.py
```

`tools/score_b.py`仅在所有在线候选保存之后离线读取真值。它生成材料/结构评分，不参与提案或接受。已发布评分可直接查看，无需重新运行它。`final_report.py`从原始追加式阶段账本生成根目录账本快照；新绘图/发布尾项可能晚于该快照，均另保留作业收据。

## 3. 有预算的真实物理入口

所有物理动作在`Book`/`SearchBook`先预约，不能绕过账本调用裸solver来规避费用。`tools/run_job.py`限制process CPU，设置单线程，并保存stdout、退出码、CPU、wall和峰值RSS。以下是本轮实际使用的入口；完整命令与次序见`results/jobs/*/RECEIPT.json`。

```sh
python tools/run_job.py replay-attribution python -m wq.experiments attribution
python tools/run_job.py replay-kernels python -m wq.experiments kernels
python tools/run_job.py replay-galerkin python -m wq.experiments galerkin-cylinder
python tools/run_job.py replay-independent python -m wq.high_precision independent
python tools/run_job.py replay-wall python -m wq.high_precision wall-qualification
python tools/run_job.py replay-resume python -m wq.high_precision wall-qualification-resume
python tools/run_job.py replay-derivatives python -m wq.high_precision fem-derivatives
python tools/run_job.py replay-pml python -m wq.high_precision wall-pml
python tools/run_job.py replay-circle python tools/circle_controls.py
```

这不是要求再次执行整个campaign。已有结果应优先复用；新复现实验必须放在独立目录和独立预算中，不能清空或重置原费用。自适应网格/内存预检会保留未执行项，原h=.13墙体配置因60k DOF预检拒绝后采用已记录的h=.145层级，没有改变物理世界或0.1门槛。

P0-Q通过后才有P0-B权限。P0-B入口：

```sh
WQ_STAGE=p0b python tools/run_job.py replay-preflight python -m wq.search preflight
WQ_STAGE=p0b python tools/run_job.py replay-search python -m wq.search search-corrected
WQ_STAGE=p0b python tools/run_job.py replay-task python -m wq.search task-endpoints
WQ_STAGE=p0b python tools/run_job.py replay-finer python tools/validate_b_endpoints.py
WQ_STAGE=p0b python tools/run_job.py replay-fitted python tools/validate_fitted_tasks.py
```

`search`是留存的错误历史版本，**不得用于算法比较或新实验**；仅`search-corrected`是修正版。原失败搜索消耗760 RHS不退回，因而本轮修正版只完成三个context，第四个明确NOT_RUN。干净复现实验若不重演该故障，其剩余额度和可执行context数会不同，不能把不同调度的结果说成原历史逐字重演。严格的原结果复核采用Release数组和保存的修正版轨迹，而不是重新制造一次失败。

端点脚本按保存的在线cal距离和稳定方法ID选定三个端点；不以truth选择。固定目标8仅覆盖1/24任务类。原未贴合目标结果和后续贴合圆周256/512段结果同时保留；joint、全部24任务、P1及三维未执行。

## 4. 数值和证据边界

全部单位、频率、材料、噪声和源/接收器位于`docs/P0_FROZEN_CONFIG.json`，物理资格和搜索配置分别见`P0_Q_PROTOCOL.md`、`P0_B_PHYSICS_PROTOCOL.md`。白化为完整实数据范数，sigma=0.01，不除以观测维数平方根。没有伪逆、隐藏jitter、自由复增益拟合或true wall减法。

相邻网格和解析圆柱支持有限的离散全波资格，不是连续PDE认证；新材料/目标需要独立审核。所有四方法搜索仍远未达到cal距离1；没有joint或posterior结果，也没有全局不存在/唯一性结论。
