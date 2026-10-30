# FBG-R2 复现与成本口径

## 固定资产

父 B1 commit：`37c9425a3bc95c702440111afb185b4cfccddee2`。A17 查询时远端 HEAD 为 `2ac81b989a2f13d62df9d5a7abd690a33c2189c7`；本地参考的 closure TeX 是非 Git 文件，不能用该 SHA 冒充其逐字版本。一次性 BLAKE2b 来源记录在 audit/BASELINE_MANIFEST.json，各阶段另有算法/source freeze；这些是来源记录，不是密码学盲测证明。

保持原 R4 的 24×24 disk-pulse scalar-TM、k=2π、实材料 box [.02,20]、48 接收器、完整多次散射。主实验12物理源，Stage3另有3/6/24源。粗材料参数始终送入完整576阶物理。

## 运行环境

- 数值 GPU：已有 i9-14900HX / RTX4060 Laptop，Windows，PyTorch2.6.0+cu124，float64/complex128。
- CPU/BLAS/Torch线程1；自己的宿主进程使用原验证过的8个P核逻辑位0/2/4/6/8/10/12/14，不修改其他进程。
- 同一 CUDA 进程4个private stream；复用已有 resident pipeline 与混合 SVD。没有安装CUDA、驱动、CuPy、服务、守护进程或新调度框架。
- 本机数据生成/离线诊断使用项目已有 Gaussian/.venv_nn/bin/python（NumPy/SciPy/Pillow）。系统Python缺NumPy，不能用它生成公共随机起点。

## 精确任务清单

实际任务在 data_contract/*_remote.json；不要把旧 stage1_plan.json 的未执行占位任务当作正式任务。Stage1为64项开发筛选，45秒wall-only；Stage2为360项，另有原物理LU/RHS上限；Stage3为64项照明消融；Stage4为132项新鲜预测（22观测×6方法，K1）。任务执行状态与科学成功由不同字段表示。

在远端研究目录配置已存在的 pipeline 路径：

```powershell
$env:CUDA_PIPELINE_ROOT='D:\AI\ZALG_CUDA_PIPELINE_20261030'
D:\python\python.exe -m src.adaptive_solver stage2_remote.json
D:\python\python.exe -m src.adaptive_solver stage4_remote.json
D:\python\python.exe -m src.adaptive_solver stage3_remote.json
```

每份 manifest 中 operators_path/public_path/output_dir/vendor_path 可改成本机已有文件位置。正式原始观测归档在旧 R4 public_data 与本研究 data_contract；使用归档观测复现同一问题，重新生成 secret realization 是新的实验。

本地采用既有 A9 SSH transport；experiments/remote_dispatch.py 负责有限前台过程与归档。audit/queued_confirmations.ps1 等待**本研究的**Stage2进程结束后立刻接续Stage4与Stage3，不运行两个竞争的CUDA任务进程，不设置永久服务。

## 开发与照明数据

```sh
python experiments/stage0_algebra.py
python experiments/prepare_runs.py stage2
python experiments/prepare_runs.py stage3
python experiments/schur_geometry.py
python experiments/illumination_diagnostics.py
```

Stage0只需读已经保存的审计收据即可理解数学结果；不必为复现主成像重复运行精度检查。Schur/current/projector诊断不作为隐藏GT校准。

## 新鲜盲测边界

STAGE4_FREEZE 在生成前固定六个方法、45秒/2000LU/60000RHS、公共初值.1、共同正则和阈值。公开目录只含匿名Y和采集算子。私有材料及生成状态在仓库外；全部132个预测关闭后，评分器才能读取真值。主solver没有真值路径参数。

生成命令使用 `experiments/stage4_blind_imaging.py generate`，明确传入 public/private 分离目录和 `data_contract/STAGE4_METHOD_STARTS.json`。不要重新生成来替换失败对象。精确评分使用封存预测与随后公开的无种子truth_after_freeze。

## 成本

所有原物理 LU、source/receiver/adjoint/current-tangent RHS 均单列；小Gram、Schur、QR、反馈拟合、试探及失败另记。共享准备记在批次 setup，不伪造为免费。旧Green构造的历史耗时未知，本轮复用缓存，其新构造费用为未重做。

单任务wall是该worker的实测历时，包含与其他worker共享GPU时的等待；不能将所有workerwall相加称并发批次历时。批次SUMMARY保存executorwall，含共享setup与结果写入，排除Python导入、SSH传输。固定总预算的45秒是分配给各起点的wall上限之和，不能冒称独占GPU下实测的45秒端到端时延。

固定参考缓存额外一次LU是新适配器显式付费的准备；记录的setup总计包含原pipeline自身的一次参考LU。返回incumbent可能跨过一次预算边界，实测费用保留。失败不会按0收费；未恢复的早期异常费用标UNKNOWN。

SM occupancy与真实kerneloverlap没有Nsight实测，仍为NOT_MEASURED。普通CUDA stream不是SM配额。本轮不把旧pipeline的43%工程提速当作新的成像速度主张。
