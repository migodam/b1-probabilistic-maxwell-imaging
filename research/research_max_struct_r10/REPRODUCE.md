# R10 复现说明

代码复用完整仓库的R1/R2/R4/R6/R7/R8/R9、FBG-R3及Z-ALG物理模块；只复制R10目录不足以运行。旧研究不覆盖。完整科学数组将按同仓库R10 Release分包，解压到仓库根目录；不要创建原机绝对路径。发布前后的状态以CHAT_HANDOFF和PUBLICATION_MANIFEST为准。

实际环境：Mac arm64 Python3.13.13、NumPy2.5.3、SciPy1.18.1、Torch2.14.0、threadpoolctl3.7.0；远端Windows Python环境为NumPy2.2.6、SciPy1.15.3、Torch2.6.0+cu124。使用既有虚拟环境入口，不解析其Python软链接为缺依赖的基础解释器。主二维反演每进程单BLAS线程，未用NN或训练。三维CPU/CUDA实际选择见vector回执，不把CUDA可用等同于已使用或更快。

以下命令从仓库根目录执行，`python`表示相应已有兼容环境。新路径包装只做路径/参数接线与语法检查，没有为验证包装重新跑一遍科学实验。

## 新二维单例重放

无需重新生成标签；使用随附公开观测：

```sh
python research/research_max_struct_r10/experiments/twod_run.py \
  --task-id r10_twod_0000__cie_single__noise_ball_99 \
  --output-dir research/research_max_struct_r10/RESULTS/twod/repeat01/cie_single_0000
```

新矩阵分母候选传入`--manifest research/research_max_struct_r10/configs/denominator_confirmation_tasks.json`，task-id取该公开清单的真实ID，并给新的output-dir。不要将canonical arm=a_plus误当普通A+；候选由material_readout.label区分。任何重放输出都须是未存在目录。

完整168次的公开重放清单可以准备为：

```sh
python research/research_max_struct_r10/experiments/prepare_reproduction.py --name repeat01
```

准备只写新输出路径，没有读取GT或运行反演。Windows有限队列例子：

```powershell
python research/research_max_struct_r10/experiments/twod_queue.py --execute `
  --manifest research/research_max_struct_r10/reproduction/repeat01/twod_tasks.json `
  --receipt-dir research/research_max_struct_r10/reproduction/repeat01/baseline_queue `
  --cpus 0 2 4 6

python research/research_max_struct_r10/experiments/twod_queue.py --execute `
  --manifest research/research_max_struct_r10/reproduction/repeat01/denominator_confirmation_tasks.json `
  --receipt-dir research/research_max_struct_r10/reproduction/repeat01/candidate_queue `
  --cpus 8 10
```

这些CPU编号来自本次i9布局；换主机必须选择该主机真实可用的P核，不宣称固定编号普遍适用。不要占用别的任务资源。单例脚本不要求Windows affinity，Mac也可运行。详细打包、远端启动、终态fetch及GetProcessTimes成本入口见[二维执行记录](delegated/twod/REPRODUCE.md)。

两队列全部终态且预测封存后，运行`experiments/score_twod_comparison.py`，明确传入新baseline/candidate manifest、两个新queue RECEIPT.json、独立output和figure-dir，以及`--authorize-postseal ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION --authorization-note`。不能先评分部分对象来调整剩余任务。评分复用付费预算点，不为GT优选重新求解。阈值合同固定在`delegated/twod/FREEZE_SCORING.json`。

## 历史开发矩阵分母重放

先恢复R9 Release内的PUBLIC_TASK及公开数据，用当前目录而不是历史工作树位置：

```sh
python research/research_max_struct_r10/experiments/denominator_run.py \
  --historical-repo "$PWD" --case obs_40ac5b79e8a74ac3 \
  --basis green_strong --projection galerkin \
  --public-files research/research_max_struct_r9/data_public/development/operators.npz \
    research/research_max_struct_r9/data_public/development/obs_40ac5b79e8a74ac3.npz \
  --output-dir research/research_max_struct_r10/reproduction/dev40ac_green19
```

`--basis twofold`、`green_equal_rank`或`--retain-self`对应已运行的机制控制，不能在确认集按GT选其中最佳版本。`onehop_solver_control.py`具有相同historical-repo/public-files/output-dir接口，用原one-hop公式和8步L-BFGS。首轮错误传参的两次中断结果保存在`RESULTS/onehop_solver_control/FAILED_ATTEMPTS.json`，不删除、不当成正式8步结果。

`denominator_fixed_state.py`、`delegated/denominator_spectrum/analyze_saved.py`、`analyze_historical_energy.py`是已保存状态的诊断入口。它们中的历史来源保持执行时路径，迁移时应指向已恢复的相应R9/R10文件，不能猜测同名目录。完整结果数组可直接复算表格，不需要重复Atlas恒等式或全部历史反演。

## 三维

实际合同与任务见`configs/vector_actual_contract.json`和`configs/vector_all18.json`；早期proposed计划不等于执行清单。原始运行入口：

```sh
python research/research_max_struct_r10/experiments/vector_run.py invert \
  --backend cpu-dense --queue configs/vector_first_queue.json \
  --max-tasks 3 --task-offset 0 --output RESULTS/vector/reproduction_v02
```

所有三维source/current/material仍是全矢量物理与完整1728材料变量。若使用CUDA，必须复用已记录的兼容环境及实际通过的精度/线性残差设置。`--backend cuda`与CPU之间不承诺逐位相同轨迹；同一个配对case三臂必须同backend。三维新输出的评分、跨网格映射与实际远端命令见vector最终交接；不可把旧绝对路径或默认拟议任务数量当作已执行数量。

## 成本及归档

R10中顶层进程费用与内部FFT/材料/corrector计时是包含关系，不能相加两次。576物理LU、5184矢量LU、19/28分母小LU分列；不同source RHS数量也不能只按“一个forward”计价。多个并行进程wall相加不是历时。没有硬件费率/工具账单的金额记null。

`experiments/package_evidence.py`在所有科学写入结束后无损归档，保留失败、部分状态、配置和原始指标；不发布第三方全文缓存或重复传输ZIP。清单只记路径/字节数，不做SHA验证。发布成本和最终文本以Git版本为准，避免为了自指的发布回执重新压缩全部科学数组。
