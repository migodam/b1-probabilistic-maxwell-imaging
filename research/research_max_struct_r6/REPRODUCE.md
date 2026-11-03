# R6 复现

从主仓库根目录运行。R6是可直接调用现有研究代码的增量包，不是自包含的新物理后端；须保留以下既有目录：

- `research_fbg_r2_adaptive_feedback`、`research_fbg_r3_experiments`：公开算子/观测、成熟CIE及物理API。
- `research_z_alg_r3_digits`、`research_z_alg_r4_diversity`：共用公开算子、原完整VIE与prior。
- `research_max_struct_r2`、`research_max_struct_r4`、`research_max_struct_r5`：global CIE、one-hop与材料块拆分实现、历史对象。
- 新的 `research_max_struct_r6`：固定公共任务清单、独立新观测、四臂与FFT运行器、评分。

环境实测版本在`audit/ENVIRONMENT.json`。依赖Python、NumPy、SciPy、Matplotlib、threadpoolctl；现有项目虚拟环境可直接使用，不需要CUDA、新服务或NN。各worker自己固定BLAS单线程。

## 只读结果

`RESULTS/scoring_primary/SCORES.csv`是全部84任务的固定最终预测评分；`RESULTS/scoring_fft/`单列4个重复观测的工程任务。原始数组/历史在`RESULTS/primary/`、首例元数据修复后的`RESULTS/primary_retry/`和`RESULTS/fft/`；图不替代这些数据。`configs/primary_final_manifest.json`是包含这次路径修正的完整84任务索引；原manifest和失败目录保留。预测文件不得按离线评分替换。

## 从公共观测重跑

默认已交付的输出目录存在。运行器采用不覆盖写入；单例必须指定新的输出目录：

```sh
python research/research_max_struct_r6/experiments/run_case.py \
  --manifest research/research_max_struct_r6/configs/primary_manifest.json \
  --task-id obs_b20d5a55677d4eca__a_material \
  --output-dir /absolute/path/to/new_reproduction/a_material

python research/research_max_struct_r6/experiments/run_fft_case.py \
  --manifest research/research_max_struct_r6/configs/FFT_manifest.json \
  --task-id obs_b20d5a55677d4eca__a_material__fft \
  --output-dir /absolute/path/to/new_reproduction/a_material_fft
```

全量复跑时复制`primary_final_manifest.json`，只更改各task的output_dir到新的目录；数值config保持不变。`run_queue.py --manifest NEW_MANIFEST --name NEW_NAME --workers 4`运行主四臂；FFT使用`run_fft_queue.py`。队列name须新建，避免覆盖冻结回执。不要用已发布GT筛选要跑的对象或修改参数。

所有公共路径相对仓库根；旧回执中的绝对工作树路径只是当时环境记录。程序从自己的路径定位仓库，不要求旧磁盘路径。实验源码快照位于`audit/source_freeze/`，没有SHA验证步骤。

## 独立评分

反演只打开manifest所列公共文件，不导入生成器或评分器。新预测全部完成后，先创建新的public seal，再在独立评分端显式授权读已公开的truth。重跑manifest中的task_id/arm/case_id须保持一致，output_dir可改变。

```sh
python research/research_max_struct_r6/experiments/score_independent.py seal \
  --manifest NEW_MANIFEST --output NEW_SEAL.json
python research/research_max_struct_r6/experiments/score_independent.py score \
  --seal NEW_SEAL.json --output NEW_SCORES --batch primary --authorize-after-seal
```

封存记录文件大小/修改时刻，作用是本地交接，不是密码学认证；旧绝对路径seal不可拿到另一机器直接验证，须给新的重跑结果创建seal。评分仅消费已付费保存的F(m)，没有额外物理LU。数据噪声σ表示复幅度RMS，每个实虚分量σ/√2；新数据discrepancy使用1152实自由度chi-square的99%半径。旧clean/旧历史门分开列，不冒充同一噪声资格。

生成器`generate_confirmation.py`属于真值保管端。无需重新生成观测来复现算法；公开NPZ是本轮实际使用数据。现在真值随交付公开，这六例之后只能称公开回归对象，不能再次称新盲测。

## 费用与解释

每个任务RECEIPT包含imports/setup、完整辅助日程、corrector、保存成本；其子计时有重叠，不能相加为总时间。物理LU和RHS以task receipt为唯一计账；辅助FFT/Green和材料读出另列。四worker任务wall之和不是队列历时。首次生成、接线、失败、评分/作图/打包的成本在COST_LEDGER及audit保留。

只重画/分析可使用`analyse_benchmark.py`，不需要正演。HTML与证据ZIP由`export_report.py`生成；该步骤无科学重跑。已有解析恒等式、Mie和历史全库tests不随复现重复运行。
