# 复现R7

Python3、NumPy、SciPy、matplotlib、threadpoolctl。实际环境沿用本地Gaussian/.venv_nn；每进程BLAS1。无CUDA、NN、后台服务或新基础设施。源码从仓库根启动。R1–R6依赖须同时存在，R7没有修改它们。

## 获取原始证据

Git仓库保存新方法、配置、逐例预测/评分、报告、图、审查和回执。大型逐步history/paid-state/budget数组在本轮GitHub Release复现档案；下载后覆盖解压到仓库根，不覆盖R1–R6。实际Release链接在CHAT_HANDOFF和发布回执中给出。第三方原始测量/全文不随档案分发；按官方脚本重新取得。

## 生成与反演

所有写入目录为write-once；默认复用已发布观测。单任务用--output-dir指定新目录；整队列先复制manifest，只修改output_dir为新目录，再运行。若要从seed重新生成，请在独立副本将同名生成产物移到备份目录后执行生成器，不能在证据原目录覆盖。保留原结果。

```sh
python research/research_max_struct_r7/experiments/generate_benchmark.py
python research/research_max_struct_r7/experiments/run_case.py --manifest research/research_max_struct_r7/configs/primary_manifest.json --task-id r7_obj_00__cie_dual --output-dir /absolute/new/run
python research/research_max_struct_r7/experiments/run_queue.py --manifest research/research_max_struct_r7/configs/first_four_manifest.json --name reproduction --workers 4
```

最终A+专用冻结config/manifest、初始化manifest及队列回执与原primary文件一起保留；不能把开发途中配置当最后算法。使用最终配置复现A+，不用默认未冻结参数。有限队列只在执行异常时停止继续派发，科学失败不会被抹除。

## A+开发、有理搜索和独立物理

A+完整开发顺序/命令在delegated/aplus/DEVELOPMENT_AUDIT.md和该目录manifest/notes。它只用四个已公开对象，不能计新盲测。每个版本保留独立source/config/result快照。

```sh
python research/research_max_struct_r7/experiments/run_rational.py
python research/research_max_struct_r7/experiments/run_rational_coherent.py
```

这两条命令固定R6三个公开失败锚点，原输出已存在时会拒绝覆盖。数值费用和source路径见各RECEIPT。前者的逐端点budget文件记录串行累计费用，不是独立matched-method曲线。

三维实际命令、初始与长校正的区别、Mie几何检查和矩阵维度见BENCHMARK_3D。实测下载、缺失通道选择、target-absent校准、训练/留出通道及调用入口见FRESNEL_AUDIT；无精确像素GT评分。

## 封存后独立评分

评分必须在180个主预测全部完成并形成PRIMARY_PREDICTION_SEAL.json之后。评分代码不会被反演导入。重新生成的GT只给独立scorer，所有算法先封存，再评分；不得在线读真值。

```sh
python research/research_max_struct_r7/experiments/score_primary.py --help
python research/research_max_struct_r7/experiments/analyse_primary.py --help
python research/research_max_struct_r7/experiments/score_rational.py --help
```

最终实际调用参数由audit评分回执与运行日志记录。匹配预算评分复用已支付原Phi incumbent的float32材料，选last checkpoint <=CPU/LU/RHS预算，不按GT选迭代；最终材料float64。无新增物理调用，bootstrap按30个对象。

没有SHA核验工作流。git_base是来源记录，不等于已验证历史科学结论；执行成本和科学资格以真实回执与报告为准。
