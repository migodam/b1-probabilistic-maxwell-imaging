# 一个真实64²公开状态的完整CUDA回调测量

状态：**COMPLETED；GPU已释放（compute_apps空、0 MiB）**。输入是在GT评分前预定的Mac已封存εr=3.5 single CIE连续80步检查点，不是另选真值更好的材料。没有读GT、训练、优化或重复有效callback；原16个Windows CPU任务没有修改或重启。GPU独占时使用host P12，complex128/float64。

| 工作 | wall秒 | CPU秒 | 物理工作 |
|---|---:|---:|---|
| CPU FFT-model setup | 0.06498 | 0.06250 | Green kernel一次 |
| **一次CPU完整value+gradient** | **1.39984** | **1.37500** | 2240 GMRES迭代，40 RHS，无LU |
| GPU冷setup合计 | **34.64241** | **31.14063** | 包含以下三个分项，不能再相加到合计 |
| 其中Torch/context/library | 3.07429 | 1.35938 | Torch2.6.0+cu124，局部cuSOLVER选择 |
| 其中原模型dense G组装 | 28.09880 | 27.90625 | 精确相同等面积圆盘自项；G为268435456bytes |
| 其中GPU上传／resident setup | 3.46025 | 1.87500 | 未改SharedOperators：额外reference.1 LU一次、20源RHS；包含未用B/DB分配 |
| **一次GPU完整value+gradient** | **1.10259** | **1.09375** | 同一LU供20source+20adjoint RHS；一次host返回 |
| GPU缓存state诊断 | 0.02678 | 0.03125 | 无额外LU/RHS；检查真实物理残差 |

整callback比值仅约**1.27×**（CPU wall／GPU wall），不是全优化器加速比。GPU callback事件为1102.65039 ms；GPU setup事件时间未测，不能把这1.10265秒叫总GPU使用时间。峰值allocated 1777947136bytes、reserved 2831155200bytes。benchmark整个进程内部CPU为33.953125秒、wall37.69335秒；表内组件都已包含其中，不能重复加。传输与外部helper另见 `CUDA_TRANSPORT_RECEIPT.json`（总wall41.26150秒，包含远端工作，不是另一份可相加elapsed）。

CPU／GPU objective为0.00435296887814／0.00435296887995；相对差4.1490e−10。梯度相对L2差3.3619e−8，绝对L∞差5.2248e−13，与CPU GMRES约1e−9状态容差相容。CPU／GPU真实state相对残差为6.4302e−10／2.9286e−15；数据残差同为约0.09330044。没有设置额外材料恢复成功阈值。

这是**4096标量材料、20平面波、dense CUDA LU对CPU FFT/GMRES**的一个已有状态。父线程另测的三维GPU FFT/GMRES约17×属于不同算子尺寸、向量场和工作量，不能搬到本二维dense-LU callback。此处的小幅resident收益与冷setup费用不足以支持重启或自动迁移正在运行的完整16项CPU cohort；按父线程决定继续四P核有限队列。它也不证明任意后续材料点、整段corrector或多lane吞吐只有同一比值。

原始来源：`RESULTS/classic_extension_cuda_callback/RECEIPT.json`、`CPU_VALUE_GRADIENT.npz`、`VALUE_GRADIENT.npz`；实现复用 `research/research_z_alg_cuda_pipeline/src/device.py::SharedOperators/ObjectiveLane`，新wrapper为 `cuda_callback_benchmark.py`。实际费用：GPU setup1LU20source RHS＋有效callback1LU20source20adjoint RHS；CPU比较另40 RHS。完整冷setup、有效回调、缓存诊断、退出及失败记录均保留，没有手写LU/FFT内核或改变默认生产库。
