# R9远端执行与费用口径

用户在本轮执行中明确要求尽量使用i9大核及优化后的4060 CUDA算子。主线程实际连接XINAN后，确认远端i9-14900HX、32逻辑核、约22GB空闲内存、RTX4060 Laptop约8GB显存；当时GPU零占用且没有其他Python研究任务。只调整本任务进程，不修改其他会话、系统驱动或常驻服务。

Windows `GetSystemCpuSetInformation` 实测逻辑0–15为EfficiencyClass1的P核，16–31为Class0的E核。每个P核仅选一个SMT线程，分配如下：

|工作|逻辑P核|设备及同硬件比较范围|
|---|---|---|
|经典Austria|0、2、4、6|四CPU lane，完整16臂Windows对照；先前Mac部分结果单列|
|可分辨空心立方体|8、10|两个CPU lane，新增五臂全部Windows；原24小3D仍是完整Mac队列|
|原3λ三维|12、14|一个专属CUDA上下文及其CPU支持；完整RHS资格已执行，随后为cross-grid五臂及matched-grid四臂|

原始探测回执：`audit/REMOTE_RESOURCE_PROBE.json`、`audit/REMOTE_CPU_TOPOLOGY.json`。远端现有NumPy2.2.6、SciPy1.15.3、PyTorch2.6.0+cu124、psutil7.2.2；没有安装CuPy、CUDA驱动或新调度框架。CPU数值库线程数在导入前设为1。远端使用有限、持续连接的SSH进程，避免Windows断开连接后作业被终止。

## 复用与CUDA范围

已读取仓库既有CUDA算子目录。576单元scalar-TM的驻留LU接口不是64²或3D矢量FFT的直接替代；A18已有Torch/cuFFT传播形式，但原接口每次matvec回传SciPy GMRES。原尺寸扩展因此复用其传播形式和成熟Torch/cuBLAS/cuFFT运算，将Krylov向量留在显存。新实现资格以一个完整RHS的真实残差、预测差和费用为准，不凭独立FFT计时宣称加速。

已采用的GPU后端使用complex128、restart60、最多20周期的左预条件GMRES，CGS2正交化，最终必须满足原未预条件残差要求。每八步检查可能多做少量迭代，属于明确记录的数值后端差异；不改Maxwell矩阵、原反演预算、材料集合或物理精度。原尺寸生成的剩余RHS及科学反演已实际使用它；全部反演是否完成、GPU计时及显存峰值仍以分支终止回执为准。

数据生成中，已知零极化电流的外部网格可精确消去：60³物理网格及h=.05、绝对坐标、材料和观测不变，只有FFT零填充区域按非零电流索引包围盒裁剪。该优化只在生成端使用；反演始终保留全部27000未知材料。CPU完整/裁剪算子接线与已付源计时均保留，未将它记为新反演算法。

## 成本不能混算

### 已执行的工作量比较

1. 原3λ fine源59、同complex128矩阵与2e−6原残差要求：i9单P核SciPy MGS为99.9054 wall秒，resident CUDA CGS2为5.8905秒，约16.96倍。GPU解回到CPU原算子核验残差1.8652e−6，端口相对差6.59e−6；Krylov/check cadence不同，不能称纯GPU硬件比。峰值分配约566MB；整体资格过程106.922 CPU秒、109.924 wall秒，含CPU对照。
2. 固定公开状态的一次完整普通CIE三维current＋material更新（β=6）：strided4为2.97547秒，packed16为1.02073秒，约2.915倍，所查材料/current/field/系数/步长结果一致。完整benchmark成本8.859 CPU秒。旧strided运行616 CPU秒仍单列保留；所有正式九臂统一新布局，不能用新快A+对比旧慢CIE。
3. 经典二维4096未知量、20源、40接收，复用已付Mac第80步公开材料的一次完整目标/梯度：CPU FFT/GMRES为1.39984秒，既有SharedOperators/ObjectiveLane CUDA dense LU为1.10259秒。梯度相对差3.36e−8，GPU原状态残差2.93e−15；cold setup34.6424秒全部收费，包括原接口多建的reference LU与未用full basis。单GPU的1.27倍回调收益不足以说明比四个CPU lane吞吐更好，因此不重启16臂或更换其后端。

回执分别在`RESULTS/original_vector3d/CUDA_WHOLE_RHS_BENCHMARK.json`、`CUDA_LAYOUT_CURRENT_STEP_PROFILE.json`及`RESULTS/classic_extension_cuda_callback/RECEIPT.json`。第三项只有一个固定状态，不代表整个优化过程的吞吐或材料恢复。三类工作量的速度比不能互相套用。

- 原Mac已付计算、资源迁移中断、打包传输、GPU资格检查和失败尝试全部保留。迁移中断是工程事件，不计为材料恢复失败。
- 同一物理对象的跨平台重跑不是新的独立材料样本。经典16臂Windows结果单独配对，不能把Mac部分CPU秒混进其平均。
- Windows完整进程CPU使用GetProcessTimes或psutil的user+kernel；Mac使用原完整进程口径。GPU执行同时报告wall、CUDA事件时间和显存，不能只看CPU秒。
- 多lane的任务wall之和不是并行历时。共享预计算只计一次，并注明摊销。物理尺寸不同的LU/RHS/FFT不能假装等FLOP。
- CPU与GPU共享笔记本功耗和散热预算，批处理中的实际时间并非隔离硬件微基准。保留作业并发与实际计数；不将一次观测比例宣传为固定加速倍数。

GT只用于生成或封存后的离线评分。普通远端反演包只含必要源代码、公开观测和配置，不随包传测试标签。原尺寸观测准备从sealed命名容器中只读取Y_clean、noise和coarse prediction，不索引材料/support；实际反演只读取公开观测，不能宣称整个准备进程从未打开该容器。

最终所有队列均已结束，GPU已释放。原尺寸九臂均使用优化后的packed16布局，后续三个八步物理校正同样使用单驻留CUDA上下文，新增1740正向、1680伴随RHS。该组主进程1557.95 CPU秒／1583.52 wall秒；逐源GMRES，没有共享LU。经典16项Windows四P核队列wall7952.01秒；任务CPU总量30037.86秒不是并行历时。完整成本见COST_LEDGER_EXTENSION.md。
