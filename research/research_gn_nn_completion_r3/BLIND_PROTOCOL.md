# GN–NN R3：数据隔离与证据范围

本研究使用新生成的匹配全波数据。历史数字只用于代码与机制背景，不进入新的确认集。任何开发联调、验证集选择和确认集结果分别报告。

## 三种信息入口

1. **生成器与监督训练。** 生成器创建材料、完整 Maxwell 观测和评分真值。训练模块只接受 TRAIN/VAL，训练归一化仅由 TRAIN 公共输入计算。材料幅度可用于离线验证指标，不作为网络输入。
2. **公开反演。** GN、推理与物理校正只读取 G/E/S、观测、公开噪声和 GN 端点。在线清单不含形状、数字、幅度、背景、字体、支撑或生成随机种子。远端 GN 作业上传包仅含公开观测和源代码，不含 truth 文件。
3. **离线评分。** 候选与在线排名写入后封存，再由评分器读取 GT。oracle 电流只能出现在离线诊断，不能加入候选集合。确认集在模型、候选尺度、校正规则与在线选择冻结后才评分。

## 分组与划分

初始数据为 256 个训练母体、64 个验证母体、64 个确认母体，另有 8 个仅用于联调的 DEV 母体。每个母体均有 clean 与相对复噪声 1e-3 观测。所有同一母体的噪声、幅度及形变衍生物始终留在同一 split。

确认集包含 32 个同分布与 32 个生成规则/字体外推母体。预指定 32 个确认母体另有 A=2.5 配对外推观测；选择按母体序号 `index % 8 < 4`，使四类形状及 ID/OOD 均衡，不按反演表现选择。低对比主实验共有 128 份确认观测，幅度外推另有 64 份，不能将其合并为 192 个独立材料对象。

主要统计单位是母体。clean/noisy 配对和 A=2.5 配对属于相关重复测量。图像保留失败与缺失条目，不按结果挑选展示对象。

## 训练选择与扩充

五个核心实验臂共享训练母体、骨干、最大训练轮数与种子 7。每臂依声明的训练目标选择验证集 checkpoint。J-null-current 的选择目标是电流误差，不能在确认集用材料 GT 替换 checkpoint。

扩充训练集的主门：在验证集完成共同物理后处理后，在线 top1 的平均幅度归一化 NRMSE 相对 GN 至少改善 5%，且完整物理合格率下降不超过 5 个百分点。这里的“5%”是相对改善，不是 NRMSE 绝对减少 0.05；表中同时保留后者以避免歧义。扩充仍受六小时总执行上限约束。

种子 17/27 与两级 Cascade 是否运行，须在确认集揭晓前根据验证证据和剩余预算冻结。未执行项标记 NOT_RUN，不用推测结果补齐表格。

## 数据合格与材料成功

完整原始 VIE 相对数据残差门为 clean 1e-3、noisy 3e-3。联合成功还要求材料 RMSE/A 不超过 0.15，IoU 不低于 0.7。真值前景由生成器的 alpha>=0.5 定义；重建分割使用重建自身的 Otsu 阈值，不逐图参考 GT 调阈值。

网络只作为学习先验。合格候选中在线选择最接近完整网络材料提案者，另保留旧数据+TV 排名。qualified-best@K 与 all-best@K 都是离线 oracle 指标，不能解释为在线恢复成功。

## 成本口径

基础 GN 构造、共享算子/参考 Jacobian、训练、模型加载、特征/推理、逐候选校正及最终验收分开记录。10 秒是热模型的增量物理后处理预算；它不包含被隐去的免费训练或免费 GN。原始继续 TRF 使用严格数值停止及同额预算，不因初始材料已经满足物理数据门就立即退出。

本机最多两个单线程 BLAS CPU worker；远端 i9 最多四个本研究 P 核 worker，不中断其他 session。GPU 训练在可用设备上执行并记录设备；本轮没有因提速改变材料自由度、正演精度或 GN 预算。

## 尚未获得的结论

本协议不保证 NN 有收益、接收零空间训练更好、确认集唯一恢复、或存在可学习的低维材料流形。完整接收零空间是辅助电流空间；经过共享实材料读出后必须重新验收非线性 Maxwell 数据。
# Transport and archive access clarification

The blind boundary is the numerical information flow into training selection, public inference and offline TEST scoring; it is not an operating-system access-control claim. Generation necessarily creates the test truth files. After the inputs and GN anchors were frozen, the release packager copied their raw file bytes into an archive while the last confirmation jobs were still running. It does not load NPZ arrays, compute GT metrics, rank candidates, or give truth to a solver. Numerical TEST inspection/scoring remains deferred until both the main and secondary prediction cohorts are sealed. Receipt wording such as “GT opened after sealing” refers to this numerical evaluation/disclosure boundary, not to absence of every backup or file-copy operation.

