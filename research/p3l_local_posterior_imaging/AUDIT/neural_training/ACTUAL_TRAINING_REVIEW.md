# 实际神经训练审计：选择了初始化，未证明学得改进

本审计只读已保存的 15 份训练日志、512 train / 128 validation 的已允许 packed scalar/image/labels/parent、固定 frame 和已保存 BEST checkpoint。没有训练更新、物理、locked 数据或新的 GT 文件读取；checkpoint 前向仅用于重现现有 NLL。最终推断与算法变更由主线程决定。

## 已经发生了什么

五种架构、各三个冻结 seed，全部训练到 update 1000，因为十次完整 validation 均未改善而触发 `VALIDATION_PATIENCE`。全部选择 update 0。独立恢复 15 个 checkpoint 后，其 base weights **逐项等于相同 seed 的原始初始化**。因此实际被保留的是 `UNTRAINED_SEEDED_INITIALIZATION_SELECTED`，不能称为训练成功的 learned improvement。

下表给出三个 seed 的范围。训练列是 update 1000 当次、参数更新前的随机 minibatch NLL；它不是完整训练集 NLL，也不是被恢复的 update-0 模型损失。连续密度 NLL 可以为负，负号本身不是实现错误。

| 架构 | 参数量 | update0 validation NLL | update100 validation NLL | update1000 train minibatch NLL | update1000 validation NLL |
| --- | --- | --- | --- | --- | --- |
| mlp_diag | 461844 | 16.8938 到 17.3002 | 45.8170 到 60.9038 | -13.7869 到 -11.5723 | 906.9516 到 1220.1354 |
| mlp_full | 467649 | 17.3396 到 17.6570 | 32.3414 到 47.5924 | -3.4073 到 -2.7674 | 293.0337 到 1639.6373 |
| mlp_mdn | 493320 | 16.8782 到 17.0873 | 24.0306 到 26.6590 | -0.0351 到 0.1996 | 57.1860 到 73.4882 |
| mlp_flow | 606800 | 15.2620 到 15.4779 | 57.9103 到 69.9384 | -13.6664 到 -12.9379 | 1181.9048 到 1240.7002 |
| cnn_full | 467441 | 17.4662 到 17.5198 | 36.9892 到 40.3057 | -4.6549 到 -2.4542 | 262.3289 到 325.0035 |

训练 minibatch 拟合变好、validation 大幅变差，是明显的样本内拟合与泛化失败模式，与过拟合相容。约 46–61 万参数面对 512 个独立训练 parent，足以解释为何存在记忆样本的可能性；参数量本身不是因果证明。不能仅靠这些日志判定条件信息不存在、decoder 不可识别或模型族必然无效。validation 每 100 次检查，1–99 次之间是否有短暂改进没有测量。

## 数据和损失合同核对

- packed arrays 全部有限，尺寸为 train `[512,1154] / [512,4,24,24] / [512,10]` 和 validation `[128,1154] / [128,4,24,24] / [128,10]`。parent 顺序精确等于冻结 train/validation 名单，各自唯一，两组无交集。
- 训练标签十个维度 SD 为 0.9342–1.2073。采用训练尺度后，validation 逐维均值绝对值最大约 0.1806，SD 范围 0.9140–1.1485，单点最大绝对值约 4.2162。没有发现 gross label shift 或尺度爆炸。
- scalar 四个块分别为 whitened observation、whitened residual、anchor q、log1p(Jz singular values²)。train/validation 标准化后的最大绝对值分别约为 observation 4.513/4.704、residual 7.099/6.209、q 7.603/5.063、spectrum 8.449/7.811。存在 28 个 exact-constant scalar 分量，代码使用 scale=1；最小非零存储尺度约 7.068e-11，但对应 validation 标准化统计未显示爆炸。图像按合同不标准化，四通道值有限；这可能影响归纳偏好，不能据此直接判为 bug。
- 重新由 packed train 拟合的四组 affine mean/scale 与全部报告逐项相同，最大误差 0。代码没有拿 validation 拟合尺度。raw-label NLL 为 standardized NLL 加 `sum(log(label_scale))=0.5086662`，符号正确；这个常数不会改变同一次训练的最佳 checkpoint。
- `batch(indices,0)` 读 train 的 scalar/image/labels 三块，`batch(slice,3)` 读 validation 三块，组内共用相同索引。随机训练索引有放回抽样，是已实现的 minibatch 策略，没有发现错位或跨组取数。
- 标签源代码为 `[Vw.T@(ztruth-zcontext), ndtri(utruth)]`；恢复 chart 为 `z=zcontext+Vw@a+Vs@b`。没有发现反号或把 eta 错当相对 anchor 属性的转换。正交 frame 残差 2.653e-15，Q 与 `[Vw,Vs]` 精确一致，绝对 determinant 约 1。没有读取 raw truth 重新验证每条标签，所以这是 source/packed-contract 审计，不是完整端到端 GT 重算。
- checkpoint 逻辑把完整 validation 的最低 NLL 作为候选，明确包括 update0；使用 clone 的 best state，最终恢复同一 state。独立前向重现所有已选 validation NLL，最大绝对差 1.513e-6，符合 float32 求和方式差异。没有发现选错 epoch、损失符号错误或保存了最后一步权重的证据。

## 一个必要的量级参照

用已经存在的 train label mean/SD 定义无条件对角 Gaussian，不加入 context、不优化网络，得到 raw-label NLL：train 14.6981，validation **15.1162**（validation 逐样本均值标准误约 0.2368）。所有被保留随机初始化模型的 validation NLL 为 15.2620–17.6570；因此“最优 validation checkpoint”在本次只表示该 run 内的最低值，并不表示优于这个简单尺度参照。这里未据此进行调参、重选算法或新训练；也未做统计显著性或条件信息下界声明。

这组标签本身接近 unit-scale Gaussian，说明随机初始化和简单尺度模型已提供非平凡基准。它不能证明 context 没有信息。当前无法把参数过多、学习率/训练目标、conditional density 方差过度收缩、特征设计、有限样本与信息弱等候选机制定量分解。最终训练权重未单独保存，本审计也不能独立重算 update1000 的完整训练集 NLL、预测 covariance 或误差归因。

## 当前可支持的结论

本次冻结神经 pilot **未取得验证集学习改进**，其保存产物全部是初始化权重。观察到的泛化失败与过拟合相容，基础 batch/scaling/frame/checkpoint 路径没有发现已证实的实现错误。信息不足、统计不可学习性或根本算法失败均尚未被识别。后续 MH/IS 即使使用这些正规化密度，也必须清楚记录它们是 update0 选中模型，并由真实 full likelihood、完整费用和有效样本诊断判断效果。

[逐 run 指标与分布 JSON](ACTUAL_TRAINING_REVIEW.json)。本审计不修改源码、checkpoint、冻结配置、训练报告或其他研究结论。
