# 冻结公共2×2对照：墙改变与强目标插入的固定chart失配

**COMPLETE。实耗12 FE /432 RHS；重FE队列已释放给ROOT。** 本试验按事前冻结的四个公共状态运行，rank16 chart只在“均匀墙＋无目标”构建。没有重训、刷新chart、改prior、改rank或新增优化；不是blind成像准确性比较。

## 四状态读数

全部误差采用相同σ=.01白化、公开参数尺度，表中为Frobenius相对误差。所有full/fixed J在同一状态比较，并共享完整物理响应。数据Schur使用公开墙prior、目标精度0、无算法damping。

| 公共墙/目标状态 | whole J误差 | wall J误差 | target J误差 | 数据目标Schur误差 | 阻尼目标步相对差 |
|---|---:|---:|---:|---:|---:|
| 均匀墙、无目标（reference） |≈0|≈0|≈0|≈0|≈0|
| 改变墙、无目标 |.338382|.343166|.062593|.006296|.584740|
| 均匀墙、强目标 |1.421331|.242723|1.979217|4.396285|.519146|
| 改变墙、强目标 |1.203040|.386262|1.918435|4.277271|.558501|

“阻尼目标步”仅为以同一已测dev y形成的**无约束局部QP诊断**：所有19个scaled坐标加1e-3I；不是实际有界优化器接受步、posterior precision或真实目标误差。无目标状态是前向reference/control，目标contrast=0，不是正contrast重建prior内的合法估计。

## 能支持的归因

1. **不能把该固定chart的全部失配归因于未知墙。** 即使墙不变，只插入预声明强目标，whole J误差已达1.421、目标J误差1.979、目标Schur误差4.396。相比之下，仅改墙、目标仍不存在时，这三个数为.338、.0626、.00630。这个受控例中，目标插入相关的状态/场/载荷变化足以破坏无目标参考chart。
2. 这里的“目标插入相关”包含目标自身改变完整L、局部场和合法材料微分载荷的共同后果。**本试验没有进一步分离其中每一个机制**，也没有证明任意目标/任意墙下总是由目标自作用主导。
3. 墙、目标对chart误差的作用不是独立可加。whole J相对误差的四个匹配差值：墙改变在无目标时为+.338382，在强目标时为−.218292；目标插入在均匀墙时为+1.421331，在改变墙时为+.864657。因此差中差为 **−.556674**。
4. 这个负号是**误差指标上的交互**，不能叫负的物理散射能量或普遍的保护作用。各状态的相对误差分母会改变；因此同时保存了绝对误差及矩阵差分，避免只靠相对数叙事。

## 绝对量与矩阵交互，避免分母误导

whole J白化scaled误差绝对范数依次为约0、21.768193、110.300331、102.224481。故在强目标状态改墙，其绝对误差范数也减少8.075850；这次相对误差下降并不只是分母变化。绝对误差范数差中差为−29.844043。

令E_ij=Jfixed_ij−Jfull_ij，i为墙、j为目标。真正的误差矩阵交互
`E_11−E_10−E_01+E_00` 的Frobenius范数为 **48.387783**。这是矩阵差的范数，不等于上面误差范数的差。其全部原始数组已保存；不能把−.556674和48.387783混成一个量。

本试验形成的是同状态Jacobian与Schur归因，没有形成完整Green矩阵；因此不是对所有载荷的Green算子分离证明，也不是有限成像可辨识性定理。它提供的是足以修正“都是墙导致固定chart失败”这一解释的预声明有限对照。

## 文件

- `SUMMARY.json`、四个状态NPZ/JSON：完整/固定complex J、相同full response、公开尺度、白化scaled J、data/step Schur、H/g/步、F27界与F29残差后果。
- `REFERENCE_f*.npz`：公开reference U、谱、质量对角及构建audit。
- `FIGURES/FOUR_STATE_TABLE.csv`：表中未四舍五入数值和绝对误差。
- `FIGURES/RAW_MATCHED_DIFFERENCES.npz`：每状态J/response/Schur及wall、target、interaction的全部矩阵差。
- `FIGURES/MATCHED_EFFECTS.json`：相对/绝对误差范数差与矩阵差范数，命名分开。
- `FIGURES/public_2x2_mechanism.png`、同名PDF：三种误差2×2图与匹配效应图；PNG目视检查，修复面板间距/标签后可读无裁切。

## 实际费用与失败保留

成功physics wrapper：`results/jobs/chart_mechanism_public_rootrepair/RECEIPT.json`，退出0；**19.386038 CPU秒、19.972133墙钟秒、轮询峰值RSS1,601,978,368字节**。物理本体receipt为12 FE /432 RHS，未触及18/600上限；计费以预动作ledger与wrapper为准，不能把本体较小CPU与wrapper相加重复计费。

第一次launch在任何FE前因frozen代码的BudgetBook默认root找不到RUN_CLOCK而失败，0新增FE/RHS；wrapper `chart_mechanism_public` 保留，CPU .531051秒。只修复诊断脚本传入`BudgetBook(root=ROOT)`，未改frozen科学代码。成功run原receipt的error字段被后续state-error字典变量覆盖，已将字段修复为null，并保留`RECEIPT_ORIGINAL_SHADOWED_ERROR.json`；没有物理失败或重跑。

两次离线绘图wrapper共2.550875 CPU秒，0 FE/RHS。首次与成功physics尝试合计19.917089 CPU秒；连同本批两次离线图约22.467964 CPU秒。输入guard日志拒绝项0；本脚本无private输入路径，未读blindtruth或scores。完整response只用于该预声明诊断，未改变冻结成像方法。
