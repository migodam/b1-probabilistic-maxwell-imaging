# MAX-STRUCT-R3 独立最终审阅

2026-11-03。结论：**可作为有条件机制分析和负面成像结果的技术报告交付；没有新的盲测成功率、成像优势或全局最优证书。** 本审阅核对实际代码、封存预测、收据、`THEORY_DERIVATIONS.md` 和主报告；只执行保存数组的算术、评分与绘图，没有正演、反演、fresh truth 生成或哈希检查。核心定理和最终科学裁决仍由主线程承担。

## 1. 封存与评分

两条新增 `paired_neighbor_rays`、`dressed_neighbor_rays` 在本次评分前均为 `SEALED_DEVELOPMENT_PREDICTION`。原六臂及九臂评分和封存快照保留，现有评分器追加两臂后统一重算11臂。`FINAL_ARITHMETIC_READBACK.json` 对保存材料、预测观测、原噪声门、同一 prior、NRMSE、继承 Otsu 和收据目标逐项回读，全部可计算项一致。最大目标收据差约1.74e−18。

|控制|NRMSE/12|继承Otsu IoU|目标相对残差|原噪声联合门|
|---|---:|---:|---:|---|
|noisy_strong_release|0.490477|0.294118|0.100661937|FAIL|
|clean_native|0.389673|0.422680|0.107821854（clean）|FAIL|
|passive_path|0.817495|0.204545|0.051017057|FAIL|
|extra_original / globaldual_restart|0.510695|0.245902|0.100734862|FAIL|
|scaled_original_restart|0.510717|0.254098|0.100733687|FAIL|
|onecell_roots|0.510717|0.254098|0.100733687|FAIL|
|neighbor_swaps / all_pair_swaps|0.510717|0.254098|0.100733687|FAIL|
|paired_neighbor_rays|0.510717|0.254098|0.100733687|FAIL|
|dressed_neighbor_rays|0.510717|0.254098|0.100733687|FAIL|

这是一个已公开材料对象的11个方法控制，不能作11个独立样本。Clean控制的噪声门在自身目标下为不适用；上表联合门只表示与原 noisy Y 比较，不能解释为clean数据拥有继承CN噪声分布。旧R2 native参考为N=0.510695、IoU=0.245902、残差0.100734862，没有重跑或改写旧盲表。

**Otsu小敏感性必须保留。** 原合同以256-bin直方图中心为阈值，并用累计量及首个 `argmax` 决定并列区间。代数等价的均值形式可能在浮点并列区间选到后一个bin；本轮强模态释放IoU变化约0.00840，缩放端点也有小变化。没有按GT选阈值，也没有修改旧实现或分数。所有N≥0.38967，远高于0.15，因此该敏感性不改变任一材料/联合FAIL。初次独立复算的差异和处理过程保存在 `FINAL_ARITHMETIC_INITIAL_FAILURE.json`。

已保存gradient的9臂，其投影梯度可由数组独立复算；`passive_path`和`extra_original`没有保存gradient，本次对应检查为 `NOT_RUN_SAVED_GRADIENT_ABSENT`，只保留其原收据记录，不补做伴随。

## 2. 机制主张与公式检查

- **硬锁代价：PASS under stated assumptions。** 满行秩A和可逆H下，最小状态能量为 `b*(AA*)^-1 b`；它是固定材料、固定模式、固定源权重的代价，不是物理corrector最终误差下界。公开GT处观测0.0289064，噪声期望0.0226618，固定模型99%上界0.0424996，数值互相一致。GT下界与原CIE阶段末0.03317位于不同材料锚点，只能比较量级，不能说解释了阶段损失的某一百分比。去噪和释放都失败，支持“不是唯一失败原因”，不能从两次控制判定所有可能释放策略无效。
- **单元圆弧：PASS under stated finite-model assumptions。** Sherman–Morrison完整响应、二次数据导数和被动Gram约束的代数相容。内部严格数据局部最低点的结论限于一条合法无极点圆弧；box边界、全prior和多单元耦合必须分开。正虚部论证针对所列有限disk-pulse Hankel核，不能推广到任意离散化/PML或作为浮点特征值认证。
- **rank-r阶数：DERIVED conditionally。** 所列有理传递形式的分子/分母次数及导数最高项抵消给出至多4r−2次根问题；不可见或不可激发方向可降阶。一个参数不等于rank1。加入smooth/TV prior后候选根不穷尽总目标驻点，数值根后向残差也不是严格根完备性证书。
- **人工损耗fiber：理论条件成立，数值为EMPIRICAL。** 商映射的fiber恒定条件和局部Lipschitz必要条件正确地分开。1.70e−17的原端口差并非精确零，不能据此否定任意不连续变换；加入i后的7.42e−5差和切向9.61e−5支持有限数值模型的稳定性障碍。没有连续PDE不存在定理。
- **一阶与候选失败：有限结论。** 缩放目标达到PGinf=9.915e−9，是声明容差内的一阶证据，不是精确驻点、严格局部极小或全局认证。单元根、全部交换和两条rank2射线字典均没有实质逃逸；这不证明所有二维方向无效，更不证明至少3个单元才可能成功。

保存数组的额外算术回读见 `FINAL_MECHANISM_ARITHMETIC_READBACK.json`：噪声谱的期望/上界、fiber端口/切向差、一单元最低已存候选与两个矢量响应谱均一致；只重复范数、标量运算和两次54×3响应SVD，零物理LU/RHS。

## 3. 新增射线与正控制

最近邻2208条方向中，985条合法区间零宽，实际1223个根问题、2470个候选，无接受。完整传播交互图选567对，其中338非邻接；1134条方向中513零宽，实际621个根问题、1250个候选。后者接受一项5.57e−5的材料L2变化，Φ下降约1.06e−13，原残差下降约2.25e−12，N反而增加约2.6e−9。称“无实质恢复”准确，称“接受成功”或“新的逃逸”会误导。

全部165600交换是576个位置的全部无序对；119881个是恒等交换，45719个需要2×2求解。它是公式规定的有限字典，不是随机初值或参数调优扫描；失败只约束这个单次交换家族。完整逆、Gram构建和候选prior均收费，没有把候选数量当成完整Maxwell求解次数。

公开一单元Hankel控制是有效的机制正例：上界m=20处导数为−2.10312，局部合法方向使目标升高；原L-BFGS-B停留20。完整圆弧候选恢复1.2，full forward残差2.71e−15。但接收S满列秩，数据已确定完整j，成熟CIE读出也恢复1.2。它解释边界停滞及current lifting优势，不是新方法战胜CIE，更不是576单元盲成像。

三维8体素DDA范围控制保留真实横向源及矢量后端。dyad响应差的相对奇异值约 `[1,3.84e−15,2.61e−15]`，isotropic为 `[1,.0110254,6.38e−5]`；足以展示后一族拥有额外可见响应方向。没有任设rank接受门，也没有连续精度、3D成像或普遍高阶多极小结论。

## 4. 数据、源度量与prior

审阅的两条射线程序只读取原公开Y、算子与已封存scaled端点。有限候选按同一原始归一化Φ和完整历史prior选择，最优提案再由原Maxwell计算核验。交互图由C/O/Q与源接收Gram决定，没有GT方向或阈值。单元/交换程序及known-controls沿用同一选择原则。

CIE clean反事实仅读公开Y_clean，未给材料或support；wd/ws仍固定于原noisy Y/E。它改变信息输入，不能宣称同数据算法改善。噪声floor和fixed-material诊断明确读取公开GT，属于离线机制证据；onecell/vector/fiber采用公开存在锚点，不能混作未知材料恢复。没有看到新fresh truth或online GT排序路径。所有这些程序的开发选择本身发生在R2曝光后，所以“封存后评分”不能恢复盲性。

原prior宽度19.5、smooth/TV各1e−4保留；box为[.02,20]。不同辅助路径的源度量、阶段数、早停和实际LU/RHS并非相同计算预算，报告没有资格作等计算成像优势比较。

## 5. 成本与可视审阅

旧 `DEVELOPMENT_COST_AUDIT.json` 是六臂时期快照，遗漏后来的单元/交换/射线以及noise/fiber/onecell/vector诊断，不能作最终总账。新增 `final_cost_consolidation.py` 保留旧文件，输出 `FINAL_DEVELOPMENT_COST_AUDIT.json` 和 `FINAL_DEVELOPMENT_COST_LEDGER.csv`。17个外层物理收据合计 **4092 LU / 101104 complex RHS**：source48922、adjoint48656、receiver-transpose70、full-resolvent3456，分项恰好相加。

各系统维数不等（576、16、6、1、24），上述是工作计数而非FLOPs。current/Gram、小2×2系统、嵌套物理计数和timers按收据保留，不再加到外层总数；历史R2参考不收费。包括各次已保存评分、只读审计和本次绘图后，已知过程wall之和约37.48s、CPU约36.30s；研究/工具/导入/写作与部分序列化未测，金额为null，没有声称全科研成本仅37秒。零物理成本的启动错误和中断算术复算也保留。

`figures/DEVELOPMENT_RECONSTRUCTIONS.png` 已从保存数组绘制并人工视读：六面板统一[.02,20]色阶、nearest显示，文字与目标差异可读，没有裁切。它直接显示损耗路径的材料恶化及ray终点与失败native相近。没有用自动缩放或平滑掩盖差异。

## 6. 主报告交付限定

主报告的主要数值和科学范围与证据相符。交付前须保留以下措辞：将“真正的一阶驻点”表述为“达到声明一阶容差/近似一阶驻点证据”；GT floor和原stage loss注明不同材料锚点；保留Otsu并列敏感性；最终成本链接指向本次完整总账。这些属于范围和可复核性修正，不要求重跑实验。

新鲜生成/新鲜反演/新鲜盲评分均为 **NOT_RUN**。在尚无具体材料收益算法时，不需要通过新样本或放宽合同救结论。可以交付已完成的机制研究、负面结果、明确物理正控制与适用边界；不能交付“新盲成像突破”。

## Parent publication disposition

Parent applied the requested approximate-stationarity wording, explicit different-anchor floor comparison and Otsu tie sensitivity. `audit/FINAL_COST_TOTALS.json` is the independent final ledger plus the separately timed mechanism plot, with unchanged 4092 LU/101104 RHS. Publication/export timers remain separate. Frozen R1/R2 files are unchanged.
