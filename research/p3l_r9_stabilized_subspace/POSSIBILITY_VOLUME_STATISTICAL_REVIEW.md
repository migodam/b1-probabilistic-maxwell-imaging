# R9 fresh32可能材料集合：独立统计审查

**主要结论：固定weak-frame方法确实以有限预算生成了比随机baseline更丰富、更宽的数据相容材料集合，并在本轮比dynamic略省wall；但主epsilon=.075下覆盖仅6/32，较共享anchor新增3个parent，尚不能支持“大概率含GT”的集合主张。固定与dynamic/random的覆盖优势都没有得到严格为正的主区间，近似ensemble动作也没有改善anchor的平均材料MSE。** 这是材料集合覆盖–宽度–费用的结论，不以posterior准确率缺失否决本实验，也不把有限cloud误叫完整posterior。

审查只使用封存后root明确授权的现存结果与必要GT数组，不生成候选、不做Maxwell计算、不重跑bootstrap、不改变门槛/权重/成员/模型/config。原始评分与旧R9负结果保持不变。

## 1. 独立性、资格和费用审计

公开元数据审计见 `audit/POSSIBILITY_INDEPENDENT_AUDIT.md/.json`：32个预声明parent×3arm×16非anchor方向=1536记录完整；8项冻结source/operator/config身份一致；同窗口轮换和最低合格dataPhi选择无违例；33个job全部关闭。封存包含parent及root/jobs的visited状态，私有truth不进入公开seal。首次日志前的duplicate-role参数修复有单独receipt，未发生任何physics，不涉及看到GT后的修改。

获root授权后，本审查独立读取保存的selected材料和必要truth，直接用norm/24重算半径、unique、全部K前缀nearest、epsilon命中、uniform与exp(-Phi)动作MSE；与canonical结果逐项一致，unique计数差0、最大nearest/action MSE差0。原在线输入仅有observation/sigma/parent，公开与私有保存的observation完全一致；saved noise重构误差最大1.43e-17。未重算forward来制造新资格。

元数据/source隔离证据未发现在线GT泄漏；它不是操作系统文件访问的形式化审计。生成器按原Gaussianized66 prior创造GT是必要数据生成，和在线优化读取GT不同。所有源代码/统计冻结都早于首次generation，all32 seal早于root评分批准与TRUTH_ACCESS_BEGIN。

## 2. anchor与完整cloud的主表

K16是16个非anchor方向预算，qualified shared anchor额外纳入，最多17个unique材料。此次32个anchor全部通过绝对noise/legal门，主cloud均非空。去重不按GT；nearest只给离线coverage capacity，部署时不能从GT选中该成员。

|方法|epsilon .03|.05|**.075 主**|.10|.15|平均离线best MSE|平均非anchor unique|平均pair diameter|平均完整方法wall/s|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|共享anchor|0/32|1/32|3/32|13/32|29/32|0.01277981|0|0|anchor优化2.038，另加runtime/setup|
|FIXED_Q8_RAY|0/32|1/32|6/32|14/32|31/32|0.01099138|11.0625|0.181499|5.9676|
|DYNAMIC_Q8_RAY|0/32|1/32|5/32|14/32|31/32|0.01084891|11.3125|0.190746|6.1364|
|RANDOM66_FULL|0/32|1/32|5/32|13/32|29/32|0.01187843|3.2188|0.038341|6.0490|

epsilon球的中心都是实际合格材料，主覆盖是这些完整材料RMS球的并集，不是逐像素min/max框，也没有估计高维Lebesgue体积。不同epsilon给不同材料误差尺度，不能看到.15覆盖31/32后把主epsilon从.075改成.15。即便.15，31/32的单侧95%exact lower≈.860，仍没有证明至少90%parent覆盖。

主epsilon=.075：fixed为18.75%，Wilson95%区间[8.89%,35.31%]；dynamic/random各15.625%，Wilson区间[6.86%,31.75%]。fixed的单侧95%exact lower仅8.50%，远小于预声明90%目标；因此“大概率含GT”不成立。相对anchor的新增命中为fixed3、dynamic2、random2个parent，这些是实质但很有限的coverage增量。

## 3. 四个预注册主contrast

|对比（fixed减reference）|点估计|98.75% parent-paired percentile区间|可作何解释|
|---|---:|---|---|
|覆盖率 fixed−dynamic|0.03125|[0,0.125]|差异仅1个parent；下限触0，不支持严格正优势|
|覆盖率 fixed−random|0.03125|[-0.09375,0.1875]|区间跨0，无法据此确认覆盖优势|
|离线nearest MSE fixed−dynamic|0.000142469141|[-4.40224412e-05,0.000494018438]|正值表示fixed平均误差略大；区间跨0|
|完整方法wall fixed−dynamic /s|-0.168759612|[-0.254402866,-0.0928114794]|本轮同窗口配对区间全负，支持该实现窗口内小幅wall节省|

这些是4项预声明contrast的Bonferroni对应98.75%区间。20000次bootstrap重采样单位始终是32个parent，方向/antithetic pair/seed/source没有增加n。二元覆盖极离散：fixed与dynamic在5个parent同时命中，fixed独有1个，双方均失26个；fixed与random共同命中4、fixed独有2、random独有1、共同失25。仅1个净命中的差异不能解释成可靠的总体提升。Percentile bootstrap在n=32下仍是近似，不能把“Bonferroni”字样误当有限样本精确保证。

secondary qualification-yield fixed−dynamic为−.015625，95%区间[−.0371094,.00390625]，高于预声明−.10 margin。它支持“这组有限方向budget里固定法没有明显损失大量合格替代数”的诊断；不是小样本nominal非劣定理，也不是GT覆盖非劣证明。

## 4. K、半径和宽度说明收益来自哪里

|方法|K0|K1|K4|K8|K16|
|---|---:|---:|---:|---:|---:|
|FIXED_Q8_RAY|3/32|3/32|5/32|5/32|6/32|
|DYNAMIC_Q8_RAY|3/32|3/32|4/32|4/32|5/32|
|RANDOM66_FULL|3/32|4/32|4/32|4/32|5/32|

K不是oracle重排后的成功候选数：它是冻结方向slot前缀，失败照占预算。K4已包括第一pair的.02/.05；K8再加.10/.15。fixed在K4为5/32、K8仍5/32、K16到6/32，新增覆盖主要来自第二pair的补充方向，不能把它说成每增加方向都稳定增加质量。

|只用该半径的4个方向，加同一anchor|fixed 命中/32|dynamic|random|
|---|---:|---:|---:|
|0.02|3|3|5|
|0.05|6|5|3|
|0.1|3|3|3|
|0.15|3|3|3|

这些radius表是预声明secondary分层，不能据其结果事后删除半径并把新算法当fresh确认。值得注意的是：本组.05 subcloud已给fixed全部6个主epsilon命中；.10/.15分层均只有anchor的3个命中。独立GT-nearest成员统计中，fixed在.02/.05/.10/.15/anchor上分别为14/15/2/0/1个parent；dynamic为13/15/3/0/1。大半径目前主要扩大集合宽度，没有提供更多epsilon=.075命中。

fixed/dynamic平均pair diameter约.1815/.1907，random仅.03834；epsilon=.075时相应tube外径均值约.3315/.3407/.1883。weak-frame方法生成更宽、更多的数据相容材料确有实证，但宽度增大本身不能代替GT覆盖。相反，random覆盖5/32虽只略少，通常只在小半径取得资格：512方向中资格数fixed354、dynamic362、random103；random .10/.15没有合格非anchor终点。

random baseline是native Gaussian66方向加full66 signed corrector，trial/J上限与成本记录匹配，但它不是维数/方向分布完全相同的随机rank8子空间对照。因此不能把所有差异纯归因于weak basis的发现，更不能把random解释成毫无校正的弱baseline。

raw/visited variableK诊断保留了更多已花费用的中间状态：fixed2386 visited、1472 noise-qualified unique（含每parent anchor汇总），dynamic2401/1502；主epsilon命中仍分别6/32、5/32，random也仍5/32。保存更多中间点未改善本轮主覆盖，这支持把后续诊断优先放在方向与有限半径覆盖，而不是先把相同路径的中间节点当成新样本。

## 5. 可部署动作和材料任务

|方法|uniform ensemble MSE|likelihood-weighted ensemble MSE|anchor MSE|
|---|---:|---:|---:|
|FIXED_Q8_RAY|0.01316807|0.01305884|0.01277981|
|DYNAMIC_Q8_RAY|0.01298733|0.01288901|0.01277981|
|RANDOM66_FULL|0.01271029|0.01277991|0.01277981|

这些是每parent可部署动作的全材料损失后求均值，和offline GT-best-of-K完全不同。fixed的两个动作都比anchor点估计差；dynamic也略差。random uniform点估计略优，但这是secondary结果，不能从三个动作中用GT挑一个后宣称在线成功。已冻结secondary paired区间显示fixed的uniform/weighted动作相对dynamic均为正差，解释应保留这个负面证据。

exp(-Phi)权重是有限cloud上的近似ensemble分数，没有prior/proposal/Jacobian修正，不是posterior。独立审查发现random的anchor权重均值≈.9915、中位≈.99999968，解释了其likelihood-weighted动作几乎完全退回anchor；fixed/dynamic anchor平均权重约.156/.154。允许approximate权重是合理研究目标，但这里的实际简单权重没有把扩大的集合转化为更好的平均成像动作。

|材料任务的集合表示/范围覆盖|fixed /32|dynamic /32|random /32|anchor基线|
|---|---:|---:|---:|---|
|ROI均值落入候选任务min/max|22|22|25|anchor是单值，不能当范围覆盖|
|硬事件真值在至少一个候选中出现|31|31|32|anchor本身31/32正确|
|soft occupancy ROI均值范围|11|11|3|仅scalar/marginal诊断|

硬事件高覆盖不能替代材料tube覆盖：本组GT事件true为27/32，anchor分类已31/32，weak cloud事件31/32没有增加表示能力。scalar均值/soft任务范围也不等于完整材料image覆盖；逐像素范围不能拼成一个不存在的候选材料。

## 6. 原噪声、prior与冻结过滤审查

GT由原未截断Gaussianized66 prior生成，两个属性经原CDF回到Uniform，原quintic/passive decoder不变；真实观测为保存的原forward加sigma_c=.00119211798514519圆对称复噪声。独立检查的GT噪声Phi均值256.89、范围209.60–289.79，28/32在原95%noise-ball内。95%是抽样模型下的概率，不能要求32个实现值必有约30个，更不能看到4例域外后改噪声或门槛。所有32仍留在主覆盖分母。

主epsilon的6个fixed命中全部位于28个noise-eligible GT中；作为**事后条件诊断**为6/28≈21.4%，仍远低于90%。这不替代预声明6/32。四例GT本身在95%noise-ball之外，说明绝对noise-qualified材料集合也不能保证逐例包含真值；epsilon球可以邻近GT，所以该现象并不自动使这些parent的epsilon覆盖不可能。

旧相对门会产生更严重的设计问题：32/32 GT都不满足DeltaPhi≤2，32/32也都不满足DeltaPrior≤10；实际中位GT−anchor DeltaPhi≈13.73、DeltaPrior≈31.42。anchor优化可比生成GT更贴噪声，Gaussian高维典型抽样的point prior penalty也远大于低范数anchor。这里DeltaPrior计数针对生成GT的latent表示；decoder可能多对一，它不证明同一GT材料的全部原像都被该prior门排除。DeltaPhi依赖材料forward，故其排除对该GT材料本身成立。故本轮采用绝对95%noise-ball与原decoder合法性、把旧相对门留作secondary是必要且正确的合同澄清。这不证明新集合覆盖充分，但避免了由旧门按设计排斥真实prior抽样。

没有GT引导的在线成员/权重/模型选择，所有后述orientation/radius/legacy数值均是seal后诊断。所谓近似probabilistic material volume应理解为有限合法材料集合及经验覆盖–宽度–费用，不应把point prior penalty、candidate计数或exp(-Phi)权重当posterior mass。

## 7. 完整费用与可保留的效率结论

当前campaign实耗293088 RHS256、8380 LU256，新增完整576阶RHS/LU为0；这不等于物理工作免费。generation=256 RHS256，全部anchor=53472 RHS256，三法proposal分别78832/78600/81928 RHS256，精确相加为293088。source/ledger各项reserve、actual、completed一致，失败action0；未达到资格的方向与全部拒绝trial也已计费。wall累计449.965s、CPU396.337s、physics span177.179s、save span154.731s。span可能有嵌套/压缩CPU重叠，不能把它们当互斥百分比相加。

单方法部署视图各自包含同一anchor与frame/runtime分配，因此standalone RHS合计fixed132304、dynamic132072、random135400不能再直接相加当research总账。历史一次operator setup264完整RHS按每parent8.25分配；本轮没有重新支付该setup，整个历史setup研究2184完整RHS仍单独标注。

fixed平均完整方法wall5.9676s，dynamic6.1364s，节省约2.75%；prototype同窗口paired wall差的98.75%区间全负。但fixed实际RHS比dynamic多232（合计约+.18%），不能说固定方法减少了物理调用；其观测wall节省来自这次实现路径/geometry/I/O等组合。原理级、硬件无关或生产端到端加速未证明。完整wall包括证据I/O，当前几GB状态保存也影响比值。

## 8. 不加新实验，值得做的事后诊断

1. **全32 parent的覆盖–宽度–成本atlas。** 直接画已有WIDTH_COVERAGE_COST_POINTS及K/radius表，保留失败与anchor基线，避免只挑6个命中例。把root的关键问题放在“多大材料tube、多少费用，覆盖了多少新parent”。
2. **方向对齐/半径误差分解。** 对已封存候选有标准恒等式error²=anchor_error²+r²−2·anchor_error·r·cos(theta)。本审查事后maximum候选–GT对齐cos均值约fixed.385、dynamic.408、random.324，与.02/.05更常成为GT-nearest的事实一致。可用该分解解释大半径为何主要增加宽度，不需forward；但不能用GT方向生成新候选后继续称fresh确认。
3. **dataPhi排序与GT-distance/动作偏差。** 仅对原封存成员画离线rank/scatter，解释为什么低dataPhi/简单加权没有改善anchor动作；保留原uniform和likelihood权重结果，不在这32上训练/选择新权重。
4. **已付费raw路径的边际贡献。** 现有variableK结果已经显示更多visited unique没有增加主命中，可展示其费用与coverage曲线，区分“增加同一路径节点”与“增加独立方向”的潜在价值。

以上都是事后解释，不能新增到四个confirmatory contrast或产生新的“胜出”门。若据此修改radii、候选数、proposal或权重，需要另行冻结后再评价；本32已暴露GT，不能再次充当未暴露确认集。

## 最终判断边界

本轮支持“固定局部weak-frame可生成多个数据相容、相距不小的材料候选，候选数量与几何宽度明显超过这组随机控制，并以略低wall接近dynamic的集合能力”。它尚不支持“epsilon=.075下大概率含GT”、可靠的覆盖优势、或简单ensemble成像优于anchor。需要改进的核心是有限集合方向/任务信息如何转化为覆盖与部署决策的价值；不是要求先证明完整posterior，也不是用posterior未合格来抹掉已有集合生成证据。

证据入口：`results/possibility_volume/offline_statistics/SUMMARY.json`及各parent scores/cost；本审查可复核补充为`audit/POSSIBILITY_SAVED_DATA_STATISTICAL_AUDIT.json`，元数据封存审计为`audit/POSSIBILITY_INDEPENDENT_AUDIT.md/.json`。新Maxwell/forward调用0，冻结文件改动0，bootstrap重跑0。
