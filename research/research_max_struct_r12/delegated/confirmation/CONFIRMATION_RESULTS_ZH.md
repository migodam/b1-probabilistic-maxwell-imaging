# R12 独立材料：封存后离线评分

本报告只给已冻结实验的数值、图像与范围；科学结论由主线程审查。确认24与开发6分开，主规则与新增次规则分开。所有材料均保留；unavailable不删对象、不补主规则结果。0次新正演、0次新优化/反演。

确认主规则在24生成前冻结；次规则在24数据生成后、首个确认逆解/GT评分前登记。次规则是新的候选选择估计器，不是纯DR机制增益。开发6未实施次规则，全部明确标NA。

指标：RMSE=||m−GT||₂/24；amplitude NRMSE=RMSE/(maxGT−.1)；relativeL2=||m−GT||₂/||GT||₂；IoU使用固定χ>1超水平集，连续材料中它不代表真实物理support。旧描述成功为继承合法box(1e−12数值容差)、NRMSE≤.15且IoU≥.7，不能替代连续指标。高/低相误差只在结构化二相生成区计算，连续随机族不强行分相。

配对bootstrap：10000次、seed20261112，以完整独立材料ID作cluster、无分层重抽；所有方法使用同一重抽索引。次规则NA留在全对象重抽册中；连续差异仅在该次重抽的双侧可用对象上计算，必须与paired_n/availability一起阅读；没有把缺失误差填0。95%为percentile区间，未作多重比较校正。

路径修正：compact manifest中stage字段含完整Windows相对路径；评分忽略该字段，以compact_primary先将反斜杠换成斜杠再取第一段。原manifest及冻结源未改。

## confirm24

|规则|方法|可用/全体|旧描述成功/全体|均值NRMSE(可用)|中位NRMSE|均值relativeL2|球内但材料差|球外|均值整child CPU秒|
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
|primary|CIE|24/24|9/24|0.202932|0.124013|0.583839|9|6|136.609|
|primary|Dual|24/24|9/24|0.185456|0.117101|0.530731|10|5|136.185|
|primary|Global|24/24|9/24|0.184007|0.117186|0.526496|10|5|135.710|
|primary|A+|24/24|9/24|0.171692|0.118582|0.492981|10|5|152.024|
|primary|G19|24/24|9/24|0.143891|0.119294|0.415800|11|4|178.758|
|primary|DR|24/24|7/24|0.184631|0.125773|0.550124|11|6|128.448|
|primary|Multi|24/24|4/24|0.619361|0.594493|1.878789|7|13|219.280|
|primary|NN|24/24|9/24|0.149514|0.115584|0.425488|11|4|85.811|
|secondary|CIE|18/24|9/24|0.091957|0.106627|0.293033|9|0|136.609|
|secondary|Dual|19/24|9/24|0.093472|0.107925|0.297084|10|0|136.185|
|secondary|Global|19/24|9/24|0.093689|0.108035|0.298298|10|0|135.710|
|secondary|A+|19/24|9/24|0.093897|0.108609|0.298587|10|0|152.024|
|secondary|G19|20/24|9/24|0.096584|0.111002|0.305528|11|0|178.758|
|secondary|DR|18/24|7/24|0.098258|0.111214|0.309753|11|0|128.448|
|secondary|Multi|11/24|4/24|0.501848|0.149844|1.218293|7|0|219.280|
|secondary|NN|20/24|9/24|0.093049|0.107290|0.298462|11|0|85.811|

次规则均值条件于可用子集，不能拿它与主规则全体均值直接声称改进；不同方法的可用集合也可能不同。完整缺失册、逐对象同物理数据判别和共同可用配对差异均在RESULTS中。

## dev6

|规则|方法|可用/全体|旧描述成功/全体|均值NRMSE(可用)|中位NRMSE|均值relativeL2|球内但材料差|球外|均值整child CPU秒|
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
|primary|CIE|6/6|2/6|0.179528|0.106746|0.562475|3|1|138.073|
|primary|Dual|6/6|2/6|0.173895|0.106342|0.545600|3|1|137.755|
|primary|Global|6/6|2/6|0.178465|0.106847|0.559479|3|1|139.880|
|primary|A+|6/6|2/6|0.174175|0.106906|0.546585|3|1|159.227|
|primary|G19|6/6|2/6|0.093441|0.107410|0.303939|4|0|186.898|
|primary|DR|6/6|1/6|0.185526|0.130211|0.630226|3|2|133.919|
|primary|Multi|6/6|0/6|0.662166|0.770476|2.260874|3|3|209.133|
|primary|NN|6/6|2/6|0.136173|0.108957|0.432402|3|1|93.396|
|secondary|CIE|0/6|0/6|NA|NA|NA|0|0|138.073|
|secondary|Dual|0/6|0/6|NA|NA|NA|0|0|137.755|
|secondary|Global|0/6|0/6|NA|NA|NA|0|0|139.880|
|secondary|A+|0/6|0/6|NA|NA|NA|0|0|159.227|
|secondary|G19|0/6|0/6|NA|NA|NA|0|0|186.898|
|secondary|DR|0/6|0/6|NA|NA|NA|0|0|133.919|
|secondary|Multi|0/6|0/6|NA|NA|NA|0|0|209.133|
|secondary|NN|0/6|0/6|NA|NA|NA|0|0|93.396|

次规则均值条件于可用子集，不能拿它与主规则全体均值直接声称改进；不同方法的可用集合也可能不同。完整缺失册、逐对象同物理数据判别和共同可用配对差异均在RESULTS中。

## DR逐对象物理诊断

物理分类区分：CN99球内且旧描述材料标准未通过、球外、球内且旧标准通过、unavailable。球内“材料差”是描述标签，连续NRMSE/relativeL2/IoU原值都保留；球外也不自动等于未成像。

- primary DR：可用24/24；旧描述未通过17个中，球内11、球外6；另有NA0个。球内未过对象：r12_confirm24_0005, r12_confirm24_0006, r12_confirm24_0008, r12_confirm24_0009, r12_confirm24_0011, r12_confirm24_0012, r12_confirm24_0013, r12_confirm24_0014, r12_confirm24_0015, r12_confirm24_0016, r12_confirm24_0017。球外对象：r12_confirm24_0001, r12_confirm24_0003, r12_confirm24_0007, r12_confirm24_0010, r12_confirm24_0019, r12_confirm24_0022。
- secondary DR：可用18/24；旧描述未通过11个中，球内11、球外0；另有NA6个。球内未过对象：r12_confirm24_0005, r12_confirm24_0006, r12_confirm24_0008, r12_confirm24_0009, r12_confirm24_0011, r12_confirm24_0012, r12_confirm24_0013, r12_confirm24_0014, r12_confirm24_0015, r12_confirm24_0016, r12_confirm24_0017。球外对象：无。

DR相对Dual/G19的全部对象差异在DR_vs_Dual_G19_per_object；没有另设事后“实质改善”门槛。下面按绝对NRMSE差排序列出双方前三个最大改善/退化，完整24册仍全部保留。

- primary DR−Dual：最大改善 r12_confirm24_0010 ΔN=-0.187071；r12_confirm24_0003 ΔN=-0.181357；r12_confirm24_0007 ΔN=-0.146175；最大退化 r12_confirm24_0019 ΔN=+0.482390；r12_confirm24_0002 ΔN=+0.005032；r12_confirm24_0004 ΔN=+0.004847。
- primary DR−G19：最大改善 r12_confirm24_0003 ΔN=-0.018400；r12_confirm24_0000 ΔN=-0.002958；r12_confirm24_0021 ΔN=-0.001980；最大退化 r12_confirm24_0019 ΔN=+0.482880；r12_confirm24_0007 ΔN=+0.361999；r12_confirm24_0010 ΔN=+0.083614。
- secondary DR−Dual：最大改善 r12_confirm24_0009 ΔN=-0.023941；r12_confirm24_0011 ΔN=-0.001909；r12_confirm24_0017 ΔN=-0.001335；最大退化 r12_confirm24_0020 ΔN=+0.009490；r12_confirm24_0018 ΔN=+0.007093；r12_confirm24_0002 ΔN=+0.005031。
- secondary DR−G19：最大改善 r12_confirm24_0009 ΔN=-0.021048；r12_confirm24_0000 ΔN=-0.002962；r12_confirm24_0021 ΔN=-0.001980；最大退化 r12_confirm24_0020 ΔN=+0.009579；r12_confirm24_0002 ΔN=+0.004995；r12_confirm24_0018 ΔN=+0.004100。

原始DR报告采用两种非GT选择：proposal阶段cap内已物理审计材料的最小原Phi，以及最后一个cap内已审计C-shadow；不是再跑一次DR，也未重建未审计的next-step shadow。源projection audit材料已从远端原记录抽取，全部audit评分在raw_DR_all_audits；与primary/secondary分开的固定raw选择在raw_DR_fixed_selections。

- raw_proposal_minPhi：24/24可用；旧描述成功7/24；均值NRMSE=0.185814；球内材料差11、球外6。
- raw_last_cap_eligible_audited_shadow：24/24可用；旧描述成功7/24；均值NRMSE=0.193849；球内材料差11、球外6。

## 成本、环境和文件

确认192个有效child外部GetProcessTimes CPU为28147.796875秒；开发48为7189.6875秒。13次零物理启动失败另计.90625秒，失败源与receipt保留。主/次规则共享同一次已付费运行，不应将二者成本相加。576LU与rank19smallLU、FFT等不同维度操作不能合并；外部CPU已经包含这些工作。

所有实验保留240CPU总上限、proposal144、corrector96、单BLAS线程和8个固定P核；不同方法/自然停止的实际消耗不同。30/60/120/240图为已付费前缀，不是四次新实验，也不是可相加的耗时。原minPhi主规则没有被本次评分改变。

Windows五个baseline wrapper源码显式设自身BelowNormal；search继承controller。EXECUTION_PRIORITY_ENVIRONMENT.json只是一张中途actual nice快照（当时baseline BelowNormal、所见controller/search Normal），不能证明整队列priority历史。开始/近结束/空闲结束频率快照同样只是成本环境采样，不证明某个硬件限频原因或严格同吞吐。没有改电源、频率、并发或线程。

NN使用冻结R8 physics权重，含既有监督训练信息；训练成本不属于本次inference CPU。multistart源码先用uniform.1到proposal deadline(或自然停止)，之后Gaussian扰动RMS1/4/10循环、每次最多24CPU/1000iter；最后再同corrector。NN-direct proposal只有feature/inference，随后最多96CPU corrector，未用闲余proposal补新优化。实际算法定义、耗时与费用须共同阅读。

详细表：RESULTS/offline_scoring_v1；科学图：FIGURES/offline_scoring_v1。每个JSON对应CSV；per_object含每个失败与NA，prefix_per_object含四个固定前缀，aggregate/pairs包含全体分母及条件bootstrap说明，per_object_pairs与DR专表给全部逐对象差异。

全24网格：confirm24_primary_ALL_MATERIALS与confirm24_secondary_ALL_MATERIALS；另有四家族逐页大图。灰格为NA、红caption为旧描述标准未过；所有对象都画出。DR_RAW_PRIMARY_SECONDARY_ALL24并列GT、两个raw定义及两个最终选择。

此处不作新颖性、普适优越性或连续Maxwell精度判断；主线程据公式、完整数据、可用性与科学边界审查。

## 追加：保存的GT噪声资格与连续材料误差（无新FW）

确认24的Yclean均已保存；直接核对Y−Yclean与同一CN99半径：24/24个GT本次也在球内，其余保留资格差异。逐对象表见GT_saved_noise_qualification。

DR主规则11个球内旧描述未过对象中，NRMSE>.15为0个、IoU<.7为11个；均值NRMSE=0.117711，均值relativeL2=0.403463，relativeL2范围0.324088–0.600206。这主要是超水平集/细节精度未过旧标准，不能把它们统称完全没有成像。

其中11/11个对象同时满足：GT合法且用已保存Yclean核验在球内，DR主材料也合法、在同球内、与GT不同。这只支持该次有限噪声半径下两个候选都与观测相容的诊断，不是exact nonuniqueness、连续模型不可辨识定理，也不证明所有先验都无法选择GT。

固定raw proposal minPhi→primary的noise状态转移计数：{'True->True': 18, 'False->False': 6}；NRMSE下降12、上升12，最大下降-0.052514、最大上升+0.050786。没有按GT挑某一步rawshadow。

全部DR连续诊断/旧标准分项在DR_failure_continuous_diagnosis；raw/primary/secondary并排数值在DR_raw_primary_secondary_comparison；raw与Dual/G19的全部诊断差异在DR_all_views_vs_Dual_G19，注明rawproposal与完整control的预算不同，不能偷换为新等成本比较。
