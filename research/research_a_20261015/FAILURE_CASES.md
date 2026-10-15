# 失败与未建立的主张

本轮的负结果和失败费用与正面材料一起保留。没有改写历史Gate、重置预算或重跑旧campaign来填补缺口。

## 没有从低残差获得更准的单张图

静态方法离线Best4=0.112949，说明候选集中包含更接近真值的材料；在线Selected=0.152732却未胜base=0.151770。排序只用full residual，未使用真值。问题是材料歧义与选择信息不足，不是可以用Best4偷偷替代在线误差。

## 强联合算法条件未通过

Best4相对multi-start改善17.1558%、paired区间排除0，但coverage只增加6.25个百分点，平均wall约降10.86%，full calls相同。冻结强联合条件要求的coverage/同质量成本增益未满足。不能因为主误差指标较好，就将整体写成GO或TAP级新算法。

## adaptive未表现出历史优势

开发版前8步仍处于warm-up，forecast未启动；第二版相同长度也没有触发新增forecast逻辑。第三版改为2步后启用，并在确认前冻结。确认结果Best4比静态差0.7436%，平均wall更高、coverage更低。它是完整保留的负结果，不以“学习策略还需继续调参”无限延长。

独立同protected cell重启能在8个development端点复现有限编辑结果，最大材料距离2.84e-5<0.03。当前有限编辑的正面增益不能归给新反馈或历史依赖。该小对照不足以证明所有搜索历史都无用。

## MDN没有产生本轮实际收益

旧C1 full-distill proposal的Best4=0.137451，multi-start=0.136339，未建立改善。没有新增训练；旧teacher成本应披露。没有因为此负结果自动训练更大网络，也没有把模型输出当full-Maxwell后验。

## 可补偿不等于材料空间宽裕

第一个预定witness使用ε≈2.2+0.16i，分别接近实部下界和损耗上界，gain≈1.276也靠近上限1.3。它的full相容性真实，但支持边界压力必须显示；不能把这样的profile-fit解释为大prior邻域或高posterior质量。

8例no-comp控制中3例FAIL→PASS全部来自无protected geometry box的global-prior步骤。保护编辑的4个小步骤本来就相容。该实验未单独证成“保持结构的非局部编辑因compensation而发现新分支”。

## E2五个条件没有合法解

|条件|base白化残差|finite archive数|结构代表数|Best4复材料MAE|
|---|---:|---:|---:|---:|
|p250065 weak|24.1987|0|0|NOT_DEFINED|
|p250065 ring|27.0654|0|0|NOT_DEFINED|
|p250069 weak|20.9292|2|1|0.0435144|
|p250069 ring|26.3769|0|0|NOT_DEFINED|
|p250073 weak|20.8919|1|1|0.0347768|
|p250073 ring|25.1159|0|0|NOT_DEFINED|
|p250077 weak|22.9950|1|1|0.0484529|
|p250077 ring|26.3348|0|0|NOT_DEFINED|

阈值23.785449固定；无gain补偿、无重新加噪。所有ring条件以及一个weak条件失败。不能剔除5个无解条件后把三个eligible均值0.0422480称全组MAE；四个parent也没有一个完成weak/ring双条件有效endpoint。

p250069 weak中唯一额外成员仅1格结构mask改变，低于6格阈值；同一材料已在no-comp原编辑中存在，MAE比base更差3.36%。因此“材料RMS去重多一个”不能记为新结构解，也不归因于补偿。

修正器不是空转：200个QP中191个报告收敛，finite残差配对下降中位数11.03%，保护单元最大缺陷0；但32条finite轨迹invalid→valid为0。缺少可行粗解、16/32实参数chart、有限修正步数和物理不可补偿尚未完全分离。不能把本轮失败写成完整TM逆问题唯一性证明。

## 首次真实audit失败被保留

E2首次dual检验对极小二阶量归一化，产生1.91668e-8的相对误差，超过冻结1e-8容差。该attempt的1state/91RHS完整保留。随后改为可分辨的1%扰动，没有放宽容差，成功audit的dual误差3.21550e-13。

为保持额外严格RHS预算，省略最初两个条件各一个原拟no-comp端点；实际no-comp30state而非32。原协议不改，偏差写在 `results/E2/AUDIT_REVISION.json`，旧失败存于 `results/E2/attempt_1_failed_audit/`。

E2累计9214/9216复LU RHS，仅余2，低于下一状态8源最小费用；即使尚余56state也停止。不是降低计数或另起run以继续。

## 缓存不等于提速

固定4对象、5次轮换回放完全重跑proxy搜索，只复用原合法endpoint；新物理0。warm wall中位数都高于原对应cold首次搜索，没有建立solver、ROM或总流程加速。旧bank/模型加载和共享coarse费用不能隐去。

## 接口限制与NOT_RUN

原A schema强制asset摘要；执行overlay如实用null，原文件保持。Mie没有暴露内部本构缺陷，null+明确原因，不能用边界残差或0填补。原jsonschema依赖未安装，使用明确范围的dependency-free schema/semantic验证，不假称通用schema引擎。

TSOM/SOM、ROM/反馈专用corrector、学习式next-edit训练和新非同心3D均NOT_RUN。没有新后验校准、全支持覆盖、全局模式发现或一般3D声明。Research B原生schema未读取，其兼容性不作假设；本轮交接采用显式版本common exchange。

## 正式交接元数据的decoder修复

主线程复核发现：optimizer请求的中心半径约2.9255e-11，经过原prior坐标变换再decode后会得到精确0。初版导出按原请求`p_r>0`数空腔，导致例如confirm_00/world_05把实际一个ring空腔误报为两个空腔；材料数组和full预测本身一直来自正确decoder，未受影响。

修复后结构半径、空腔数、连通材料区间和feature scale统一从canonical decoder和已保存径向复ε数组推导；原请求坐标继续保存在history/provenance，不被静默改写。所有193个材料的真空连续区间数与导出空腔数一致。32份包在独立`export-actual-packages-semantic-repair` job重建，新物理0；旧导出receipt和原metadata保留，未改ledger、材料数组、预测或评分。具体影响条目与回执见 `results/ACTUAL_PACKAGE_DECODER_REPAIR.json`。


## 后处理错误与评分进程边界

最终汇总第一次把stage字段重复传给dict导致TypeError；该0物理失败job保留，独立repair job修复，没有重算物理。E2原执行在全部搜索冻结后于同一进程打开离线真值；未发现反馈至在线搜索，但不符合严格进程隔离。交付代码改为独立e2_offline进程，并实际重放0物理评分，32行指标和684数组完全一致。旧结果和执行偏差在legacy_offline_inprocess保留；不假称原执行本来就已隔离。
