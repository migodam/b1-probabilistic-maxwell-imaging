# Fusion强对照：廉价ordinary矩提供器

本次只改接口和元数据校验，不生成实际特征、不训练、不做物理；旧Spline与Fusion已训练检查点不改。风险改善只能支持具体风险特征，不自动说明Flow概率质量或后验优势。新增对照的目的，是区分Full66 Flow machinery和更便宜的全局y矩特征是否带来独立增量。

同一个tools/r4_fusion_moments.py现在支持显式provider：

|provider|固定K|来源|四通道|
|---|---:|---|---|
|spline|32|既有rawy full66 Spline|原decoder pushforward的经验meanRe/meanIm/VarRe/VarIm|
|diag_gaussian|32|models/density/extended_diag_gaussian_s2026102501/BEST.pt|相同32个固定seed规则及原decoder经验矩|
|field_moments|0|models/field/extended_mlp_small_moments_s2026102501/BEST.pt|确定性学习meanRe/meanIm与joint ReIm点矩的PSD协方差边际VarRe/VarIm|

Field不是完整conditional density，没有latent sampler；K0如实声明确定性，绝不造x_samples或伪称32个世界。其原raw variance与PSD correction warning另存。使用joint2×2矩的边际不等于声明Re/Im独立，四通道也不提供跨像素协方差。所有提供器都是冻结model特征，不能叫true posterior mean/variance。

生成命令保持原SplineCLI兼容，新增`--provider diag_gaussian`或`--provider field_moments`；checkpoint仍显式指定。例：

`PYTHONPATH=src <python> tools/r4_fusion_moments.py --provider diag_gaussian --checkpoint MODELS/density/extended_diag_gaussian_s2026102501/BEST.pt --observations ORIGINAL_PACKED_Y --parents FIXED_PARENTS.json --out NEW_DIAG_MOMENTS --device cpu`

Field同入口换provider/checkpoint。提供器只读原y/parent，不读取True x/material/s标签、anchor或J；保存provider、K、source identity及sealed每parent文件。密度提供器只接收rawy512/full66匹配模型，Field只接收真实mlp_small+second_moments检查点；错类型先拒绝。Field不使用sample_seed，元数据为deterministic，其他两者沿用同一parent seed函数、sample-seed20261025121及K32。

risk_fusion.read_moments、matching_moment_gates及FusionPredictor新增provider/K校验。TRAIN/validation必须同provider、K、checkpoint和sample-seed并truth_read=false；错误provider/K拒绝，Field文件出现x_samples拒绝。新训练payload记录moment_provider、moment_K、contract_version2；未来predictor调用必须显式传对应provider/K，评价loader=fusion自动传入并冻结这些信息。旧payload/旧gate未含provider时仅按既有Spline/K32兼容，保持原数值路径；旧检查点不重写，FORMAT名称保留是架构兼容标识而非将Field冒称Spline。

训练继续用完全相同FusionRisk6418参数、1024TRAIN/256validation父样本、原8在线通道+4矩通道、lr5e−4、batch64、weight_decay1e−3、clip5、max4000更新、同3seeds。新两个provider对照强制`--refit-feature-normalization`且contexts必须1024/256：原8通道TRAIN refit过程及数据与strongbaseline相同；新增4通道各自依同一TRAIN-only均值/尺度程序归一化，validation不参与。材料GT只用于材料MSE监督。每provider的moments控制3seeds，共父调度器计划6fits；没有自动启动或扩大到更多controls。若development无增益停止，不按locked结果调参。

训练入口不需要额外provider参数，直接从sealed特征gate取得显式provider/K：

`PYTHONPATH=src <python> -m p3r4.train_risk_fusion ... --moments-train NEW_PROVIDER_TRAIN --moments-validation NEW_PROVIDER_VALIDATION --control moments --refit-feature-normalization --seed 2026102501 --updates 4000 --out NEW_RUN`

PROVIDER_COMPATIBILITY_CHECK.json：3个provider精确加载重放、6418参数、legacySpline无新字段兼容、错provider/K、TRAIN/validation不匹配、Field伪K32及伪latent samples拒绝、Field mock从未调用sample，均PASS。测试无实际已训练模型featuregeneration、fit或Maxwell。两指定checkpoint的真实header已核对：Diag kind=diag_gaussian/dim66/context512；Field architecture=mlp_small/second_moments=true。
