# 协议与实际资产的修订记录

本次直接用户任务优先于压缩包中的建议。原包及历史冻结文件保留，不修改旧 Gate、数据、参考或账本。

| 差异 | 本轮处理 |
|---|---|
| 原包只授权 E0，E1/E2 execution_authorized=false | 用户本次明确授权合计4096 attempts、7200受控wall；独立 AUTHORIZATION_OVERLAY 生效，不修改原配置、不借P2剩余541 |
| 原包首推参考后验/EIG | 优先部署有限候选集合；候选无posterior。EIG/material-risk仍作为独立参考辅助比较，不能把它们选action后的profile结果称全程candidate-only |
| 有相同schema文件名 | A原生1.0、execution profile及CommonExchange与B schema不同；用显式version adapter，不覆盖A |
| A manifest要求digest字段、内部constitutive residual | 按源execution profile保留digest=null与NOT_EXPOSED；未计算任何SHA。边界方程残差不冒充内部本构缺陷0 |
| A只声明new2.4 | 原A菜单保留；B独立已授权六频菜单写入自己的版本化overlay，使用自身已付费资格核对与噪声锚，不修改A |
| A prediction是g_source*F | 独立保留源mean，另导出source_fitted_gains与unscaled_old_fields，禁止二次乘gain |
| A union-grid、physical coordinates与P2系数不同 | 材料数组为准，显式decoder转换，记录tiny-h浮点inactive atom与坐标roundtrip误差；没有强迫13-cell模型 |
| accepted NPZ中的partition_* | 三条refined摘要已更新，但partition_*仍是历史诊断；不用旧字段充当refined evidence |
| P2 full correction有full字段 | 只在1.6有192/384个paidfields，不能称双频full posterior；E0候选字段读取已有37个双频代表，Mode2单列连续proxy |
| 固定K更新与发现新材料 | 固定K的Best-of-K不变；E2仅扩decision catalogue、不扩2048 posterior支持；真实A独立corrector只有full端点通过后才记新增材料 |
| 旧数据曝光 | E0七例为暴露开发；E1方法家族先冻结、新32parent后生成，完成实现后评估；E2复用同32parent为明确暴露跟进；A四例为独立小规模交接、非32parent泛化证明 |
| 原包数值门槛 | GATE_DECISION只汇总真实科学结论，不引入额外审批流程或为了发布重跑实验 |
| E2_FROZEN的A导入预估200次 vs E3实际240次 | 运行前以E3_IMPORT_FROZEN明确240硬子上限，增加双频corrector及参考端点验证；仍计入本次4096总上限，原200估算保留为历史预估，不覆盖旧配置 |
| 新history参考数值不稳 | 保留旧表，增加一次限定独立审计，不放宽准则；Mode2不称完整posterior或upper bound |
