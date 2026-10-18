## 2026-10-18更新：真实P1资格已运行

[真实物理报告](PHYSICS_QUALIFICATION_REPORT_ZH.md) · [网页](physics_report.html) · [实际重建图](physics_gallery.html)

26states / 976RHS / 5.477s封账：两个固定对象取得1个合法且数据相容的physics-only基础解，另1个未收敛。没有新发现替代world；A17-X比较未运行。下面保留原零物理材料构造记录，不能把历史状态当最新结果。

# P1 已实现：相关被动材料先验与共同 Maxwell target

2026-10-18。**P1现在有可运行生成器、明确prior law、192份实际材料样本、固定采集/noise、材料运动代码及计费的观测/反演入口。真实Maxwell观测和合格物理anchor尚未运行。**

[中文报告](RESEARCH_REPORT_ZH.md) · [网页版](report.html) · [材料图](gallery.html) · [模型与测度](MODEL_AND_PRIOR.md) · [下一次真实计算](PHYSICS_NEXT_RUN.md) · [实际指标](METRICS.json) · [复现](REPRODUCE.md)

- 独立P1：64维相关Gaussian场+2个uniform被动材料参数，不输入拓扑label，不作peak normalization。
- 原24×24 scalar TM、8sources/32receivers、固定gain1；仪器sigma=.00119211798514519，沿用已有尺度，不随新物体调噪声。
- 已生成128 calibration、32 development、32 material-check样本；这些是材料侧样本，不是192次Maxwell states或同测量多解。
- 6个固定材料点、两种实际RMS距离.02/.04：冻结步与material-native各12次，均达到距离；没有验证likelihood相容性。
- native优化器仅2/12正常收敛，10/12迭代上限；两方法各1/12达到预设结构变化尺度。距离匹配不等于探索成功。
- 接口当前0 states/0 RHS，实际调用会在导入物理backend前拒绝；没有训练NN。

**状态：MATERIAL_MODEL_READY；OBSERVATIONS_NOT_RUN；PHYSICS_ANCHORS_NOT_RUN。** 两个物理parent事前固定为dev000/dev001。下一资格阶段建议独立32states/1024RHS/600s CPU，尚未授权；不是复用旧F3或A17-X余额。
