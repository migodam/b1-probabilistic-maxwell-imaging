# 本轮原任务和现状

本轮的唯一研究对象是从多照明的带噪端口读数直接证明无损介质某处或指定区域的介电对比度下限。用户原文见 ORIGINAL_REQUEST.txt。旧证据入口为 ../z_alg_secant_tasks/highdim/sector_inertia。按独立材料与实物计数。

- Astra 保管证明和结论。三组 Sol 分别负责 core、synthetic/roi/multifrequency、fresnel。Neon 从工作区外独立的 Z-ALG-Secants-Neon 项目通过 CLI 自行读研究证据、做创意建议。
- 原五份存档，六个预定的新形状（全/半环、2%/10%配对）、三个 Fresnel 实物（4/6/8 GHz）构成证据。用 public 资料选乘子，离线评价代理误差及优化器输出。离线已知点估计可用于构造反驳松弛的存在性见证，并明确署名。
- results/theory/hidden_radiation.json 与 all_scalar_multipliers.json 含两个全环上对全部标量乘子的严格限制。core 保存 32² 扫描和证书，roi 保存 16 区图、严格负惯性下限；fine 直接在 96² 正演算子上证实峰值。fresnel 按实际接收交集比较；单列空场误差传播和无标定硬上界的限制。
- 恢复时先看本目录 RESEARCH_REPORT_ZH.md、MATHEMATICAL_THEORY.md、EXPERIMENTAL_RESULTS.md、NEXT_PRO_HANDOFF.md，再根据具体问题读有数值的源码/JSON。发布沿用项目研究网页。

最终结果：在 96² 原生网格上，六件半环介质的峰值大于 0.5，2% 噪声 6/6、10% 噪声 3/6。两件 2% 全环物体在 0.75 处严格超出所有全局单源标量乘子的能力。区域 0/112，Fresnel 0.5 为 0/9。报告、复现、英文论文草稿和下一轮交接均在本目录。
