<!-- P3LR10_CURRENT -->
## P3L-R10-BI：双边传递、条件物理复用与材料覆盖边界

[阅读完整中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r10_bilateral_structure/) · [140条公式集](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/maxwell_formula_atlas/) · [无损证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr10-bilateral-20261030)

固定形状的320次新幅度查询全部通过误差门；实际32节点冷形状费用节省51.82%。新材料对象的覆盖优势未成立，Full66 posterior仍未认证，强TAP贡献HOLD。完整源码、失败、独立QA与论文草稿已保存。
<!-- /P3LR10_CURRENT -->

<!-- Z_ALG_R2_CURRENT -->
## Research Z-ALG-R2：行列式商代数、材料解枚举与盲成像

四材料六维商与六材料20维商求解器已实现。六材料新盲测返回2/1/1个合法根，固定64起点对照为1/1/1；双方包含GT，额外材料相距1.15879 RMS。三目标单批次0.1044秒对1.5180秒，仅支持本协议。48单元成像及SVD噪声失败独立保留。53,951完整复RHS；全球新颖性与严格全根认证仍HOLD。

[中文报告](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r2/) · [15组推导公式](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r2/FORMULA_COMPENDIUM.html) · [原始Markdown](research/research_z_alg_r2/RESEARCH_REPORT_ZH.md) · [阅读包与完整代码数据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r2-20261030)
<!-- /Z_ALG_R2_CURRENT -->

<!-- W_GQ_R2_PUBLICATION -->
## Research W-GQ-R2：墙体未分辨，目标分辨率与合法补偿

微结构材料 RMS 差 1.13137，held-out 目标族最大传播白化距离 0.159971；受控有限参考定位 0.079224 λ₀，宽墙材料仍可补偿 0.05 λ₀ 目标位移。一般未知墙在线与净加速未成立。722 FE / 5212 RHS，全部失败保留。

[完整中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gq_r2/) · [中文 Markdown](research/research_w_gq_r2/REPORT_ZH.md) · [GPT Pro 交接](research/research_w_gq_r2/NEXT_PRO_HANDOFF.md) · [完整原始证据与复现包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-gq-r2-20261029)
<!-- /W_GQ_R2_PUBLICATION -->

<!-- P3LR2_CURRENT -->
## P3L-R2 已完成：监督材料风险改善，完整后验仍未合格

[完整中文报告](research/p3l_r2_risk_posterior_exploration/RESEARCH_REPORT_ZH.md) · [网页与固定材料图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r2_risk_posterior_exploration/) · [GPT Pro研究交接](research/p3l_r2_risk_posterior_exploration/GPT_PRO_HANDOFF/START_HERE.md)

32,768独立训练材料、2,048验证材料、两批分别冻结的256全新测试材料。首批相对强Probit点的材料MSE降低25.27%（调整97.5%区间22.34%–27.99%）。第二批确认Maxwell训练谱相对固定随机decoder包络增量约0.70%。FULL posterior仍UNQUALIFIED；轨迹风险平均未改善；不把风险图称为新材料世界或后验样本。
<!-- /P3LR2_CURRENT -->

<!-- P3L_CURRENT -->
## 已完成：P3-L 非线性材料概率成像 · REDIRECT

[中文报告](research/p3l_local_posterior_imaging/RESEARCH_REPORT_ZH.md) · [实验网页与图像](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_local_posterior_imaging/)

128独立对象：当前近似材料均值MSE .0615021，对应点基线 .0110309，显著恶化。11组四链Full参考均未合格，15次神经训练无验证增益；8,374,472 complex RHS封账，GPU0。保留负结果，停止当前Gaussian/NN/A17增量路线；不把近似失败等同于真正Bayes成像无益。
<!-- /P3L_CURRENT -->

## 历史：P2材料探索完整Campaign

[中文报告](research/p2_material_fiber_20261019/RESEARCH_REPORT_ZH.md) · [网页与材料图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p2_material_fiber_20261019/)

24个新物体，15个合格anchor；有限噪声替代存在，强基线削弱软谱优势。整体REDIRECT，概率不估计；2236states/46728RHS封账。

## 历史：P1零新增物理审计与下一最小实验

[中文报告](research/p1_zero_physics_audit_20261018/RESEARCH_REPORT_ZH.md) · [网页与冻结协议](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p1_zero_physics_audit_20261018/)

新物理0：材料尺度弱方向仍在；A/GT不是重要模式见证；B支持标准trust-region修复。A 12states/288RHS与B 16states/640RHS仅为建议，未授权、未运行。

## 历史：P1材料先验与Maxwell target实现

[中文报告](research/p1_material_target_20261018/RESEARCH_REPORT_ZH.md) · [网页与先验材料](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p1_material_target_20261018/)

材料构造已完成；新增真实物理资格26states/976RHS：1/2基础材料合格，B未收敛。无新替代world，不训练NN。[真实物理报告](research/p1_material_target_20261018/PHYSICS_QUALIFICATION_REPORT_ZH.md)。

## 历史：A17-X × F3 竞争审计

[中文报告](research/a17x_f3_competition_audit/RESEARCH_REPORT_ZH.md) · [网页与数值图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/a17x_f3_competition_audit/)

零新增物理：17份旧F3增量均rank575；小端口假设不能由七维latent推出。合成机制修订完成，真实P1及新预算未就绪，比较NOT_RUN。

## 历史：Research F3 商坐标与有限噪声探索

[中文报告](research/research_f3_quotient_exploration/RESEARCH_REPORT_ZH.md) · [网页与实际材料图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_f3_quotient_exploration/)

55 states / 712 RHS封账：两个新资格物理驻点，weak小幅竞争端点，无新拓扑、高阶方法无同成本增量。结束原P0局部Explorer路线，不训练NN。

## 历史：Research F 材料几何与独立物理反演

[中文报告](research/research_f_intrinsic_transport/RESEARCH_REPORT_ZH.md) · [网页与材料图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_f_intrinsic_transport/)

24 states / 960 RHS封账：分离坐标与真实prior曲率；独立物理反演合法且噪声相容，但驻点资格未过，F2未运行。不训练NN。

## 历史诊断：有限Maxwell反馈与GT差量

[中文报告](research/finite_feedback_fiber_diagnostic_20261016/DIAGNOSTIC_REPORT_ZH.md) · [网页与数值图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/finite_feedback_fiber_diagnostic_20261016/)

45个二维states、1880个complex RHS；高阶响应已核验，GT合法局部路径未建立。材料decoder曲率是当前关键限制。

## 历史研究：LXZ 局部–全局多解验证

[中文报告](research/research_lxz_20261016/RESEARCH_REPORT_ZH.md) · [网页与材料图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_lxz_20261016/)

96个二维states封账，未发现显著结构新world；已有C1提案经正确MH为本次proxy基线最佳。N3合法anchor与分离mode仍未证实。

> 2026-10-15 原先验竞争解释：[本轮入口](research/prior_constrained_20261015/START_HERE.md) · [中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/prior_constrained_20261015/)。E0/E1已执行，原先验内0组合格竞争见证；45states/680复RHS封账，E2未启动，N3保留。

> 2026-10-15 Research G：[下一解搜索中文报告](research/research_g_20261015/START_HERE.md) · [材料图与网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_g_20261015/)。三个观测新增一个合法边界结构；Best-of-K未改善，历史引导无稳定优势。291states/6936复RHS，暂不扩大。

> 2026-10-15 Research T：[物理、基础解与补偿成核](research/research_t_20261015/START_HERE.md) · [公开成像页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_t_20261015/)。8/8合格默认base，15/16首轮补偿成核；材料MAE没有改善，成本192states封口。

> 2026-10-15 独立并行研究：[Research B 主动Maxwell测量](research/research_b/START_HERE.md) · [中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_b/) · [原始证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-b-20261015)。真实E0/E1/E2与四例A交接已执行；复杂选测未稳定胜过强固定基线，候选补全更值得继续，reference限制完整披露。

> 2026-10-15 新研究：[Research A 多解材料成像](research/research_a_20261015/START_HERE.md) · [实际图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_a_20261015/)。每观测3–10个径向候选；Best-of-4有改善，A30反证与机制限制完整保留。

# 历史研究：B1 × C1 Joint P2

[中文报告](research/joint_p2/RESEARCH_REPORT_ZH.md) · [交互Top-K图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/joint_p2/) · [科学裁决](research/joint_p2/GATE_DECISION.md) · [完整证据下载](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/joint-p2-20261014)

同一连续材料prior保留实心/中心/环形空腔解释，提供事件质量、条件不确定性、full-Mie补偿与独立传感器更新。新增3459/4000次物理attempt；C1未达到完整连续posterior目标。下面原有内容保留为历史研究。

# 历史入口：B1 E0

已完成零新增物理的连续posterior摘要修复：[E0入口](research/b1_next/e0/START_HERE.md)、[独立中文报告](research/b1_next/e0/RESEARCH_REPORT_ZH.md)、[资源与原始数组](research/b1_next/e0/results/SUMMARY.json)。旧报告/图/384次账本保留；新阶段公开更新到同一仓库。

以下为历史固定模板baseline入口。

# B1：多假设 Maxwell 概率成像

这是一个实际运行的3D矢量Mie概率成像研究原型。输入16个复散射通道，在三个固定同心材料模板和连续复介电常数先验下，计算各拓扑事件的后验质量、代表性实/虚部材料切面与不确定性。

- [完整中文研究报告](RESEARCH_REPORT_ZH.md)：背景、数学目标、两轮实验、失败归因、结论。
- [在线重建图与概率面板](https://migodam.github.io/b1-probabilistic-maxwell-imaging/)。
- [离线面板](dashboard.html)；[主要重建图](figures/fig01_competitive_representatives.png)。
- [理论](THEORY.md)、[结果](RESULTS.md)、[科研裁决](DECISION.md)、[研究日志](RESEARCH_LOG.md)。
- [复现方法与可执行入口](README.md)；[原始概率表](results/CASE_PROBABILITIES.csv)。
- [逐次物理资源账本](results/PHYSICS_LEDGER.jsonl)、[运行账本](results/JOB_LEDGER.jsonl)。

竞争案例：中央空腔50.64%，中央空腔＋环形空腔40.43%，实心剩余质量8.93%。Top-2不重新归一化。百分比是连续材料拓扑分支的近似后验质量，单个精确材料图像没有非零质量。模板外几何的质量尚未估计。

结论：受限Top-K原型已实现；一般几何恢复、严格全域概率认证、真实仪器统计校准与强TAP原创贡献尚未建立。原始参考失败、错误重放、历史负结果与访问缺口均保留。384/384次新物理评估全部记录；没有大型NN训练或旧locked数据访问。

## 给网页版GPT的阅读顺序

先读本入口和中文报告，再读THEORY.md与RESULTS.md；涉及概率可靠性时查阅DECISION.md、results/refinement/SUMMARY.json及research/delegated/probability_review/REVIEW_NATIVE.md。不能把公开发布状态、软件测试或development重复观察当作独立科学确认。
