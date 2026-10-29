# Research W-GQ-R2 — Target Resolution through Unresolved Dielectric Walls

实际二维 scalar-TM 全波开发研究，2026-10-29。物理机制与算法经济性分别裁决。

- 两宏观墙、两已解析微纹理：墙内材料 RMS 差 1.13137；最细层 12 held-out 目标、四发三频最大白化响应差 0.159971。
- 4 墙×12 目标 proper finite reference：目标均值位置 RMSE 0.079224 λ₀，墙纹理熵仍 0.998388 bit。受控有限证据，不是一般未知墙在线识别。
- 合法墙材料重新补偿 0.05 λ₀ 目标位移：三层离散白化距离 0.889152 / 0.912855 / 0.913910。
- 返回-only rank 8 最坏白化误差约 2.05e-4；普通 Schur/T-matrix 同享结构。未知宏观墙在线定位、可靠未查概率包络和冷部署净加速尚未成立。
- 完整科研费用：722 FE / 5,212 complex RHS / 665.087712 wrapper CPU 秒；峰值 RSS 3,348,283,392 bytes。失败、拒绝与修正保留。

网页：https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gq_r2/

中文原文：https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_w_gq_r2/REPORT_ZH.md

## 完整复现包

`Research_W_GQ_R2_Results.zip`：400,963,935 bytes，4,427 entries。解压入口 `Research_W_GQ_R2/`。包括报告、附件理论、实际复场/响应、材料、Jacobian、完整有限 posterior、查询与物理账本、源码、失败和作业回执。没有 locked 数据。

最终科研费用以包内 `delivery/FINAL_ACCOUNTING.json` 为准；早期 PACKAGE_RECEIPT 是历史快照。科学归档保留原完成时“本地未自动发布”状态，本次发布由用户另行授权，记录于 Git 的 publication/。

这是有限二维开发合同下的离散全波证据，不是连续 PDE 认证、任意墙目标唯一性、一般 posterior 或新算法加速证明。
