# MAX-STRUCT R15：twofold分母与周期DR

[中文完整报告](RESEARCH_REPORT_ZH.md) · [轨迹与反例](TRAJECTORY_ANALYSIS.md) · [精确公式](THEORY_DERIVATIONS.md) · [下一轮Pro入口](PRO_HANDOFF.md) · [费用](COST_LEDGER.md) · [复现](REPRODUCE.md)

本轮已实际完成264主反演、56同状态干预、16读出反事实和44终点验证。24对象均是已公开R14材料，不能称新盲确认。没有新增NN训练、三维或实测。

同120秒整进程CPU，原primary规则的过线数：CIE9、Dual9、G1911、S19/U11/邻域/路由DR各12、独立DR13，均/24。S19相对G19救回4个、丢失3个；融合DR相对邻域版救回0020、丢失0021。正面恢复和严重反例均保留，不以一个成功数代替结构评价。

U11并非全局G/S共同SVD，本批S强G弱空间为空。current弱空间不是独立材料集合；周期融合还须重新注入原Fourier/DPCS表示。早期反射材料变化大多落在实U01外，注入误差及错误活动面解释了需要继续研究的具体接口，尚未证明完整失败因果。

公开网页由本目录报告和真实数组生成。大数组通过同仓库release分卷提供，下载清单以发布后实际匿名回读的DOWNLOADS.json为准。冻结算法源提交ba6d68a87；报告、分析和发布版本随后另记。旧R4–R14原始结果不改写。
