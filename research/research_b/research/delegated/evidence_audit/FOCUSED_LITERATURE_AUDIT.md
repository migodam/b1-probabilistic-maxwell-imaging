# Research B 最近邻：五篇定向核对

日期2026-10-15。使用 `scholarqa-research` Evidence QA；primary transport为已安装`scholarqa-cli verify`，对包给出的五个DOI批量返回canonical paperID、title、authors、year、externalIds、abstract；5/5 resolved、unresolved 0。没有宽泛新检索，不替主会话作新颖性裁决。原始网页仅核对包已列链接及同篇全文，不引入新文献。

## 主张账本

| 包ID | 核对后的论文身份 | 阅读层级、实际支持 | 本轮不得声称 |
|---|---|---|---|
| R01 | A. C. Atkinson, V. V. Fedorov (1975), *The design of experiments for discriminating between two rival models*. DOI[10.1093/biomet/62.1.57](https://doi.org/10.1093/biomet/62.1.57)；S2 `7f52876bbcc17547f7718719af9ff3ec341792b0` | **Tier C**：身份完全解析；S2 abstract被publisher移除，DOI网页读取失败 | 本次未独立核对“对rival parameters重新拟合”的全文细节；该细节仍是源包的阅读声明，不能包装为本次全文验证 |
| R02 | A. E. Fouda, F. L. Teixeira, *Bayesian compressive sensing for ultrawideband inverse scattering in random media*. Journal **2014**，preprint 2013；DOI[10.1088/0266-5611/30/11/114017](https://doi.org/10.1088/0266-5611/30/11/114017)；S2 `f903578812c7508ce3949fa1e68ec558c4761a21` | **Tier A**：[作者全文§IV-B、Eqs30–31](https://arxiv.org/html/1401.1092)，按posterior covariance选下一传感器、myopic placement；另有Born与distorted-Born框架 | 主动EM测量、相同传感器数量下比较uniform均有先例；不说明其已解决本轮多世界共享校准问题 |
| R03 | A. I. Sandhu, B. M. Dia, O. Dorn, P. Soupios (2023), *Bayesian Experimental Design for Efficient Sensor Placement in Two-Dimensional Electromagnetic Imaging*. IEEE Access11,65649–65662；DOI[10.1109/ACCESS.2023.3288997](https://doi.org/10.1109/ACCESS.2023.3288997)；S2 `3e49504d3f091aacf07c6961516c74d376f8200b` | **Tier B**：[机构摘要](https://pure.kfupm.edu.sa/en/publications/bayesian-experimental-design-for-efficient-sensor-placement-in-tw/)明确EIG/KL、Laplace与sensor-network geometry设计；本次IEEE PDF访问失败 | 本次未独立重读其VIE/GMRES细节；源包称先前已读全文需保留为source-package attribution。机构页copyright2013与Published2023冲突，以正式bibliographic year2023，不误抄copyright |
| R05 | A. Bartuska, L. Espath, R. Tempone, *Small-noise approximation for Bayesian optimal experimental design with nuisance uncertainty*. Journal **2022**，preprint2021，arXiv v4修订2023；DOI[10.1016/j.cma.2022.115320](https://doi.org/10.1016/j.cma.2022.115320)；S2 `0149f7f6036cb296e2d41375162469900f9bbf7b` | **Tier B**：[作者摘要与journal reference](https://arxiv.org/abs/2112.06794)支持nuisance内层积分的small-noise近似、Laplace降本及EIT案例 | nuisance-aware design不是B独有；近似不认证多峰Maxwell任意参数/补偿范围 |
| R13 | M. Naghshvar, T. Javidi, *Active sequential hypothesis testing*. Journal2013，preprint2012；DOI[10.1214/13-AOS1144](https://doi.org/10.1214/13-AOS1144)；S2 `3855abd3c3ddac3dd171e82f17a7ca66a0ce022d` | **Tier B**：[作者摘要](https://arxiv.org/abs/1203.4626)支持adaptive actions、wrong-declaration penalty、total cost与有条件渐近policy分析 | sequential停止/误判代价不是新概念；有限候选hypothesis的结论不保证开放材料世界coverage |

这里每条实质支持来自明确abstract或实际全文段落；标题/元数据不足以支持方法细节。表中“不得声称”是审计推论，不是论文作者原话。

## 重要metadata冲突

S2返回的year分别1975/2013/2023/2021/2012，后三项与preprint年份混用。Fouda和Bartuska的arXiv `Journal reference`明确给2014和2022；Naghshvar保持Pro包的journal2013并分列preprint2012。不要机械把数据库year写成journal发表年。全部五篇的title/author/DOI身份吻合。

## 对本轮主张的约束

已有证据足以要求把ordinary Fisher/EIG、uniform、nuisance-aware以及停止/误判代价作为真实竞争基线。它没有证明Research B的新准则会胜出，也没有完成Maxwell-specific novelty认证。可检验机会仍是：在明确有限候选、共同物理与calibration合同中，补偿冲突是否带来更好的材料选择和完整成本；这必须由本轮真实实验回答。

核对命令、稳定标识与取得层级另存 `LITERATURE_VERIFICATION.json`。未执行Semantic Scholar snippet search或新graph扩展；本任务限定为既有五篇核对。本次未重读其余R04/R06–R12/R14–R18，不能把源包的阅读层级转记为本次已完成。
