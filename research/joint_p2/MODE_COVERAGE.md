# 材料事件覆盖与未显示质量

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

这里的“hypothesis”是一个明确material event及代表材料；事件质量不是代表图像的点概率。“mode coverage”按预注册事件阈值评价，不把十二个bins自动称为十二个后验modes，也不按观测接近程度把材料不同的世界合并。

事件由四个central-radius bins `{0,(0,.10],(.10,.20],(.20,rmax]}` 与三个annular-halfwidth bins `{0,(0,.06],(.06,hmax]}` 形成，ID=`centerbin*3+annularbin`。生成器没有收到这些标签。rmax=.292551825R、hmax=.117020730R；材料图明确为球对称中心截面。定义、prior mass与材料尺度见 [protocol](configs/protocol.json)。

实际metric规则在模型冻结、fresh evaluation前保留于 [EVALUATION_RULES.json](configs/EVALUATION_RULES.json)：reference事件质量≥.05为重要事件；count discovery要求8192模型draw至少一条，mass-aware discovery要求model mass≥.01，false exclusion要求重要事件model mass<.005；false confidence要求model maxmass≥.9且reference有至少两个重要事件。它们是有限scope的统计指标，不是“全模态发现”的证明。

## 可信reference输出

large-cavity C package的Top4选event9、0、6、1，质量分别.677226、.162540、.071828、.041259，represented mass **.952852**，unshown modeled mass **.047148**。solid C package的representatives累计质量 **.964808**，仍显式保留 **.035192** remainder。K依据全十二事件质量自适应，不预先强制三张图，不将展示项重新归一化到1。

代表材料是事件内部在体积加权复材料距离下选出的有限medoid，并通过full-Mie endpoint验算。gallery显示ε实部、损耗与几何；forward prediction、bounded conditional gain、noise-distance关系和conditional covariance随JSON一并交付。每个连续event可包含许多结构/ε解释，medoid不是全部条件不确定性。

主要数据：[large-cavity package](results/packages/large_cavity.json)、[solid package](results/packages/solid.json)、[all packages](results/packages/)、[full reference](results/FULL_MIE_CORRECTED_REFERENCE.json)。它们分别声明full residual或proxy target。模型外质量保持`null/unquantified`；不能用“已建模事件和为1”来证明现实世界没有其他材料支持。

径向位置概率与材料不确定性另保存于 [MATERIAL_FIELDS.npz](results/MATERIAL_FIELDS.npz)：七case的`P(vacuum at r|y)`、材料实部/损耗均值与SD。它们由实际wholeweightedpackage在同心支持内计算，不是生成多张外观合理图，也不能把径向概率场误称一般三维位置概率。fullcorrectedcase采用已披露的corrected事件质量与proxy事件内样本，场的不确定性仍是有限数值近似。

## 生成器fresh覆盖

64 fresh parents、128observations、三个fixed seeds上，full与diagonal distillation按本事件划分达到重要事件mass recall1、false exclusion0、false confidence0。full平均eventMAE约.013215，diagonal约.013028；平均Top-K覆盖的**reference**质量分别约.97764、.97481。Gaussian distillation虽然重要事件recall也为1，eventMAE约.027550，说明粗事件发现不等于连续密度准确。

full/diagonal distillation的高覆盖只对该support、noise、12事件分辨率与fresh panel成立；没有测试更稀有、更窄或一般未知形状的现实材料。其rawdraw概率属于learned proposal，对比的是共同proxy reference，不自动获得full-Mie posterior身份。

最终eventMAE与重要eventrecall采用 [REFINED_FRESH_EVENT_COMPARISON.json](results/REFINED_FRESH_EVENT_COMPARISON.json) 的acceptedreference。全128event/mean/SD已核对，原13阶三例升级17后对照21通过。原Top-Kcoverage、latentKL、truthcoverage与其他旧13统计继续标历史近似，不把未经重算的列合成一个“全部accepted”的概率摘要。

joint-training density存在明显退化：diagonal joint的三个seed false-exclusion观测率约2.34%–4.69%，false-confidence约.78%–3.13%；full joint虽然重要事件count/mass recall高，但eventMAE .045–.053与latentKL15.6–25.4明显不可靠。用“都命中过事件”不能掩盖错误质量。

原fresh node-resampling CDF值不作为continuous-material准确性证据。最终positivepartition comparison的full/diagonal重要eventconditionalgridCDF .141761/.161768，明显未达.01；marginal误差 .098906/.101614。原截断CDF的.00265042非单调失败保留，最终同密度正区间累计CDF单调性缺陷0；不同partition的evidence/event与原13阶估计误差及subset-onlyrefinement完整披露，仍只作empiricalfixed-grid比较。事件MAE达标不能替代conditional CDF/SD/tail正确性。

## 高ESS与合法correction的负证据

两个generic tempered SMC runs均到达beta1，分别295/301full calls，但event0质量约.296875/.086544，而full reference为.162540；event6约.046875/.182959，而reference为.071828。第一个末端ESS约64，仍有显著质量误差。合法MH/SMC只说明声明的target修正，有限64particles仍有seed变化、rare-event遗漏和normalizing-constant误差。

三个matched64-call的MDN full-likelihood correction均发现了三个≥.05reference事件，但不同事件仍可零访问，并出现错误质量/较大MC区间。未访问事件的bootstrap区间可能退化为零，绝不能据此认证其posterior被排除。

本项目故意保留这两类失败：重要材料世界遗漏（false exclusion），以及多解释被压成单一自信输出（false confidence）。数值参考及界面能够揭示它们；系统没有普遍消除它们。
