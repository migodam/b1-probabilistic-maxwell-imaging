# 最终文字与数字一致性只读核验

核验范围：`RESEARCH_REPORT_ZH.md`、`SCIENTIFIC_DECISION.md`、根目录四份 `P1_B_*.md`，对照已冻结 `offline_scores/SUMMARY.json`、`PAIRED_CI.json`、`R2/RESULTS.json`、第一 joint 见证的 `RESULTS.json` / `GRID_COMPARISON.json`、C 的 `PHASE_TABLE.csv`。父线程随后明确许可只读已冻结 `inference/RESULTS.json` 提取 branch/CDF 各项最大误差；没有打开其中引用的任何数组路径。行号为本次读取的一基行号。只写本文，没有修改报告、结果、源码或账本，没有 GT/生成计划读取、Book、PDE 或新物理。

结论：**没有发现本次允许来源中主要结果数字转写错误。两项轻微措辞需明确；两项已登记 tail 的结果待最终插入。** 这不是科学裁决，不改变父线程的 REDIRECT 或各项条件结论。

## 1. 应明确比较对象：branch ensemble 误差是相对 finite

位置：`P1_B_TASK_RISK.md:7`。

原文：“branch ensemble最大误差direct0.012936、progressive0.017026”。前两句刚讨论相对 continuous 的 CDF 误差，接着省略 branch 的比较对象，读者可能误以为这也是相对 continuous。

实际独立核对如下：

|字段/比较对象|direct 最大值|对应 observation|progressive 最大值|对应 observation|
|---|---:|---:|---:|---:|
|`ensemble.*.branch_mass_max_error_vs_finite`|0.012935831202173786|15|0.01702564516972327|3|
|`max(abs(ensemble.*.average_branch_mass-continuous_ref.branch_mass))`|0.01311247398510551|15|0.017093751414709257|3|
|`ensemble.*.CDF_max_error_vs_finite`|0.012935831202173786|15|0.01702564516972327|3|
|`ensemble.*.CDF_max_error_vs_continuous`|0.01745799848616547|24|0.020374264044473928|21|

来源为 `results/P1B/inference/RESULTS.json`，按上述 observation 与 ensemble.direct / ensemble.progressive 字段定位。

**原数字正确，先前疑似把CDF数值误称branch的怀疑撤回。** 本例 gap 处 CDF 与 branch 质量对应，使二者相对 finite 的全体最大值恰好相同。建议仅补为：“branch ensemble**相对finite参考**的最大误差……”；不改任何数字。父线程已表示会补此限定。

## 2. 轻微单位措辞：独立的是噪声实例，观测共享 parent

位置：`P1_B_2X2_ABLATION.md:3` 与 `P1_B_TASK_RISK.md:9`。

- 前者“32份独立新噪声观测”可能被理解为32个独立总体观测单位；建议“16个parent各生成两份独立噪声实例，共32份新观测”。
- 后者“32个相关噪声副本”又可能被理解为噪声向量本身相关；建议“32份共享parent的观测，不是32个独立校准parent”。

这不是 bootstrap 或算术错误。`PAIRED_CI.json` 明确以16个 parent 配对，先平均两次noise；`RESEARCH_REPORT_ZH.md:55,73`、`P1_B_2X2_ABLATION.md:3` 的后续句也已经正确说明。仅建议去掉两种相反的字面歧义。

## 3. 已登记圆腔/argmax：待最终插入，不是缺失实验

位置：`RESEARCH_REPORT_ZH.md:39,41`。

- 39行“独立argmax网格核验在最终验证报告另列”：现有第一 joint `RESULTS.json` 确已保存24项任务最大值 `d_T_sup_24=4.81819338103861`、`task_argmax=6`。这不等于另外的独立 argmax 网格核验已经完成。父线程说明后者正在已登记 tail 中，最终应插入其真实结果/来源；完成前可明确写“已登记，待 tail 完成后插入独立网格核验”。
- 41行“其最终资格见B2正式表”：圆腔 fine 验证也属于已登记队列，最终资格和数字尚待插入。本轮不读其新输出、不猜测成功或失败，不将这个正常待完成项列作“缺失实验”，也不要求追加实验或预算。

父线程最新说明为Green全部112行已完成、tail将开启圆腔。本文只记录待合并位置，不据此评价Green结果或圆腔资格。

## 4. 已确认正确的关键转写

### R2 中心与单侧 FD 已分开

`RESEARCH_REPORT_ZH.md:17` 没有把单侧 FD 写成约3e−9。允许源 `results/R2/RESULTS.json` 的逐状态最优记录为：

|case|scheme|三个step中的最佳 relative_error|
|---|---|---:|
|air|one-sided|2.4154926393553342e−6|
|weak|central|1.069975687984216e−9|
|high|central|3.080534886241988e−9|
|box_boundary|one-sided|1.0257637132352606e−5|
|area_boundary|central|1.381330851201064e−9|
|cavity|central|2.122417884162132e−9|

因此“中心差分各状态最佳误差的最大值约3.08e−9”正确；单侧2.42e−6与1.03e−5正确。它是跨三个step择优后的摘要，不是所有step都达到此误差。adjoint最大缺陷为3.999978156336922e−14；“约1e−14”是数量级描述，没有发现单侧精度偷换。

### B2 joint / cal 与几何资格没有混同

`RESEARCH_REPORT_ZH.md:37–39,96–98` 与 `SCIENTIFIC_DECISION.md:5–6` 的核心数字一致：fine joint `0.7987161019689071`、cal `3.5605119112114565`、target shift `.1`、mass change `0`。joint world0/1网格差 `0.005674623106771554` / `0.015201453253063586`；cal分别 `0.008109629255612308` / `0.017857310015052252`。

原验证字段为 `joint_physics_threshold=true`、`cal_physics_threshold=false`、`strict_original_prior_GO=false`、`conditional_fixed_support_geometry=true`、`posterior_probability=null`。报告没有把joint见证写成cal近等价，也没有从其存在推断posterior mass。

### C 固定对的孔径/带宽数值正确

`RESEARCH_REPORT_ZH.md:45` 对应 `PHASE_TABLE.csv` 中 loss scale=1、joint_contrast3：narrow full `0.7987161019689071`、central `0.3609093994299449`；wide full `7.086518859632612`、central `3.599735443053257`。表中全孔径240、central60复通道，30行固定对全部满足central范数不超过full；文本明确只有孔径严格嵌套，没有把非嵌套频点当数据追加定理。`49`行contrast7的full narrow `5.0935`也与表值 `5.093473559730694`相符。

### 七臂风险、CI、prior 与 KKT 范围正确

- `RESEARCH_REPORT_ZH.md:63–75`、`P1_B_2X2_ABLATION.md:5–19` 的全部7臂MSE/MAE按保存值正确舍入。
- direct point→direct SMC 的MSE相对下降为 `42.54034958713797%`；强连续point→direct SMC为 `37.03053721205589%`。对应主区间分别 `[0.00043417855174534414,0.0034984395145395077]` 与 `[0.00029833710995545183,0.0027941375609500966]`，均转写正确。
- 强point→SMC的MAE差为 `−0.0017165884033550766`，区间 `[−0.012960217246638497,0.010142021197895233]`；“MAE没有确定改善”没有把负的点估计改成正收益。
- direct减progressive SMC 的MSE差 `−1.5016154188923992e−6`、区间 `[−1.7841531637676967e−5,1.592213763656574e−5]` 正确，没有把调度跨0区间写成明确改善。
- prior MSE `0.0023661589336290852` 不仅出现在主表71行，也在75行、消融表12行及19行、科学裁决9行明确讨论。prior→direct SMC差 `−0.00017894859036529293`、区间 `[−0.000739362424696961,0.00038532434189392505]` 正确；没有隐藏零观测prior或声称已建立观测增益。
- `RESEARCH_REPORT_ZH.md:73,79` 与 `P1_B_PROGRESSIVE_ALGORITHM.md:12,18` 正确区分 finite target、连续参考、profile point与平方损失后验均值动作；明确没有Full KKT/global MAP认证，没有把gap内均值称第三个物理世界。
- `P1_B_TASK_RISK.md:7–9` 的每schedule128 runs、continuous CDF超差36/39、progressive ensemble observation21的0.020374例外及prior/SMC经验区间宽度均与已审核输出一致；没有以高经验覆盖宣称精确连续校准。

## 5. 原先验范围与本次没有重复认证的数字

`RESEARCH_REPORT_ZH.md:7,25,53,57,79,94,100`、`SCIENTIFIC_DECISION.md:5,7`、`P1_B_PROGRESSIVE_ALGORITHM.md:3,12,24` 均保留原minimum_feature未定义/未PASS、Oracle非online、新proper fixed-support prior非原宽墙posterior、家族外质量null等限制。没有发现将原minimum_feature回写PASS或将两个固定墙场夸成16种独立墙体几何的语句。

本轮允许源没有覆盖R1/R3全表、Born/压力测试单独表、warm cost replay、银行/ledger最终收费、全部原始字段inventory。因此本文不对这些外部数字（例如主报告83行warm wall、110行约224GB）作“已独立再次验证”的声明，也不据来源不在本轮清单而判其错误。它们由父线程相应来源审计及最终封口负责；本文没有扩大读取范围。

最终科学裁决仍完全归父线程；本核验只交付上述定位、数值一致性和有限措辞建议。
