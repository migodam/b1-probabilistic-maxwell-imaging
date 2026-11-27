# Austria 与 FoamTwinDiel 的原端口复算

**两个材料格上的精确有限修正可以略降训练端口残差。三个 Austria 观测中，等预算继续跑原物理目标降得更多，但六个配对端点的材料误差都升高。FoamTwinDiel 实测的两个起点给出相反的留出排序，无法建立有限修正的稳定预测收益。**

## 已保存的原件和旧费用

新试验直接读取 2026-11-20 top_venue 的已保存训练场、观测及 30 秒口径 light 端点。Paper-NIE 使用论文 NIE 的积分方程适配，Dual-CIE 使用 sourcewise 双 \(\beta\)。两条都是历史真实运行。新实验自各自保存端点分叉，后续共用原离散 Maxwell 数据目标和 smooth/TV。表中历史 CPU 和墙钟按实际消费，包含各自原求解及已付轻校正。观测和原输出未重造，同一 Austria 几何的三观测是配对物理条件。 [top_venue 原报告与合同](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/top_venue_killtest_20261120/VERDICT_REPORT_ZH.md)

|冻结观测|原端点|训练残差|材料相对 L2 / 半峰 IoU|未见源误差|旧 CPU / 墙钟 s|
|---|---|---:|---:|---:|---:|
|Austria εr=4，400 MHz 全环|Paper-NIE|.12476|.77718 / .64118|—|6.74 / 6.97|
|同观测|Dual-CIE|.12100|.78678 / .63768|—|8.95 / 9.06|
|Austria εr=3.5，600 MHz 全环|Paper-NIE|.15009|.76368 / .56403|—|7.51 / 7.78|
|同观测|Dual-CIE|.15860|.73092 / .60504|—|9.66 / 9.76|
|Austria εr=3.5，400 MHz，180°|Paper-NIE|.18251|1.58140 / .19936|—|3.93 / 3.98|
|同观测|Dual-CIE|.21181|1.51018 / .17932|—|6.71 / 6.81|
|FoamTwinDielTM 8 GHz，13 个训练源|Paper-NIE|.40771|—|.59880|29.63 / 30.21|
|同观测|Dual-CIE|.33569|—|.50145|31.18 / 31.81|

三个合成观测用原 2 m 区域、150² 生成和 32² 反演、20×40 全环或 10×20 半环，以及各自保存的 10% 定范数复噪声。600 MHz 条件约为 4 个波长宽，400 MHz 约为 2.67 个。用历史解析参考材料直接代入 32² 模型时，三观测原数据残差分别 .26731、.54169、.20663，跨网格误差显著。相对 L2 用真材质范数归一，幅度 NRMSE 用真对比幅度，IoU 用固定半峰；细网原观测全沿历史设计。 [原始指标与历史/新增费用表](results/classical/score_v2/OUTCOMES.csv)

R8 的作者式 CIE 一并作为历史参照。R8 原 εr=4 400 MHz 独立生成观测上，CIE 单 \(\beta\)、Dual-CIE、Paper-NIE 的材料相对 L2 为 .8316/.8257/.8350，CPU 为 16.85/17.22/10.60 秒。另一个 R8 εr=2 750 MHz 上，CIE / Dual 相对 L2 为 .3351/.3189，CPU 15.29/13.17 秒。R8 与 top_venue 的噪声实现、先验适配、辅助日程和轻校正上限不同，数值只在各自观测及历史合同内配对。 [R8 Austria 论文对照](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_max_struct_r8/PAPER_BENCHMARKS.md)

## 同状态最多两格的有限修正

从端点重算真实原场 \(U\)，按原物理盒投影梯度的两项最大绝对值选两个像素。逐格以及两格一起测试 .5、1.5、3.0 三个幅度，方向取负梯度符号，截到原 \([0,20]\) 范围并去重。令 \(L=I-D_m G\)，\(V=L^{-1}I_J\)，沿候选盒内改动 \(d\) 以 2 阶系统求端口差。

\[
\Delta Y_s = P_s S V (I-D_d (GV)_{J,:})^{-1}D_d U_{J,s}.
\]

\(P_s\) 是第 \(s\) 个源实际接收器选择。合成观测为给定接收阵，Foam 使用公开的逐源掩码。排名用原 \(\Phi=.5\|Y(m)-Y_{\rm obs}\|^2/\|Y_{\rm obs}\|^2+R(m)\)，只用训练复场。最好的两格提名再完整做一次原波动方程重算，并与保存端点同池选最小 \(\Phi\)。原物理续跑用原 top_venue 未缩放目标的 L-BFGS-B，按同准则保留已审计候选。冻结上限合成 1 CPU 秒、实测 4 CPU 秒，并另以每对实际两格计算完成的 CPU 前缀作精确配额比较。V1 沿用了 R13 的目标缩放，改变了线搜索。V2 在读取标签及留出以前修正为原 top_venue 的尺度，V1 完整保留为敏感性记录。 [冻结配置](configs/classical_v2.json) · [两版实际封存](results/classical/v2/PREDICTION_SEAL.json)

|Austria 观测 / 起点|两格训练残差|续跑训练残差|两格材料 L2|续跑材料 L2|
|---|---:|---:|---:|---:|
|εr=4 400 MHz / Paper-NIE|.12465|.11050|.77745|.78671|
|εr=4 400 MHz / Dual|.12100|.10938|.78678|.79187|
|εr=3.5 600 MHz / Paper-NIE|.14973|.12155|.76368|.76806|
|εr=3.5 600 MHz / Dual|.15860|.12865|.73092|.75134|
|εr=3.5 400 MHz 180° / Paper-NIE|.18184|.13742|1.58073|1.62454|
|εr=3.5 400 MHz 180° / Dual|.21116|.15235|1.50937|1.54184|

上述为共同 1 秒预算的当时最小原目标。两格臂在 0.115–0.158 秒内完成，用 2 次全阶 LU（含起点），最快前缀上的两格端口残差 4 胜、1 平、1 负，幅度很小。六端点在续跑中都出现数据残差显著下降与全图材料误差上升。两格未恢复新的 Austria 结构。两版独立重现八个旧端点复端口，重放最大相对差为零。V2 八个新两格候选的小方程端口预测与全波复算的最大相对差为 \(9.20\times10^{-15}\)。 [候选选择、真实端口与两种费用口径](results/classical/score_v2/OUTCOMES.csv)

## FoamTwinDiel 真正留出端口

复核 June 2004 Institut Fresnel 原七列 TM 文件，在 8 GHz 取散射端口 \(\overline{Y_{\rm total}-Y_{\rm empty}}\)。原始 4338 个端口，18 个源，每源 241 个实测接收器。13 个训练源的 13 组掩码各异，共 3133 个端口。源 1、5、9、13、17 整列留出共 1205 个端口。训练、留出与原文件逐元素一致，最大差为零。未测端口不进入目标或评分。空场多极照明的独立校准误差中位 .11490，原逐源照明沿用且在选择前冻结。真实复散射误差只在实际端口计算，原模型的 15% 参考球并未被解释为噪声置信区间或材料正确性。 [原数据来源](https://doi.org/10.1088/0266-5611/21/6/S09) · [原始行与掩码复核](results/classical/score_v2/RAW_PORT_ALIGNMENT.json)

|FoamTwinDiel 8 GHz|旧训练 / 留出|两格训练 / 留出|续跑同实际前缀训练 / 留出|续跑 4 秒训练 / 留出|
|---|---:|---:|---:|---:|
|Paper-NIE 起点|.407706 / .598796|.406112 / .597976|.407051 / .598051|.407051 / .598051|
|Dual-CIE 起点|.335688 / .501446|.335469 / .501237|.334771 / .500680|.298588 / .457328|

NIE 起点的两格留出值比续跑低 .000075。Dual 起点在与两格完全相同的 2.575 CPU 秒前缀已优于两格，完整 4 秒时留出值进一步降至 .45733。NIE 的实际两格前缀为 3.407 秒。八个观测和起点的配对中，Foam 的两组均来自同一件实物。留出考查同一试样在新照明方向上的复端口预测。这里不存在逐像素实材质真值。 [实测训练与留出误差图](figures/classical_foam_ports.png) · [留出复预测数组及有效掩码](results/classical/score_v2/HELDOUT_PREDICTIONS.npz)

SRAP 的历史匹配掩码对同一 Foam 实物提供另外一组真留接收参照。两种训练布局各 733 个端口，5 个同一留出源 1205 个端口及已见源的 390 个留出接收端口。原 Dual 在 30 秒 light 下，连续/交错布局的未见源误差为 1.0007/.2571，未见接收为 .9820/.1715，实际 CPU 22.39/23.96 秒。这个对照变更采样几何，成本和预测值与上表独立。R8 按另一套成熟 8 GHz air0 的 Dual 计算，训练/留出为 .08134/.13724，CPU 741.18 秒。它保留了相同实物上更长优化和不同接收合同的历史正例。 [SRAP 实测掩码与成本](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/srap_cie_array_study_20261121/RESEARCH_REPORT_ZH.md) · [R8 原实测报告](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_max_struct_r8/FRESNEL_REPORT.md)

两进程各单线程，V2 整批墙钟 13.81 秒，单任务 CPU 合计 22.09 秒，133 次 1024 或 3600 阶全 LU。超过 CPU 前缀才完成的全解计入实际成本，超前缀候选未参加选择。V1 另耗 22.51 CPU 秒及 122 个全阶 LU，未用于数值结论。评分端复用相同已付材料，并在此基础上独立重算 7 个留出端口场，7 个 3600 阶全 LU 和 35 复 RHS。原训练包、空场校准和旧端点费用按历史回执保留。 [评分端逐项新增费用](results/classical/score_v2/POSTSEAL_PHYSICS_COSTS.csv)

## 复现和源文件

通过 [输入路径清单](configs/classical_history_inputs.json) 恢复公开合成数据、历史源码和各自保存端点。原 Foam 测量以及含原测量的训练/留出缓存保留于用户本地。单一解压根可设置 HC_CIE_R2_INPUT_ROOT，历史 top_venue 的所在目录可单独设置 HC_CIE_R2_TOP_ROOT。如果原始 Foam 测量放在其他本地路径，设置 FRESNEL_2005_DIR。结果输出为合法候选、预测端口、数据聚合误差和掩码。

~~~bash
PY=Gaussian/.venv_nn/bin/python
ROOT=research/hc_cie_r2
REPLAY=.artifacts/hc_cie_classical_replay
$PY "$ROOT/src/classical_run.py" --config "$ROOT/configs/classical_v2.json" --out "$REPLAY/v2" --workers 2
$PY "$ROOT/src/classical_score.py" --sealed "$REPLAY/v2" --out "$REPLAY/score_v2" --figures "$REPLAY/figures"
~~~
