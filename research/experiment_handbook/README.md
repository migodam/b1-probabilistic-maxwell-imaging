# Maxwell 逆散射实验手册

**下一轮比较 A+、G19、DR 或 CIE，先确定要解释哪一组材料上的差别，再取那一组已保存的观测。** 本手册登记了 56 个数据集、子集及派生视图和 43 个可调用的算子入口，逐项列出生成者、物理合同、训练或曝光史、观测文件和评分入口。另有 28 个可调用的算子入口，以及 15 条可追到原论文、代码和逐例结果的基线勘误。登记更新于 2026-11-22。

[有损材料与有限孔径新比较](../z_alg_optical_inertia/ip_unified/)

## 按研究问题取数据

|想回答的问题|先取哪组已有证据|可以做的比较|
|---|---|---|
|A+ 的材料反馈能否救回 CIE 失败的材料？| `r13_fresh48` 的同一 48 个材料；再查 R14 独立 24 件|同对象、相同付费时间前缀、相同选根；比较多救回和丢掉的对象。R13 这批已评分|
|DR 能否离开已有坏终点，NN 的作用能否由 AP 解释？| `r12_paired_restart` 的原始材料及 A/B/C 组父运行|同端点或同 NN 初值，完整物理双场接 AP/DR，G19 从合法辅助状态重启|
|半环优势来自恢复目标还是抑制背景？换阵列是否保持？| `srap_main24` 的七个阵列；`srap_followup8` 的细桥和连续介质|同一材料配对，全图、内部、背景和拓扑分别判读。三维另取 `srap_vector12` 平面阵和半球阵|
|能否复现已有论文，能否预测真实测量？| `r8_austria`、`r8_paper_vector`、`fresnel_2005` 及冻结训练/留出掩码|文献仿真查尺寸、频率、源数与停机；实测查复场校准，按未见源、未见接收分别评分|
|训练一个低成本初值是否划算？| `r8_nn_train`、`r8_nn_validation`、`r8_warm42` 与冻结权重|按材料家族和物理条件分层。把推理、状态转换、后续全波更新与训练摊销分别计费|
|能否不重建整图，直接排除介电对比度全在某一上限内？| `z_alg_secant_sector_inertia` 的五份原有二维合成观测|同上限、同乘子网格比较正定与不定残差下界。R10 旧环的 576 个像素若均不超过 5，严格误差下界约 0.014695，超过原噪声半径 0.010517。[证书与对照](../z_alg_secant_tasks/multiple-illuminations.html)|
|光学定理能否从有限孔径读数约束材料峰值？|[六件新强散射材料，24 份配对观测](../z_alg_optical_inertia/synthetic/public/manifest.json)|同一材料比较全环、半环及 2% / 10% 噪声。先复核 96² 有限模型数据球，再比较实部、光学定理和空间乘子。[实验与复算](../z_alg_optical_inertia/experiments.html)|
|全环噪声变动时还能排除什么峰值？| `z_alg_optical_fullring_paired30` 的 6 件材料和 30 份配对复观测|逐件比对标量、空间乘子、噪声球和高阈值限制。[配对曲线](../z_alg_optical_inertia/submission/noise_curve.html)|
|独立实验室的材料峰值结论如何受校准约束？|`manitoba_eigor2010_dataset1` 的 4 件实物与 28 份频率观测|[校准误差分解与复算](../z_alg_optical_inertia/validation/manitoba.html)|
|能否用较少发射排除镜像材料？| `z_alg_secant_reflection` 中按真实端口闭合的已存带噪视图|预先给定误差上限，按镜像发射组读前缀。152 份合规 SRAP 视图中，同样判定 150 份，逐组 428 次、全批 1,517 次、事后同布局固定前缀 926 次。[中文报告](../z_alg_secant_tasks/index.html)|

[数据集与覆盖缺口](DATASETS_ZH.md) · [其他研究系列及其用途](PROJECT_FAMILIES_ZH.md) · [算子与设备](OPERATORS_ZH.md)

## 三个可以立即查清的问题

在项目根或 B1 完整检出目录运行。查询只定位已有文件。只有网页与代码时，查询会说明取得相应冻结数据的位置。

```bash
python3 research/experiment_handbook/lookup.py r13_fresh48 r13_fresh_0017
python3 research/experiment_handbook/lookup.py srap_vector12 --list-views
python3 research/experiment_handbook/lookup.py --operator measured_ports
```

同一数字尾号不足以标识数据。`r11_independent48` 与 `r13_fresh48` 是两批不同材料。R13 的同批 48 件，在历史联合标准下 CIE 为 30、Dual-CIE 为 32、A+-MM 和 A+-Arc 各为 36、G19 为 30、DR 为 29。R13 中的 12 件跨网格对象全用全环接收，16 件半环全部是同网格。因此这批数据支持 A+ 的配对收益，也给“半环跨网格下仍然有效”留下了明确缺口。[基线与错题本](BASELINES_AND_ERRATA_ZH.md) · [按实际物理条件统计的设计单元](results/design_cells.csv)

SRAP 的 24 件标准半环在 30 秒 raw 口径下有 21 胜、1 平、2 负，配对全图材料相对 L2 中位降幅为 39.1%。其中 18 个胜例的目标内部误差反而变大。逐对象原图显示，轻校正能让 `s2d_010` 出现两个孔，而 `s2d_007` 在同一半环上留下严重高估的粗亮环。选择指标会影响人对这项正结果的理解。[封存预测复算与原图入口](WALKTHROUGHS_ZH.md)

## 从开题到交接

- [实验设计、历史指令与中断恢复](DESIGN_AND_CONTEXT_ZH.md)说明怎样冻结比较、验收子任务、追究失败并续跑。
- [基线与错题本](BASELINES_AND_ERRATA_ZH.md)按“误差如何影响主张”解释 15 个实际问题，附原始来源与可复算结果。
- [四条历史路径和一次配对复算](WALKTHROUGHS_ZH.md)实际走通坏终点、Austria、Fresnel 和三维端口，并从材料数组复算半环收益。
- [磁盘的保存、归档和回收](STORAGE_ZH.md)给出受保护证据、运行暂存和历史副本各自的处理办法。本地历史盘点从 `communication/EXPERIMENT_STORAGE_INDEX.md` 进入。
- [给下一轮 Chat 的入口](CHAT_HANDOFF.md)列出当前可引用判断及继续研究所需文件。

可以在[网页目录](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/experiment_handbook/catalog.html)按名称、来源与物理类别查找。机器可读的登记表为[数据族和物理合同](COHORTS.csv)、[观测与材料对应表](OBSERVATIONS.csv)和[算子入口表](OPERATORS.csv)。新增数据、算子或改变论文解释的勘误与该轮结果一同登记。
