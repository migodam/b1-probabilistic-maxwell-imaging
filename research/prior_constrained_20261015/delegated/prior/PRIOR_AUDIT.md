# 原 A25/A30 材料 prior：E0 零物理核对

本文件只给出先验来源、确定性重放与无物理材料生成证据；竞争解释与下一阶段科研裁决由主线程作出。

## 已核实的实际材料定义

源代码来自本地 A26 checkout `2603f3bd16c48e48587fedc4363975779def17f5` 的 `src/a25/dataset.py` (`generate_manifest`, `parent_material`, `generate_dataset`)、`src/a24/cases.py` (`backgrounds`, `noisy_observation`) 和 `configs/a25/campaign.json`。A30 `717d3d26a290248bd2b9aeb55f3d03694b407b84` 中公开 A29 输入副本供应实际 TRAIN 与几何缓存。

定义为

`chi = chi_b + a * p / max(p) * v`, `v = (1 + .2i)/sqrt(1.04)`。

`a1` 是下标 1，实际幅度 **1.2**，不是 1。另一个原始幅度是 .3。`rho = Re(conj(v)*(chi-chi_b))` 是非负实增量，不是绝对介电材料相。`rho=0` 回到公开背景；不代表空气。

背景 weak 为 `.1+.04i`；ring 的 `.2<=r<=.35` 内为 `2+.1i`，其他位置为 weak。因此所有 P0 材料的绝对 Re/Im chi 均不低于背景，且不能有 chi=0 的空气孔。原 phase 规则限制 `chi-chi_b`，不是绝对 chi 的复相位。

实际物理几何来自公开 `geometry.json`，24×24 cell-centered points，side=1，native ordering；8 source 与32 receiver 未重排。这里只读几何和材料，不构造 Maxwell 对象、不调用求解器。

## 原 TRAIN 生成测度

四类等权混合：shell、asymmetric_piecewise、smooth、mixed_scale，每类原16个独立 TRAIN parent，总64个。每 parent 两背景×两幅度×两噪声条件，共512 rows；不能当512个独立几何。

几何参数从独立 RNG draw 产生，原顺序完整保留：rotation U[-π,π]，center两个坐标U[-.12,.12]，scale U[.78,1.04]，axis_ratio U[.82,1.18]；shell radius U[.14,.19]、width U[.035,.055]；disk radius U[.065,.10]、weight U[.45,.8]；second-shell weight U[.45,.75]；rectangle width U[.17,.25]、height U[.19,.29]；subtype为等概率0/1；三个Gaussian centers各坐标U[-.13,.13]，widths U[.05,.095]，weights U[.4,.9]；segment offset U[-.06,.06]。TRAIN为ID，未启用 structural OOD 分支。

原 rasterizer 在 transformed x/y 坐标下硬栅格化 shell/rectangle/disk，smooth为两个Gaussian和，mixed为shell、Gaussian、segment的逐元max；最后统一除以实际栅格最大值。不存在自由的全材料复补偿，也没有新增的hole-count标签或生成桶。

## 精确重放与 P0_EXACT 的含义

`src/prior_constrained/prior.py` 只提取原 `parent_material` AST函数并提供 NumPy 环境；没有导入原 dataset模块的物理入口。`draw_train_row` 精确重现 TRAIN RNG 分支，不调用会同时生成 locked geometry 的 `generate_manifest`。

以原 master seed251008、domain250、split-index0、family-index0..3、local0..15，重放64 TRAIN几何。与实际 `train.npz` 512条chi逐元比较，**最大绝对差0.0，512/512精确一致**。见 `TRAIN_REPLAY.csv`。没有打开 locked manifest、locked材料、locked几何或clean arrays。

`phase_audit` 只给 BOX/PHASE必要条件，绝不把 PHASE写作P0_EXACT。`verify_exact_endpoint` 先检验所有原参数support，再用原 rasterizer重建端点，要求逐元误差<=1e-12；只有同时具备合法参数和实际端点重放才标P0_EXACT。未重建成功是RETRACTION_NOT_FOUND/PRIOR_UNKNOWN，不证明不在P0。参数域外反例已被拒绝；见 `SUPPORT_VALIDATION.json`。

## 新 material-only reservoir

已用全新 master seed20261015、原TRAIN domain/branch、每类64个，总**256个独立新材料parent**。每个同时生成weak/ring两种公开背景端点，a1固定1.2，保存`rho/weak/ring/points/parent_id`，完整seed与参数保存于 `MATERIAL_ONLY_PROVENANCE.json`。新seed与TRAIN seed集合无交，未使用development或locked geometry。

256/256样本在两背景均通过固定phase和幅度检查。它们是新增 material-only 样本，**不是新 paired Maxwell training data**。没有forward、label、teacher、likelihood或训练。

## 预冻结结构尺度建议

归一化rho按a1=1.2，excursion阈值固定 `.1/.25/.5/.75/.9`；Gaussian尺度0/1/2 cells。报告连续shape、occupancy、boundary和4-foreground/8-background连通分量/孔洞；组件小于4cells只作为小结构诊断剔除，不把绝对chi幅值识别为空气。

连续shape RMS先移除每尺度的DC。DC独立报告，不能单独过结构novelty。TRAIN64独立parent最近邻的低十分位尺度与±.5 cell插值扰动95%分位取最大，得到：

| 指标 | TRAIN nearest p10 | half-cell jitter p95 | 建议固定下限 |
|---|---:|---:|---:|
| centered multiscale field RMS | .06316233 | .08796801 | .08796801 |
| multithreshold occupancy difference | .01766204 | .02549190 | .02549190 |
| multithreshold boundary difference | .01971065 | .03399884 | .03399884 |

结构还要求至少三条excursion阈值具有平均跨尺度occupancy disagreement≥.05。这只是预冻结尺度；主线程需明确联合判定语义，不能因后续最小预测pair结果调低阈值。全部校准不读y、clean、development truth或MAE。

256个TRAIN图像插值变体只用于栅格稳健性，不属于P0_EXACT候选、不放入reservoir。保守将它们也计为material-only计算数组时，本次共512，小于2048授权上限；真实新generator parent仍为256。

## 成本与可复用缓存

主生成/重放脚本wall .4823905s、processCPU .474753s；参数support复核有独立receipt。新Maxwell state0、RHS0、GPU0、paired labels0。读入缓存不对旧physics重复计费。

实际缓存根：

`Gaussian/A30/a30-maxwell-material-manifold/inputs/a29_release/inputs/a25/`

- `train.npz`：TRAIN-only chi/metadata；本任务没有加载data/current/clean字段的数值。
- `geometry.json`：公开物理网格与acquisition。
- `manifest.json`：公开train/validation角色元数据；没有原locked geometry。
- `background/{weak,ring}.npz`：已付费背景；本任务从公开config重构chi_b，未调用其J或physics。
- A26 N3 `training_manifest.json`：train parent 250001–250064；validation 250065–250080，locked_truth_opened=false；checkpoint路径见原manifest，未反序列化checkpoint。

本分工没有评价任何development材料的MAE，没有调prior以获得多解，也未作全局唯一性/后验质量判断。
