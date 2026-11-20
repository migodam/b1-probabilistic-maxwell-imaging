# R16 使用的历史二维散射数据

本目录读取 MAX-STRUCT 原先保存的带噪二维标量 TM **仿真**观测和真实材料数组。共 130 份不同材料、138 份观测。R11 跨网格实验给同一组八份材料保存了两种观测。各组的逆问题使用逐元素相同的 $G,E_0,S$ 和网格坐标。历史确认组名称保留，现有存档在 R16 属于回顾性评估。

|接口名称|观测数|原始来源和形状|
|---|---:|---|
|`r10_twod`|24|R10 主组。同网格。环、桥、包覆包体、梳齿、弯条、多岛、穿孔星形、开缝壳各含低、高、混合对比。|
|`r10_conditional4`|4|R10 条件策略确认。同网格。环、桥、包覆包体、开缝壳各一案。|
|`r11_independent48`|48|R11 主组。同网格。沿用八个形状类别，含高对比及混合对比。|
|`r11_bothblocks16`|16|R11 两块分母确认。同网格。八个形状类别，每类两案。|
|`r11_crossgrid8`|8|R11 细网格生成，粗网格反演。环、开缝壳、桥、多岛各两案。评分材料是 $48^2$ 网格的 $2\times2$ 体积均值。|
|`r11_crossgrid8_matched`|8|R11 配对对照。与上一行共享八份材料，使用 $24^2$ 体积均值材料和 $24^2$ 算子生成观测。|
|`r12_dev6`|6|R12 当轮开发材料。同网格。包含枝状、连续随机、水平集和混合岛。|
|`r12_confirm24`|24|R12 当轮确认材料。同网格。四类形状各六案。|

R12 历史盆地逃逸试验重用了 R10 的 `0001`、`0004`、`0022`，对象仍归 R10。表中各独立材料组的 $24^2$ 真值数组没有逐元素重复。形状和对比标签来自原始生成记录，供离线分析使用。[逐对象标签](results/baseline_case_labels.csv)

## 正演和输入

令 $D=\operatorname{diag}(m)$，$m$ 是实对比度。原正演以电流为未知量。

```text
(I − D G) J = D E0
E = E0 + G J
Y_pred = S J
```

等价的总场形式为 $(I-GD)E=E_0$，$Y_\mathrm{pred}=SDE$。$G$ 是 $576\times576$，$E_0$ 是 $576\times12$，$S$ 是 $48\times576$，$Y$ 是 $48\times12$，均为复双精度。数组没有转置或共轭。$m$ 按 C-order 排列，连续元素沿 $x$ 方向，每行后沿 $y$ 增加，单元中心间距 0.04。$k=2\pi$，脉冲盘半径 0.018。材料盒为 $[0.02,20]$，已知宿主 $\chi=0.1$。[原正演实现](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_z_alg_r4_diversity/src/physics.py)

带噪读数为 $Y=Y_\mathrm{clean}+N$，每条通道的 $\mathbb E|N_{rs}|^2=\sigma^2$。存档的 `noise_radius` 是按 $\sigma$ 标准化后的 CN99 半径，约为 25.1654。求解端读取 `Y` 和 `sigma`，以及公共算子、网格和材料盒。形状标签、生成配置、干净观测和材料真值用于离线工作。

```python
from cohorts import load_cohort

model, cases = load_cohort("r11_crossgrid8")
# model 包含 G、E0、S、shape、coords、bounds、sources、receivers、k、radius
# cases 的每项包含 id、Y、sigma
```

## 原始文件位置

当前工作区直接读取 `research/delegated/max_struct_r1/r{10,11,12}_worktree`。外部复现时可用一个固定版本的公开历史仓库。将 `R16_MAX_STRUCT_ROOT` 指向仓库根目录或存放三个历史工作树的目录。[固定版本](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research)

```sh
git clone --filter=blob:none --no-checkout https://github.com/migodam/b1-probabilistic-maxwell-imaging.git max-struct-history
git -C max-struct-history sparse-checkout set --no-cone \
  '/research/research_max_struct_r10/data_public/twod/' \
  '/research/research_max_struct_r10/data_public/conditional_confirmation/' \
  '/research/research_max_struct_r10/sealed_truth/twod/' \
  '/research/research_max_struct_r10/sealed_truth/conditional_confirmation/' \
  '/research/research_max_struct_r11/data_public/' \
  '/research/research_max_struct_r11/sealed_truth/' \
  '/research/research_max_struct_r12/delegated/confirmation/data_public/' \
  '/research/research_max_struct_r12/delegated/confirmation/sealed_truth/'
git -C max-struct-history checkout --detach ec3a37da6f3ee90beb61e753067d4f639f4118f5
export R16_MAX_STRUCT_ROOT="$(cd max-struct-history && pwd)"
```

## MAX-STRUCT R3 弯条压力案例

另保留一个有独立观测几何的历史例子。MAX-STRUCT R3 重用了 R2 存档 `obs_b20d5a55677d4eca`。它不计入上表的 130 份材料，不能拿上表公共 $S,E_0$ 跑该例。`load_r3_curved_bar()` 读取它的原算子与带噪读数。`load_r3_curved_bar_truth()` 用于离线评估，读回材料和原支撑掩模。原 CIE 终点保存在 R2 的 `predictions/obs_b20d5a55677d4eca/native.npz`。旧 NRMSE 和 Otsu IoU 定义及失败数值见附录。[历史 CIE 更新与失败案例](HISTORICAL_CIE.md)

使用精简检出时，可将原观测、真值及旧端点添加到同一历史检出目录。

```sh
git -C max-struct-history sparse-checkout add \
  '/research/research_max_struct_r2/data_public/obs_b20d5a55677d4eca.npz' \
  '/research/research_max_struct_r2/sealed_truth/obs_b20d5a55677d4eca.npz' \
  '/research/research_max_struct_r2/predictions/obs_b20d5a55677d4eca/native.npz'
```

## 离线评分和历史基线

保存预测后，用 `load_truth(cohort, case_id)` 取材料，再调用 `score_materials.material_metrics`。幅度 NRMSE 的分子为 $\lVert m-m_\mathrm{true}\rVert_2/\sqrt{576}$，分母为 $\max(m_\mathrm{true})-0.1$。R11 跨网格原评分固定用 11.9，等于八案的这个分母。固定支撑 IoU 使用两幅图各自的严格超水平集 $m>1$。历史恢复合格需要材料在盒内、幅度 NRMSE ≤0.15、IoU ≥0.7。相对 $L^2$ 的分母为 $\lVert m_\mathrm{true}\rVert_2$，单独列示。

R16 的 `SSIM_fixed_range` 是新增评分。使用完整 $7\times7$ 窗口、样本协方差、$K_1=0.01$、$K_2=0.03$ 和统一值域 19.98。历史主基线表没有 SSIM 分数。

[逐对象历史基线](results/baselines.csv) 汇集 864 条原先记录的结果。每条保留原方法名、历史候选视图、原表位置和已报告费用范围。R10 主组含 CIE、Dual-CIE、A+ 全局与逐源、G19 等七臂。R10 条件组记录 CIE 及两条回退策略的部署终点。R11 主组、跨网格组和两块分母组用各自终点表，含 CIE、Dual-CIE、Global Dual、G19 及原组实际报告的控制。R12 开发和确认组使用历史 primary endpoint，含 DR、A+、CIE 和其他实际执行的基线。附加的 R12 盆地逃逸结果为三份 R10 旧观测的四臂续跑，从公开旧坏终点出发，沿用独立的历史视图标记。

表内的相对数据残差为 $\lVert Y_\mathrm{pred}-Y\rVert_F/\lVert Y\rVert_F$。R11 跨网格组和两块分母组由原表的 $\lVert Y_\mathrm{pred}-Y\rVert_F/\sigma$ 与已存读数换算。空缺表示原输出没有对应指标。历史 CPU 口径按 `CPU_scope` 解释，费用来自当时的软硬件和预算。已存基线用于材料效果的对照，R16 的同机计时另行统计。
