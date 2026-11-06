# 原3λ四立方体：CM九臂实际结果

**PUBLIC_DEVELOPMENT。cross-grid五臂与matched finite-grid四臂全部完成、封存后统一评分。18个阶段均达到冻结Eq16系数停止，无CAP、无数值失败；没有把停止条件当成原物理驻点或材料恢复证书。**

本例两组数据上，R9 vector A+的两项材料误差及原物理Φ均略高于Global Dual，wall也更长。CIE系保留四块轮廓但明显模糊；LSIE缺失较高contrast块。matched数据降低了CIE系的物理失配，却仍远高于已知GT可行点。以下只给该固定场景/噪声的实测边界，不作统计推广或作者软件/连续物理复现声明。

## 合同与尺度

原示例与Eq16/17参照[Zhong与Xu，2019原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)。这里采用独立CM+RR DDA实现及明确披露的DPCS/约束适配；积分物理资格差异已经单列，不冒称作者程序的数值复现。

原3λ立方域、60电偶极源、同位置60接收点×3分量、four cubes ε=2/3/4/5。反演网格30³，**27000实材料未知**；180×60复观测至多提供**21600实数据分量**。这是信息维数上界，不能宣称全局唯一；GT含box边界时，也不能仅凭计数断言该特定GT点多解。small3D的512未知与同数量端口的数据/未知量比例不同，不能仅按“三维”合并结论。

所有九臂采用冻结CM+RR点DDA、实χ∈[0,4]、air χ0=0/APCS γ0=0、数据驱动DPCS Jd、R=0、ρ=1、MF6→10、最多500/阶段、原Eq16根号RMS<.001、无额外物理corrector。八角低频mask对应每源5184/24000个复系数，原生正交归一IFFT坐标，不做SVD whitening。DPCS固定.1σmax得到rank162；此L规则是明确适配，不冒称作者唯一规则。

原尺寸物理Φ = mean_s[||prediction_s−Y_s||²/||Y_s||²]，**没有1/2且R=0**；不能数值横比PhaseI Φ。Auxiliary data/state是当前系数模型的两项，并非完整Maxwell Φ；Dual/A+不同current/material β时也没有共同辅助目标。主终点始终为stage2，stage1只保留诊断，不按GT挑选。

cross Y逐元素保持冻结CM60³→CM30³数据；matched复用已付CM30³ clean响应，没有新增forward generation。两个cohort不同data_id、分表报告，matched是有限离散匹配控制，不是同数据算法收益。CN10%为**期望**相对范数，实际两组均.0997754141612836；matched复用同seed2026110619及噪声方向，仅按clean范数缩放sigma，未强制样本norm=10%。sigma分别.0172764795207761/.0173157557648113。

data prep打开sealed容器，但仅提取Y_clean/noise/已付coarse prediction，没有读取fine_material/inverse_material。反演接口仅Y、公共算子、box/grid；材料GT首次由全部九seal之后的独立评分程序读取。

## 主终点材料与物理数据指标

Air参照：paper Err=.1896620343、χ-relativeL2=1、Φ=1。大背景占比能让论文平均ε误差显得较小，必须并列χ-relativeL2及图像。

### Cross-grid：CM60³数据/CM30³反演

已知同一GT材料的可行物理Φ=0.294687185（由已付coarse响应计算，无新正演；不是最优值或唯一性证明）。

| 方法 | paper Err | χ-relativeL2 | 原物理Φ | 阶段步数 |
|---|---:|---:|---:|---:|
| Source CIE | 0.216272 | 0.414601 | 0.228611 | 113/161 |
| Source Dual | 0.210070 | 0.410882 | 0.237815 | 113/133 |
| Global Dual | 0.209837 | 0.410724 | 0.237800 | 138/139 |
| R9 vector A+ | 0.211773 | 0.411815 | 0.238453 | 134/141 |
| Source LSIE | 0.179943 | 0.910057 | 0.454097 | 195/165 |

| 方法 | final aux data | final aux state | Eq16 RMS |
|---|---:|---:|---:|
| Source CIE | 0.00027423 | 0.00099829 | 0.00098774 |
| Source Dual | 0.00074181 | 0.00275300 | 0.00096646 |
| Global Dual | 0.00074172 | 0.00275334 | 0.00098624 |
| R9 vector A+ | 0.00074099 | 0.00275285 | 0.00098132 |
| Source LSIE | 0.00032615 | 0.00353912 | 0.00090069 |

### Matched finite control：CM30³数据/CM30³反演

已知同一GT材料的可行物理Φ=0.010203349（由已付coarse响应计算，无新正演；不是最优值或唯一性证明）。

| 方法 | paper Err | χ-relativeL2 | 原物理Φ | 阶段步数 |
|---|---:|---:|---:|---:|
| Source CIE | 0.175899 | 0.427272 | 0.134780 | 101/143 |
| Global Dual | 0.168225 | 0.422750 | 0.141197 | 119/124 |
| R9 vector A+ | 0.169388 | 0.423472 | 0.142384 | 110/128 |
| Source LSIE | 0.176669 | 0.929532 | 0.494142 | 171/145 |

| 方法 | final aux data | final aux state | Eq16 RMS |
|---|---:|---:|---:|
| Source CIE | 0.00028101 | 0.00083664 | 0.00098510 |
| Global Dual | 0.00075233 | 0.00249976 | 0.00099555 |
| R9 vector A+ | 0.00075181 | 0.00250014 | 0.00098838 |
| Source LSIE | 0.00032647 | 0.00375632 | 0.00096198 |

matched仍有明显的auxiliary/full-physics间隙，因此纯跨网格失配不能解释全部现象；现有结果尚不能区分全部representation、约束与优化路径贡献。CM coarse真值forward没有保存完整J，只有`sealed_truth/coarse_truth_response.npz`的prediction(180,60)；没有为找J重新启动正演。

父线程另用已付norm作零新solve的[残差传播诊断](../../RESULTS/original_residual_transport/DIAGNOSTIC.json)，见[THEORY §9](../../THEORY_DERIVATIONS.md)。既有恒等式`r_phys=d−S L(m)^−1 e`及三角不等式给当前残差方向、日志ordinary-current `p e`规范下的gain区间：matched CIE约12.11–13.27、Global约6.97–8.06、A+约7.00–8.09。它解释aux小而端口物理残差大的现象，Global/A+该项未明显改善。区间上端不是operator-norm上界，没有严格线性求解/离散化误差包络，不是新定理或严格误差证书，也不按不同β的gain排序。该metadata分析由root在union另计，不吞入原九臂总账。

独立只读接线审计未发现baseline削弱或A-current混淆。Source CIE/Source LSIE全部八个阶段的已付辅助目标最大增加为0；Source CIE与Source Dual首阶段材料差=0，packed16与旧strided4已付Source CIE首阶段材料差=0；27个端点从stored prediction重算Φ与记录差=0。2473个新辅助步，1020个实际非零完整RHS，最大原线性状态残差1.73742631e-06。这些是接线/状态求解证据，不是全局凸性、GT盆或唯一性证书。

## 物理诊断边界

冻结CM跨网格GT响应差global L2=.535254、mean-source RMS=.535244，高于噪声。单独预声明source59积分诊断：CM30/60差.508374→IGT.359841；同grid IGT/CM变化30=.276261、60=.072778。IGT两RHS代回原CPU矩阵的残差1.56813e−6/1.88546e−6，CPU/GPU接收预测差<5e−16。采用固定near8/far3/self16、同一source-cell S与self/kernel、中心incident、同一LUT local Ablock。它有局部改善，仍未达到10%噪声尺度；不生成IGT60源包、不扫网格、不做IGT反演。

R9 vector A+保留两次2³块坐标扫掠及按完整固定J的A损失接受普通CIE种子提案，属于已披露R9变体；ordinary global current更新与Global Dual相同，不把readout变化暗换为另一current方程。

## 实际费用与复用归属

下表CPU是host资源时间，**不含GPU算术时间**。A+的CPU较低但wall更高，不能用CPU列宣称更快。每臂actual为真实共享queue interval；归属估计在声明setup摊销后补入cache命中的首次evaluation CPU/wall，不是新独立冷启动计时，也不再次计费。

| 数据/方法 | actual CPU/wall秒 | 含setup和复用归属估计 CPU/wall秒 | 实际完整RHS |
|---|---:|---:|---:|
| cross_grid/Source CIE | 336.781/338.594 | 337.447/339.493 | 120 |
| cross_grid/Source Dual | 295.672/297.914 | 313.884/316.516 | 60 |
| cross_grid/Global Dual | 345.016/347.080 | 345.681/347.979 | 120 |
| cross_grid/R9 vector A+ | 149.922/445.483 | 150.588/446.382 | 120 |
| cross_grid/Source LSIE | 382.609/385.055 | 383.275/385.954 | 120 |
| matched_grid/Source CIE | 294.828/296.803 | 295.473/297.679 | 120 |
| matched_grid/Global Dual | 296.156/298.346 | 296.801/299.223 | 120 |
| matched_grid/R9 vector A+ | 125.062/380.331 | 125.707/381.208 | 120 |
| matched_grid/Source LSIE | 338.078/340.164 | 338.723/341.040 | 120 |

九臂whole-process=2570.593750CPU/3138.316540wall秒，只计一次；1次SGram/EVD/Vr、2次Y-fit、27物理value calls、9cache hits、1次air解析零预测构造（air实际RHS=0）。峰值GPU allocated/reserved=2.693/3.077GB，host RSS=1.486GB。

历史工程、IGT两个进程、九臂及两次离线绘图/评分共6221.913316个已记录CPU秒：Mac 2709.225816，XINAN 3512.687500。16个additive行明确分离，不能将各进程wall相加当makespan；也不能把混合平台CPU和作为单平台速度。旧Mac475.441/2224.278CPU秒、旧strided4部分616CPU秒、所有基准/失败均保留。首次CPU轴图也留存并计费，主成本图用wall。

已知计数缺口：第一次uncropped中断源的GMRES callback数未知，但其全部Green/CPU保留；旧strided4部分没有导出末步完整Green counter，230完整辅助步与120物理RHS已观察，可能的in-flight剩余计数未知，不填0。补入缓存的归属估计及shared/每臂子计时均为nonadditive。

整体CPU/GPU source59比16.96倍含SciPy-MGS与resident-CGS2/检查节奏差异，不能称纯硬件提速；一次同状态完整current/material步骤的packed16收益2.915倍，输出差全部0。没有为包装层额外重算整套数值实验。

## 交付与复现

- [完整原始评分与aux/停止量](../../RESULTS/original_vector3d/NINE_OFFLINE_SCORES.json)、[主终点CSV](../../RESULTS/original_vector3d/NINE_PRIMARY_SCORES.csv)。
- [最终去重账](../../RESULTS/original_vector3d/ACCOUNTING_FINAL.json)：additive进程行与缓存归属视图分开。
- [跨网格材料PNG](../../FIGURES/original_CM_cross_grid_materials.png)、[PDF](../../FIGURES/original_CM_cross_grid_materials.pdf)；[matched材料PNG](../../FIGURES/original_CM_matched_grid_materials.png)、[PDF](../../FIGURES/original_CM_matched_grid_materials.pdf)。
- [wall成本/材料/物理/系数轨迹PNG](../../FIGURES/original_CM_nine_convergence.png)、[PDF](../../FIGURES/original_CM_nine_convergence.pdf)。固定共同χ0–4色标、z≈±.7切片与max-z投影；已实际查看渲染，标题/标签清晰。
- [冻结增补协议](CM_NINE_PROTOCOL_AMENDMENT_ZH.md)、[实际复现命令与环境](REPRODUCE_ORIGINAL_CM.md)。完整repo+已有PyTorch CUDA；薄路径入口要求显式新run-dir，默认不继承历史绝对目录、不覆盖封存输出。
- 原9冻结不改。父线程另授权matched air/Global/A+三起点各8accepted原物理L-BFGS-B控制，作为独立DEV补充另账另报告；不扩原CIE预算、不替换原stage2，也不吞入本ACCOUNTING_FINAL。
