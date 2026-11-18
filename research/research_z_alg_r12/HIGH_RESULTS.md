# 高对比度历史端点与 R11 折叠提案

从历史实际保存的预测出发，使用原圆盘脉冲标量 TM 体积分方程和原观测。反演网格统一为 $24\times24$，材料为实对比度 χ 的 576 个逐单元参数，盒为 `[0.02,20]`。每例照明与接收均来自保存的物理算子和坐标，分别为 6 或 12 个照明及 24 或 48 个接收点。历史形状标签与端点已公开，因此这批是开发诊断。新优化与提案只读取保存的端点、物理算子与观测。新材料保存后，另行读取离线真值打分。

误差列的 $r/\delta$ 是原观测的全复 Frobenius 残差除以原 CN99 半径。原噪声为独立圆复高斯，每个复样本总方差为 $\sigma^2$。全部重新计算的物理目标为 $\Phi=\tfrac12\|Y(m)-Y\|_F^2/\|Y\|_F^2+R(m)$。$R$ 是相邻格点差除以 19.5 后的平均平方和平均平滑 TV，二者权重均为 $10^{-4}$，TV 平滑尺度为 0.01。NRMSE 用全域 RMSE 除以真值峰值减 0.1，R10 三例使用其归档的同口径幅度。IoU 固定比较大于 1 的对比度支持，R10 三例使用原归档阈值支持。表中历史参考按原物理 $\Phi$ 最低选择，材料真值没有参与挑选。

## 11 个端点

| 对象与方法 | 历史 ID | 初始 $r/\delta$ | 初始 NRMSE / IoU | 同起点 Full-J 后 $r/\delta$ | 历史参考 $r/\delta$ / NRMSE / IoU |
| --- | --- | ---: | ---: | ---: | ---: |
| 环，CIE-single | r10_twod_0001 | 13.044 | .423 / .353 | 13.044 | G19 .915 / .099 / .756 |
| 壳，矩阵 G19 | r10_twod_0004 | 56.255 | .515 / .087 | 56.255 | CIE .924 / .111 / .759 |
| 弯桥，矩阵 G19 | r10_twod_0022 | 16.062 | .337 / .330 | 16.054 | CIE .929 / .097 / .718 |
| 高对比分离块，pair00 | r14_new_0007 | 7.032 | .305 / .186 | 7.032 | CIE .897 / .065 / .628 |
| 高对比岛，G19 | r14_new_0019 | .886 | .051 / .788 | .885 | DR .885 / .060 / .635 |
| 高对比开孔体，DR | r14_new_0023 | 12.743 | .531 / .138 | 12.743 | CIE 2.656 / .569 / .439 |
| 高对比片区，DR | r14_new_0021 | 19.344 | .516 / .251 | 19.339 | CIE 2.055 / .515 / .485 |
| 高对比桥网，A+-Arc | r14_new_0020 | 10.493 | .454 / .213 | 10.471 | DR .901 / .069 / .783 |
| 高对比连续场，DR | r13_fresh_0031 | 61.817 | .386 / .039 | 61.803 | G19 .884 / .027 / .937 |
| 高对比连续场，DR | r13_fresh_0032 | 2.920 | .564 / .345 | 2.917 | CIE .912 / .047 / .885 |
| 高对比连续场，A+-Arc，跨网格 | r14_new_0004 | 4.105 | .521 / .142 | 4.105 | DR 4.105 / .522 / .132 |

R14 跨网格连续场用 36×36 离散生成，24×24 反演的模型误差未包含在 CN99 球内。R13 两个连续场为匹配网格，保留其历史 12 发 24 收的半圆几何。岛样本已在噪声球内且材料较好，属于对照。它继续 Full-J 后 NRMSE 从 .0512 增至 .0580，IoU 从 .7875 降至 .7778，尽管残差从 .8863 降到 .8853 倍噪声半径。

MAX-STRUCT R12 的同坏起点历史 A+、DR、G19 和纯物理续跑另行评分。环的 A+ 留在原坏起点。壳和桥的 DR 归档端点在原 CN99 球内，分别为 .924× 和 .929×，NRMSE/IoU 分别为 .113/.756 和 .097/.706。它们来自原先实际运行。上表 R10 的替代 CIE/G19 来自另一条历史初始化轨迹。完整记录见 [材料与成本表](results/high_scores/scores.csv)。[统一色标材料图](figures/high_materials.pdf) 同时显示实际残差和 IoU。

## 盒面曲率与续跑

环、壳、分离块及开孔体的活动边界乘子符号在保存位置均正确，原目标在当前自由盒面上的 Hessian 最小特征值分别为 (6.686\times10^{-9})、(3.586\times10^{-7})、(8.547\times10^{-9})、(7.721\times10^{-9})。对应投影梯度无穷范数分别为 (9.68\times10^{-9})、(8.91\times10^{-9})、(4.93\times10^{-8})、(9.25\times10^{-9})。壳的可行中心差分检查相对误差为 (9.43\times10^{-10})。这组有限精度观测支持局部盆地解释。环的原归档坐标并未精确满足一阶条件，原自由面 Newton 校正会改变活动约束。

弯桥的盒面 Hessian 有 24 个负特征值，8 个活动乘子符号违例。高对比桥网、片区及两个匹配网格连续场也有负盒面曲率，说明这些轨迹截停位置仍有优化空间。对四例延长 L-BFGS-B，各从原端点最多评估 1800 次 (\Phi)，弯桥、片区、桥网、连续场 0032 的 $r/\delta$ 分别降至 14.877、19.232、9.263、2.874。对应 NRMSE 变为 .402、.517、.460、.575，IoU 为 .168、.238、.213、.317。残差下降没有带来材料改善，四例仍在噪声球外。纯数据长续跑的残差比为 15.180、19.246、10.072、2.874，也保存在结果中。详细特征值诊断见 [盒面 Hessian](results/high_hessian.json)。

## R11 二次提案与对照

在当前自由盒面上用原观测实数 Jacobian 求两条最弱奇异方向，另取一条奇异值最接近 $\delta/\sqrt{576}$ 的方向。记 $Jv=\sigma_v\ell$ 与 $q=D^2F(m)[v,v]$，用同一保留 LU 计算二阶方向，需 2 组照明 RHS。测试 R11 同响应跳步 $\tau=-2\sigma_v/(\ell^Tq)$，并测试带初始残差的二次标量方程实根。强像截断为奇异值大于 $\max(10\sigma_v,10^{-7}\sigma_1)$ 的模式。必要时按原 R11 式补偿强分量，只在材料满足原盒约束时使用补偿点。预先限定未补偿的方向步 RMS 在 0.01 到 2 之间，以测试有限可行跳跃。

33 个方向的法向曲率相对值 $\|P_\perp q\|/\|q\|$ 中位数为 .00177。法向分量与所选左奇异方向的绝对余弦中位数为 .138，最大 .726。去除该一维方向后，剩余法向分量占比中位数 .990。这些值同时反映法向信号通常很小和单标量消去覆盖有限。六个可作双侧差分的方向中，二阶实现检查的最大相对误差为 $2.30\times10^{-9}$。全部同响应非零候选都超出可行步长或低于有限跳跃下限。唯一合格的是壳的残差平移根，跳步 RMS .0253，补偿后的预测和实际旧数据残差均约 56.255×。完整 Full-J/TRF 校正后的最低 $\Phi$ 仍高于旧端点，故未替换旧端点。

在相同的 4,800 物理复 RHS 和 150 次完整正演上限下，另跑普通弱方向、固定随机序列平滑多起点及从原点直接续跑 Full-J/TRF。弱方向使用同尺度两符号启动，折叠无合法有限提案时，用合法射线的一半长度，射线太窄时退回原点。普通多起点保留匹配位移的扰动，折叠无有限点时使用预设三种尺度并含一个独立均匀附近起点。每臂最多三个启动且按原 $\Phi$ 保留最佳候选。高对比分离块的第三方控制有效，普通多起点将残差从 7.032 降至 1.128 倍，NRMSE 从 .305 降至 .107，IoU 从 .186 升至 .551。已存 CIE 可进一步达到 .897 倍、NRMSE .065、IoU .628。其他 10 端点未产生类似改善，0 个端点从 R11 提案得到更好的候选。这组结果限制当前 3 个弱模、原旧照明和盒面迁移的适用范围。

## 实际成本与复现

每次完整正演使用一个 576 阶复 LU，RHS 另计照明、伴随、接收转置和二阶方向求解。11 对象的短续跑连同各自的首尾诊断共 6,046 次完整正演、124,254 个复 RHS、被计时数值段 wall 之和 56.42 秒。四例长续跑含两个目标与诊断，共 14,400 次完整正演、345,888 个复 RHS，被计时数值段 wall 之和 98.99 秒。11 对象的四臂折叠比较及实际首尾诊断、共享设置合计 4,621 次完整正演、212,460 个复 RHS、被计时数值段 wall 之和 234.21 秒。RHS/正演上限针对搜索，首尾诊断和二阶差分单独计数。新实验 wall 为数值段计时之和，文件输入、写出和进程启动在计时段外，原档 wall 沿用其回执口径。例，高对比离散块的多起点使用 94 次正演、4,788 个复 RHS，wall 6.51 秒。原 CIE 档案是 1,369 次正演、16,284 个复 RHS、wall 122.83 秒。原档由其工作机测量，当前 wall 随并发与硬件变化，按完整运行范围比较物理计数。

输入压缩包保留本轮实际使用的物理算子、档案端点、离线评分真值及已报告费用对应的历史回执。历史目录、NPZ 键和来自 ZIP 的原始成员可从每例 JSON 和评分表追溯。原始归档可查 [MAX-STRUCT R12 配对附件](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r12-paired-20261112)、[MAX-STRUCT R13 全结果及输入](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r13-20261112)、[MAX-STRUCT R14 结果与离线真值](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r14-dualmetric-20261115)。

在本目录，用带 NumPy、SciPy、Matplotlib 的 Python 执行。`high_inputs` 对应解压后的目录，其内直接是 `research/`。历史输入已恢复到完整旧仓库时，可以直接将两个环境变量指向对应根目录。

```bash
unzip results/high_inputs.zip -d /tmp/zalg-r12-high-inputs
export ZALG_HISTORY=/tmp/zalg-r12-high-inputs
export ZALG_R13_HISTORY=/tmp/zalg-r12-high-inputs
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
cases=(ring shell bridge disconnected islands porous patchwork bridge_ap continuous_dr31 continuous_dr32 continuous_crossgrid)
long_cases=(bridge patchwork bridge_ap continuous_dr32)
extras=(ring_ap shell_dr bridge_dr continuous_mixed)
python src/high_endpoints.py "${cases[@]}" "${extras[@]}"
python src/high_endpoints.py "${long_cases[@]}" --limit 1800 --output results/high_long
python src/high_fold.py "${cases[@]}"
python src/high_hessian.py shell ring bridge disconnected porous patchwork bridge_ap continuous_dr31 continuous_dr32 islands continuous_crossgrid
python src/high_score.py "${cases[@]}" "${extras[@]}"
```

[本轮必要输入与历史材料](results/high_inputs.zip) · [离线评分与历史费用](results/high_scores/scores.csv) · [材料对照图](figures/high_materials.pdf)

## Good deflated Gauss–Newton 对照

Bloor Riley、Webb 与 Baker 的 good deflated GN 用原问题的 GN 步 $p$ 计算 $s=\langle\nabla\log\mu,p\rangle$。当 $s>0.01$ 且 $|1-s|\ge10^{-8}$ 时，走 $p/(1-s)$，保留反向或上坡跳步，原目标线搜索只用于普通步。这里在 $m/\sqrt{576}$ 坐标使用移位逆平方 $\mu=1+\mathrm{RMS}(m-m_\star)^{-2}$。消去中心 $m_\star$ 取原始归档坏端点。论文的无约束不收敛保证要求满足其假设，盒面和活动约束会改变迭代，本适配对此不作保证。 [论文第 3.2 算法](https://doi.org/10.1137/24M1695191) · [作者实现](https://github.com/AlbanBloorRiley/DeflatedGaussNewton/blob/615c185f3d61f044eab7cfd27f97f698ff341d8f/Methods/Good_Deflated_GaussNewton.m)

两臂复用同一组三个非真值多起点，使用相同旧观测、盒界与物理正则。平滑 TV 用精确平方残差表示，完整 Jacobian 同时包含数据和正则项。两臂的 $p$ 均来自同一个有界、带 $10^{-10}$ 阻尼的线性化最小二乘子问题。子问题以对角尺度变换后的凸二次式求解，从可行截断牛顿点、Cauchy 点和零步中取最低者初始化。消去步按论文做标量修正，超盒部分逐坐标截回。每臂合计至多 4,800 个物理复 RHS 与 150 次完整正演，包含共同原端点评估。保留整个已评估集合的最低原 $\Phi$，归档端点也参与两臂比较。真值仅在保存所有新材料后评分。

| 对象 | 方法 | 旧数据 $r/\delta$ | NRMSE / IoU | 消去跳步 / 上坡 | 完整正演 / 复 RHS / wall 秒 |
| --- | --- | ---: | ---: | ---: | ---: |
| 环 | 普通 GN | 13.044 | .423 / .353 | 0 / 0 | 121 / 4764 / 2.68 |
| 环 | Good deflated GN | 13.044 | .423 / .353 | 12 / 11 | 109 / 4764 / 2.77 |
| 壳 | 普通 GN | 56.255 | .515 / .087 | 0 / 0 | 150 / 4392 / 2.33 |
| 壳 | Good deflated GN | 56.255 | .515 / .087 | 77 / 44 | 83 / 4788 / 2.70 |
| 分离块 | 普通 GN | 7.032 | .305 / .185 | 0 / 0 | 149 / 4782 / 3.26 |
| 分离块 | Good deflated GN | 4.961 | .547 / .140 | 5 / 3 | 128 / 4800 / 3.34 |

两法在环壳上均保留归档端点。分离块的消去候选距原端点 RMS 9.742，原 $\Phi$ 从 .01132 降至 .00566，仍在噪声球外，材料误差变大。此前同三个起点的有界 Full-J/TRF 多起点对照用 94 次完整正演、4,788 复 RHS，得到 $r/\delta=1.128$、NRMSE .107、IoU .551。较早的未缩放子问题版本中，若干普通 GN 轨迹在迭代子问题还未找到下降方向时提前停住。该版的分离块消去结果曾达到 $1.086$ 倍噪声，但 NRMSE .525。两版都保留，说明此对象的跃迁与有界子问题数值实现有关，分离块的物理目标改善没有转化为材料改善。

除去各臂重复分摊的三次共同原端点评估，正式比较实际完成 737 次完整正演、28,260 个复 RHS、数值 wall 17.05 秒、CPU 16.62 秒。未缩放尝试另计 559 次完整正演、24,012 复 RHS、wall 12.11 秒、CPU 11.81 秒。复现时先生成上节的 `results/high_fold` 三组公共起点，再运行 `python src/high_deflated_gn.py`。 [评分和费用](results/high_deflated_gn.csv) · [完整轨迹与代价](results/high_deflated_gn.json) · [输入、起点及新材料](results/high_deflated_gn.npz) · [未缩放尝试](results/high_deflated_gn_first/high_deflated_gn.csv)
