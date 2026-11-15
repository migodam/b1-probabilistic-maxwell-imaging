# Shared L-BFGS-B60：固定四对象机制控制离线评分

父线程已验收source/solver/cost并另行授权attempt02 seal后的GT评分。所有saved-view、raw min(initial,proposal)、secondaryfallback和paidCPU-prefix选择计划都在读取GT前按公开Phi/eligibility写盘。原seal/预测保持不变；没有forward/optimizer/新参数搜索。

本组为0002/0017/0032/0037四个已曝光对象×四metric。Material统一原chi坐标、同current-m起点、同box/prior、L-BFGS-B60/ftol1e-12/gtol1e-8/maxls25、完整metric导数、无scalarseed/Arc/MM，旧m非增保障。它控制materialtarget与有限innerprogram同时改变的混杂；不是新主算法、全48恢复结果或独立确认。

## 归一化与选择

复用原R13 truth_spec/load_truth/material_metrics。Amplitude NRMSE=全域RMSE/(max(coarseGT)-.1)，另报relativeL2与boxrange_NRMSE；固定prediction chi>1与truth chi>1支持。Phi=.5全Frobenius相对residual平方+原prior；GT不选择candidate、阶段、prefix或threshold。

保存view共80行（含16initial），其中raw/corrected/primary/secondary64行；paidprefix15/30/60/90/120共80行。Raw按公开Phi回退initial的条数为0。缺失corrected保留missing；缺失secondary明确fallbackprimary。

## 每对象Primary完整四组合

### r13_fresh_0002

| pair | NRMSE | IoU | Phi |
|---|---:|---:|---:|
| 00 | 0.5736061 | 0.382263 | 0.000659042552 |
| 01 | 0.4860269 | 0.477966 | 0.000825366823 |
| 10 | 0.6305576 | 0.295890 | 0.000495086414 |
| 11 | 0.5768162 | 0.366667 | 0.000652316988 |

### r13_fresh_0017

| pair | NRMSE | IoU | Phi |
|---|---:|---:|---:|
| 00 | 0.4684052 | 0.505085 | 0.0003376735 |
| 01 | 0.0695107 | 0.878906 | 0.000182464901 |
| 10 | 0.4450119 | 0.534884 | 0.000320111723 |
| 11 | 0.4754782 | 0.493289 | 0.000341164577 |

### r13_fresh_0032

| pair | NRMSE | IoU | Phi |
|---|---:|---:|---:|
| 00 | 0.0602181 | 0.881789 | 0.000194084156 |
| 01 | 0.0595544 | 0.878594 | 0.000194083459 |
| 10 | 0.0587976 | 0.881410 | 0.000194085924 |
| 11 | 0.0593286 | 0.881410 | 0.000194085121 |

### r13_fresh_0037

| pair | NRMSE | IoU | Phi |
|---|---:|---:|---:|
| 00 | 0.6931308 | 0.283286 | 0.000711697992 |
| 01 | 0.6258711 | 0.336478 | 0.00140142014 |
| 10 | 0.7034746 | 0.261838 | 0.00125279858 |
| 11 | 0.5997602 | 0.380368 | 0.00133904804 |

## 描述性汇总与材料数组差

| pair | primary mean NRMSE | primary mean IoU |
|---|---:|---:|
| 00 | 0.4488401 | 0.513106 |
| 01 | 0.3102408 | 0.642986 |
| 10 | 0.4594604 | 0.493506 |
| 11 | 0.4278458 | 0.530433 |

MATERIAL_ARRAY_DIFFERENCES保留每case、每raw/corrected/primary/secondary、全部六pair-pair比较的L2、原amplitude归一化RMS、relativeL2、Linf、supportflips及共同Phi。NPZ保存完整材料/GT/固定mask及数组差；相似NRMSE不代表相同材料数组。

Corrector后NRMSE增加的条数为14；所有D/R/Phi和材料变化逐条保留。描述性统计不作最终新颖性判断。

## 同host费用和失败

成功16children完整OSCPU=952.745000s；旧psutil-bootstrap失败CPU=0.514267s，累计=953.259267s。Controller CPU=0.832362s。旧16次setupfail均在物理前，0LU/RHS，完整日志/原失败seal保留。

实际physical LU=43330，complex RHS=886170。Inherited M1 whitening计数的名义offset单列于已完成costreconciliation；CPU/LU/RHS不改。Offline评分/图形CPU单列SCORING_RECEIPT，额外LU/RHS/optimizer/CUDA=0。

## 解释边界

图MATERIALS_RAW与MATERIALS_PRIMARY均为4×(GT+4pairs)，共享绝对chi[0,20]色标，不插值。只有本地同16组内CPU具有共同host边界；不比较Mac与旧Windows绝对速度。

本地算力使本组各臂都完成3000nativeupdates；原Arc对照在旧host/72proposal cap下的迭代能力不同。即使与原Arc结果有差异，也不能把全部变化归因于materialsolver。此次共同innerprogram只提供这四个曝光对象上的机制辨别，主综合与最终新颖性/科学解释由父线程负责。

补充机械验收：公开selector一致性不匹配0；80条可报告view中有11条为已授权且显式标注的secondary→primary回退，实际归档saved-view文件69个。两张raw/primary图均已查看标签、完整数组和共享色标。
