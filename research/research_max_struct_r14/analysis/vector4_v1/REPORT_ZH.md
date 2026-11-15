# Vector4：封存后离线材料评分

24/24原任务完成后，主线程审阅固定公式/目标/费用/预测seal并另行授权GT评分。原seal和预测数组均未修改；本程序没有forward、optimizer、新调参或CUDA。

## 合同和对象数

四组two-metric maps为已有CM+RR vector-DDA的3D几何对应，四组合都用原full analytic material derivative与L-BFGS-B inner12；不是R7 Arc或authorcode复现。强对照为原R13 CIE_sourcewise/Dual_sourcewise。每voxel一个实chi，216材料变量，共享xyz及60源；[0,8]box，原prior和per-source物理目标无1/2。

四个条件包含三个已曝光独立几何；crossgrid8→6 frame重复旧frame，不算第四独立材料。matched/crossgrid分开汇总。strict48CPU corrector allowance和absolute120整进程cap保留；实际converged trajectories远早于很多公共prefix，后续prefix只保持已有paid incumbent，不伪造继续计算。

## 指标与选择

完全复用R13 truth_spec/load_truth/material_metrics：amplitude NRMSE = whole-box RMSE / (max(coarse truth)-.1)，另报boxrange_NRMSE和relativeL2；固定prediction chi>1对truth chi>1支持；support内signed amplitude bias另报。没有新jointsuccess阈值。GT从未选择candidate或prefix。

raw incumbent审计24/24：受latent proposal-only输出问题影响的预测数为0。本次全部raw都已经等于min(initial,proposal)的paid Phi incumbent；未改数组，latent来源与future一行report修复单独记录。

## Primary：每个条件的完整六方法比较

### mixed_crossed_bracket

| method | NRMSE | IoU | Phi | physical residual |
|---|---:|---:|---:|---:|
| pair00 | 0.070940 | 0.675000 | 0.000253870984 | 0.011247 |
| pair01 | 0.070937 | 0.675000 | 0.000253870986 | 0.011247 |
| pair10 | 0.070935 | 0.675000 | 0.000253870982 | 0.011247 |
| pair11 | 0.070895 | 0.769231 | 0.000254055517 | 0.011236 |
| CIE_sourcewise | 0.070895 | 0.750000 | 0.000254055528 | 0.011236 |
| Dual_sourcewise | 0.070941 | 0.675000 | 0.000253870985 | 0.011247 |

### slotted_high_frame

| method | NRMSE | IoU | Phi | physical residual |
|---|---:|---:|---:|---:|
| pair00 | 0.256946 | 0.388350 | 0.000410638322 | 0.015375 |
| pair01 | 0.256922 | 0.388350 | 0.000410638279 | 0.015374 |
| pair10 | 0.256921 | 0.388350 | 0.000410638259 | 0.015373 |
| pair11 | 0.256921 | 0.388350 | 0.000410638242 | 0.015373 |
| CIE_sourcewise | 0.283961 | 0.263636 | 0.000412538695 | 0.015579 |
| Dual_sourcewise | 0.256937 | 0.388350 | 0.00041063829 | 0.015375 |

### slotted_high_frame_crossgrid8

| method | NRMSE | IoU | Phi | physical residual |
|---|---:|---:|---:|---:|
| pair00 | 0.229973 | 0.515152 | 0.000380984442 | 0.015342 |
| pair01 | 0.229982 | 0.515152 | 0.000380984422 | 0.015342 |
| pair10 | 0.229984 | 0.515152 | 0.00038098446 | 0.015342 |
| pair11 | 0.229985 | 0.515152 | 0.000380984427 | 0.015342 |
| CIE_sourcewise | 0.229989 | 0.515152 | 0.00038098446 | 0.015343 |
| Dual_sourcewise | 0.229982 | 0.515152 | 0.000380984451 | 0.015342 |

### unequal_coupled_steps

| method | NRMSE | IoU | Phi | physical residual |
|---|---:|---:|---:|---:|
| pair00 | 0.091832 | 0.777778 | 0.000198935049 | 0.010937 |
| pair01 | 0.091862 | 0.777778 | 0.000198935087 | 0.010938 |
| pair10 | 0.091844 | 0.777778 | 0.000198935056 | 0.010938 |
| pair11 | 0.091841 | 0.777778 | 0.000198935049 | 0.010938 |
| CIE_sourcewise | 0.091837 | 0.777778 | 0.000198935047 | 0.010938 |
| Dual_sourcewise | 0.091858 | 0.777778 | 0.000198935076 | 0.010938 |

## Corrector与负面结果

24条中，原物理corrector后material NRMSE增加的条数为1；D/R/Phi变化与material变化逐条保留在CORRECTOR_DECOMPOSITION。不得把较低Phi直接称为材料恢复提升。各方法全部raw/corrected/primary/secondary和secondaryfallback均有完整行，未按GT挑选阶段。

## 实际费用

24个solver children完整OS CPU合计 297.930588秒；controller CPU 0.060425秒，含sourcefreeze的controller wall 156.623597秒。子receipt累计CPU 296.787076秒是另一个body/serialization边界，不能与OS总数重复相加。

实际physical LU=14331，source RHS=859860，adjoint RHS=854940，合计disjoint solver RHS=1714800。basis/densecurrent/materialHop/setup/failed/output记录在原run COST_LEDGER和本表nestedcounts；不要再次加nested delta。

离线评分/图形CPU、wall、RSS单列于SCORING_RECEIPT；离线physical LU/RHS=0。没有GPU/CUDA；旧Windows/其他物理域时间不并入本地速度结论。

## 图形与范围

每个条件的raw/corrected/primary/secondary均显示GT+全部六方法的完整XY/XZ/YZ正中央cellplane(index3, +.02lambda)，共享绝对chi[0,8]色标与物理坐标；不插值。图像来自真正3D648-current/216-material数组，centralplane并不展示每个voxel，也不代替全volume/continuumreview。

PAID_COST_CURVES仅画按原Phi选择的paid incumbent；CORRECTOR_PHI_AND_MATERIAL_SCATTER分离fit与material；PRIMARY_FIT_VS_MATERIAL_SCATTER保留同一case不同物理解。跨格噪声球不覆盖模型失配；没有新Mie/PDE证明。

最终科学解释、所有图形审阅和publication验收仍由主线程负责。
