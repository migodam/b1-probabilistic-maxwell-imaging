# 实际科研迭代与失败归因

本记录保留改变方法的原因。统计目标始终是原 P1；新坐标、初始化、几何与提案不能改变 prior、噪声或 likelihood。各轮费用来自独立作业回执及追加账本，失败不清零。

## 初始实验：先检验完整 posterior，而非训练 teacher

假设是：经 full-physics 反演得到合法基础材料后，固定或更新局部 GN 几何可以支持完整后验的有效移动。32 个独立 development 基础反演全部进入原 noise-ball，27 个达到预注册驻点标准；随后运行四链 MALA、pCN 和变化 decoder mobility。

结果没有支持参考资格。接受率尚可的链仍出现材料函数跨链不一致。边界 native MAP 还对应无穷 probit 属性，导致原始初始化与生产混合问题混在一起。修订采用精确 probit 测度变换下的有限坐标 mode 初始化，原生 MAP 和旧失败都保留。初始化改善后，完整链的 Rhat、ESS 和 MCSE 仍未通过；因此不能把问题全部归因于初值。

实际费用与原始链分别在 `A/jobs/`、`C/jobs/` 和 `RESOURCE_LEDGER.jsonl`。每条链均保留拒绝重复状态；没有依据 GT 选择链或 production 区间。

## 第一次实质修订：显式处理属性与强方向非线性

假设是：强方向或属性被冻结，导致了错误的材料平均和采样几何。

两个独立检验已经完成。其一，在同一概率测度下联合处理属性，对 supplied Gaussian reference 作 GH5/9/17 积分；解析 decoder 矩只消除该 Gaussian 内层积分误差。32 parent 的 joint Gaussian 材料均值 MSE 为 0.0185700，仍比对应坐标 mode 图的 0.0120131 差。GH9 到 GH17 已稳定，而四个 Full 重要性检查的 ESS 均约为 1，因此不能再把失败解释为积分阶数不足。

其二，使用非线性三角 strong-response spine，但保留全部 64 个横向 Gaussian 坐标和两个属性，以真正 Full Hamiltonian 作 MH。四条链的材料函数最大 Rhat 为 1.6054、均值 RMS MCSE 为 0.01082，未优于普通 Full-gradient HMC pilot。这个结果否定了当前单一二次 property-spine 的实际增量，并未否定所有非线性坐标变换。

费用全部在 `B/jobs/`；`B/JOINT_GAUSSIAN_IMAGE_RISK.json` 是独立离线评分，`B/references/DEV_030/spine_hmc_pilot/` 是采样诊断。二者均不是合格 teacher。

## 独立核验与学习负结果

原生 Gaussian field 的 Elliptical Slice 与 bounded property reflection 作为独立机制继续运行，未限制 slice bracket 来伪造低成本。最终资格按同一门槛由保存的四链判断。

512 train、128 validation 的 joint-supervision 后，Gaussian、MDN、Flow 与 CNN 共 15 次训练全部由 validation 选择 update 0。没有把这些检查点称作学习增益，也没有为它们支付无收益的神经 Full-MH 扩展。随后只拟合已预先允许的简单归一化密度：PCA32/ridge1000/shrink0.5 的 validation NLL 为 14.8138，无条件 Gaussian 为 15.1007。这是预测密度对照，不是后验参考。

新增 128 个锁定对象在任何观测生成之前冻结全部估计器、失败 fallback、主比较及费用。预观察阶段修正的仅是报告合同：补上已登记的主比较 bootstrap 和 Full-IS 完成检查，见 `AUDIT/PRE_OBSERVATION_METADATA_CORRECTION.json`。不在读取 locked GT 后调整模型。

## 第二次实质修订：保留非线性 decoder 的代理力，而不是用错误 Gaussian 直接成像

新可证伪问题是：在相同起点、随机数、轨迹长度和物理预算下，保留真实 decoder 的 material-response 代理力，是否比 latent-linear 代理力更接近 Full posterior 的有效移动？

两臂都使用标准对称 leapfrog，只在最终端点使用完整 Maxwell Hamiltonian 作接受判定。质量矩阵与代理仅在 warmup 刷新，production 冻结。这里不要求代理力是完整目标梯度；其作用是提出可逆、保体积的相空间运动。这个通用思想已有先例，本轮只检验 Maxwell 专用的 decoder/response 分离，不主张发明 HMC。

纯代数、保存数组的 gradient/packing 检查已完成，未新增物理。付费两臂将在冻结的 locked-risk 序列完成后串行运行。当前尚无两臂结果，不预写成功。普通 Full-HMC pilot 的 warmup 与轨迹长度不同，只作实际成本背景，不能冒称仅替换力函数的配对消融。

## 内存中断的处理

中断前 109,488 RHS、失败和缺失尾部 wall 的保守上界均保留。修复的是 campaign 适配器的 transfer/LU 强引用与并发；没有修改历史 solver、物理目标或旧账本。恢复后的 context 和采样进程均在受控内存内运行。网络、代理或缓存速度不构成概率正确性证据。

## 最后实际结果与停止

两臂实际完成：每臂164384 RHS。decoder材料代理的平均绝对Phi预测误差从8.27591降至.545205，但material MCSE .0168920未优于linear .0166675，wall393.749s高于245.671s，仍未合格。128独立对象的Gaussian近似均值MSE .0615021，点基线 .0110309，负区间明确。DEV32合法投影MSE .0208401也未改善。没有据此调大noise/prior、改变locked模型或继续训练RL。

首次代理入口在Book前因重复metadata关键字报错，0RHS；原ERROR/stdout/SOURCE保留。最小parser修复后沿用同一paid job ID/ledger，Cfix1不是新campaign或新预算。

最终REDIRECT，所有未用额度封账，不为寻找正结果追加实验。
