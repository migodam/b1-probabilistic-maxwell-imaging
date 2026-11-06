# CM九臂协议增补：先冻结，后运行与统一评分

本增补在新九臂启动前写入。原3λ four-cubes、60源、180接收行、30³全27000实材料未知量、CM+RR算子、χ∈[0,4]、air0、R=0、ρ=1、MF6→10、每阶段最多500、Eq16根号RMS<.001及无物理corrector保持不变。全部采用已核验的complex128/float64连续FFT packed16；旧strided4部分轨迹和所有费用单独保留，不计为新九臂完成。

两个不同数据cohort分别报告：

- **cross_grid：五臂**sourceCIE、sourceDual、GlobalDual、AplusGlobal、LSIEsource。复用已经冻结的CM60³生成/CM30³反演观测，Y逐元素不变。独立网格真值响应差global0.535254；不能称作者积分物理复现或连续Maxwell精度验证。
- **matched_grid：四臂**sourceCIE、GlobalDual、AplusGlobal、LSIEsource。复用已经付费的CM30³全60源clean预测构造新的独立公开data_id。没有新增forward generation。此cohort是有限离散模型匹配诊断，不作为同数据的算法收益、独立网格或连续物理复现证据。

matched噪声严格沿用原冻结事实：原数据为iid circular CN，sigma=.1·||Yclean||/sqrt(n)，期望相对范数10%，未强制样本范数精确10%。matched复用原seed2026110619的同一噪声方向，并按||coarse_clean||/||fine_clean||缩放sigma及噪声。因此两cohort的实际相对noise norm相同，预期均10%；分别保存sigma、实际范数、噪声方向复用事实，不归一化或修改旧fine观测。父线程原“固定范数”用语已明确更正，不引入新噪声规则。

一个GPU进程驻留共享G/E/S及一次SGram/EVD/Vr，rank阈值固定.1；两个Y各自构造Jd/Qd/Dd及相对数据权重并计费。sourceCIE的current/material beta6→1；sourceDual/GlobalDual/AplusGlobal current6→3、material6→1；LSIE为普通p=1。AplusGlobal保留R9向量A+固定J完整A损失接受CIE种子提案的已披露变体，不称R8原实现完全相同。

只在air与两阶段末做完整物理验证，最终阶段为主要端点，不选GT最佳迭代。允许以材料数组严格相等为条件复用确定的物理预测，必须记录cache来源、命中和首次付费；air解析零预测实际完整RHS=0，物理value调用与求解RHS分开。实际queue成本与共享/缓存摊销分别记录，不能把后运行方法少付一次共享工作解释为独立算法速度优势。

原尺寸物理指标明确为Φ=mean_s[||prediction_s−Y_s||²/||Y_s||²]，没有1/2且R=0；不能把其数值与PhaseI的Φ横向比较。费用表并列实际共享queue费用、声明setup摊销口径及按每个cache命中首次evaluation CPU/wall补入的归属估计；该估计不是新独立冷启动计时，不再次计入全进程总账。

data prep虽打开sealed容器，但仅提取Y_clean/noise/已经付费的coarse prediction，没有读取材料字段。材料真值仅由全部九seal之后的独立offline评分程序读取；反演接口仅接收公开Y、公共算子和box/grid。

阶段内sourceCIE与LSIEsource的精确线搜索+精确R0材料readout应使其各自辅助目标不增；只对已付轨迹报告最大增加，不再做新数值探针或人为成功门槛。Dual/A+没有同一个共同阶段目标，不施加此要求。

**九个最终端点全部封存后，才统一读取材料GT评分和作图。** 两cohort分表，分别关注GlobalDual/A+、CIE/LSIE及cross-vs-matched。新增数据诊断不构成A+优势；不做全局凸性、唯一性、强新颖性或论文软件复现声明。

IGT仅完成预声明source59的30/60两个clean RHS：相对差CM.508374→IGT.359841，仍高于10%噪声。保存负结果；不生成IGT60源包、不扫网格、不启动IGT反演。该诊断不改变当前CM九臂。

最终去重账`RESULTS/original_vector3d/ACCOUNTING_FINAL.json`：历史十个CM数值/工程进程、IGT prepare/solve两进程、新九臂whole-process分别为不重叠additive行；每臂、共享预计算、每Y准备、cache与GPUevent为nonadditive明细。所有失败、中断与未知计数范围保留。
