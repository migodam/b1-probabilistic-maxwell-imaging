# R4正在执行：恢复入口（非最终报告）

工作区：/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r4-conditional-transport；Python ../../.venv_nn/bin/python；目标研究目录 research/p3l_r4_conditional_transport。分支research/p3l-r4-conditional-transport-20261025。初始base c7594d4f，首次R4源码冻结03c4b6b7。

已完整读取R3正式报告与执行方案：主线程主报告，density代理附录B，physics代理附录A/C；原P1/p3r2与固定6eadb469无源代码差异。授权以用户本轮R4任务为准，包内静态未授权文字不阻断运行。只用native6.1-sol任务代理，medium/high，不DeepSeek、不SHA、不新建repo。

独立授权：0/I/II/III/IV RHS0/24M/6M/8M/20M，总58M；作业wall2/24/12/24/32h，总94h；GPU0/48/16/8/24h，总96h。原P0/P1/P2/P3L/R2/W全部历史封账不改。R4账本未封账，必须继续本campaign，不能误认为W已完成就停止。

已完成：Stage0probability/cache合同PASS(0新physics)；首个新anchorFD相对误差1.67e-9，真实source-adjoint与receiver-J梯度差2.15e-13(实付56RHS)。恢复32768train/2048val完整原数据。16项matched512/128风险拟合/经典控制完成，3seed新localCNN旧validation平均clippedMSE.00722639，R2冻结CNN.00752941，探索性约4.02%改善，非新locked结果。阶段I首轮risk日志logs/risk_first_sequence.log。

新R4DEV128已生成并完成native96+probit64上下文，128/128noiseball(非MAP/后验资格)。所有GT只在inputs/offline_only，线上只读inputs/online与results/contexts。

StageII相事件有限反馈8旧DEV：Full/reuse误差1.31e-15，Full7.427svsreuse+setup8.374s，无wall优势；pCN星形审计54接受/202拒绝/0误判。原full-gridRHS8003加reducedcomplexRHS6144保守计费，合14147，closed原始回执不改，根账本追加reconciliation。此支线STOP_REUSE_NO_WALL_GAIN。

StageIII几何两轮真实完成：459712RHS、16060states/LU、449.28swall。首次32对象192config四段Euler含h²self缺陷，不能拿它否定bracket理论。Heun修订8固定DEV48config小h门48/48经验通过；任务尺度Heun仅1局部teacher，direct94，Heun288RHSvsdirect8，未建立路线增量。全部报告、失败与源码在results/III及AUDIT/geometry，切勿重跑。

StageII新决策精算prototype已完成：results/II/decision_refresh_star_v1，8新DEV×128proposal，R0一次576RHS+Richardson/Neumann残差势区间，1024个Phi区间全包含Full，0fallback/0decision错误。这里只是fixed-star、exact-arithmetic条件normbound，float64非区间认证，不是后验链。GMRES较慢，所有Fullaudit都实付；计费ADDENDUM补充GMRES自身8source/proposal，原实际总62770RHS未改。

当前后台：
- logs/density_runtime_fix1_sequence.log：raw-y全66条件模型；原logs/density_first_sequence.log首启动因MPS默认low1.4>显式high.5失败，已保留成本/日志，新job_runtimefix1仅修runtime，runner现在low.4/high.5。LR3次+4families×3seeds+4096/16384数据曲线=19fits。每fit6000updates(除LR1000)、MPS、max1950sec，逐个串行。原first失败不得重置。
- logs/context_expansion.log：原TRAIN前1024、VAL前256重新建立相同96+64强上下文，64对象一块；已预登记CONTEXT_EXPANSION_01，整个上界8.192M属于I。每个上下文自己Book/端点，不能冒称新独立材料truth。
- configs/INTERVAL_CHAIN_RELEASE.json已登记priority；上述context循环下一块前调用tools/run_r4_interval_chain.py默认pair，分别interval/full/compare三个进程，DEV16..23×1024 pCN66 beta.002，sharedGaussian/U，II总cap200k，1200s/arm。必须等PairCOMPLETE/比较结果后继续，任何mismatch NO_GO。当前/候选双Phi区间，accept后状态仍approx，不能假装新exactcache。阶段II仅floating条件包络+Full重放，非认证posterior。

需要继续执行，不等用户微批准：
1. 完成当前raw66训练，检查实际validation/NLL+更新曲线；补diag_gaussian容量匹配3seed(682943vsSpline682860)，谱triangular vs固定随机frame对照。CFM完整trace实际NFE16/32/64只作数值density审计，不作exactMH。
2. 运行exact-bin phase revision1(phase8/16三seed)，单独risk_bin.py/train_risk_bin.py；pointwise标签恒等误差~1e-16(double)，不改变P1，移除点阈值积分偏差。原16训练不改。
3. strongcontext1024/256完成后pack_r4_context_training.py，训练新anchorfull66 Spline/diagGaussian/CFM三seed；GN/随机坐标新变体由density代理正在实现(所有66维、完整logdet)，缺省translation变体已做代数检查。普通raw-y不能声称使用了anchor。
4. 全66生成器与风险动作输出分开；新512locked+256确认尚未生成。必须先冻结checkpoint/超参/normalizer/K128/样本index0/物理预算/主比较，再生成、线上推断完成后独立PID读GT评分。旧R2两批256均EXPOSED。
5. tools/evaluate_r4_risk.py及evaluate_r4_density.py提供freeze/infer/score；概率checker必须另实际paidFull检查Gamma256/sample0与likelihoodSBC。CFM样本不是exactdensity，不能直接以近似divergence做exactMH。Fullposterior失败不能关风险和几何。
6. II还需strong learnedJacobian/DINO-like或refresh对照；可以复用已付III的FullJ/field标签，但必须披露teacher获取成本，禁止oraclefallback假称可部署。
7. 完成报告、实际数组图、各GO-A/B/C/D/E独立裁决、GPTProHandoff，发布到现有B1Pages+Release并匿名回读，给网页和rawMarkdown链接；现在尚未发布R4。

资源与内存：1physical+1training由不同锁控制，共享budget.lock原子计费；receiver强引用cache<=2。所有失败、验证、拒绝成本都保留。后续source_snapshot在Bookjobdir记录实际新src，gitrevision外补workingtree内容。早期代理数值神经审计1.8853s已保守计入I的GPU allocation，并非Stage0GPU。
