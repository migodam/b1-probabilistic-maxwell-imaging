# R14 冻结独立24对象与 self-only 对照合同

状态：准备实现；所有规则在新预测前固定。仅准备公开输入任务，不读密封真值、不评分、不启动实验。

## 独立24对象

24个新对象来自原8种几何家族×3种材料regime；不是4种新材料家族。四种算法度量组合00/01/10/11，与成熟CIE、Dual-CIE、G19、DR加一个预先冻结的01→11stage候选，共10臂；全部24对象跑全部9臂，共240任务。对象ID以生成worker实际公开MANIFEST为准（预期r14_new_0000–r14_new_0023），不得从密封真值补字段。

生成器冻结方案由主线程给定：继承原几何、采集和噪声；原H~U[12,18]、L~U[2,5]，定义uH=(H−12)/6、uL=(L−2)/3。low: high=4+uH,low=2+uL；high: high=15+3uH,low=10+2uL；mixed:原H/L。所有regime显式h=high，host=.1。生成器负责记录实际公开范围/数组元数据；本准备器只读取生成器公开MANIFEST，不猜缺失字段，也不打开密封材料或评分。源输入、生成实现和公开收据在预测前冻结。

全部算法使用同一对象观察与m_init，同120秒whole-process CPU cap；proposal72秒，剩余48秒物理corrector；导入、准备、失败操作收费，不为某臂重置进程CPU。公共prior保持smooth=tv=1e−4、epsilon=.01，合法box[.02,20]。仍用原物理Phi选择已付费合法候选，corrector为既有物理L-BFGS-B，maxiter3000/maxfun12000。失败不删对象、不调整case-specific参数、不延长预算；全部臂即使别臂失败仍执行。

四个metric-pair使用R14原StatefulDualMetric/run_experiment，不改核心。M0/M1及stage、globalPR、Fourier/DPCS与R7Arc2sweeps/继承gtol按既有dualmetric48_v1公开protocol冻结。基线使用逐字节复制的原R13 frozen_C2 runner与其历史运行模块，保留各基线既定初始化/表示及sourcewise规则，不将其强改为pair配置。

run_independent通过同进程runpy分派，保留全部命令行参数，单独记录wrapper/历史runner来源；不导入优化器或新backend。queue_independent是run_queue的runner路径替换，case分组在同一核顺序跑10臂，12个固定逻辑核0,2,4,6,8,10,12,14,16,18,20,22；保存host资源审计提供的core-type标签，未提供则显式UNVERIFIED，禁止跨cohort旧绝对CPU值比较，NORMAL、单BLAS线程；无服务、无自动重试，不动其他任务。

## 暴露集 self-only16

固定对象索引顺序[2,5,7,16,17,24,29,37,38,42,44,45,8,20,32,46]，只跑pair01与pair11，并设置onehop_self_only=True，48任务（其中32个self-only、16个原global-Dual）。输入仍是data_public/exposed48原公开NPZ。其作用是区分M1对角自作用与原四cell非局部项，不是独立确认集；费用、prior、停止、corrector及读出版本与原pair合同一致。此对照保留M1公式并将实际blocks.goff清零，不增加传播阶数。

预测完成且封存后由主线程另行评分、科学审查和发布实际结果。本文件和准备manifest本身不构成完成实验或新恢复证据。

## 冻结stage候选与公共选择回放（预测前补充）

pair01_to11固定metric_schedule=[01,01,01,11]：前3个原stage保持01，进入原mf12stage后使用11。只利用既有stage切换已有PRreset，不加阈值、先验、观察、selector或预算。动机来自暴露集同状态11对物理consolidation的观察；本候选尚无独立恢复证据，必须计新增与丢失，不能预称改进。由run_control执行；普通pair仍用原run_experiment。

另冻结Pro公共选择规则为独立集候选回放基线：pair01对应原a_plus_arc_r7；顺序pair01→G19→DR。严格采用bool(当前incumbent.noise_feasible)停止；否则调用下一独立uniform子结果、先计其external_child_CPU_seconds，仅当trial.Phi严格小于incumbent.Phi才替换。只用公开列，选择后才评分；不读取或调节NRMSE/IoU/family。所有任务仍为预先全跑，回放总费用为实际被规则调用子程序CPU之和，不是实测条件控制器部署，也不是总120秒组合。stage候选不参与该冻结回放。

预测前补充：全部新24加原global-Dual；曝光self-only16各加一条原global-Dual，区分material beta变化、Gii自作用与Goff非局部传播。其目的为归因，不新增自由参数。

## 评分/独立预测前的预算接线修正

参见FAILURES第3项。以上较早“48秒corrector”和“逐字C2runner”描述为初始计划，实际独立部署统一改为proposal最多72、总120剩余时间可用于corrector；C2仅改一行deadline并完整保存原文件。四组合与基线、stage候选、self-only均同口径；没有改数值优化器或GT选根。另报原48口径前缀。主曝光48四组合已经按总120执行，将按实际合同解释，不覆写原预测。三维六臂保持其独立严格48合同。

最终预测前补充：加入同m0、同数据/box/prior、总120CPU的原物理Direct L-BFGS-B（冻结C2 maxiter6000/maxfun24000），避免只比较辅助算法而遗漏强物理控制。独立24最终11臂264任务。Direct本身是全程物理优化，无proposal/corrector分段；严格48corrector前缀对它不适用，其全部实际CPU/物理RHS仍完整计费。其余方法/冻结选择不改；不再扩增方法名单。
