# 冻结生成样本的付费概率审计：实现合同

新增tools/run_r4_probability_audit.py及纯统计src/p3r4/probability_physics.py。没有生成新样本、contexts或训练，没有启动物理求解；父调度器决定所有数学解释及实际付费运行。

spec示例：

```json
{
  "likelihood_samples":64,
  "importance_weights":true,
  "generated":[{
    "id":"spline_frozen",
    "directory":"/absolute/ONLINE/spline_variant",
    "online_gate":"/absolute/ONLINE/ONLINE_COMPLETE.json",
    "observations":"/absolute/stage/inputs/online",
    "parents":["R4DEV_00000","R4DEV_00001"],
    "extended_parents":["R4DEV_00000"],
    "density":"explicit"
  }]
}
```

每个group的directory内为parent.npz；online_gate必须是既有evaluate_r4_density.py产生的观测限定ONLINE_COMPLETE，样本文件的byte/mtime身份必须在该gate中已封存，并沿用同一个原target。density取explicit(保存原始dx66下log_q_x)或none(CFM等不提供合格显式密度的生成器)。parents及extended_parents在物理读取前固定，不能按残差重排。没有extended_parents的父样本只评价index0；扩展父样本评价原固定index0…63，index0包含于64中而非再重复一次。likelihood_samples只能0或64。

freeze命令：`<python> tools/run_r4_probability_audit.py freeze --spec SPEC.json --out FROZEN.json`。冻结样本目录/父样本/原观测路径、样本文件封存身份、source gate、原target、完整预定样本索引和最坏RHS计划。生成失败或缺失样本保留为冻结失败记录，不能因没有文件而把该父样本从分母移除。冻结不读GT，不付费。

run命令：`<python> tools/run_r4_probability_audit.py run --freeze FROZEN.json --stage III --job UNIQUE_JOB --max-rhs N --max-wall SECONDS --out NEW_PHYSICAL_DIRECTORY`。stage由父调度器选择，示例III不是自动授权。显式hard RHS/wall必须足够覆盖冻结最坏样本名单；预算不够时先拒绝，不选择性评价好样本。既有Book拥有CPU physics.lock，3GiB RSS上限；Target对每次8源Full收费8复RHS及实际LU，失败/拒绝不会消失。生成样本q_samples直接用于原decoder，x_samples用于原Gaussianized目标及logq权重；不截断属性、重抽或换sample0。相同或拒绝重复状态按原索引保留，仍逐个实际求解及收费，不去重冒充更多独立父样本。

物理driver仅加载生成输出parent/sample_index/q_samples/x_samples/log_q_x/observation_features和原在线observation/target_id/sigma；不接受truth参数，不读offline_only。原始复观测必须与生成器保存的条件y512逐项一致，不能用错一组噪声。每个Full点的q/χ/current/field/prediction由Target立即持久化；审计也保存固定索引、x、Phi和原观测。每sample结果及时追加，逐parent流式处理。异常样本保留NaN和明确错误，不按成功样本重归一化IS。PHYSICS_COMPLETE在Book关闭后写入，封存逐样本及parent向量的byte/mtime身份；预算/墙钟中断为STOPPED，不能作完整概率评分。

所有固定index0的Φ用于Gamma(256,scale1) prior-predictive基准：在正确prior-predictive Y和正确posterior conditional draw下，单个固定draw与Y的联合law恢复原噪声关系。该基准不需要物理driver读取GT。摘要保留计划父样本数、失败数、PIT描述及全父样本均值；有任何失败时全父样本均值不报告，有限样本均值只标描述性，不悄悄改分母。不同方法共享Y，不能把各方法的父样本当作额外独立观测。

对explicit密度和完整固定64样本，保存logw=−Φ−||x||²/2−logq_x；省略的是所有样本相同的prior/likelihood常数。先减最大logw再做logsumexp，保存所有64个原索引的logw和归一权重，报告1/sum(w²)及最大权重。它是权重集中度诊断，不能变成独立MCMC ESS、MCSE或FULL资格。若输入样本已经来自纠正链而不是声明的显式生成q，不能把proposal logq冒充该链的实际sampling law；此时该代数比值须单独解释。所有重复/拒绝状态仍保留，不能去重人为抬高ESS。

score命令在另一个进程执行：`<python> tools/run_r4_probability_audit.py score --freeze FROZEN.json --online PHYSICAL_DIRECTORY --truth ISOLATED/offline_only --out NEW_SCORE_DIRECTORY`。先检查物理gate、不同PID、Book关闭及所有封存物理文件，再允许孤立评分器读已有GT prediction；不读取q/chi并转交任何物理Target。使用已付真值prediction和原观测计算Φ_GT，无新增真值Full。固定64个Φ样本对数据依赖统计量Φ(theta;y)作含ties的随机inclusive rank/SBC诊断；没有真值prediction或物理失败时明确NOT_RUN/FAILED，不补求解。该score可省略--truth，此时只汇总Gamma及已保存IS诊断。

PURE_DIAGNOSTICS_CHECK.json通过极端共同log权重1e300、均匀重复64样本ESS64、单点集中ESS1、含失败不归一化、重复likelihood rank和Gamma失败分母检查。全部新物理RHS/状态及训练更新为0；尚无实际概率审计结果。
