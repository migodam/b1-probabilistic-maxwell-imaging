# 可能材料集合独立元数据审计

**状态：PASS_METADATA_AUDIT_WITH_SCOPE_LIMITS。** 本项只读公开封存元数据、冻结source/config身份、ledger与parent完成/visited摘要；不读取私有GT或材料数组，不执行physics/统计评分，不更改冻结方法。后续root已另行授权统计结果审查，二者分开记录。

## 完成与冻结

- 32个预声明fresh parent，3个arm，每法16个方向slot；合计1536条方向记录，slot0…15完整。无parent按结果剔除。
- 封存文件23,016个，共3,226,644,701 bytes，路径/存在性/size元数据核对错误0；没有private/truth文件进入公开seal。没有新增SHA检查。
- 8个冻结source/operator/config身份与记录一致；统计冻结先于首次generation ledger。运行期修复仅为首次ledger前duplicate-role日志参数错误，receipt显示当时无physics/online文件/private truth；旧源码/空失败目录/receipt保留。
- 所有parent声明data_only_selection=true、truth_read=false；候选GT_used=false。visited摘要中的lowest-qualified-dataPhi选择核对无误，同窗口方法轮换顺序无违例。

## 资格与对照解释

|方法|执行status计数|primary-qualified方向数/512|
|---|---|---:|
|FIXED_Q8_RAY|{'EXECUTED': 494, 'RADIUS_BRACKET_NOT_FOUND_WITHIN_LATENT_CAP': 11, 'NO_RADIUS_PROGRESS_ENDPOINT': 7}|354/512|
|DYNAMIC_Q8_RAY|{'EXECUTED': 494, 'RADIUS_BRACKET_NOT_FOUND_WITHIN_LATENT_CAP': 11, 'NO_RADIUS_PROGRESS_ENDPOINT': 7}|362/512|
|RANDOM66_FULL|{'EXECUTED': 512}|103/512|

资格采用原绝对95%噪声球Phi≤282.8737921622988、原decoder legal、目标材料半径2%及signedprogress≥.5r；旧DeltaPhi/DeltaPrior仅secondary。此处只核对记录内部一致，独立材料RMS重算留给已获GT授权的统计审查。qualified方向数不等于unique材料数、经验GT覆盖率或posterior质量。

RANDOM66_FULL是独立native Gaussian66方向经材料线性范数归一化，并使用普通full66 signed-row corrector。它与weak-rank8方法共用trial/J上限、观测、anchor和target radii，但方向分布及允许校正几何类不同；因此是有意义的matched-budget随机探索基线，不是严格的随机rank8子空间消融。其资格率较低不能单独证明实际GT覆盖更差。

## 费用闭合

- ledger共124,379行，sequence连续，最后快照一致，33个job全部关闭，无active jobs。角色：{'DATA_GENERATION': 1, 'ONLINE_ANCHOR_AND_CLOUD': 32}。
- 新增293,088 RHS256、8,380 LU256；完整RHS576/LU576均为0。forward实际/完成67040，J实际/完成226048，合计与reserve一致；失败action=0。
- 费用独立分解：generation={'reduced_RHS256': 256, 'reduced_LU256': 32}；anchor={'reduced_RHS256': 53472, 'reduced_LU256': 1512}；全部proposal（含未达门/拒绝工作）={'reduced_RHS256': 239360, 'reduced_LU256': 6836}。总和等于snapshot，没有成功样本免账。
- 冻结上限RHS256400000/LU25612000；实耗低于上限。controlled wall=449.965s，CPU=396.337s，单进程RSS峰值301.8MiB，低于wall7200s/RSS8GB上限。
- 精确ROI setup历史264个完整RHS及整个setup研究2184个完整RHS保留来源；本轮cache hit没有重付setup，也不称reduced工作免费。单方法部署anchor/setup分配视图不可相加当research total。

## GT隔离证据与边界

- 原online_observation仅接受observation/sigma/parent三字段；在线代码不import独立truth generator。统计模块的私有truth读取仅在score入口、全部32 seal与root后置批准之后。
- UTC顺序：统计冻结 2026-10-29T16:11:21.700592+00:00；首次ledger 2026-10-29T16:19:35.901207+00:00；seal 2026-10-29T16:27:07.817939+00:00；评分批准 2026-10-29T16:28:01.968035+00:00；TRUTH_ACCESS_BEGIN 2026-10-29T16:28:07.109453+00:00。
- 当前元数据/source隔离证据支持“未发现在线GT泄漏”，但不是操作系统文件访问日志的形式化证明。生成器创建truth本身是必要的数据生成，不等于在线优化读truth。
- fresh seeds为[20261030,900901,role,parent]，role101/102/103/104区分prior/noise/weak/random；antithetic正负共享pair基础draw是设计要求。统计独立单位始终parent，不能把方向/源/seed当额外n。
- 本审计不判断覆盖效果好坏、不主张posterior准确率；有用性需看已封存集合的经验GT覆盖、材料宽度和完整费用。
