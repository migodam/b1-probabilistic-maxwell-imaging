# Research B 必要数值合同与 E0 代码审阅

2026-10-15。审阅对象为 `research/research_b/code/{sensing.py,cost.py,run_e0.py}`；没有修改 `sensing.py`、`run_e0.py` 或原 Research B package。仅修复新 `cost.py` 的执行限制，并添加 `research/research_b/tests/test_scientific_contracts.py`。

**已执行的是软件数值合同，新增真实 Maxwell attempts = 0。** 没有运行 solver、参考积分 campaign 或训练，没有重新运行主线程完成的 E0。本文不能作为物理有效性、主动测量收益或未知拓扑成像证据。

## 执行回执

12项 unittest 全部通过，其中共享 gain 更新测试包含27种组合：old/new通道数分别遍历32、16、1；模拟无约束 gain 为1、1.8、.2，覆盖内部及上下截断极端情形。额外测试无响应的新通道不能改变材料权重。

- 归一化 posterior 与独立原 Joint P2 `log_likelihood(joint)-log_likelihood(old)` 最大绝对差：`4.443667656062189e-12`。
- 扣除每个draw的常数后，conditional log increment 最大差：`2.7284841053187847e-11`。
- 被测未截断gain均值离合法支持最远约`80.36`个条件SD；精确体积L1的529对最大绝对差`2.220446049250313e-16`。
- 测试 Python CPU：`.538068s`；完整运行 wall：`.656741000013426s`。含导入和测试，不是物理计算费用。
- 完整逐组合指标：`NUMERICAL_CONTRACT_RESULTS.json`。
- 复现：从 `b1-research-b` 运行 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ../../.venv_nn/bin/python -m unittest discover -s research/research_b/tests -p test_scientific_contracts.py -v`。

所有计费边界测试使用自动清理的独立临时目录；合成4096行账本不是当前实验已消耗4096物理次数。实际 Research B `PHYSICS_LEDGER`、`COMPUTE_LEDGER`、旧 Joint P2 账本均未被测试写入。

## 共享 gain 条件增量

`conditional_log_increment(..., model='shared')` 从旧数据条件 gain 的截断高斯积分新观测，保留新旧 precision比和截断区间概率项。它有意省略仅依赖新观测draw的Gaussian常数；该常数在材料后验归一化时抵消。因此不能单独拿它的未补常数原值称完整 predictive log density，但用于归一化 posterior、event posterior与本代码的 risk/EIG 计算是正确的。

独立对照调用未更改的 `src/b1_joint/likelihood.py`，先拼接全部 old/new fields、observations、noise，再做一次共享 gain marginalization。测试同时保留非均匀材料先验，避免只有等权情形才通过。极端截断测试仍得到有限值和一致 normalized posterior。此结果验证实现等价，不代表截断gain prior就是现实仪器的正确模型，也没有为任意无限浮点范围证明误差界。

`model='independent'` 是故意改变 calibration 假设的负控制，不是同target替代优化；`model='fixed'` 相当于gain=1条件模型。二者若用于比较应明确标为mismatch/不同target，不应混称公平共享gain方法。

## 材料 MAE 的定义与独立核对

`field_l1` 正确计算 **R=1球体归一化体积上的 relative epsilon 实部/虚部分量 L1 和**，s_ref=1。它不是 complex modulus L1，也不是旧 `material_distance` 的平方根体积L2。

独立实现分别调用现有 `decode`，合并两材料所有真实层界，在每个完整球壳中点取各自epsilon，并按 `r_hi³-r_lo³` 精确积分分量差。比较23个材料全部529对，包括solid、中心空腔、环形空腔、combined、随机合法材料以及不同inactive coefficient对应同材料。归一化球壳体积中的4π/3抵消。

解析检查：solid epsilon2.5+.1i 与半径.25的同epsilon cavity，损失为 `.25³ × (1.5+.1)=.025`；两个solid epsilon2.5+.1i与2.7+.15i，损失为`.25`。这些与vectorized实现一致。该软件检查只验证声明同心支持上的损失，不表示一般3D材料或未知拓扑已成像成功。

## 计费实现发现与已修复内容

原实现有两个真实执行漏洞：只拒绝名字恰为E0的stage，任意其他stage字符串可以避开E1/E2 wall检查；overlay被读取但更低cap和`execution_authorized=false`未执行。

现已修复为：

1. 仅E1/E2允许付费callback，且必须有job identity和显式`execution_authorized=true`；E0、未知stage和无job调用均在callback前拒绝。
2. attempt cap取用户硬上限4096与overlay值的较小者；wall cap取7200秒与overlay值的较小者。无效、负、非有限上限拒绝加载。
3. `BEGIN`在callback前追加，失败保留`END FAIL`且消耗预约；重建Cost对象不会重置计数；重试同样必须占剩余额度。
4. 第4096次允许，第4097次在callback前拒绝；E1/E2共用同一上限，不因换stage取得余额。
5. 既有失败job wall也计入共享wall上限，更低overlay wall达到时在job body前拒绝。

没有改变历史授权、已有config、原 Gate 或实验结果。当前4096/7200 overlay的科学预算数值不变。

## E0 结论必须保留的范围

以下已告知主线程，主线程负责实验解释及后续算法，未在本子任务改动E0：

- **realized MAE是reference-derived decisions的选择实验。** 当前 `score_realized` 的selection/update使用2048 reference particles、96个reference-derived decision materials；candidate-only输入的3–7个full-validated代表仅用于action评分。因此不能将该MAE直接称为原候选集合上的candidate-only闭环恢复成绩。主线程将分开这两种评价。
- `structure_correct` 当前为 posterior event argmax的event决策是否正确，未必等于MAE最优所选材料的event。它应独立命名，并另报selected-material structure。
- `uniform_one_block` 与`predetermined_second_frequency`当前完全相同，均选2.4，不能计作两个有独立证据的竞争策略或多频uniform sweep。
- `random` 每case使用相同seed，对同menu抽同一个action。它是一个冻结随机选择，对随机策略期望的证据很有限；后续独立parent比较宜使用预先冻结的不同parent随机流。
- `set_scores` 只比较不同预定义event且materialMAE≥.005的pair；其`robust`是该有限**跨event**pair集合的heuristic，不包含同event不同介电/几何的所有材料区分，更不是全材料支持的minimax保证。
- `conflict` 是在旧/新响应间比较分开或共享 **relative gain scale** 的有界secant分数，未做远分支epsilon/geometry的full nonlinear bounded compensation。它可以作为screening heuristic；不能仅凭E0排名就证明完整H1机制。
- calibration action为冻结的unit-reference/gain-only Gaussian模拟；它不是新硬件实测，也不是不付费的外部信息。其测量block成本和误差模型必须保留，未来真实Maxwell calibration需重新匹配实际标准和单位。
- 同一batch的common random numbers有助于paired utility比较，但256个predictive draws与有限inner/decision池仍有MC/离散decision误差。排序接近时需要已预算内的固定refinement或明确保留不确定性，不能把排名第一自动当显著优越。

本次没有发现需要修改 `sensing.py` 的shared likelihood或volume loss数值错误。主策略是否改善材料MAE、是否值得付费及H1是否成立，仍需主线程的真实E0/E1结果裁决。

## E1/E2 运行源代码只读审阅（追加）

同日新增只读审阅 `active_experiment.py`、`prepare_e1.py`、`configs/E1_FROZEN.json`、`configs/E1_IMPLEMENTATION_DETAILS.json` 及其调用的 `e0_candidate_controls.profile`。主线程正在执行E1评价时，**未修改上述任何文件、没有启动campaign、没有调用solver，也没有新增不相关测试**。下列是执行合同与解释范围，不能替代主线程的材料MAE结果。

### 真实材料和未来观测隔离

`prepare_e1` 独立生成32个prior parents，各2个noise/gain episodes，并保存隐藏truth和各action的潜在观测。在线 `recommend` / `candidate_scores` 接收当前history、候选材料、相应响应和noise；不接收truth latent、隐藏gain、未来noise或未来观测。真实材料与未来潜在结果由 `Simulator` 的初始观测、测量和offline scoring接口读取。基于所读代码未发现用truth挑action或对未来观测预看。

`Simulator` 与在线函数位于同一个本地程序，属于接口/代码隔离而不是操作系统级sandbox；不能声称在线过程对隐藏NPZ具有访问权限层面的强隔离。当前函数调用链保留了约定的数据方向。

### gain、noise槽与 likelihood

- 每个parent/noise episode只有一个proper truncated real gain，生成公式将其广播到全部频率和全部step。
- 初始1.6观测使用`outcomes[parent,noise,0,index(1.6)]`；第1/2/3新增散射block使用step槽1/2/3。新repeat1.6也取新槽，未复用初始noise。不同策略共享相同parent/action/step潜在结果，用于paired比较，不能把这些共享结果当更多独立物体。
- calibration与external geometry各有三个独立noise槽。已知solid参考的充分统计量为unit-gain观测，SD为原1.6 real-channel sigma除以固定已知响应范数；与准备文件和在线实现一致。其隐藏均值是同一episode gain，未单独为calibration重抽gain。
- 对新scattering/calibration先以完整旧history计算conditional increment、更新当前材料权重，再追加新block到history。它等价于合法joint-history更新，未另乘一份旧likelihood，也未给不同频率重新独立gain。
- external geometry使用独立void-volume likelihood一次，candidate-only增加对应profile residual penalty。它不进入散射gain history，因其不依赖gain；这与所声明独立传感器模型一致，不代表现实中必然独立或免费。

### Mode 1 确实使用625个已验证材料

`active_experiment.run` 从E1 bank的四轴网格显式构造625个latent；六频`CF`直接读取对应原始`fields_k`，不是把proxy预测冒称全full-Mie候选。1.6/2.4复用既有1250节点，四个新频率使用本阶段已计费2500节点。

Mode1 update调用 `profile` 在**全部625材料**上最小化joint-history residual；同一候选跨所有观测保持同epsilon/geometry，只允许各世界自身一个合法共享gain。它没有归一化材料权重或residual softmax。采集screening另取best+9内的最多48个当前compatible candidates、材料L1≥.02的pair，包含同event pair；该set heuristic不具有全材料空间保证。screening没永久删除其余625候选。

Mode2保持原coefficient prior，初始order13 conditional reference抽取2048 particles、96 decisions；随后使用相同shared-gain增量更新这些有限样本。两个模式的selected-MAE、Best-of-K及结构选择已分开字段。**reference-based EIG/risk选action、再用625材料profile选最终材料的组合，仍依赖Mode2 acquisition资产；不能因为最终选择字段名为`mode1_candidate_only`就宣称整条策略无需reference。** raw/robust/Fisher等仅使用candidate/history的action才对应完全candidate-only部署。

### 已保存数值资格（只读）

64独立prior probes对每个新频率的最大response relative/whitened误差为：

| kR | relative | whitened |
|---|---:|---:|
|1.2|.0000496825|.00191499|
|2.0|.0001404304|.00641981|
|2.8|.0006127812|.02812274|
|3.2|.0013910752|.07114990|

每频前8个probes的n14/n18白化差最大约7.22e-17，独立recursion与matching白化差最大约2.07e-14。full-Mie数值收敛和proxy插值误差是两类误差，后者随频率增大，不能以multipole收敛认证proxy的全支持action-ranking精度。

读取时 `E1_INITIAL_REFERENCE_AUDIT.json` 已包含全部64 episodes；13/17 event质量最大差`.00027739751598679363`，保存的各坐标mean最大绝对差`.00009860684103735196`。它只核对初始观测参考，不是全部新history的continuous posterior certification。后续reference weights仍是有限粒子近似；高ESS也不能证明没有遗漏支持，低ESS时须按实际误差保留限制。

### 需要主线程处理的诊断与成本限制

1. **geometry residual sentinel自由度漏计。** external geometry每次增加一个残差项到`cq`，但阈值仍使用`len(scattering/calibration_y)-1`。这不影响选action、选材料和MAE，影响sentinel标志的名义阈值。可在独立派生分析中按已获取geometry数修正；考虑625择优、gain边界与有限候选，仍应标作启发式，不宣称精确1%误警控制。已即时告知主线程，未改运行源文件。
2. **pre-kernel非法latent未统一收费的入口。** `prepare_e1.forward` 先decode再进入metered backend；decode异常会在预约前抛出。所读625网格和fresh truth均合法，没有据此发现本次已执行attempt漏收，但未来非法setup路径需要单独失败预约或明确记录拒绝发生在physical attempt之前。没有重置或补造已执行账本。
3. **E2没有实装数据驱动early-stop。** 除none外，policy当前会走满设定blocks；pair为空时零分也会按menu顺序选首项1.2。`stop=False`变量未使用。当前只能报告固定最大block性能/成本曲线，不能宣称已学会在无价值测量时自动停止。主线程可根据E1决定不启动后续campaign，或在新的已记录协议中实现stop；不能静默改变正在运行实现。
4. **refresh只扩充展示/decision catalogue。** E2每次把2048已有reference支持中的64个样本追加到初始96决策池；它不扩展这2048个支持材料，也不是Research A的新候选发现。固定625 Mode1 Best-of-K严格不变；Mode2 Best-of-K改善只能归于输出catalogue扩充。
5. **共享portfolio计时是保守上界。** 初始recommend为多个方法一次计算整组scores/utilities，单策略字段可能承担整组portfolio时间；不是每种优化实现的孤立推理计时。总job wall包含实际工作，可报告真实inclusive费用；不能从单策略字段直接声称某criterion有精确速度优势。
6. **parent统计仍需在汇总中落实。** raw rows保留parent0..31、noise0/1、policy和step；独立单位是32个parent，须先平均其两noise再做policy paired差与2000次bootstrap。不能把64 episodes当64独立物体，也不能把重复policy step0基线当独立数据。此次审阅的源码尚无最终统计汇总函数，故未认证最终区间实现。

E1初始refinement和forward audit均为本阶段已支付/已完成资产，本子任务仅读取；本次追加审阅新增Maxwell为0，没有将软件审阅提升为物理验证结论。

## E2覆盖诊断执行与审阅状态更新

主线程后续E2实现已经增加 `risk_cost_stop`、`robust_abstain`、two-block bundle，并在新源码sentinel自由度计入geometry观测。上节“不支持early-stop”和geometry自由度发现描述的是**当时读取的E1实现快照**，不应用于否定后续已修改E2实现；本子任务没有因此重跑E1或修改主源码。

按主线程独立任务，新增并直接执行 `code/coverage_stress.py`；冻结定义在 `configs/E2_COVERAGE_STRESS_FROZEN.json`，结果 `results/E2_COVERAGE_STRESS.json`。这是同阶段E2诊断，非新科学迭代。指定计费job `e2-coverage-stress-v1`：PASS，CPU `.742883s`，wall `.7777362910564989s`，新增真实Maxwell attempts **0**。

固定625材料先做full、删central、删annular、删combined四组（625/125/125/225）。每个parent/noise在各删除组**重新计算raw/robust的全部七action分数**，只接收当前观测和保留候选。真实材料仅用于测量模拟与离线score。robust所选action相对full组的改变次数分别33/64、13/64、9/64；没有沿用full集合的推荐冒充deletion后的策略。

全32 parents、每个2noise、4组、6策略、2个独立列出的场景，共3072逐episode行；噪声先在parent内平均，形成1536 parent行、48 summaries，并保留2000次固定parent bootstrap的paired差。

| in-model支持组 | none Selected-MAE | fixed3.2 | robust | 固定Best-of-K MAE | Coverage@.05 |
|---|---:|---:|---:|---:|---:|
|full625|.162279|.101982|.103714|.062224|.40625|
|删central|.173102|.108934|.123012|.068414|.375|
|删annular|.291167|.282930|.282914|.206102|.3125|
|删combined|.167504|.106912|.105510|.065289|.375|

本有限面板不支持robust普遍优于简单fixed3.2。删annular后的robust-minus-none paired均值−.008253，95% parent bootstrap区间[−.020611,.003944]跨零；测量无法补回已经不在固定候选集中的ring材料。Best-of-K和coverage对固定集合的policy与新增measurement保持不变。

H3的直接负证据：删central使15/32 parent的真实结构family不在支持；robust更新后缺失的30个noise episodes中29个没有residual警报。删combined缺失的16 episodes中15个没有警报。**低residual不能认证未遗漏材料family。** 这是启发式sentinel漏检，不是有posterior的错误校准结论。

`sentinel_first3.2`与fixed3.2采用同一次测量和相同profile结果，只把警报作为诊断abstention；始终另报无条件proposed-material MAE和接受率，不用只剩容易case的selective MAE制造优势。gain-drift场景以相同已付费Maxwell field和原noise改变未来instrument gain，保持错误共享gain拟合，独立列出、未与in-model池化。它不是新物理计算，也没有改写E1潜在观测。
