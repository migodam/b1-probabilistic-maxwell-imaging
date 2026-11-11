# R12 冻结 NN 接口四对象比较

四个 ID 在新推断/新反演前冻结，全部 16 个预测及四个 NN 初值封存后才打开 truth。无训练、无在线 GT、无物理 backend 替换。

实际 R8 ID 为 0001、0004、0022、0027；0001 为 6 源，其余 12 源。0027 是既有暴露的 heldout-family 反例。它们不是 R11 同编号对象，也不是独立确认集。

原 physics NN（epoch86）及原十通道 features 的四个输出与历史初值完全一致。所有比较使用原 [.02,20] box、原 Maxwell Phi 和原 smooth/TV prior。

DR 使用完整物理 NN 双场开始，scale6/omega1，3000 步，每50步物理审计，global CN99 数据球及 grouped_ports=True。CIE 为原 R8 cie_single/原DPCS、Fourier和3000步调度；不初始化物理 current，保留原接口。所有四臂使用同一个最多3000迭代、120CPU秒的原 corrector；proposal 费用另外实测，不声称总费用相同。

|对象|方法|幅度NRMSE|IoU|原Phi|solverCPU秒|corrector停止|
|---|---|---:|---:|---:|---:|---|
|r8_test_0001|uniform_direct|0.668228|0.1895|0.004114395|18.33|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|r8_test_0001|nn_direct|0.076345|0.8824|5.20408e-05|13.62|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0001|nn_dr|0.077289|0.8736|5.204397e-05|16.88|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0001|nn_cie_single|0.076404|0.8824|5.204381e-05|20.03|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0004|uniform_direct|0.368955|0.0000|0.5227501|0.06|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0004|nn_direct|0.080315|0.6230|0.000418066|10.27|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0004|nn_dr|0.080656|0.6179|0.0004180661|16.88|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0004|nn_cie_single|0.080315|0.6230|0.000418066|21.12|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0022|uniform_direct|0.289828|0.0364|0.1004618|4.18|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0022|nn_direct|0.289830|0.0364|0.1004618|4.64|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0022|nn_dr|0.066429|0.7843|5.351844e-06|21.12|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0022|nn_cie_single|0.065756|0.7692|5.350793e-06|29.98|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0027|uniform_direct|0.460948|0.3243|0.02242874|16.70|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0027|nn_direct|0.120932|0.8562|0.0004286581|14.90|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0027|nn_dr|0.120458|0.8431|0.0004286547|23.55|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r8_test_0027|nn_cie_single|1.012975|0.2345|0.0004626687|33.33|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|

DR 相对 NN→direct 的材料幅度NRMSE：改善 2/4，变差 2/4。四臂均保留原物理 Phi 的 incumbent，材料标签从未参与选择。

|对象|DR−NNdirect NRMSE|DR−NNdirect Phi|DR/NNdirect CPU|DR−初始NN NRMSE|
|---|---:|---:|---:|---:|
|r8_test_0001|0.000944|3.172093e-09|1.239|-0.196625|
|r8_test_0004|0.000341|7.964805e-11|1.644|-0.050079|
|r8_test_0022|-0.223400|-0.1004565|4.549|-0.177459|
|r8_test_0027|-0.000475|-3.396401e-09|1.580|-0.143544|

投影本身的 NN 保留诊断（共同 corrector 之前）：

|对象|初始 NN NRMSE|DR最低Phi候选 NRMSE|CIE最低Phi候选 NRMSE|
|---|---:|---:|---:|
|r8_test_0001|0.273914|0.107522|0.205175|
|r8_test_0004|0.130735|0.074037|0.130735|
|r8_test_0022|0.243888|0.055806|0.073886|
|r8_test_0027|0.264002|0.166300|1.049422|

成本与解释边界：PREPROCESS_RECEIPT 给出实际 feature、NN加载和推断费用；各 RECEIPT 的 physical_cost 与 projection/auxiliary_cost 分开保留。OFFLINE_CPU_PREFIXES 按已付费原Phi候选给10/30/60/120CPU秒前缀，不按GT挑选。材料误差下降与物理残差下降分别报告。

该诊断只回答四个暴露对象上 NN 信息通过不同接口后的保留/损伤。未运行人为错误 prior 的新臂，故不能推出任意错误 prior 的可靠性。科学主张与 novelty 由主线程独立评审。

复现：`run_nn_interface.py prepare`（会重写协议/初值，仅新目录复现），随后 `run_nn_interface.py run`，整组 seal 后运行 `score_nn_interface.py`。原历史来源只读；无 push。
