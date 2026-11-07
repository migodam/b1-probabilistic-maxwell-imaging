# 新24例的离线评分合同

在未读取任何新真值数组、未进行评分时，由root要求冻结 [FREEZE_SCORING.json](FREEZE_SCORING.json)。此更正只处理三水平材料的评价失配，不改算法、先验、原Phi候选选择、预算或正在运行的队列。所有方法与全部24例使用共同规则；没有按新预测表现或真值结果选择阈值。

主连续指标保持原定幅度NRMSE：全576像素RMSE除以`max(chi_true)−.1`，使用已封存的原生成幅度metadata。主support指标统一为预测与GT都按严格`chi>1.0`分割的IoU。这定义要检测的对象材料，不是在线先验。主material_success为NRMSE≤.15、fixed_support_IoU≥.7且材料合法；合法性沿用实数、有限、576维及[.02,20] box，float64容差1e−12，确实量化的float32预算快照保留原单ULP容差。收敛、噪声相容及物理Phi另列，不并入这个材料success门。

历史Otsu必须保持独立标签：`legacy_Otsu_IoU`沿用预测自身256bin Otsu及历史GT半幅度support（`chi_true>.1+.5*amplitude`）；`legacy_success`沿用NRMSE/.7/合法的旧组合。它在mixed对象中可能只度量高材料部分，不能替代主support success。另记录`Otsu_full_support_IoU`：同一预测Otsu分割对全部`GT chi>1`support，用来明确诊断低材料被Otsu漏掉的情况。二者的GT分割定义不同，不能静默互换或合并。常量图像的旧Otsu分割为空，按既有定义处理。

真值区域RMSE只在所有终点封存且root明确授权后计算：true low区域为生成预定chi2.5，true high区域为生成预定chi12/18，使用固定1e−12等值容差；分别给原chi单位的RMSE与区域像素数，空区域记null。它们不能选择在线迭代、候选或方法参数。

原因来自三水平定义本身：背景.1、低材料2.5与高材料12/18不应由单一Otsu/半幅度分割解释为完整support，即便输入是精确GT也可能漏掉低水平。连续误差及共同固定support阈值可以明确评价全部对象。冻结本合同没有读取新GT、调用评分或新增正演。
