# R6 成本账本

每次物理调用按实际task receipt计一次，失败也收费；嵌套timing不重复加。完整机器表为[AGGREGATE_COSTS](RESULTS/analysis/AGGREGATE_COSTS.json)、TASK_LEDGER.csv及ALL_CHARGED_ARM_STRATUM_COSTS.json。

| 项目 | LU | complex RHS | CPU秒 | wall口径 |
|---|---:|---:|---:|---|
| 主84个终态任务 | 157418 | 3773976 | 4290.925623 | task wall之和4372.222497秒，四worker重叠 |
| 保留的public loader失败 | 0 | 0 | .373810 | task .382559秒，外层process .497679秒；两者重叠 |
| FFT四个终态任务 | 6607 | 158376 | 176.356921 | task wall之和178.884822秒 |
| 六新观测生成 | 6 | 72 | .070935（内部窗） | 内部.119363秒；含import成功process .535279秒 |
| 生成前inventory失败 | 0 | 0 | 未记录 | process .692822秒 |
| FFT新action及单步接线 | 8 | 96 | .278204 | .275911秒，576维有限物理；不是完整成像比较 |
| 观察器不改变求解的有限接线 | 12 | 36 | .015157 | .019605秒；四个独立16维模型，不与576维LU成本等同 |
| 主独立评分含初次全图 | 0 | 0 | 1.300602 | 1.475121秒 |
| FFT独立评分含初次图 | 0 | 0 | .233892 | .236399秒 |
| 汇总分析及费用/轨迹图 | 0 | 0 | 2.596025 | 2.791970秒 |
| 最终plot-only排版修复 | 0 | 0 | 1.507805 | 1.528111秒；前两布局尝试1.233759及1.119171秒另保留 |

**新增已知物理总调用164051 LU、3932556 complex RHS**，其中12 LU/36 RHS是16维观察器小模型，其余为576维；不是把不同维度的LU数作为等算量。封存后GT可行目标见证只读取已存Yclean/材料，额外物理为0，其单独CPU窗未记录。静态审查/CLI、一次错误使用系统Python读取numpy、网络429/404与远端只读盘点均无新物理；未记录的时间写未知，不填零。

主队列两段elapsed 789.191126+351.179639＝1140.370765秒；两段之间修复间隔41.261817秒单列。FFT队列elapsed55.359759秒。它们不是整轮研究从阅读到发布的墙钟时间，也不能与各task wall之和相加。所有运行本地Apple Silicon CPU，四worker各单BLAS线程；GPU/CUDA执行0，远端只读观察已有作业，没有占用其他会话资源。

## 成本归因

主21对象A总CPU1504.322188秒，Dual-CIE932.812807秒，比例1.6127；A的辅助CPU1049.64秒、corrector444.91秒，后者减少不足以抵消前者增加。FFT两例A各节省15.7%/22.9%总CPU，但都增加了对应物理LU，且端点不同；只有阶段/总实测成本可据此比较，不能称逐步bitwise轨迹的纯速度测试。

auxiliary与corrector是顺序独立计时窗，task总量还含imports/setup/序列化。每个物理模型中的LU/source RHS/adjoint RHS、FFT各2D变换、block乘法/除法及LBFGS内评估都保留在原回执。未记录的Global/A逐步辅助operation_counts为null，不能伪造四臂完全同粒度的算术曲线。

所有拒绝的辅助端点、物理line-search评估都包含在task总成本。DPCS/Green公共数组复用的是历史产物，本轮不重复收费其历史组装，但本轮导入、每任务receiver/APCS分解、FFT kernel setup及重复工程运行均计入。没有训练成本，未启动NN。

主任务未量测进程峰值内存；FFT runner记录自身进程RSS，但窗口含dense物理模型和快照，不用它宣称全局内存优势。未知peak字段保持未知。真实货币/API模型费用没有账单接口，金额未知；本地CPU与0GPU不表示整轮科研免费。

HTML/ZIP、Git传输与网页读回的时间属于publication，不混入科学运行。它们的回执由`audit/EXPORT_RECEIPT.json`及`audit/PUBLICATION_RECEIPT.json`记录；发布完成不代表科学主张通过。
