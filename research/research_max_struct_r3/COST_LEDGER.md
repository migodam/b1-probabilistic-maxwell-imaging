# R3成本账本

物理工作合计 **4,092次LU、101,104个复RHS**，来自17个实际外层数值回执，矩阵维数不同，不是FLOP总量。11个高维开发控制之外包括固定材料诊断、1次新被动adapter接线、噪声谱、fiber、单圆盘及3D范围实验。

|类别|复RHS|
|---|---:|
|Forward source|48,922|
|Adjoint|48,656|
|Receiver transpose|70|
|完整inverse作用|3,456|

额外小系统分开：H2有36个576维current Cholesky/36RHS；H1有3个17维Gram Cholesky/36RHS，噪声谱另1个17维Cholesky/29RHS和一次17维eigenproblem；GlobalDual有100个12维Gram factor/4800RHS。邻接swap195个2×2solve/2340小RHS，全部swap45719个2×2solve/91438小RHS。解析配对射线直接算2×2transfer，没有把每个根算成完整LU。CIE的receiver/APCS SVD、显式Fourier列、PR/readout/内层求解在原回执中单列；不是FFT复杂度声明。

物理数值进程wall之和 **36.410166s**，CPU之和 **35.495631s**；加已计读证、评分及绘图后为 **38.087146s wall / 36.827485s已知CPU**。这是进程时间总和，不是并行历时或全部科研耗时。导入、未计时失败、工具/模型推理、写作和部分序列化开销未知；费用金额null。没有GPU/云/NN训练。最终export、发布访问回执另计，避免递归打包时间重复。

独立逐行账本：[FINAL_DEVELOPMENT_COST_LEDGER.csv](audit/FINAL_DEVELOPMENT_COST_LEDGER.csv)；主线程补计机制图后的汇总：[FINAL_COST_TOTALS](audit/FINAL_COST_TOTALS.json)。旧six-arm/中途ledger作为快照保留，不与final相加。原R1/R2历史成本没有再次收费。

source/setup/online/validation/scoring按实际回执范围而非人为分摊：开发inverse每run包含其载入/准备/候选验证/corrector和序列化；setup内嵌timer不重复相加；diagnostic与新接线单列；所有评分/图只读缓存、不增加物理RHS。setup和online逐操作精确互斥拆分并非所有runner都有，缺失项保留unknown。完整峰值内存未统一测量；各记录的process RSS仅作局部数值。
