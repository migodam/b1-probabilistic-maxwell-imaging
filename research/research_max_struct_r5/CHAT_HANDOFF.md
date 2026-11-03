# MAX-STRUCT-R5 给 Chat 的入口

先读 [中文报告](RESEARCH_REPORT_ZH.md)，再读 [数学推导](THEORY_DERIVATIONS.md)、[轨迹分析](TRAJECTORY_ANALYSIS.md) 和 [研究循环](RESEARCH_LOOP.md)。所有R5二维对象已公开，非新盲测。

## 最值得继承的结论

1. 旧弯条成功可由**普通global current + one-hop material readout**复制：NRMSE=.041480、IoU=.939759；仅改current失败（.558514）。本例收益定位于非局部材料读出，不需current传播。
2. 自然mf3阶段后，两种完整m/J/A交接均失败：global→material-only .388613；material-only→global .558738。指定协议中不是只需一个早期好初始化；也未证明全局盆地/所有切换点。
3. 简化仍未恢复42cb/63ab：NRMSE=.766323/.404823，环=.072525。63ab材料有所改善但物理残差更差，禁止GT选根。
4. 168个base-plus全域cross均正；几何固定块谱说明同相残差被加重、交替残差被放松。成功与失败共享此结构，因此没有可靠预测器。
5. 非根处 dR=−P C dm+(dP)r，metric导数不可省略。根处加权LS transport不是实际prior/box/非线性readout的完整模型。
6. 两体素真实三维vector-DDA给出极化反例：β=.2同相纵向weight1.34284，横向.85619。是离散6维全波机制例，不是三维成像。

## 证据与限制

8重放/168检查点+6历史态/120反事实；7次带同原物理corrector干预；1个矢量小模型。只有4个二维材料对象。原始数组、source snapshots、失败和各run receipts全保留。新增7717LU/183686RHS（不同系统维度分开），全部本地CPU。没有新照明收益、额外prior收益、NN或参数扫描。

R4的六对象4/6共同结果、32网格比较优势消失和进一步corrector负结果全部继承。R5九数字/新盲测/真实测量/3D成像NOT_RUN。历史九数字完整single/dual表见[附录](RESULTS/HISTORICAL_DIGITS_REFERENCE.md)，不是R5运行。

原Chat候选A/B/C：R4的[A/B/C完整答复](../research_max_struct_r4/ABC_VALIDATION_ZH.md)仍有效。R5深化A并定位作用位置；未把B/C重命名成新成功。

## 下一步最有价值的问题

不要再用cross正负、单独||transport||或相同零集给成像成功作保证。优先研究完整原物理residual-feedback在**合法材料切锥**中的协同方向，结合接收数据方向；先解释R4候选C的box投影破坏，再做一个失败例加控制。只在有限原Phi下降及独立确认有新证据后扩大。现有结果是机制定位/反例的论文素材，不是成熟“稳定胜过CIE”的算法。

标准乘积法则、Walsh对角化、dyadic纵横分解均不主张原创；原创性仍需围绕具体可检验机制及有效算法单独建立。

