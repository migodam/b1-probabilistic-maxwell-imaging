# 定向历史审计与开放问题

审计开始 2026-10-29 15:40:44 UTC，首轮于15:49 UTC完成，未重算历史实验。仓库 `migodam/b1-probabilistic-maxwell-imaging` 的本地与远端 main 均为 `9cd323cbc88af89c5dcbdef8b39049218f8dfcdf`。本轮仅新增 `research/research_z_alg_r2/`；同一工作树已有未跟踪 R9 工作，保持不动。

## 已建立与未建立

- **Z-ALG E1**：6个独立实材料、2源×1复接收，真实 Hankel–Green 标量TM圆盘pulse VIE、固定非对称标签。直接材料与完整电流图册通过相同147/147局部网格目标；统一保存成本窗口节时13.601%，未过20%门。合法二维局部纤维存在；没有全域覆盖、跨连通分支或一般物理对称群结论。历史参考材料公开，不能当新盲成像。
- **E2**：4材料、4实观测，33启动31合格，去重1根；2次失败。往返路径不是非平凡monodromy。未找到第二合法根不等于全域唯一。
- **早期Z-ALG**：两cell交换C2及商(p,q)为真实的特定对称正例；不可推广到非对称模型。
- **P3L/R5-B**：原24×24、Full66、8×32复通道。512新对象中蒸馏器相对Strong CNN MSE降低3.11%，是有统计支持的成像风险改善；至少80%增量保留仍HOLD。它没有建立材料解遍历；模块压缩与全系统速度必须分开。
- **R6**：低成本物理锚点在两噪声条件减少约33%原窗口时间，强普通对照也很强；更强完整基线比较未全部确认。真实同材料QD补偿属于gauge，不是新材料解。
- **R7**：37456 RHS；Full-J在材料RMS=.02达到8/8有限噪声竞争端点、7/8 pair≤1，高阶路线NO-GO。材料Jacobian和响应rank的数值现象不能变成统一精确纤维。
- **F/F3**：native材料坐标可减轻参数曲率障碍；有限高阶Explorer没有胜过普通材料回缩，少量接近端点不是新拓扑。改变先验、噪声或固定几何来救正结果不可接受。
- **A17X/F3审计**：7个参数可引起575个变化像素，精确端口rank可达575；7端口遗漏约58–64%算子谱范数。相同试验空间的A17与ordinary Galerkin等价。固定算子、Schur/SMW本身不新；普通基线也应获同样缓存。

## 复用而不重跑

1. E1/E2圆盘积分Green公式、精确材料Jacobian与接收伴随计费。
2. 固定材料变化支持的精确低秩有理表示；秩等于变化单元数，不等于decoder参数数。
3. E1已冻结几何用于新测量目标，旧参考仅用于公开代数根结构实验。
4. 公共先例审计入口：E1 `audit/prior_art/PUBLIC_CITATION_REVIEW.md`。Everett1996已有Helmholtz多项式化与homotopy多根搜索；CSI/SOM/TSOM已有电流与SVD分解。

## 真正开放问题

- 能否在不读取真值的情况下用有限代数根求解器找到切片全部合法候选，并保留失败路径？
- 能否利用一个接收、多个照明的共同材料结构，将一般multiaffine的阶乘根界降低？
- 固定Green和采集SVD何时直接给出精确非线性材料逆，何时因秩/小奇异值而失败？
- 代数准备、根生成、完整验收全部计价后，是否胜过动态Full-J或标准同伦？

本轮不把连续纤维采样点计成不同连通分支，不重开posterior积分，不按GT选择弱坐标、活动支持或初始化。

## 源文件

- `../research_z_alg_asymmetric_validation/RESEARCH_REPORT_ZH.md` 与 `src/physics.py`
- `../p3l_r5_b_independent_confirmation/RESEARCH_REPORT_ZH.md`
- `../p3l_r6_task_optimal_maxwell/RESEARCH_REPORT_ZH.md`
- `../p3l_r7_material_solution_fibers/RESEARCH_REPORT_ZH.md`
- `../research_f_intrinsic_transport/RESEARCH_REPORT_ZH.md`
- `../research_f3_quotient_exploration/RESEARCH_REPORT_ZH.md`
- `../a17x_f3_competition_audit/RESEARCH_REPORT_ZH.md`
