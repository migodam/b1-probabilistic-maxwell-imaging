# Maxwell有限反馈、数值势区间与实际费用

## 有限状态更新的正确作用

所有候选使用同一材料与8照明。本构增量为实际ΔX，不是错误复相位或任意current。精确关系 Δj=(I−R0ΔXG)^{-1}R0ΔXE0 与原Full端点核验。Jacobian变化分为传播反馈、内部总场注入和decoder导数三个有序望远镜项；单独冻结后两项不正确。高阶恒等式本身不是多解存在性或加速证明。

active-row exact T-matrix 支持集通常180–256/256，含setup复用8.374秒、Full7.427秒，净成本失败。原记录保存reduced RHS、缓存内存与fallback，不将低latent维数误写为低端口秩。

## 决定驱动误差区间的真实迭代

首轮完整链与Full相同决定却更慢（355.79比229.78秒），已封存。修订只收紧可验证的逆范数包络，使用min(Frobenius,sqrt(1范数×∞范数))，并优化残差动作，不改变统计目标。

独立锁定8对象×1024步结果：Full65,536RHS/71.397秒；Richardson4,608/24.153；同缓存globalGMRES14,691/37.605；固定精度GMRES20,759/46.832。总105,594RHS全部计费；同一innovation、U、current-cache权利与落盘合同。所有接受、拒绝、状态路径及势包络重放零差异。

这支持约66.17% wall节省，范围是小步长beta=.002、当前TM24和缓存策略。普通Krylov也是强正对照，不能称A17独占。更大beta的边界实验：.005/.01/.02的3072次提案中接受302/95/6，实际移动快速变小，fallback增多。相同坏混合核的加速不等于得到合格posterior更快。

包络在普通float64及已验证离散后端条件下成立，没有定向舍入或连续PDE稳定常数证书。因此称“核验后的数值决定区间”，不称连续Maxwell certified posterior。

现有缓存取得需每parent576列inverse actions；setup/rebase/低阶线性求解/残差/fallback全部保留。未观测到稀疏active-set净优势，不隐藏failed runs。详见evidence_summary/DECISION_CHAIN_COMPARISON.json及原closed收据。
