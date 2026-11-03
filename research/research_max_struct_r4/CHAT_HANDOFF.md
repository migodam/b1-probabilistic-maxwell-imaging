# 给Chat / 后续研究Agent的R4入口

先读 [中文报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/RESEARCH_REPORT_ZH.md)，再读 [A/B/C答复](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/ABC_VALIDATION_ZH.md) 和 [完整推导](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/THEORY_DERIVATIONS.md)。网页与全部图在 [这里](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r4/)。

## 必须继承的结论

1. A是实际四cell Green的本构缺陷传播，不是近似物理forward。精确block可缩成一次交换 `R=(I+L)u`，整个公开box有解析可逆界。旧开发弯条RMSE/12约.04167、IoU.93976；Windows独立栈复验约.04160/.95181。
2. self-only、同幅π反相和成熟native8000在这个旧例失败，支持具体局部机制；但32网格的普通global也恢复，优势不是跨网格普遍事实。
3. 六个新鲜对象、18个冻结端点、同prior/同corrector：onehop/global/native8000全是4/6，相同两个新弯条都失败。onehop相对global NRMSE仅1胜5负；三候选最低Phi选择仍4/6，两个失败弯条都选到最差材料臂。
4. 失败后用四种分块平移平均构造8邻接传播，未调权重：两个已暴露开发例均失败，旧positive丢失；42cb NRMSE稍改善但数据残差更差。这不改原18和4/6。
5. 一例成熟native再给4000步物理额度，实际3361步达到原PG精度，仍没有恢复，NRMSE还变差。A在251步FTOL停滞，PG未合格。原1000步确认和追加费用分开，不称局部极小证书。
6. B纯port-only全current消元是原目标缩放；有限eta版本确实改变metric，但lambda1更低residual更差材料。C完整有限协同射线未逃逸，原因涉及box可行性和方向，而不是公式不精确。
7. 有用的读出判别：`D_M−C†=(CᵀMC)⁻¹CᵀM(I−Π)`。局部权重抵消取决于材料切空间不变性；非局部耦合一般改变读出，但改善没有保证。它是已有加权LS/R3取消关系的精化，非首创新基础定理。

## 证据入口

- [六对象原始分数](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/results/CONFIRMATION_SCORES.md)
- [独立确认与成本审计](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/audit/confirmation/CONFIRMATION_COMPARISON_AUDIT.md)
- [最后科学审查](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/audit/FINAL_SCIENTIFIC_REVIEW.md)
- [费用与维度口径](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/COST_LEDGER.md)
- [失败后如何推导/转向](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r4/RESEARCH_LOOP.md)
- [全部代码和数组目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_max_struct_r4)

不要把单个公开开发正例、最低residual、额外cost、更多观测或GT事后指标改写成新盲测优势。preconditioned CSI、NIE/CIE、T重整化、Neumann和variable projection都有直接先例；具体算法的普遍收益和论文原创性没有成立。3D新反演NOT_RUN，R1单调性证书缺口未被解决。

所有已评分GT现在都是公开开发资料；不得再次称独立blind。新的算法版本允许分析这些失败，但必须另立开发版本，并用新独立数据确认。
