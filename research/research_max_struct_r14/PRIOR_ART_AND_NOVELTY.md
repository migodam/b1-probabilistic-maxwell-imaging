# 先例、新颖性与论文定位

本轮是对既有A+真实恢复增量的结构归因，不把“两个块换成不同度量”本身当作已有文献之外的新发明。机械文献审查、精确条目与检索失败见`audit/prior_art/`。ScholarQA两次检索遇429/无候选，随后核对原始论文与公开原文；这不构成穷尽全球文献的证明。下载的第三方全文不随本项目发布。

## 最接近的既有机制

1. **CIE/NIE/CSI与source归一化。** Zhong等NIE（2016，DOI10.1109/TAP.2016.2535492）和Zhong/Xu CIE-I（2019，DOI10.3390/jimaging5020027）已有材料相关状态残差重写、内部场及辅助current搜索。本轮保留这些成熟骨架，Möbius/Arc属于R7已做成果，不重计原创。
2. **CC-CSI及非单一辅助下降。** Sun等Cross-Correlated CSI（IEEE TAP2017，DOI10.1109/TAP.2017.2673758；公开稿arXiv1906.10864）已讨论传播/观测与构成误差、材料依赖归一化、CG记忆。不能声称首次发现辅助误差和物理恢复不同，或首次借跨块耦合改路径。
3. **Relaxed physics、variable projection及Gram。** van Leeuwen/Yang的Relaxed PDE-Constrained Inverse Problems（2024预印本arXiv2403.15292；2025期刊DOI10.1088/1361-6420/adab89）相关于状态消元、Gram加权及材料依赖结构；固定Gram与variable Gram的导数不应混同。本轮没有把省略metric导数包装成同一目标。
4. **交替最小化/非精确下降。** Attouch–Bolte–Svaiter（Math Programming2013，DOI10.1007/s10107-011-0484-9）等收敛理论需要共同下降、相对误差、正则性等假设。01/10的两块不共享同一个势，不能直接援引这些定理；00/11也还需证明真实有限读出满足其条件。
5. **非势动力学与CG。** Differential games的对称/反对称结构和Hager–Zhang的非线性CG先例解释了为何应审查向量场与历史；同时梯度games的结论不自动适用于顺序PR+材料读出。一步精确线搜索的方向非线性和浮点停止deadband并非新数学发现。

## 本轮可检验的增量

- 用真实非周期FFT/DPCS程序建立了字面可重入的完整state接口，保持原A+更新而不是另造弱基线。
- 明确β与material相关metric的函数身份，使四组合真正分离一致/不一致；额外self-only、global-Dual、共享材料求解器和次序实验排除常见混杂。
- 将理想根谱、有限current解误差、真实PR记忆/guard、同state短期响应、长程材料恢复作为五层不同证据。一个层次不支持结论时保留反例，不强行用根谱解释恢复。
- 通过同对象/同数据/公开选根/完整计费的实验判定不一致何时有利、何时有害；新增恢复和丢失恢复分别列出。新24独立预测与曝光48分开，三维另做几何对应而非改变量名。

这些是否构成充分算法新颖性，须由最终独立结果、机制可预测程度与文献进一步审查决定。01→11是无新增阈值的简单待验stage控制，若不能新增恢复就只保留负结果，不将其命名成“统一全局框架”。

## 论文应如何定位

无需所有case全胜才有论文价值。一个明确困难子类上的可重复增量、正确机制与公平成本可以支持算法论文；但应避免只展示开发正例、用更低Φ替代材料质量、或忽略失败和物理误差。

本轮最自然的候选叙事是“auxiliary metric consistency versus finite recovery dynamics in Maxwell inversion”。强SIAM JIS需要比理想块GS重述更有内容的真实有限程序命题与有辨别力的独立证据；强TAP还需要更可靠的跨离散/三维/实测物理与实际恢复优势；Inverse Problems通常需要新的可辨识性、稳定性或连续理论，而本轮有限矩阵/程序观察不足以自动满足。最终资格判断见中文报告，不能由期望期刊反推实验结论。

原始来源链接：[CIE-I](https://doi.org/10.3390/jimaging5020027)、[NIE](https://doi.org/10.1109/TAP.2016.2535492)、[CC-CSI](https://doi.org/10.1109/TAP.2017.2673758)、[relaxed inverse problems](https://doi.org/10.1088/1361-6420/adab89)、[非精确下降理论](https://doi.org/10.1007/s10107-011-0484-9)。所有精确作者/年份/页码核对记录在本轮citation ledger。
