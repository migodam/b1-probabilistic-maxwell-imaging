# R3基线与公平性

基线复用实际 `research_fbg_r3_experiments/src/cie_tsom_baseline.py::run_cie` 和 `full_maxwell_corrector.py::corrector`，不另造弱CIE。继承R2 native mf3/6/10各500，βcurrent1、后两段βmaterial.2，DPCS相对截断1e−3，共同原物理L-BFGS最多1000。历史基线结果直接引用，不重复收费。

强模态释放只在末阶段增17个strong-current方向，零系数起步，保留其余配置；这增加current容量/路径，Atlas已有先例，不叫新理论。clean对照固定原noisy源权重，仅用clean目标；原corrector归一化按clean目标，明确是信息反事实。

损耗路径从R2 native已校正端点起步，γ4、1、.25各最多250，γ0原corrector1000。阶段使用不同辅助物理介质，最终仅原γ0目标选候选。普通续跑1750、GlobalDual1000+corrector1000、同原目标×1e6/gtol同步缩放的L-BFGS1750分别控制额外工作、已知feedback及停止尺度。相同迭代允许量不是相同LU/RHS/wall；实际提前停机和全部费用均保留。

有限单元/交换/配对路径从同一scaled-original一阶驻点出发。单cell最多10个接受，实际1个极小移动；swap最多10，实际0；配对射线仅1个最好proposal，再同corrector1000（只有接受时执行）。每项选候选用相同完整Phi；current/root候选不读GT。无额外support/字体/decoder/材料prior。rank2swap暂时保留histogram，但最终corrector不受此限制。

CIE原文和近邻方法已有文献依据。只有有限路径的数值候选枚举新接入；没有证明其精确根覆盖带prior全部驻点。onecell公开控制中普通CIE readout本来已能恢复，不能拿它写成beyond-CIE胜利。

本轮无足够开发材料收益，未冻结新盲测算法或新生成确认对象。R2历史2/3没有变成本轮0/11或11对象统计；全部开发结果逐臂解释。
