# R7 A+ 独立开发总结

状态：20/20完整运行已完成；4个既定公开开发案例全部保留；GT与R7 primary结果访问均为NONE。最终工程候选为Möbius状态驻点二次方程提案，加原先验MM接受与继承梯度阈值。原A仍保留为独立控制。

固定参数：`aplus_mode=mobius, aplus_arc_candidates=True, aplus_inherited_gtol=True, aplus_sweeps=2, aplus_step_tol=1e-8`；继承`material_inner_gtol=1e-8`。每次运行沿用3000步辅助日程与原3000次迭代corrector额度，一CPU进程、单BLAS线程。

## 方法与边界

精确保留原固定J/Q的onehop A目标、固定源度量、实材料box及原smooth/TV先验；实际差分算子除以19.5，先验按全局边数归一化。最终候选解精确条件状态驻点二次方程，把合法根、端点、旧点代入包含原先验的完整MM surrogate，选最低者；这是一种近似正则化条件提案，不是精确quintic最小化。原A精确loss最后复核并保留fallback。辅助循环与物理corrector仍使用原Maxwell模型。

## 完整开发成本

| 候选 | 运行数 | 完整CPU s | 辅助CPU s | 物理LU | 物理复RHS | 留存材料块loss上升 |
|---|---:|---:|---:|---:|---:|---:|
| original_a | 4 | 284.209923 | 193.722288 | 6967 | 167016 | 0 |
| mobius | 4 | 347.436909 | 210.205427 | 7021 | 168312 | 0 |
| projected | 4 | 323.025623 | 188.671858 | 7511 | 180072 | 0 |
| mobius_guarded | 4 | 268.020790 | 169.392407 | 7014 | 168144 | 0 |
| mobius_arc_candidates | 4 | 197.781556 | 104.944029 | 6709 | 160824 | 0 |

20次完整运行合计付费CPU 1420.474801s。固定current检查与报告整理是额外开发开销，未混入这些完整运行比值。

最终arc四例合计CPU 197.781556s，对原A 284.209923s减少30.41%；7b8d单例完整CPU仍增加约4.8%。这是本四例的实测付费成本，非通用速度定理。共享机器状态会影响CPU时间；quintic特征值求解被实际移除，才是结构性工作量变化。

## 全部20次结果与负证据

详表见[FULL_TRAJECTORY_COMPARISON.md](FULL_TRAJECTORY_COMPARISON.md)，机器可读行见[ANALYSIS.json](ANALYSIS.json)，每行含完整原Phi、物理相对残差、CPU/LU/RHS、迭代数、fallback与轨迹路径。

- projected8在b20d的原Phi从6.4736e-6恶化至0.0058324，约901倍；即使材料块每步留存loss非增，也不能保持全程物理状态，因此拒绝该替代。
- unguarded quintic两sweep保留四例公开物理状态，但总体更慢；加入继承gtol后，实际减少的quintic根很少，不能把所有CPU变化归因于早停。
- 最终arc在两旧bar保持同一公开Phi数量级与残差状态；easy也维持原状态；slit-shell仍在约0.087物理残差的失败/困难状态，没有改称成功。
- 是否真正保留两旧bar材料恢复，须由主线程预测封存后的独立GT评分判断。这里没有以GT选择候选，也没有读R7 primary结果。

## 检查与更正

新实现接线检查PASS：30个固定current readout均有限、在box内、留存原loss非增；6个stored-state梯度比较最大绝对误差6.05e-18，24个逐坐标残差更新最大相对误差3.90e-15；这些检查物理LU/RHS均为0。完整20次共60000个辅助材料readout，留存loss上升计数为0。这不构成原物理目标单调或恢复保证。

更正记录：[DEVELOPMENT_AUDIT.md](DEVELOPMENT_AUDIT.md)。最初7b8d文件ID后缀写错导致首次检查只覆盖b20d，随后修正并补齐；worker曾擅加1%相对Phi自动gate，主线程指出无物理依据后已从脚本/协议移除，四个既定案例全部完成。无样本因小幅结果差异省略。

原12运行、guard4、最终arc4分别保留manifest、结果和源快照。没有覆盖原R1-R6科学源码或历史结果。

工程冻结记录：[APLUS_FREEZE.json](APLUS_FREEZE.json)。最终科学解释、噪声合同资格与材料恢复结论保留给主线程。
