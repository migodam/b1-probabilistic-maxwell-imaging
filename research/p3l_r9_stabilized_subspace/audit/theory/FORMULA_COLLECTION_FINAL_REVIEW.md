# 公式集独立理论一致性复核

2026-10-30。范围：本审阅者负责的02/06/07章与论文核心说明，第05章关键误差公式，以及各章对历史实验/最终E5合同的身份区分。本文件不表示对所有历史章节逐行重新证明，不替主线程作最终科学结论。

## 数学结论

1. **固定载体与非线性状态分开。** 完整host/ROI消元保留非均匀Full66材料支持；field/current矩阵顺序、非零host观测项、8源加32receiver导数RHS均在02章一致。精确support消元与低秩近似明确分开。
2. **投影与材料合法性分开。** SOM/TSOM的投影像不是一般的空间交集；接收暗空间不必在非线性内部散射下封闭；多源电流必须对应共同decoder材料。没有将receiver-null电流直接当作合法材料样本。
3. **stable solve与受限探索分开。** 全列QR GN、prior metric、仿射约束与profile消元的条件完整。固定Q本身不减少完整Jacobian物理费用；profile不是概率积分。
4. **误差公式条件明确。** 已向root核对05章真实residual符号、resolvent误差、缺陷校正、dual remainder和likelihood区间；需要可逆性/稳定性上界的地方不能降格成浮点残差即可认证。复投影与实512维residual的维数记号已提醒并由root处理。
5. **E5保持weak line，不保证正向weak ray。** 最终协议允许signed弱幅度；B固定弱直线、C逐步当前frame弱直线，signed material progress不等于弱系数正号。历史`*_RAY`ID保留，但数学描述已更正。历史REV1只在Qs步后作正向radial缩放，与E5不同，不改写历史算法。
6. **候选覆盖与posterior分开。** 合法云、ε管状覆盖、范围、oracle容量、部署动作和启发式权重使用不同定义。准确posterior是可选升级，不阻塞有限可能性云的主检验，也不由云成功自动推出。
7. **新增F-VOL-015条件完整。** 非空有限云、非负归一化权重、PSD任务metric、任务Lipschitz条件和近GT权重a>0给出加权动作误差界。证明明确写出组内/组间方差，处理a=1与a=0。a依赖未知GT，该式是分析关系，不是部署在线证书；硬事件和空云另行处理。

## 最终协议身份

当前E5依据 `config/possibility_volume_protocol.json` (`approved_by_root=true`) 与 `config/possibility_statistics_frozen.json`；DRAFT仅保存历史。此次只读最终配置，不修改任何配置或代码。

- 32个fresh parent、三臂、每臂K16非anchor槽位、最多17个主云点、四个半径、绝对95%噪声门及RHS/LU上限均已冻结，不再称待批准设计。
- anchor只在合法且噪声合格时进入主合格云；无条件anchor baseline和空云动作fallback另记。存储allselected+flags不等于已经过滤的合格云。
- 空云保留parent，coverage=false，材料最小误差记未定义/无穷右删失；不作complete-case主对比，不追加有限失败罚值。
- radial bracketing失败仅表示既定搜索/cap下未找到，不证明非单调decoder的全局不可行。
- 平滑tether、满秩混合支持、任务导向槽位分配保留为 `NOT_IMPLEMENTED`；不属于当前E5，也不因公式出现而改变其评分合同。
- 03与workflow章节的历史NOT_RUN/未来冻结措辞由对应作者更新；05的历史未运行范围由root更新。不会用E5后续结果回填旧门。

## 保留的科学限制

B/C和R的corrector几何不同，当前对照是整体策略比较；不能声称纯proposal因果消融。三臂共享ROI，不能用其差异证明ROI本身的收益。32个parent的总体结论需要区间；候选数增加、权重为正或prior支持内输出均不能替代覆盖证明。RTO/rMAP/RML/EKI具有直接先例，因此“多个可能性”不是独立新颖性依据。

**独立审阅结论：** 在已列假设与限定主张下，未发现阻止该有限候选云研究被解释为一项清楚实验的致命代数矛盾。它仍可能在GT/task覆盖、合格率、范围或总费用上失败；这种失败必须由正式封存结果报告，不能被理论叙事覆盖。

## 本次操作边界

文档编辑与静态读取；新增物理RHS=0，新增数值测试=0；未读取E5 GT、候选结果或评分；未改变冻结代码、配置、预算或历史结果。最终完成状态与科学主张只能由主线程依据封存E5证据确定。
