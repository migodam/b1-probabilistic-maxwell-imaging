# Codex → Pro：FBG-R3-E 证据交接

## 请先看实际结果，再推下一步理论

1. [九数字真实材料图](figures/fullwave_digits_1to9.png)：同一历史高对比576材料任务，成熟dual-beta CIE + 原数据corrector 10/12联合成功，数字1–9全部通过，ring/curved_bar仍失败。固定beta8/12，Native TRF0/12。这个正面结果属于成熟基线，不归为N1。
2. [native-affine开发新增案例](figures/e1_native_tsom_finite_gain_case.png)：同容量192实current DOF，四种finite current控制均恢复一个原强基线未达到的高对比双岛；7/12对6/12。不是task选模专属。
3. [第二组新确认结果](results/e1_native_confirmation_scores.json)：648次独立封存预测。task local多恢复一个clean curved16的材料与数据；primary joint仍6/12，noisy没有新增。finite task反而3/12，普通Fourier/random finite6/12。
4. [N2](results/e2_mechanism_summary.json)：33个合法subset roots，全部观测不相容；2次代数成功来自anchor fallback。没有全问题新根或576优越性。
5. [独立科学审阅](audit/SCIENTIFIC_REVIEW.md) 与 [费用](results/COST_LEDGER.csv) 包括失败、truth访问偏差和未测量成本。5748次唯一推断执行不等于5748个独立对象。

## 新鲜task的有限正面证据应精确理解

对象 obs_7ad98676df61dbb6，6×6、幅度16、uniform、clean。共同起点0和1下，task选择随机字典中的 random_native_3；实际current DOF192→192，各一个共同TRF seed。最终材料NRMSE约.000310、完整数据残差约8.9e-5，而固定Fourier约r=.135。它是真材料恢复，不能因联合门未过而说成仅辅助拟合。

但IoU=.625，而该GT自身Otsu与预声明occupancy mask的IoU也为.625；原joint门仍FAIL，未改门。对应noisy观测 obs_b3b3447fe12c27c7 未恢复。后续24次披露后组件诊断表明普通normal选择在一个起点、soft在另一个起点也能进入相同正确材料区域；random/energy未进入。这个对象不能证明A17分数的独占作用或跨对象统计优势。

一个成功起点的paid physical local gain为负（约−.00703），最后joint rollout + full corrector却恢复了GT；另一起点局部gain约+.00906。故不能用简单“局部预测正收益”作为排除全部探索的证书。这里不是全局反例定理，只是付费轨迹观察。

## 本轮改变了哪些认识

- 成熟CIE早停原型不是强基线；公开DEV相同beta/mask从22/23步延至305步，物理残差从.3266/.0782降到约4.27e-4/2.66e-4。
- 同总DOF不能代替相同DPCS/APCS自由度配置。初始N1三模式有明确数据地板；补全可见方向后仍只有很少dark方向，不能据此否定在合理TSOM空间中的normal机制。
- 新native分支固定jD与原生APCS，解决了这项混淆，但整体收益仍有限。
- 真实物理anchor的0/2/4阶不能直接套到非物理DPCS起点。CIE完整材料导数必须包含DP·e。
- Hessian downdate与gradient downdate都必须保留；同span dressing和非零模式缩放不是新自由度或soft release。
- 对纯法向电流，正负seed的材料readout二阶项相同、差仅O(t³)。这来自读出恒等式，不是每个混合N1seed的实际结局保证。
- 材料固定时，把辅助j精确投回Maxwell方程会得到同一个j(m)；有用探索必须通过材料变化实现，不能只把电流投影本身当branch生成器。
- 条件小核求得的远根可精确拟合少量通道，却完全不适合全数据corrector。outside误差是关键缺口，不是Q6/Q20恒等式错误。

## 理论、代码与来源

完整推导 [MATHEMATICAL_AUDIT.md](MATHEMATICAL_AUDIT.md)，代码表 [FORMULA_AND_CODE_REUSE.md](FORMULA_AND_CODE_REUSE.md)，源manifest [BASELINE_MANIFEST.json](audit/BASELINE_MANIFEST.json)。

R3-T完整报告作为输入保留。本轮A17仅提供材料任务相关性的理论动机；没有运行其原生Galerkin vector-current exchange后端。不要用两个Schur形式相似就建立A17独立因果链。CIE2019原则已按2D共享实材料适配，非字面3D实验复现；Fourier显式列没有CG-FFT复杂度主张。

## 下一轮最值得集中攻克的一个问题

**给定非物理current seed与其有限材料读出，能否用可支付、无GT的物理量预测它经过共同corrector后能否恢复，而不只预测冻结profile模型的瞬时收益？**

建议重点研究：

- 联合有限材料/current变化中的真实state-defect，而非只用一阶normal norm。
- 原物理恒等式 F(m)−Y=(Sj−Y)−SL(m)⁻¹e 如何变成有信息的、不过度保守的proposal评价；不要忽略resolvent放大。
- 同rank current replacement选择随机字典中某个方向，究竟来自任务分数的符号/曲率、可实现性，还是更一般的材料读出路径。
- 在共享实材料、box、source whitening和共同prior下，怎样保留合理的探索但防止current替错误材料解释数据。
- 把这种有限预测先在两个以上全新对象验证，再设计source-block profile实现扩到576。不要直接把巨型实化QR照搬上去，也不应把工程优化本身当成新成像定理。

保留最强成熟CIE、固定beta、原数据多起点和普通random/normal/Fourier current controls。需要同时保持总预算、actual real DOF、候选获取、有效seed及corrector成本，而不只比较迭代数。

## 当前论文状态

[论文候选](PAPER_CANDIDATE.md)：GO代码/基线/可复现实验，HOLD广泛新算法优越性。N1有值得保留的局部、开发与新鲜材料恢复证据；N2无本轮增量。没有全球收敛、完整根枚举、posterior、唯一GT或连续PDE结论。

不要让接下来的理论叙事抹掉成熟CIE已恢复九数字这个事实，也不要让它掩盖N1/N2仍未建立的独立贡献。

