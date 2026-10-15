# 消融：区分有界修正、结构分层与历史反馈

本页只解释已实际运行的对照，不追加物理，不改变冻结确认集。所有候选均为candidate-only，概率null。

## 1. 同结构区间，改为独立重启

后确认、development-only归因固定solid与annular_cavity两个对象、步骤0/1/3/5，在同protected geometry区间从独立补偿变量初始化，然后使用相同corrector与原full-Mie verifier。8次新增物理已计费。

8个端点全部full-valid；相对亲本有限编辑端点的最大材料L2差2.83877e-5，远低于0.03重复阈值。因此，在这8个端点上，亲本历史不是取得该材料的必要条件。现有确认收益应主要归为**结构区间分层+普通有界补偿优化**，不能归为新relative feedback、材料null continuation或learned next-solution机制。该对照只覆盖有限development端点，不证明所有可能搜索问题都不需要历史。

证据：`configs/ATTRIBUTION_FROZEN.json`、`results/STRATIFICATION_ATTRIBUTION.json`。

## 2. 同proposal，拿掉材料补偿

|development对象/步骤|raw full residual|corrected full residual|raw→corrected相容性|
|---|---:|---:|---|
|solid/1|4.92287|4.89738|PASS→PASS|
|solid/4|17.15029|4.80420|FAIL→PASS|
|large_cavity/1|4.95662|4.91771|PASS→PASS|
|large_cavity/4|18.47129|4.86483|FAIL→PASS|
|annular_cavity/1|5.92349|5.89784|PASS→PASS|
|annular_cavity/4|7.35931|5.45607|FAIL→PASS|
|combined/1|5.13838|5.13059|PASS→PASS|
|combined/4|6.61117|4.45503|PASS→PASS|

两边使用同一个bounded real-gain profile，没有给补偿方法更多仪器自由度。8次raw full evaluation均计费；corrected endpoint复用历史已付费预测。

3个FAIL→PASS均出现在step4的20% global-prior步骤，其corrector没有protected geometry box。它们证实有界非线性材料修正有真实响应作用，**不能单凭此表证明“受保护结构编辑被补偿救回”**。4个step1的小保护编辑原本已经相容。full nonlinear secant数组及ε/gain变化见 `results/NO_COMPENSATION_ABLATION.json`。

## 3. 下一解历史策略

初版8步开发试验在warm-up结束前就停止，adaptive还没有真正进入forecast分支；第二次加入forecast后仍用8步warm-up，是保留的工程负结果。第三版在2步之后启用history+proxy forecast，并在确认前冻结，不能将前两版称为正面学习证据。

32parent确认中，adaptive相对静态Best4差+0.000839862，95%区间[-0.00219920,0.00378843]；coverage更低、平均搜索wall更高。利用历史成功率与失败反馈，并未显示更快找到下一种材料世界。保留20%全局探索不等于覆盖了全部材料分支。

## 4. 旧MDN proposal

冻结旧C1 full-distill模型、同一corrector、同一12个full端点budget，没有新增训练。Best4=0.137451，multi-start=0.136339；差值区间跨0。旧模型是已支付teacher上的proposal，不是当前full target的后验认证。没有证据支持加GRU、attention或训练更大模型。

## 5. 增加第二频率

固定前8个确认对象、每个base与材料最远候选；不按真值MAE挑对象。kR=2.4额外付费24次，含8次新测量生成与16次候选预测，跨频共享同一个gain。联合噪声不调大，阈值按64个实测量固定。

7/8 base与2/8最远候选保留相容；新增频率能排除一些单频候选。也有一个base在新噪声实现下被拒绝，必须保留，不能把所有候选剔除都称成功。实验没有比较主动测量选择策略，也没有给posterior odds，状态只是已运行的固定第二频率控制。

## 6. E2反证与机制边界

A30 finite-history的唯一额外去重材料在same-edit-no-comp中原样存在；结构mask只变1格，不满足≥6格预注册尺度，且材料MAE更差。有限编辑32条轨迹没有invalid→valid；multi-start和weak/null虽有2/4条转为相容，但归档后仍是base同一材料邻域。无新增结构、无Best4收益。

真实finite secant、Schur、dual和材料mass-adjoint恒等式通过，证明接口一致；它们没有自动带来更好的搜索或更快成像。当前没有ROM/TSOM/反馈专用的匹配成本消融，相关算法增量为NOT_RUN，不能把恒等式本身当原创或实证优势。
