# 搜索历史的预测信号：E0 已执行诊断

结论：旧径向历史对未来失败/材料重复有可检测的关联，但对严格“下一份结构不同的 full-valid 材料”没有稳定额外预测增益。不能据此宣称强化学习会发现更多解。当前没有训练 PPO/GRPO 的证据。

分析来自冻结 Research A `CONFIRM_CAMPAIGN.json`，32 个独立材料 parent、5 方法、各12步，共1920条实际付费端点；原确认集现已曝光。本轮只读缓存，新增物理与 RHS 均为0。这里的径向 Mie target 不与 A30 24×24 scalar-TM 混合。

## 标签与防泄漏

主标签先要求本次 corrected endpoint full-valid，再与该方法该 parent 的全部严格前缀 full-valid ledger（含历史 duplicates 和合法 base）比较：结构新颖采用物理半径 δr≥0.05R 或环半宽 δh≥0.025R；材料新颖采用真实径向壳体积积分 L2/|1.5+0.1i|≥0.03。它们不是 latent 距离，也不是孔洞类别。原 archive greedy representative 标签只保留为 secondary，不能与全部前缀 ledger 标签混称。

A 使用实际测量32通道、噪声、当前 parent material、raw连续proposal、原方法/动作/保护box、剩余端点预算，以及原运行在 correction 前已计算的 raw proxy residual。B 加严格前缀的成功/失败/重复、最近两步与同操作摘要。C 将 B 的历史在观测相近 parent 间作无自配对置换，保持 arm/step，在训练/验证分区内分别置换；D 仅用历史摘要。C 是负控制，不是可部署在线策略。

5折按 parent 分组；所有32个 parent 的全部方法/步骤同时进入一侧。归一化仅在训练折拟合。固定 L2=1 的小 logistic，无超参数搜索。当前/未来 validation、corrector nfev、convergence、corrected residual 均未作为当前特征；脚本包含当前结果扰动不改变当前特征的断言。真值/离线MAE未读取。

## Cross-fitted 预测结果

|标签|模型|log loss|Brier|ROC-AUC|ECE(10 bins)|
|---|---|---:|---:|---:|---:|
|valid|intercept|0.151492|0.031718|0.267961|0.000513|
|valid|A|0.271590|0.048222|0.773959|0.049984|
|valid|B|0.189069|0.037918|0.877018|0.034335|
|valid|C|0.254864|0.044107|0.790800|0.045776|
|valid|D|0.076275|0.020383|0.940692|0.005972|
|corrected_material_novelty|intercept|0.542320|0.178380|0.463847|0.000116|
|corrected_material_novelty|A|0.538572|0.164124|0.770596|0.102834|
|corrected_material_novelty|B|0.480053|0.136864|0.826959|0.082931|
|corrected_material_novelty|C|0.554417|0.164652|0.772272|0.104669|
|corrected_material_novelty|D|0.366241|0.112554|0.865033|0.033963|
|corrected_structural_novelty|intercept|0.457933|0.141744|0.445806|0.000370|
|corrected_structural_novelty|A|0.318084|0.097752|0.878991|0.034318|
|corrected_structural_novelty|B|0.321869|0.096589|0.897186|0.053637|
|corrected_structural_novelty|C|0.318379|0.097582|0.882302|0.025632|
|corrected_structural_novelty|D|0.275499|0.089193|0.903075|0.027092|

下表为 B−A / B−C 的 parent 平均 log-loss 差，负值更好。95%区间是32个parent的2000次描述性配对 bootstrap，固定 cross-fitted predictions，不在每次重采样中重训；不作确认性政策收益检验。

|主标签|B−A [95%]|B−C [95%]|
|---|---|---|
|valid|-0.082520 [-0.135691, -0.040901]|-0.065795 [-0.117698, -0.028359]|
|corrected_material_novelty|-0.058520 [-0.113384, -0.007399]|-0.074365 [-0.144784, -0.009125]|
|corrected_structural_novelty|+0.003786 [-0.031108, +0.040507]|+0.003490 [-0.033365, +0.040816]|

结构新颖的 B−A 为 +0.003786，区间跨0，且比匹配历史置换的C也没有优势。材料新颖与合法性虽然有历史关联，B 的合法性 log loss 仍劣于仅用训练折频率的 intercept（0.189069 vs 0.151492）。History-only D 更稳，提示高维 observation/action 线性模型在32个 parent 上有泛化/校准限制，而不是证明历史不包含信息。D中成功率与累计结构数也可直接编码当前搜索处于哪种耗竭阶段；预测阶段性重复不等于能选择更好的反事实proposal。

## 旧 adaptive 为何未超越 static

|原方法|full-valid /384|重复|新增物理材料(原archive)|全prefix材料新颖|全prefix结构新颖|原archive结构新增|proxy调用|
|---|---:|---:|---:|---:|---:|---:|---:|
|finite_static|349|188|161|161|141|128|22259|
|multistart|384|334|50|48|24|24|39717|
|population|384|343|41|39|20|20|37162|
|mdn|384|338|46|43|19|19|36878|
|adaptive|357|203|154|154|123|122|84433|

这些是每parent每method顺序中的 discovery events，不是把不同parent/方法的材料合并后得到的全局distinct world数。

旧 adaptive 仍在预设八个径向编辑box中用 UCB 和额外proxy forecast选择；每五步第5步固定global，revision3仅先探索两个box后就利用历史。相同384个endpoint下，adaptive的proxy调用84,433，是static22,259的3.793倍；搜索wall 23.83s vs16.38s。原archive新增材料154 vs161，结构新增122 vs128。全prefix几何新颖123 vs141，未显示新的结构发现优势。

320个受保护提案中，adaptive重复已用region131次，static96次；平均覆盖5.906 vs7个编辑region。这符合重复利用和可行性forecast偏向已知盆地的风险，但不能单靠日志确定因果。static原规则因global step覆盖t=4而根本未访问edit_id4；它也不是均衡、topology-agnostic的G baseline。原box控制不可直接移植为按孔数分类的新策略。

纠错后的成功/重复已经用于旧UCB；可预测失败并不会自动产生“怎样走到新结构”的动作。没有日志propensity、随机策略重叠与反事实corrector端点，因此不做off-policy policy value，也不把重排同一候选池称作生成增益。G应比较共同动作库、共同corrector和相同完整物理预算；当前仅支持轻量历史策略的有界实测，不支持预先宣称RL有益。

## A30/C1适用性与边界

A30历史E2只有4个material parent、8个weak/ring条件，新增结构族为0；32条有限编辑轨迹中 invalid→valid 为0。保护格实际未改变，不能说corrector把孔填回去了。这里不拟合“新结构”分类器、不与径向样本混合、不把steps当独立对象。新目标的历史信号必须由G/T实际full-state轨迹另行积累。旧C1 frozen-MDN作为此1920行中的mdn arm已包含；这不是把C1的四维径向网络当成A30图像生成器。

## 实施问题与成本

第一次无阻尼Newton在合法性极不平衡标签上出现数值饱和，不能作为有效比较；完整初次输出保留在 `results/history_attempt_v1_undamped_solver/`。仅增加单调backtracking修复优化器。根任务同时要求以全部prefix full-valid ledger作主标签，因此在新标签CV前记录v2冻结说明；原archive标签仍保留secondary。没有根据新标签表现挑选阈值、模型或case。

两次受控分析合计 wall 0.802152s、CPU 0.771077s；BLAS=1、nice+10，GPU=0，新增Maxwell=0、RHS=0。详见 `results/E0_HISTORY_JOB.json`，主进程统一记账。

建议：暂缓 PPO/GRPO。优先积累共同corrector下的全prefix结构发现轨迹，用轻量CEM/bandit实测是否优于连续尺度static和multistart；预测关联与实际新增结构收益分别报告。
