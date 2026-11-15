# 双度量主对照合同 v1（首次实验前）

## 物理、度量与分母

原公开 2D scalar-TM 有限 disk-pulse VIE，exp(-iωt)，J 为相对电介质 contrast current，Q=Einc+GJ，m=εr−1 实标量共享所有照明。全部576材料未知量，原box [.02,20]、原smooth/TV先验及source归一化不变。源/接收在未知域外。

设 P0(m;b)=diag(b/(1+bm))。原阶段 current b0=1，material b1=1,.2,.2,.2。P1(m;b1)=(I+diag(m/d)Goff)diag(b1/d)，d=1+(b1−Gii)m。Goff 仅原不重叠2×2四cell块、对角为0，顺序TL,TR,BL,BR；不替换成全域G。M0=P0*P0，M1=P1*P1。

|ID|current|material|严格语义|
|---|---|---|---|
|00|P0(b0)|P0(b0)|同一个材料相关scalar残差目标|
|01|P0(b0)|P1(b1)|原global A+-Arc，非一致|
|10|P1(b1)|P0(b0)|反向非一致|
|11|P1(b1)|P1(b1)|同一个onehop残差目标|

所有current PR固定m时的梯度为 U*{S*(D Ωd)+ρ Lm* Pc*Pc u Ωs}，Lm=I−diag(m)G。精确二次线搜索也用同一Pc。material readout使用Pm(m)的完整导数，原Möbius/Arc两扫描、MM先验接受；00/10普通readout使用b0。PR旧梯度/旧方向保留原规则，stage边界重置。已有sourcewise CIE/Dual、global-Dual、MM、G19、DR单列参考，不将00冒充原Dual。

11的目标只在单一固定stage内相同；阶段改b1不是一个全程固定势。01相对11的计算费用与迭代数可不同，必须计费。M1正定性只能在实际合法状态上核对，不由命名假定。

## 选择与GT

主选择：所有已支付且预算内、合法、完成原Maxwell审计的候选中最小Phi=D+R。secondary：合法CN99球内最小原R；没有则primary回退。delta继承公开CN复噪声约定；跨网格只作模型失配压力，不称噪声证书。raw、corrected、primary均保存。历史NRMSE≤.15且IoU≥.70只作既有事后恢复统计，连续误差与整体结构另报，不作停止门或在线控制。

48对象为曝光开发/机制后续；不读其GT调整规则。所有新预测先封存，再调用独立离线评分。离线成功并集不是在线算法。

## 预算与状态

主2D比较继承120秒整进程CPU cap、72秒proposal、随后到120秒共同原物理corrector、每2整进程CPU秒首次完成readout时审计，以及原stage审计。所有初始化/SVD/分解/FFT/失败和超限计费；超限完成候选不用于预算内结果。不同host实际CPU不得直接合并成速度结论。

checkpoint保存m、A、J/Q、DPCS/receiver基、source weights、旧precurrent gradient/direction、stage/iteration/global counters、所有stage flags和readout flags。同状态干预保持完全相同state再切换目标；另做PR reset控制。恢复试验与诊断重放分别收费。

## 扩展策略

先执行代表性正反例、原共同五失败和同状态诊断；四组合随后覆盖完整既有48，不以结果好坏删对象。不靠扫β/幅度/额外Green阶数取得正例。可由机制导出一个简单候选后冻结，交给新独立多样材料确认；若机制不足，继续有辨别力的反例与实际物理实验，明确未实现的主张。

## 评分前执行审计更正

2D新runner实际corrector为绝对总120上限，早结束proposal的余额会给corrector；原书面最多48与实现不同。见FAILURES第3项。主四臂都使用同一实际策略，保留原结果；与旧R13历史表分开，并加付费前缀分析。独立24全10方法统一采用该总120规则，旧C2只改deadline一行，数值算法/selector不变。不得将额外corrector时间视为metric机制收益。3D另保留统一严格48口径。
