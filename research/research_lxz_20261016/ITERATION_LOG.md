# LXZ 两轮研究与实现修正

## Round 1

冻结原 Joint P2 C proxy target、4 observations、4 seeds、6 arms，每臂4096次proxy likelihood（含64点初始化及geometry的8次差分）。新增物理0。普通local event MAE .012039、global .004364、hybrid .007243、tempered .016145、prior .011206。所有旧teacher、response-bank与训练费用保留来源。

## Round 2：nuisance-aware metric

预注册预测：未知gain可能使原response metric过度收缩可补偿方向。保持target、anchor、seed、步长规则，仅将proposal metric改成A^T[I−aa^T/(a^Ta+1/.15²)]A。likelihood始终是proper bounded-gain marginal，绝不以profile替代evidence。只比较geometry，其他基线保持原数组。

## 两轮共同的实现缺陷及修复

独立复核发现差分行按minus/plus成对排列，代码却用f[2j+1]−f[2j−1]，后者错误跨坐标。冻结R/B仍满足Gaussian prior reversibility，因此不改变MH不变target，但这两臂不能解释成Jacobian几何。原round1和round2全部数组、费用保持。修复为f[2j+1]−f[2j]后，仅各重跑geometry，目录增加_geometry_repair，不覆盖旧记录。这是同两项算法的实现修复，不引入第三种技术方案，不支付新增Maxwell。

协议原393216+65536为无失败时的proxy分配估计；用户没有给总proxy次数硬上限。本次两geometry修复额外131072次，全部明记。本轮用户的物理state/RHS/wall硬上限保持不变，不借用任何历史预算。

另一次写研究记录的零物理辅助命令因stdin编码失败，文件未写成功；不影响已冻结JSON/运行，新增physics及proxy均0。本文件保留该失败，不静默删除。

## E2 finite-chart失败、修正与最终核验

首35states/312RHS后continuation把Gaussian-CDF参数推到数值边界；原35状态数组因末尾才写盘而缺失，ledger不变。改为每个endpoint先持久化，链从10步缩为5步，并以.5 prior-coordinate trust半径限制candidate-only延拓。完成56states后仅剩5states，全部用于第一例导数h/h10独立核验；总96states/896RHS封账。没有第三轮算法、追加state或新prior。

E2两个material charts秩7/8；continuation的重复field明确去重，结构合格新world为0。Gain-aware geometry改善proxy混合但不胜简单globalMH，所以停止复杂混合与A17优越性主张；下一阶段缺的是合格inverse anchor和分离mode证据。
