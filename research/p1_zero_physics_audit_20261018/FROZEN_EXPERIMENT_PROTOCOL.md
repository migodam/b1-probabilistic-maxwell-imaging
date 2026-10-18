# 两项最小真实实验 v1：已冻结，未授权、未运行

共同target仍为P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018。prior、decoder、固定sigma=.00119211798514519、gain1、mesh、source/receiver、数据packing全部不变。未来任何改动另建协议，不能与本协议混合结果。

## A：合法finite-noise移动，不预设结构模式

从已付费合格state_0024和缓存完整J开始。方向固定为在线冻结的三个prior-aware material weak directions；只用其进入损耗上界可行切锥的一侧。实际任务RMS=.02，不拿参数距离替代；.04已见不良的decoder-only压力结果，本最小试验不支付它。

两臂：普通有界Full-J proposal与普通material-native retraction。前者最小化anchor处原posterior二次模型；后者最小化材料目标误差加1e−8 latent tie-breaker。共同使用真实decoder、同一半径/进度约束、Gaussian位移≤1与原property bounds；这只是局部搜索约束，不改变prior支持或宣称posterior截断。

每组最多两次新full-forward：先计算合法proposal及J，再允许一次共同Full-J corrector，保持实际RMS与有符号进度，不允许回到anchor算成功。首proposal8RHS、其receiver共享J32RHS、修正后full8RHS，共2states/48RHS。3方向×2臂总12states/288RHS，受控wall300s。全部失败和拒绝在此上限内，未使用额不转移。

报告三个分开的结局：原noise-ball+backward residual+prior legality；DeltaPhi≤2的likelihood-competitive启发式；原P1结构变化检查。DeltaPhi≤2不称95%区域，density不称event mass。实际RMS匹配容差2%，但task-floor仍为.02（仅1e−10浮点容差）；progress≥.01。候选按固定arm/direction顺序、pairwise材料RMS≥.02去重；GT不参与选择。

终端J不是本任务资格要求，因此不为A支付额外stationarity求解；若声称驻点则需另计费用，不能从旧缓存虚构。没有posterior采样/模式概率/credible interval目标。

## B：修复独立物理反演的收敛机制

两臂均从原独立physics-only B末端state_0048及其已付费J开始，不从GT/N3或新随机物体重启。强baseline为原有界LM+固定.5 Gaussian postscale继续；新臂为直接有界trust-region二次子问题。

目标仍为Phi+.5||z64||²。新臂初始Gaussian半径.5，范围[.03125,2]，property步长∞范数≤.2，保原property bounds；rho>.75且边界步时半径×2，rho<.25或拒绝时÷2。每个实际目标下降才接受。局部子问题直接处理半径/边界，不能先求长步再统一缩所有坐标。两臂精度相同。

每臂最多8次新full-material state attempts（含line-search拒绝/异常），至多8次共享receiver J（包括终端梯度/KKT记录），因此8×8+8×32=320RHS。总16states/640RHS，受控wall600s。状态/梯度/J/全成本立即持久化，早停节省不追加到其他实验。

资格保持原prior合法、白化残差≤23.785449004057032、状态backward residual≤1e−10。stationarity单独报告，不再混入可行anchor资格。失败时保留失败，不解释为先验中不存在解。

## 缓存、授权与A17-X

两个实验合计建议28states/928RHS/900s CPU，NN/GPU/Mie/新3D均0；**尚未获授权**。本轮新增physics预算仍0，旧26/976及余6/48均封账。共享缓存可以降低增量费用，但从零部署必须披露旧反演、J和背景setup的费用。

本协议不包含A17-X额外物理费用。A17-X未来作为独立finite-state/current enrichment轴，固定合法材料后比较full、exact T-matrix、普通Full-J、ordinary goal-ROM与exchange；同V/C应一致。完整256-cell端口setup若采用，最多256个传播RHS/anchor，不能用latent66替代。当前一份合格anchor不能自动满足原两anchor竞争协议；不得据本文件启动它。
