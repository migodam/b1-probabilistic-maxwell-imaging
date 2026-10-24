# 固定核实际链：双端势区间与独立Full对照

新增src/p3r4/interval_chain.py及tools/run_r4_interval_chain.py。本轮只实现及小矩阵核查，未启动Maxwell链。星形结果不再用于反事实加速推论；下面必须实际运行两条独立arm再比较全部路径与整作业费用。

冻结首8个新DEV_00016…00023；各1024步；full66Gaussianized pCN，beta=.002，seed20261025090。每parent的独立进程随机流使用SeedSequence([seed,parent_index])；每步先normal66再uniform1，两arm只有相同初值、固定核和随机流合同。运行过程中两arm不读取对方任何输出。初始onlinecontext/cache的字节数及mtime身份记录在各自gate中，比较时须相同；不用SHA。

若当前Φ∈[a,b]、候选Φ'∈[c,d]，则log MH ratio位于[a−d,b−c]。同一个logU若小于min(0,a−d)即接受；若不小于min(0,b−c)即拒绝；否则确定性地先精化较宽的不确定端点，同宽时先候选。当前和候选都能被精化，沿用同一个U。每个端点从最初线性近似起累计最多4次Richardson修正；不能精化或仍不确定时，显式付费Full刷新较宽端点，直到判断完成。这里使用的是决策刷新规则，不改变pCN提案本身。

accepted approximate state保留其current/field及Φ区间，下一步继续使用该区间；它没有因接受而成为精确cache。拒绝时当前端点可能已经被精化或Full刷新，该信息仍保留，同时显式保存重复x状态。η使用原paidR0的近似左逆缺陷公式，始终相对于最初同一个base；没有自动重基。η≥1的候选必须Full fallback。Full回调要求实际账本增加8个复源RHS，返回同一x/q/χ的精确势。所有接受与拒绝依据同一U、当前和候选两端势的合法包络；浮点实现仍不是定向区间证书。

interval arm每parent仅一次通过backend.solve_anchor(I576,base.state)取得R0，实际576逆列RHS必须核对。Full arm只使用共同已付初始cache，随后每步实际Full候选8RHS。R0密集应用只计矩阵动作、时间、内存，不再计为求解器RHS。两arm均无GT读取，没有新网络训练。

父调度器须写明确release：

```json
{"released":true,"experiment":"interval_chain_v1","campaign_id":"P3LR4-20261025"}
```

然后分别以不同进程运行：

`<python> tools/run_r4_interval_chain.py interval --release RELEASE_PATH`

`<python> tools/run_r4_interval_chain.py full --release RELEASE_PATH`

不得在一个进程中复用任一arm的决策。既有physics.lock保证物理作业串行；两arm可以任意顺序运行。各自StageII Book不超过1200秒、3GiB。为不让一个arm读取另一个活跃账本结果而仍保证200000 combinedRHS上限，预分配Full arm65536RHS(8×1024×8)，interval arm134464RHS(含4608逆列setup及所有真实fallback)。任何新增Full、失败及求解器RHS均通过原Book记账；达到各自硬上限即停，不拓预算。

两arm有相同的公共持久化权利：initial、每step候选及每step接受后或拒绝重复状态，均立即保存x/q/χ/current/field/prediction/residual、Φ区间及exact标志；区间arm的每次近似精化另行立即保存。Full后端的标准attempt/endpoint证据也保持。每step的写盘和计算都计入该arm实际wall，比较以Book关闭回执和整个作业费用为准，不能只比矩阵计算时间。

只有两个ARM_COMPLETE.json都存在，状态均COMPLETE、各自Book CLOSED均COMPLETE且PID不同后，才能运行：

`<python> tools/run_r4_interval_chain.py compare`

事后逐parent比较全部1024个accept标志、含初值的1025个x状态路径，并检查interval决策时记录的当前及候选Φ区间是否包含独立Full arm的对应真值。使用严格数值包含与路径逐项相等，不加GT或结果依赖容差。任何路径/判断/包络不一致均写NO_GO并以失败退出，保留全部证据；不继续重调阈值救结果。两个arm尚未关闭或任何arm未完成时，不打开对方链数组。

即使全部一致并有整作业费用改善，也只支持这批固定核、同随机路径的决策成本结论，不赋予posterior mixing、Rhat、MCSE或FULL后验资格。链很小的beta可能保持高度相关，这一问题不能被路径一致性检验解决。

INTERVAL_CHAIN_MATH_CHECK.json核查512步任意12×12复线性系统的独立Full/interval小链：x路径逐项相等、0判断不一致、0势区间遗漏，493个接受状态仍为近似状态，另一个强制η≥1例子正确调用显式Full callback。新Maxwell状态/RHS及训练更新为0；这不是本计划的实际物理运行结果。

星形已封账结果费用补记另见results/II/decision_refresh_star_v1/COST_ADDENDUM.{json,md}：普通GMRES等价费用补回每候选8个原始源方程RHS，共8192，使等价总和46402→54594；实际已计费62770RHS不变，不再次加费，不改原始结果或账本。

## 调度兼容入口（实现完成后冻结）

省略arm时，`<python> tools/run_r4_interval_chain.py --release RELEASE_PATH`默认为pair协调模式。协调器不构造Book/Target，不读GT，不执行物理或链数学；它依次启动interval、full、compare三个不同子进程。显式arm入口保持原合同。只有三者成功且比较PASS，才写results/II/interval_chain_pair/PAIR_COMPLETE.json；该文件是调度兼容完成标记，明确引用interval_chain_v1/COMPARISON.json和固定核费用范围，不是额外实验或后验资格。比较NO_GO时同名标记状态为NO_GO并非成功；子进程失败以非零状态停止数据扩展，保留PAIR_FAILED/已完成arm，不自动重跑付费作业。若目录已有任何pair/arm证据，协调入口要求父调度器先核对，不悄悄覆盖。
