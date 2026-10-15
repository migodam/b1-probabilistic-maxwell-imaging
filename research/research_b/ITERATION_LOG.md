# Research B 迭代记录

## 迭代 1：缓存上的共享补偿与材料风险

研究假设：当前数据允许的材料世界，可能在第二频率或独立 gain 校准后无法沿用同一合法补偿；raw separation 与材料风险未必一致。

先用 Joint P2 已保存 full-Mie 响应和 accepted reference 做 E0。禁止缓存缺失时自动调用 solver。固定 material loss 为球体体积加权 Re/Im ε 分量 L1、s_ref=1，候选选择风险与无约束图像中位数分别报告；事件沿用 P2 分区。E0 所有新噪声标为模拟，不能把数字组合当新增观测。候选集合不输入概率；集合策略用复合假设距离并声明其筛选性质。

预期成本：E0 新 Maxwell=0；后续物理总额 4096 attempts、受控 wall 7200 秒，先冻结菜单和子分配，再运行。所有运行回执追加到独立账本；不改历史配置、数据或 Gate。

主要比较：no measurement、random、preset 2.4、uniform、raw separation、Fisher、joint projected Fisher、event EIG、material risk、robust composite、shared compensation conflict。

## E0 结果与理论修订（在E1新truth前）

七个暴露案例的连续proxy条件预测风险均选择2.4。单独最大化共享补偿冲突却均选calibration；这表明冲突量必须与材料任务方向和条件可用信息一起解释，不能当作材料收益。旧七例单次新噪声下的实际MAE顺序不同于predictive风险，不能用它来挑选获胜策略。所有E0新Maxwell=0。

可证伪的下一轮预测：新频率菜单可能使material-risk/robust与普通Fisher/EIG分开，但在完整成本和新parent上未必超过fixed2.4或uniform。优先检验这些简单基线是否已经足够；保留原conflict作为失败诊断。冻结E1_FROZEN.json后才生成32新parent和独立噪声。E1/E2合用本轮4096 attempts、7200受控wall；不借旧余额。

## 迭代2：新频率验证结果与下一轮预测

E1共32独立parent、每parent2噪声。固定3.2、Fisher、joint nuisance Fisher和material-risk在第一block完全同选；Mode1 Selected-MAE均0.101982，对比不测量0.162279、固定2.4为0.122534。不能宣称复杂准则优于固定3.2。Robust为0.103714，并未胜出；raw-conflict几乎只选calibration，MAE0.141574。

下一轮明确修订：停止追求更复杂准则的普遍优势；测试是否能在最多3block中通过条件停止节约测量，同时诊断固定候选支持的MAE下限。保留所有原基线，新增有限参考风险停止与校准二步bundle，比较完整成本前沿。E2复用E1对象，明确为暴露后续分析，不称新blind test；新增物理只为预算内的机制/独立验证和A真实包交接。

Research A已产生真实Full-Mie候选，固定lexicographic前四包confirm00/01/02/03共27候选，不按truth或表现选择。只读snapshot，未来响应及校准更新由B独立计算，本次账本E2计费；无需改A目录。

## E2/E3实际结果与最终理论修订

三block的Mode1 fixed3.2 MAE=.090885，robust=.093842，material-risk=.093775，event-EIG=.090107；parent配对区间未证明复杂策略优势。risk-stop平均2.984375 blocks，节约很小。625池Best-of-K=.062224，只有40.625% parent达到.05 coverage；删除family不能靠残差sentinel可靠识别。

因此修订预测：候选支持扩充比继续细调acquisition更可能改善材料任务。E3在新观测前冻结前四真实A包和240次子分配，保护原候选、进行双频fullvalidated ε corrector与三个无拓扑标签的连续prior起点，检验是否真正改善Best-of-K。结果protected Selected-MAE=.064836，prior多起点=.036170；固定候选选测=.119000，原base=.125480。只有四例，缺相同容量无新测量corrector对照，不能声称提升全由选测造成。新增端点分别20/12；prior优化有2个未达收敛flag，仍保留状态及通过full/data检查的端点。

概率参考的负诊断独立记录：原13/17审计12/16唯一history未稳，有限样本事件差可达.187；一次不同坐标adaptive仅2/16达原准则，不以其输出重选policy。停止当前付费campaign，剩余812 attempts。下一阶段优先coverage/corrector与公平信息增量控制；保持简单固定频率强基线及小孔歧义。
