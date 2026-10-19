# 完整失败及解释

9/24基础解未进入noise-ball，包含4/8 held-out。其真实材料由同P1生成；不支持先验无解结论。部分合格端点也在subproblem失败处停下，不称MAP驻点。

材料distance/progress、Gaussian局部半径或bounds不满足时不付费forward，普通数值时间记录；实际forward失败/拒绝、line search均计入。B blocks保留每次decoder失败与BLOCK_PHYSICS_LIMIT。

C原始ERROR、CLOSED、BOUNDARY_ERROR_DIAGNOSIS及恢复receipt都保留。机器精度恢复阈值64 eps；不是任意物理材料clip。

高温采样没有给出可信posterior，ESS或低残差不证明mode覆盖。Flow/Gaussian训练目标是少量corrector轨迹，不能称完整条件posterior或进行概率校准。

主Full-J对照有angular-progress confound，补充sphere-only揭示其不能被当弱陪跑。补充已exposed，不能回溯升级主blind结果。

11个原archive拓扑标签变化中，只有1个同时通过结构检查；没有同时在三个阈值稳定的结构拓扑改变。主held-out无拓扑改变。不能以阈值闪烁制造多解。

Semantic Scholar focused collect返回429，保留访问缺口；后续只核对少量原始arXiv/出版方记录，不称穷尽新颖性搜索。原生subagent因已有线程槽位满未启动；核心全部由主Codex执行，无DeepSeek。

逐solve RHS、LU与总wall完整记录。cached restore矩阵audit次数可由CACHE_REUSE准确推回；dense surrogate/SLSQP callback的逐次矩阵作用计数未全面记录，只包含在受控总wall中。不声称完整FLOP账本。
