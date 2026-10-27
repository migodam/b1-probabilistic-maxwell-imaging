# 理论、实现与证据的对应

|命题|必要条件/实际测试|最终证据身份|
|---|---|---|
|精确背景重表达|同一公共FE/PML空间，材料依赖b_X，完整输入输出；Full/Schur/resolvent场及方程比较|标准恒等式；本次真实离散全波数值核验通过，不是原创理论|
|跨目标cache|固定公共ROI、合法材料积分；12个fixed-v-fitted案例|QUALIFIED_DISCRETE_FULLWAVE_FOR_TESTED_CONTRACT，非uniform-prior/continuum界|
|相对保留率授权单点化|需绝对信息、proper joint、任务；9物理parent/12fixture|C2约98%variance仍存；不能angle-only；部分rank UNRESOLVED|
|材料不确定但任务稳定|有限15 proper prior、完整响应、固定权重预测替代|C3窄任务支持；不推导新likelihood安全|
|任务传播差|C4独立argmax审核|C4_1单端口d=.729、future max3.068；主观测可区分，无目标位移|
|单墙保joint|Q(W=g)=1时TV≥1−P(W=g)|有限概率精确命题；M1最佳单墙仍11/24失败；不涉及保留随机残差的换背景|
|概率压缩|原prior/likelihood与条件lift；2950全表重放|M5未发生任何merge；M4未省查询；物理posterior未认证|
|严格数值区间|原F64 vs独立归一化|原界FAIL，正越界保留；一次修复仅限固定log-potential|
|净成本优势|减少真实FE/RHS超过setup/cert/reopen，强基线同权限|当前机制NO-GO，G2依门停止|
|旧joint见证|保留oracle/解析support/旧minimum-feature/经验网格身份|历史条件证据，本轮未重算；cal距离3.56不能作cal omega(1)下界|

PDF的总场常量源公式转为backend散射场时须增加真实b_X，不是理论创新。PML非自伴；最小实部谱、condition估计和经验网格差不自动给稳定下界。Bayesian推送、Schur、sin²/广义保留谱、ordinary smoother/ROM均属已有工具。没有开展新的穷尽文献优先权检索，不作“首次”声明。
