# 图注与来源

## 图 1 共享墙体的时序 Maxwell 因子图

静态墙体 \(W\) 同时连接全部观测因子 \(f_t=p(y_t\mid W,X_t,q_t,\gamma_t)\)。状态 \(s_t=(X_t,q_t,\gamma_t)\) 包含目标、收发器位姿和仪器变量，\(\psi_{t-1,t}\) 表示已声明的时间转移模型。圆为变量，方形为因子，灰色圆为已观测数据；每份 likelihood 只出现一次。图中使用 Markov 仪器状态作为可选建模例子；一般跨时刻仪器先验 \(p(\gamma_{1:T})\) 应添加相应高阶或长程因子。图示只说明概率分解，并不宣称消场后的 Maxwell 材料灵敏度具有局部空间稀疏性。

文件：figure_01_maxwell_factor_graph.png / figure_01_maxwell_factor_graph.pdf。

## 图 2 有限维线性 Gaussian 理论的三个例子

（a）两个单位白化噪声测量 \(y_1=w+x+n_1\)、\(y_2=w-x+n_2\)，分别的有效墙信息都是0，共享静态目标下的联合有效信息为2。这是标准共享补偿不一致现象，不单独作为原创主张。（b）报告D2的秩亏目标状态例子中，不正确地把 \(G^{-1}\) 替换为 \(G^\dagger\) 会制造 \(4/3\) 的虚假信息，而直接联合 Schur 和正确 quotient 均为0。原始矩阵和零残差补偿向量保存在 observability_toy.py 与相应JSON中。（c）两类实验 \(y^{(1)}=w+x+n_1\)、\(y^{(2)}=2w+x+n_2\)，每类独立重复 \(n\) 次，噪声标准差 \(\sigma=0.2\)，永久墙偏差 \(b=0.1\)。完整联合估计的目标均方误差为 \(5\sigma^2/n\)，错误冻结的均方误差为 \(9b^2/4+\sigma^2/(2n)\)。绿色线在预定 \(n=20\) 时重新开放墙变量并使用全部历史数据，从而恢复完整联合估计曲线。该切换时刻不是数据驱动重开策略的性能结果。

三幅面板均为有限维解析例子或矩阵核验，未运行 Maxwell PDE，不代表物理墙体可实现性、材料重建效果或计算加速。

文件：figure_02_linear_algebra_counterexamples.png / figure_02_linear_algebra_counterexamples.pdf。

## 再生成

运行 deliverables/code/make_research_figures.py。图像生成源代码、figure_metadata.json、两项代数核验脚本和JSON共同记录来源；图中没有使用未经运行的 Maxwell 实验结果。
