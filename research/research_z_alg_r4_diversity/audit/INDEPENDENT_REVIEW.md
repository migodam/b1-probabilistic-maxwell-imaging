# R4-DIVERSITY 独立审查

状态：`PRE_SCORE_STATIC_REVIEW / SCIENTIFIC_OUTCOME_HOLD`。日期：2026-10-30。独立审查只判断冻结协议、代码可见机制与推导范围；主代理保留最终科学裁决。本文件是本代理唯一拥有的产物。

没有读取 `private_truth/`、映射或任何新 GT；没有启动求解、测试、数值精度核对、重复源码比对、hash、新文献检索或其他代理。已读合同、SOURCE_FREEZE/FREEZE.json、算法说明、solver/physics/CUDA 实现、score_results；为审配对和硬件范围定向读取 generator 与 remote runner 的源代码、公开生成收据及调度说明，没有执行生成器或 scorer。正式分数尚未评阅。

## 首审结论

未发现足以要求中止冻结实验的确定数学或实现错误。公共物理问题、在线目标和材料可行域相同；Data L-BFGS 与 GlobalDual 的对照较直接。两种 TRF 允许不同 coarse 路径、变量度量和终点，构成成熟求解器基线，不是与 L-BFGS 的逐步相同迭代。

预算是带回收规则的共同上限；不能将其说成每个 arm 实际消耗相同计算量。评分的 NRMSE、practical tie 与 offline oracle 边界在代码中明确。**结果、实际执行完整性、GT 恢复与方法优势继续 HOLD。** 当前审查不为将来可能出现的失败、缺失或后端差异预先背书。

## 冻结与数据层次

`SOURCE_FREEZE/FREEZE.json` 记录 UTC 04:29:42.767，公开 generation receipt 记录 UTC 04:29:44.367，声称冻结时新数据尚未生成、正式 GT 不传到逆求解主机。合同的 72 观测为 12 形状 × 3 材料增量(.4,2.5,12) × clean/noisy；GPU 18 观测为预选 digit2/digit7/ring 三形状的同一设计。CPU 288 与 GPU 72 method tasks 合计360。审查接受这些记录的自述，没有做进程隔离或不可变性独立证明。

generator 源码对每个形状固定 alpha、一个随机背景(.08,.15)和一个复噪声方向，再跨三个增量复用；clean/noisy 共用相同 clean forward。noise 按各自 clean 数据范数缩放到相对 .001。因此72不是72个独立形状重采样；自然配对/聚类单位只有12个形状，每个形状有六个相关观测。跨 contrast 的背景和噪声方向配对有利于比较，但不能用72个独立 Bernoulli 试验的直觉解释 win counts。

GPU18是CPU72的预定子集，不能并入得到90个独立样本，更不能当独立复现实验。CPU全量与GPU子集的总体成功率差异混合了样本组成、硬件/库、预算内迭代路径与资源竞争；算法比较应在各 backend 内配对，后端差异只在重叠18条观测上作描述。没有运行数值一致性测试，不宣称CPU/CUDA输出逐位或任意容差等价。

## 方法公平性

| 检查项 | 首审判断 |
|---|---|
| 固定576物理网格、12源48接收、公共Y/G/E0/S | 同权限；8/12/24只改变材料系数基，没有换粗PDE或增加观测 |
| 边界[.02,20]、m0=.1、共同R3正则 | 相同；没有真值依赖初值、mask或online support |
| Data L-BFGS vs GlobalDual coarse | 8和12层均原目标L-BFGS；不存在GD独享粗层启动策略 |
| final24阶段45/55 | 两者都reset optimizer memory并用同资源分割规则；第一阶段目标不同，第二阶段均用各自第一阶段的最佳原目标状态 |
| 原目标最终选择 | 所有完整物理trial按同一normalized data loss+R选择；GD辅助目标不直接替代最终排名 |
| native vs reciprocal TRF | 同c可行集合；eta=1/c只是坐标变换，coarse材料仍m=Bc。并不是m=1/(B eta)。两者路径不同是合法算法差异 |
| TRF vs L-BFGS | 相同物理权限/目标/总上限；TRF全阶段解析J、dense SVD、x_scale=jac与L-BFGS变量缩放不同，不能解释为只改一个残差度量的纯消融 |

Data/GD的“同粗层”是同算法合同，并不保证在wall截止时各backend、各worker完全相同的c；真正结果应保留stage/phase事件。GD第一阶段增加Gram/伴随成本，数据L-BFGS对应阶段使用原数据伴随；这份差异正应由共同预算承担。

正则使用固定差分尺度19.5，各contrast与method一致。它使横跨材料增量的曲线同时包含归一化损失尺度、有限预算和正则相对强度变化，不能把全部contrast差异单因果归于多重散射或Gram变化。

## 预算、计费与硬件

solver的三个stage采用累计 .20/.45/1.00 wall/LU/RHS端点；早收敛余量可留到后层。24层第一phase端点位于当前剩余资源的45%，第二phase到共同global终点，首phase未用完的份额会转给第二phase。Data和GD同规则。所谓45/55是分配边界和回收策略，非实际耗时/操作数恒等。

物理operation开始前执行global及phase guard；已经开始的12源forward、48接收J或12源伴随块可越过剩余RHS边界，dense/opaque optimizer工作也可能越过wall。计费在尝试开始时登记，异常保留failed_operations；这与合同相符。预算内中止可能已有forward但未得对应梯度，best物理trial仍可被保留，不能称其为收敛终点。`all_unknowns_released`仅表示实际出现stage24 trial，未必top1来自24层或已探索完整576维。

物理计费分开了LU、源、接收转置、伴随RHS；12×12 Gram的48个小系统RHS独立计数，没有混作576维RHS。TRF内部SVD工作包含wall/CPU而未伪造调用数。CUDA每个counted block结束同步，并在solve前/最终计时前同步；complex128路径明确，没有float32替代。CUDA仍用CPU SciPy优化器和正则，因而不是端到端纯GPU算法。

remote runner读取Windows core topology，以观测i9的“同一物理核多个逻辑siblings”为P-core识别规则，每worker只绑一个逻辑sibling；四CPUworker不共享同一已选择物理核，两个CUDAhost worker用另外两核。它不独占GPU或其他程序资源，两CUDAworker共用device0，wall含彼此竞争。公开计划不是实际affinity回执；最终应据worker receipt记录执行状态。A9占用时只推迟GPU的调度说明不改变算法，但CPU/GPU不是同一时刻/相同负载环境，不能据本批宣称硬件speedup。

硬终止/worker异常收据可能没有已付费LU/RHS，scorer按known_count与UNKNOWN保留。最终成本须同时报告已知成本、缺失成本和失败任务，不能把未知当零或只比幸存成功任务的平均耗时。

## 评分与GT泄漏边界

scorer不import求解器、不做新forward；只有显式关闭planned dispositions后才读取local truth映射。在线求解路径未import generator或score。公开观测payload只有Y与匿名ID；solver的load_public并非通用安全沙箱，但当前worker只加载显式operators/public NPZ，模型只访问所需公共keys，未见读取truth路径。

`seal`会记录每个planned task为complete/failed/NOT_RUN；它允许缺失任务被明确关闭，而非只有全成功才封存，这是合理的预算结算方式。但seal保存路径，score时重新读取该结果路径；这是协议性的“关闭后不再改写”，不是不可改动快照证明。父代理应保持结果关闭后不改，不需要新增hash或重复比对。全部计划的CPU/GPU状态均应封存后再开启新真值评分；若GPU未运行，保持NOT_RUN，不能使用GT后续调参并继续将其称为fresh。

NRMSE代码确为RMSE/已知材料增量，属于offline真值指标；background不进入该归一化分母。它使三个contrast的相对成像精度可比较，但同一个NRMSE门对应不同绝对RMSE，必须并列保留RMSE/MAE。practical tie按合同 `abs(delta)<=max(.02,.10*min(NRMSEleft,NRMSEright))`，代码一致，连续delta也保存。它是预声明的实际差异阈值，不是统计显著性、等价检验、置信区间或方法同等可靠的证明。

Otsu阈值只取reconstruction自己的min/max和histogram；与GT mask计算IoU属于评分，不是online阈值泄漏。known-amplitude midpoint IoU明确OFFLINE，不可用其做online选图或更换门槛。data gate与material gate和joint gate分开，clean .001/noisy .003阈值冻结。finite_endpoint当前只检查576个有限材料数，不额外证明方程可逆、KKT或全部未知量执行；正式报告还须看原残差、state residual、状态和stage事件，不能把finite当physics/scientific PASS。

oracle best-K只对有限phase-best材料做GT后验筛选，candidate不是独立root/basin。第一名在线选择按原目标，不能用oracle指标替代它。pairwise按observation_id、backend内配对，clean/noisy也按shape/amplitude/backend/method配对，没有把CPU和CUDA揉在一起。总体pairwise可描述本数据集，但低/中/高增量的差异应看冻结strata与continuous delta；不得挑有利stratum讲普遍优越。缺失pair保留missing，不得忽略。

## 弱材料解释的必要条件

固定有限VIE中，若 `gamma=||G D_m||_2<1`，则

`Q=(I-G D_m)^(-1) E0`,

`||Q-E0||_F <= gamma/(1-gamma) ||E0||_F`。

令`deltaQ=Q-E0,H0=E0^H E0`，有

`H-H0=E0^H deltaQ+deltaQ^H E0+deltaQ^H deltaQ`，

`||H-H0||_2 <= 2||E0||_2||deltaQ||_2+||deltaQ||_2^2`。

因此在确实弱传播且稳定的范围，GlobalDual主要是在原物理残差的source维度上使用近固定的加权度量 `tr(R H0^-1 R^H)`，而不是引入新的观测或普遍消除非线性。若H0与标量单位阵近似，该度量接近原数据范数的缩放；若H0明显各向异性，则会改变残差方向的权重。共同正则和有限预算下，两种度量均无普遍成像优越性。

**H近H0不自动给逆矩阵近似。** 还需H0正定，以及 `eta=||H0^(-1/2)(H-H0)H0^(-1/2)||_2<1`。此时可得到相对inverse/weighted-loss扰动界，量级为eta/(1-eta)。若源Gram病态，微小绝对扰动也可使eta很大；仅写“低材料所以H固定”不足。thin-QR回退改善数值求解，却不会创造源独立性或使统计度量变成原范数。材料增量.4是协议分层标签，连同背景.08-.15并不单独证明上述小范数条件。

此外，目标值近似不自动证明梯度/迭代方向近似：精确GD梯度含H随材料变化的项，代码中的`-G^H(Q X^H X)`正保留它。若要省掉该项，还需邻域导数有界与相应小残差/弱传播尺度假设。不能从Q≈E0一句话推出全域GD等价于固定预条件器，更不能推出优于数据L-BFGS或TRF。

## 结果到达后的追加范围

只追加读取封存收据与实际scores，分别报告CPU72×4、CUDA18×4的planned/finite/failed/NOT_RUN及分层joint gate；审查matched GD/Data L-BFGS的paired continuous NRMSE差与practical ties、成本和stop reasons，再看两种TRF。CPU/CUDA重叠观测可描述路径差异，但不合并样本量、不作精度等同性或hardware speedup声明。

最终科学状态现在保持`HOLD`：无普遍优越保证、无新理论或后验/概率解集结论；正式GT尚未评阅，本首审不得被引用成方法有效的独立确认。
