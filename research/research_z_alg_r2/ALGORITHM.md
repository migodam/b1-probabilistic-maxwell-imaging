# AMBS-Q：观测到材料候选的可运行算法

## 输入与适用条件

输入是固定几何/背景Green所对应的采集配置、复观测Y、事前材料box。核心Q6适用于一接收、两独立照明、四个独立实材料；可以是整个四单元模型，也可以是明确固定其余材料后的四单元切片。代码不读取真值、真值支撑或生成seed。

充分接收的成像另提供SVD直接逆和SVD初始化的Full-J校正。两条算法使用相同完整非线性VIE，但数学适用条件不同。

## Q6求解步骤

```text
准备 G,E0,S；取与GT无关的m0，确认M0可逆
固定活动坐标U，形成K,L,q,Y0；若从已有大核切片，使用Schur复用
确认rank(q)=2；从wq=Y-Y0消去两个复坐标，w=w0+cN
构成四个共享实材料兼容式，独立复化共轭参数d
整理A(c)[1,d1,d2]=0，计算四个最大子式的系数
分出三次块T与低次块R；若数值条件不合格，计费六路径fallback
Q = solve(T,-R)；构造Mx,My
对(Mx+eta My).T求特征向量，联合恢复c，再回代d
检查六候选分离、回代秩、兼容残差；未决状态保留
逐候选检查d=conj(c)、节点、实材料box、resolvent分母
去重材料；对每个剩余候选用fresh full Maxwell solve验收
输出全部合格材料与全部拒绝/失败分类，不利用GT选其中一根
```

实现：`src/quotient_solver.py`、`src/branch_enumerator.py`、`src/physics.py`。固定η、阈值、最多polish步与fallback均保存在各阶段冻结合同；结果中区分数值资格和严格证书。

## 正维纤维的遍历

E1有六个材料、四实观测。首次粗网格固定后两材料于{.75,2.25,3.75}²，多数切片无box内根。该失败原样保留。

第二阶段从观测独立求一个合法anchor；后两坐标仍事前固定，未依据弱Jacobian或GT选择。在anchor自由坐标附近用冻结±.6网格取切片，每个切片独立枚举Q6/六路径全部候选。切片的已知自由坐标与未知四坐标组成完整材料，全部独立验收。这能跨同一切片上的多个离散根，不局限于从anchor连续跟踪一支。

这是材料切片遍历，不是整个二维纤维、不可约分量或连通分支认证。想覆盖更大范围须换切片或处理投影边界，不能凭当前网格推断。

## 成像与强普通对照

- `svd`：J=S†Y，E=E0+GJ，按共享实本构行拟合恢复材料；接收满列秩且无噪声时精确。噪声可严重放大。
- `hybrid`：同一SVD候选裁剪到事前box，再用原非线性有界Full-J TRF校正；不使用GT调参数。
- `fullj`：固定box中点出发，精确Jacobian，动态LU共享；较大模型使用min(材料列数,接收数)方向的导数求解。
- 枚举对照：标准24路径、多起点Full-J、相同固定核/解析Jacobian的ReducedFullJ8与64。普通方法获准使用相同SMW缓存。

48单元实验只在最终验收计算完整resolvent最小奇异值，避免每次普通试探额外SVD拖慢基线。失败、原始SVD越界、未合格残差均保存。

## 缓存和复杂度

全n单元、k活动单元的稠密准备需一次n×n分解及s+k RHS，代数求值只用k×k系统。Q6固定k=4后的主要在线工作为4个小多项式行列式、4×4六右端solve、6×6eig和最多6个4×2回代；最终每个候选另支付完整验证s RHS。

不同Y共享物理核和测量子空间，不共享整个目标依赖T逆。改变自由材料切片时，F11的Schur式复用全核，仍支付真实小矩阵运算。跨目标总成本和实测crossover见研究报告/成本账本。

## 数值失败模式

rankq不足；T近奇异；重根或η碰撞；回代条件差；有限node_floor附近不确定；极点；实共轭失败；材料越界；完整验证失败；路径跳支/碰撞。标准24路径的8个切片发生端点碰撞，故其“运行完成”不代表完整覆盖。没有通过放宽合法box或误差门修复原实验。

## 直接运行

在本目录，任意装有NumPy/SciPy的Python可运行：

```bash
python run_solver.py --geometry experiments/QUOTIENT_HOLDOUT_CONTRACT.json --observations experiments/quotient_holdout_public/QH4_03.npz --box 0.5 30 --method quotient --output /tmp/qh4_reconstruction.json
```

输出包括所有候选、完整验收和费用；文件已存在时拒绝覆盖。新观测只需NPZ中的`y[receiver,source]`及同格式几何JSON。更完整复现、冻结原始实验与评分见[REPRODUCE.md](REPRODUCE.md)。

## Q20与一般商构造：六材料、三照明全部未知

同一算法保留s个自由复坐标，构成2s×(s+1)矩阵A。计算全部最大子式，按总次数s+1与≤s划分T、R；在T合格时解Q=solve(T,−R)，生成s个D×D乘法矩阵，其中D=binom(2s,s)。固定分离线性组合的转置特征向量给c，回代d，再执行同样的物理过滤。

新实现`src/general_quotient_solver.py`用于冻结s=3研究：15个四次子式，15×15系数求解（20右端），三个20×20乘法矩阵，20个复候选。它不修改原Q6源文件；此研究遇到条件失败会明确停止并保留UNRESOLVED，不临时改参数。新目标仍重新计算系数，并计入在线成本。

```bash
python replay_general6.py --output-dir /tmp/zalg_general6_replay_new
```

该入口只读公开观测，按已冻结参数重放Q20与ReducedFullJ64，把全部候选和费用写到新目录。发布后重放不再是新盲测试。新原始观测的程序接口是`enumerate_general_quotient(kernel, Y)`、`filter_general_roots(...)`，之后必须逐个fresh full-forward验收；完整冻结驱动位于`experiments/run_six_material_holdout.py`。

一般公式的组合维数并非低复杂度保证。s=4的70维及更高未运行，六材料当前实现也没有认证Hankel系数误差、极小节点或重根。完整推导见[general_quotient](theory/general_quotient.md)。
