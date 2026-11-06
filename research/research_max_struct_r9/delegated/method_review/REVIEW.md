# R9新增方法的独立接线审查

范围：`THEORY_DERIVATIONS.md`、`src/landscape_solver.py`、`src/baseline_adapter.py`。仅为核对成本计时入口，另读 `experiments/run_landscape.py`；已读父线程保存的 `audit/gradient_wiring/RECEIPT.json`。**未修改代码、未运行测试/求解器、未打开GT。最终数学及科学判断由父线程负责。**

结论：未发现新的公式到代码接线阻断错误。发现一项影响等CPU曲线的实际记账不一致；另外两类比较合同差异须保留，不能包装成单一机制的因果结论。未见“已凸化原材料问题”或“新发明WRI/CSI”的表述。

## 需要处理的实际问题

### [P1，若用于等CPU前沿] Landscape预算轨迹缺少准备前缀，不能直接与baseline轨迹配对

- `src/landscape_solver.py:18` 才设 `self.cpu`，`:31-34` 的预算时间从此计起。
- 但 `experiments/run_landscape.py:22-32` 在创建LandscapeRun之前已完成NumPy/求解模块导入、公开数据读取和MaxwellModel准备。`:39-40` 的 `process_CPU_seconds` 才是覆盖这些工作与输出tail的总量。
- 对照的 `src/baseline_adapter.py:270-278` 明确在继承的startup-inclusive预算上再加R9准备前缀；`:285-289`另存全过程CPU与scope。

因此，同一原始预算CPU数值代表不同费用边界。直接拿Landscape `BUDGET.json` 对baseline `R9_BUDGET_STATES.npz`做等CPU选择，会将新方法的准备成本暂时免费化。最终整进程CPU仍可用于端点成本比较，公式和最终Phi不因此失效。

当前 `process_CPU_seconds − result.cost.CPU_seconds` 只能给出 **prefix＋序列化tail**，不能精确分离prefix；`run_landscape.py:33-40`表明tail发生在run计时结束后。已有轨迹可保留，选择以下一种清楚标记的处理：只比较最终process总量；将全部未分配overhead保守加到所有Landscape预算点，并称“保守上移”；或将精确startup-inclusive前沿标为未取得。不能把两个总量相减后称精确prefix，也不需要为本次审查重跑。

## 必须保留的高影响合同差异（不是新增公式bug）

### 1. 松弛方法的主要代价不能仅按physical LU/RHS比较

`src/landscape_solver.py:84-94` 每次profile评价都构造完整576²的 `L*L`，进行一次576维复Cholesky，并求解全部source RHS。这些被记入 `auxiliary.current_Cholesky/current_RHS/dense_state_Gram_builds`，没有冒充physical LU，这是正确的分类。相应成本仍是真实稠密三次方工作；只以physical LU或physical complex RHS对齐，会漏掉新方法的大项。

`THEORY_DERIVATIONS.md:79` 所说读取成本曲线，必须同时使用含上述工作的CPU或明确保留这些辅助计数。少physical LU不能推出更少总因式分解工作、等总算量或加速；不应把Cholesky与LU简单按次数相加当同成本。

### 2. 当前比较改变了路径、监测候选集合与corrector起点策略

`src/landscape_solver.py:119-127` 的profile仅每25个accepted点及stage终点支付/保留原Phi；data homotopy在 `:37-46` 每次函数评价均保留原Phi，包括line-search评价。`:129-133` 把最后辅助iterate送入原corrector，而最终 `:134-141` 从所有已支付物理候选中选择Phi incumbent。baseline adapter则继承R8的“先选initial/stage最小原Phi，再启动corrector”策略。

`THEORY_DERIVATIONS.md:73`已披露最后iterate传递与历史选择不同，故这不是隐瞒或GT挑根；所有额外物理评价也由model计数。只是不能把结果解释为**只改变λ/目标而其余选择完全相同**。此外，profile_homotopy三次调用minimize会重置L-BFGS历史，fixed profile的一次300步调用不会（`:116-124`与`:49-53`）；这里观察到的是整个冻结schedule/重启/传递合同，不能单独识别λ变化的因果效果。

## 公式与实现一致的关键处

- **完整current松弛：** `THEORY_DERIVATIONS.md:50-62`的K、rhs、包络梯度，与 `src/landscape_solver.py:82-92`一致；e是原LSIE缺陷，无材料相关CIE p(m)，因此没有漏掉一个本应存在的p导数。
- **λ0：** `src/landscape_solver.py:69`先解 `L* X=S*`，再取共轭转置，确为 `S L⁻¹`；`:70-72`求其谱范数平方并乘ed/yd，与理论`:73`一致。
- **data homotopy归一化：** `:40-43`先去除prior，再将继承函数按target能量归一化的数据项乘回 `||target||²/yd`，最后加一次原prior；value与gradient同时重标度，符合固定原yd的推导。对本次非零target合同，未见归一化接线错误。
- **噪声规则：** `src/baseline_adapter.py:129-140`使用公开原始Y、完整数据球 `sigma² chi2(2Y.size)/2`，正向累加正交尾能量，并限制为数值支持的奇异前缀；`:154-155`确实选最小合格前缀。`:239-251`只把选出的rank传给继承CIE，额外SVD/准备CPU另计。没有GT选秩或按rank反复使用不同点态99%门。
- **概率与凸性措辞：** `THEORY_DERIVATIONS.md:42-44`用同一个完整噪声事件同时控制所有投影，即使r依赖Y也不需要逐r多重检验；此论证成立。`:48,71,77,83-85`均保留已有方法归属、非凸性与物理合同限制，未见从可逆Möbius换元推出原Phi凸化的错误。

父线程已有必要导数检查：profile方向导数相对差约5.99e−10、target约4.03e−11，检查点内解normal residual约1.22e−15。本审查没有重复该检查。这支持该检查点的接线，不是全box、所有λ或所有迭代的数值精度/收敛证明。profile normal residual目前仅在stage终点显式记录（`src/landscape_solver.py:98-100,125-126`），不能把终点值写成全部内解的最坏残差。

## 任何当前结果仍不能推出的结论

1. `K(m)`对固定m正定，只给current内问题唯一解；不意味着profile对m凸、没有错误边界点或延拓必达真材料。
2. full-noise-ball保证不是rank最优、保留模式无噪声或真材料在截断Fourier族中的保证。固定r的噪声二阶矩公式（理论`:35-39`）不能不处理选择依赖就直接代入随机r(Y)，宣称是给定所选rank后的条件期望或误差证书。
3. 中间stage的 `projected_gradient_inf` 来自当前profile/target梯度（`src/landscape_solver.py:55-58`），不是原Phi梯度。最终prediction还可能来自较早已付费候选，因此某个stage端点的PGTOL声明也不自动传给最终选中的材料。
4. 更低松弛目标、原Phi或同数据残差，均不等于材料GT改善；四个公开DEV案例不能重新充当独立blind确认、普适基线胜出或新颖性证明。此审查没有重做文献优先权判定。

审查产生的新物理/测试/训练调用：0。只写本文件；不请求扩展预算或额外实验。
