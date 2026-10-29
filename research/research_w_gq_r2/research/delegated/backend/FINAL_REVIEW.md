# R2 最终实现与数据流复核

只读复核；本次未运行 PDE、测试或重跑，也未修改主代码。以下行号相对于本轮
`research/research_w_gq_r2/`。核对了列出的八个工具、冻结合同，以及其直接依赖
`src/r2/worlds.py`、`src/r2/physics.py`、`src/r2/return_model.py`、记账器和少量 JSON 顶层结构；没有读取大数组。
这些是实现与报告边界意见，不是最终科学判断。

1. **源 RHS 与返回分解的数据流一致。** `src/r2/physics.py:515–526,571–575` 实际组装
   `b_W` 和 `b_X`，每列对应一个解析线源；`src/r2/return_model.py:190–214` 保留
   `q=b_X+V u_W` 并通过 `L0-V` 做完整孤立目标自散射，未把入射场换成节点值。
   `tools/return_operator.py:93–100,119–138` 由 16 个接收伴随构造 `B_W`，再另外付费做
   真实 `G_W j` 验证。该验证未反馈给拟合或选秩。未发现漏掉 `b_W`/`b_X` 或源维度混用。

2. **打包与嵌套通道正确，但 norm/RMS 必须分开评级。** `src/r2/worlds.py:42–44` 为全部实部
   后全部虚部、每实坐标除以 `.01`；`tools/imaging.py:44–60` 只生成一次完整噪声张量，再取
   `f=1` 的前两源/四源，A0/A1 的已见噪声没有重抽。`tools/compensation.py:43–44` 也先在
   复数 source-major 行上截取，再实化。`tools/microstructure.py:49`、补偿曲线和
   `tools/return_operator.py:103–109` 的 white **norm** 都是此实化欧氏范数；返回工具同时
   存 RMS，E3 两源×16 RX 的 64 个实坐标使 `norm=8*RMS`。不得把 RMS≤1 当作
   合同中的 norm≤1。`tools/imaging.py:94` 的复介电相对误差实际在 `f=1` 定义，外频仅另有
   `loss_error`，不能称为跨频介电误差。

3. **有限微结构的界面与宏块均值在源码上匹配合同。** `src/r2/worlds.py:4–24` 使用 η/2 网格，
   `.25` 宏块、孔洞及占据边界全部对齐；冻结 n=8/12/16 对应每宏块 4/6/8 个横向子格。
   两种纹理每行横向正负交替且列数为偶数，因此各完整占据宏块的 ±.8 扰动面积均值为零；
   `mean`/loss 则在宏块内常数。`src/r2/physics.py:239–246,251–285` 将这些轴作为实际 PSLG
   约束，refine2 不改变材料数组。**限制：** `tools/microstructure.py:54–64` 的独立网格是
   每墙 parent 在 f=1/1.15 各选最强差异目标的端点复核，不覆盖全部目标、尺度或 f=.85；
   不能据此给出整个连续目标族的一致误差界。

4. **有限成像的在线/真值屏障成立于函数输入，属于有限候选库而非盲目标泛化。**
   `tools/imaging.py:50–69` 的 `infer()` 仅读观测及已付费公共宏墙库；真实微墙表只在
   `score():78–89`、在线封存后用于标明的 `KNOWN_WALL_ORACLE`。
   同文件 `16,22,41–47,51` 显示 12 个生成目标也就是公开候选目标，不能把这里的 heldout
   解释为在线不知道其候选坐标/材料。`CONTRACT.json:220–222` 已限定有限 counting prior、
   12 target parents、2 暴露开发墙 parents；24 个噪声案例不是 24 个独立墙。
   seal 是阶段化输入约束，不是操作系统隔离，也不是源代码不可修改的证明。

5. **连续 FWI 的目标轴、梯度与调用上限相符，结果必须保留未收敛。**
   `tools/continuous_fwi.py:19–37,45–56` 从公共库似然选相同 disk/ellipse 两个起点，按
   `(px,py,contrast,loss,radius)` 接上目标梯度；真墙只用于标记的 oracle，普通方法固定所选
   支持，未用真实目标初始化。`39–64` 每起点最多十次实际 full objective+gradient，重复
   坐标缓存不消耗新调用。读取的 `results/continuous_fwi/RESULTS.json` 三个所选结果均为
   `NONCONVERGED_CALL_CAP`。等调用数不是等 CPU、等网格或全局最优；oracle 无墙先验罚项，
   普通方法有 `49–54` 的同形式材料罚项，故不同方法的 objective 值也不宜直接排优劣。

6. **对称性复核使用正确的源/RX 镜像及受限未知量；修复后的空气材料合法。**
   `tools/symmetry.py:8,17–32` 同时交换 ±1.5/±2.25 源并反转 RX，墙基列成镜像对，目标仅取
   B 的 px 列；不能把该结果扩大为深度、半径或材料都被保护。`33–38` 的差分幅相 gain
   是另加的显式仪器 nuisance；本次只实际计算 f=1，不能把注释中的跨频共享表述当作已测
   跨频结果。`14` 同时把 cavity 的 eps 设 1、ell 设 0，符合 `src/r2/physics.py:83–86` 的
   严格 ambient 合同；失败的第一次 symmetry job 仍见 `FAILURE_LEDGER.jsonl:2`。

7. **E2 的参数轴和有限浮点秩标签清楚，双向信息展示有一项输出缺口。**
   `tools/compensation.py:22,40–54` 对 144 个 eps 加 144 个 ell 列分别缩放 .5/.05，五目标列
   缩放顺序与物理梯度一致。`information():14–18` 区分了数据范围投影后的目标 K 与
   proper、未截断 `N(0,I)` 切空间先验 covariance。JSON 结构记录 A0/A1/A2 分别
   64×288、128×288、384×288，浮点秩分别 64/128/288；这不是网格认证的精确物理秩。
   **新 caveat：** 第 16 行算出的墙侧 Schur `Kw` 未返回/保存；保存的是墙平均方差和
   wall-target cross covariance，不能声称已保存/审阅 `Kw` 的谱。`45–53` 的墙单方向和
   px 有限差分也不能认证全部 Jacobian 的统一误差或最小奇异值。

8. **有限补偿端点受真实非线性目标与材料盒约束，但不等于全局 profile infimum。**
   `tools/compensation.py:65–84` 和 `tools/compensation_extended.py:17–32` 的步长来自有界线性
   子问题，每个 trial 都重新计算全 Maxwell 距离，只有下降才接受；WallGrid 严格拒绝非法
   材料，没有为了保留结果而 clip 介电系数。**新 caveat：** `lsq_linear` 的 success/status
   没有写入 history；保留的是实际端点与下降，不是线性子问题收敛证据。六次/五次 cap、
   两个 ker(A) 方向及两档 secant 步长（extended `34–47`）均为局部有限构造。
   `tools/compensation.py:25–37` 的独立网格只复核被选中的合格端点，不是整个补偿曲线的认证。

9. **返回模型两次尝试的表示、度量和对照范围有明确限制。**
   `src/r2/return_model.py:122–174` 的 Q 是 32 个 Gaussian H-action 的 Euclidean FE coefficient
   QR，未做 L2 白化，也没有 oversampling；其小核/action 奇异值不能解释为物理 L2 算子谱。
   `tools/return_operator.py:67–81,212–219` 的固定 ROI L2 Gram 只用于场误差；whole-C 对照
   `src/r2/return_model.py:254–263` 使用同一个只拟合 H 的 Q，不能据其表现概括所有普通 ROM。
   `tools/return_transpose_repair.py:47–68,144–169` 只改为复对称的右转置投影，沿用既定 Q/秩。
   `91–100,117–130,180–183` 检查旧尝试、公共空间和解析 forcing，旧尝试保留。
   两个 SUMMARY 顶层均为 `COMPUTED_UNCERTIFIED`，各 64 行、0 numerical closure failure，
   cold Bayesian 均为 `NOT_RUN_SEPARATE_GATE`；“计算完毕”不能升级为证书通过。

10. **RHS/LU 与 CPU 责任没有重复收费迹象；重用、失败和验证成本仍需显式呈现。**
    `src/r2/physics.py:365–380,540–555` 在每次 LU/RHS 前独立预约，返回模块没有再次预约。
    第一轮 SUMMARY 的 708 RHS/19 LU 与 `tools/return_operator.py:48–51,334–337` 分项一致；
    修复的 256 RHS/3 LU 与 `tools/return_transpose_repair.py:40–41,131–150` 一致，其中验证
    RHS 分别 128/128，不能当成免费在线工作。CPU 以 `tools/run_job.py:19–34` 的一次
    completion 收费，局部计时不能与它相加。`tools/microstructure.py:27–35` 明确重用已付费
    benchmark，所以 RESULTS 的 `newPDE=True` 不表示每个文件都另做了一次 PDE。
    `FAILURE_LEDGER.jsonl:1–4` 保留 freeze、symmetry、初次 plot、文献调用失败；包装器
    不退还失败消耗。E3 的外层异常会写 FAILED（`tools/return_operator.py:342–348`、repair
    `240–246`）；其他脚本依赖 wrapper 失败回执，不能仅凭中间输出存在判断该 job 成功。

结论限于上述代码与已存 JSON 结构的一致性。未发现需要立即改主代码的源维度/打包错误；
第 2、5、7、8、9、10 项的报告边界应由 parent 在最终解释中保留。
