# Heun revision1：先修复数值commutator，再谈任务幅度

父线程已明确批准：首轮32对象封账后，固定前8个新DEV进行单次机制修复；不按GT或首轮成败选对象。独立job `III/dynamic_geometry_heun_r1`；物理上限168,192 RHS（低于父线程250,000额外cap），1800秒。第一轮所有数组、成本和失败保留。

## 组合展开与修复理由

以下主项由父线程独立数学复核。对C¹向量场，Euler `E_h^F(x)=x+hF(x)`，而局部流有二阶Taylor项 `.5h² DF·F+o(h²)`。四Euler组合 `E_-h^Y E_-h^X E_h^Y E_h^X` 的端点为

`x+h²([X,Y]−DX·X−DY·Y)+o(h²)`。

交换X/Y得到负bracket，但相同self缺陷不消失。因此首轮正反Euler回路都不是一致的bracket提取器，其task-RMS失配是实现诊断，不是theorem反证。

Heun单段先取 `p=x+hF(x)`，再取 `H_h^F(x)=x+h/2(F(x)+F(p))`。Taylor代入给 `x+hF+.5h² DF·F+o(h²)`，补齐Euler缺少的self项；四段组合在相应正则条件下恢复 `±h²[X,Y]+o(h²)`。只有在局部足够光滑且导数受控时才加强为O(h³)。P1 seam下场通常只有C¹，协议不承诺全局O(h³)，实际小h收敛必须检验。

Heun负步仍不是正步的精确逆；此修复不自动给可逆数值map、Jacobian determinant、MH核或posterior资格。每段predictor和corrected endpoint都是实际付费field evaluation点，有限采样也不构成连续路径证书。

## 冻结机制

固定R4DEV_00000–00007，64/66×rank2/4/8，沿用首轮同一pair的Pci材料归一化场、seed20261025030和base。参考bracket用已付h=.001 centralFD，首轮已付三个h合计480 RHS/config不隐去；所选48配置对应23,040 RHS历史共享费用，仍归属首轮ledger。

正反四段每条Heun曲线需8个Full states+7个J，共288 RHS。先运行h=.001/.002/.004×两方向。全48配置small-h预算82,944 RHS。每配置须满足：两方向fine h的材料scaled-displacement相对参考bracket误差≤.2；fine–mid差≤1.25×mid–coarse差（绝对floor1e−8）；无材料rank变化或relative-gap≤.001 crossing标记。seam单列，不借seam标签假定高阶光滑。

只有通过门的配置才运行.02/.04/.08两方向Heun任务曲线与新付费direct±bracket。direct用真实decoder material RMS匹配；Heun用 `s=sqrt(radius/||D bracket||)` 后报告实际距离，不调步长救结果。最多另85,248 RHS；总168,192。失败配置不增h、不调seed、不加refinement。

局部teacher仅在实际半径1%内、沿对应bracket材料方向进展≥半radius、Full ΔPhi≤2、prior能量增加≤2、原prior合法及backward residual≤1e−10时标eligible。这个资格只描述无GT的局部计算监督，不是posterior样本、逐例GT风险或部署收益。

两项纯合成测试已经PASS：可交换线性场的Euler h² self漂移而Heun normalized误差消失；非交换线性场的Heun signed bracket收敛及正反方向。测试不是Maxwell状态，不计物理RHS。

冻结执行合同：`configs/GEOMETRY_HEUN_PROTOCOL.json`。脚本必须验证第一轮CLOSED COMPLETE，源代码快照在付费前保存；所有gate失败、部分路径和paid endpoints均保留。后续cost归并还须把旧reuse的小系统6144 complex RHS计入总cap，首轮raw8003仍只表示full-grid RHS。
