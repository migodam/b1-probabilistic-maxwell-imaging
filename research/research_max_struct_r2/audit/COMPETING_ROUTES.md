# MAX-STRUCT-R2：竞争路线与机制失效审计

2026-11-02。证据工作而非最终数学裁决；不实现solver、不跑sweep、不改R1。已读R1完整 `THEORY_DERIVATIONS.md`、`BASELINE_AUDIT.md`，D1/D3真实代码、R3 CIE实现/审计、R1源清单及相关Atlas限定快照。以下公式为供主线程审核的既有结构在CIE上的核对，不以再次做数值恒等式测试建立新颖性。

## 路线A：精确 current profile

固定material m、β_current、current表示 `J=J_D+Q A`。CIE状态残差可以写

`r_s(m,J)=P_β(m)[(I−D_mG)J−D_mE]`，`P_β=diag[β/(1+βm_i)]`。

对同一个有限罚二次目标加共同 `R(m)` 精确最小化A，是标准variable projection。已有Atlas06 F-OPT-003、Atlas10 F-BI-022以及文献PRIOR_ART中的penalty/VP原文直接覆盖。**新数学不能只是伪逆、envelope gradient或Hessian Schur补；本轮可能新增的是CIE合同中的可审计专门化及反例。**

必须冻结同一个目标才有profile。如果β_material与β_current不同，历史双βreadout不再是该profile的材料下降步。把双β改为一个β是新对照/目标选择，不是无声“修复成熟基线”。阶段切换β或current mask也改变目标/可行集，不能把跨阶段数值当一个固定函数的下降史。

`P_β(m)`随m变化。离物理根时梯度必须含 `dP_β` 乘原状态residual，不能只对原本构算子求导。R3已有完整off-manifold导数；R1的根处β一阶抵消不能外推到离根、prior、clip或单PR迭代。精确内层常秩时envelope导数中dA*项消去；不精确内层需控制正规残差及引起的梯度误差，不能把一次PR结果当作精确profile。

**维数/成本**：若Q为独立q列，普通增广QR复杂度约O((N+n_r)q²+q³)，每个新m都更新状态块。全current q=N时并非免除大系统；普通强基线的原物理求解本来共享一个LU、多源RHS，不应显式构造完整逆作弱对照。用normal equations会平方条件数，须在精度/成本间明示选择，不能将Cholesky一概当严格最优。

原Fourier APCS保留native系数但含冗余列时，A可不唯一；物理current profile可仍唯一。此时不能直接套“系数矩阵满列秩”的证明，需换成独立current基、常秩伪逆或明确商空间。若对PR基线正交化Q却保留未运输的欧氏CG度量，会更换基线算法。

## 路线B：全current profile的有效数据度量与认证缺口

当CIE状态H_m可逆，令物理current满足 `H_m J_phys=b_m`，写 `J=J_phys+H_m⁻¹Z`。白化后data residual记e，`K_m=W S H_m⁻¹`。同一source metric兼容时profile有已知形式

`Ψ_full−R(m)=½ <e,(I+K_mK_m*)⁻¹ e>`。

这是WRI/reduced penalty的标准结构；vanLeeuwen–Yang2024v1 Cor1 Eq10直接写出对应variable metric形式。**不是新增观测或新likelihood**；若解释成对隐状态做Gaussian积分，还会涉及determinant，profile与marginal不能混淆。

它也解释可能的失败：沿K_m的大奇异方向，原物理residual e会被强烈下权。仅证明profile很小/gradient很小，不保证原数据residual小，更不保证material正确。单点有

`||e||²/[2(1+||K_m||²)] ≤ Ψ_full−R(m) ≤ ||e||²/2`。

要由小profile推出原物理门，需要**该candidate处**的K上界（或直接原物理解验收）；要给全域保证需统一resolvent稳定性。当前导数或有限几个端点不足以控制近共振全域。stationary只关gradient，小gradient不能代替这个residual条件。

把resolvent certificate转为数据域small Gram可能节省维数，但形成K或其Gram需要状态/adjoint求解及有效误差包络，不能只报n_r阶逆。R3已有物理缺陷/校正身份，本轮不为再证它追加公式误差底测试。

## 路线C：固定仿射current限制与容量

对固定m和同一个二次目标，若完整current块Hessian `B=S*W*WS+H_m*H_m` 正定，完整极小值J_full*唯一。仿射集合内极小值J_aff*满足精确二次差

`Ψ_aff=Ψ_full+½||J_aff*−J_full*||_B²`。

这是quadratic projection/Pythagoras，Osborne2007 pp2–3 Eq1.3–1.5有直接分解先例。可以作为本研究CIE响应解释，但不宜称新一般定理。

在J_D一致且集合真正嵌套时，增加current空间不能提高profile。**这不是material改善的单调性**：更大current自由度可以更容易吸收物理不一致，改变material objective的最小点。不同affine offset、different DPCS锁定、differentβ不构成这个嵌套比较。最有判别力的证据应将profile下降、原物理验收、material变化和全成本并列，而非仅比辅助loss。

若证明使用B正定但H_m接近奇异、Q冗余或测量存在完全不可见方向，须分别验证所需条件；physical forward可逆是H_m相关CIE变换可逆的条件之一，β pole-free也需保留。

## 路线D：global source-covariant更新

完整酉源变换 `E'=EU,Y'=YU,J'=JU`，两个右侧metric分别运输 `Ω_d'=U*Ω_dU`、`Ω_s'=U*Ω_sU`；全局Frobenius内积下梯度、PR scalar、搜索方向、单全局line coefficient及历史重置一起运输，才保持完整迭代等变。该机制已经是R1 D3，不是R2新物理方向。

逐源PR或逐源不同步长定义了一个依赖当前源基的坐标分解，即便objective本身协变也会不同；R1首次梯度步已出现差异，所以不能把差异只归因PR历史。global更新可作为合理稳健baseline，但sourcewise能用不同步长适应不同源的曲率；global共同步长也可能被最难源拖慢。等变性本身没有性能序关系。

一般非酉可逆source变换时，仅运输objective右metric不保证欧氏gradient按J的坐标映射运输；还要运输变量空间内积/预条件。R1/R2所证范围应先保留在完整酉变换；不要把“source-covariant”扩成任意缩放都无需改metric。

**精确current求解的关键成本**：有两个full source metrics时正规方程一般为

`(S*W*WS) J Ω_d + (H_m*H_m) J Ω_s = RHS`。

它是generalized Sylvester结构。只有Ω_d=cΩ_s等共同缩放情形可无条件共用一个空间factorization和所有源RHS。即便两者都是对角阵，wd_s/ws_s不同也产生不同加权空间矩阵。两个SPD source metrics可作generalized eigendecomposition/congruence以解耦为若干加权源系统，但这是标准线代，须支付source分解、状态更新和不同空间factorization。不能将“六源”一概宣称只多六RHS而不改变factorization；R1 current_exact事实上逐源计了Cholesky。

## 路线E：receiver nullspace与数据约束

SOM/TSOM已有DPCS/APCS机制。只有**精确**接收核N满足SN=0；截断接收奇异值后 `(I−P_ell)Q` 往往含弱辐射方向，SQ不为零，不能删除data term。把含噪强模态固定为J_D通常也改变优化问题；Atlas10 F-BI-021已提出显式修正强模态以保留同一data likelihood。

任意receiver暗电流不等于合法material。必须满足所有源共用m、state equation、real/passive/box/prior；标量某行E=0还必须J=0，minors alone漏掉这个例外。将接收暗空间直接当material fiber会重犯旧失败。

所谓低阶端口精确消元通常产生随m变化的左零空间或Schur块；其构造不免费，所需块可逆性不由整个physical matrix可逆推出。Atlas06 F-OPT-004已有明确反例。固定接收空间不能一般冻结内部传播与多重散射。

## 竞争矩阵与停止规则

|路线|已有部分|本轮值得分辨的量|主要失败/成本风险|
|---|---|---|---|
|exact current VP|标准VP；R1/Atlas已写|同目标physical结果、总wall、内层正规残差、外层梯度误差|全current昂贵；不同源metric需不同factor；过解内层|
|controlled-inexact VP|已有partial minimization文献|在同误差/同预算下是否已解释exact收益|不能只对一PR弱baseline；容差选择必须公开|
|global PR|R1已验证等变；标准矩阵优化|对同noise/metric、同stop/重置的actual material结果|等变不增信息；common step可能慢|
|receiver-nullspace constraints|已有SOM/TSOM|是否保留data与state合法性、实际capacity|弱模态当严格核；DPCS硬锁含噪量|
|effective data metric|WRI及2024 Cor1直接已有|何时小profile隐藏大physical residual；resolvent验收|近共振K大、成本漏记、profile误称likelihood|
|affine capacity gap|标准二次投影|增大空间降低辅助loss却损害material的反例|非嵌套/变β后仍套单调比较|

若只有标准结构的重写或一份公开开发小模型，应交付限定的机制审计/反例，停止“普适新反演器”“新全局保证”宣传。公开锚点的新增真实测量可以验证**声明候选族内部**的分离，却不证明GT包含于该候选族；candidate-coverage miss须独立保留。任意给定source使局部本构核被破坏，也不自动保证测量分离超过噪声。

## 读源定位与未执行事项

R1 `THEORY_DERIVATIONS.md`全文，重点§4–6；`BASELINE_AUDIT.md`全文；`src/dynamics_audit.py:40–140`；R3 `src/cie_tsom_baseline.py:52–104,139–171,185–292`、`CIE_TSOM_BASELINE_AUDIT.md`全文；R1 Atlas快照中F-OPT-003/004/005/009与F-BI-020/021/022；`audit/D_INPUT_PROVENANCE.json`。

没有打开新的隐藏GT、没有重跑R1/历史结果、没有新numerics或sweep。文献证据与失败详见 `../literature/PRIOR_ART.md`。最终核心证明、方向判断及发表强度由主线程负责。
