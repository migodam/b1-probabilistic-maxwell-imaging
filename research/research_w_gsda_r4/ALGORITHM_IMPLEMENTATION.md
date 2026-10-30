# 实现范围与必要公式

本轮使用静态二维 scalar-TM，exp(-iωt)，epsilon(f)=epsilon_infinity+i loss/f。源为iH0^(1)(kr)/4，FEM仅离散散射场。物理目标是全波测量与真实材料图，不将小残差当作材料成像质量。

## 全波主基线 B0

在公共3×2墙格和7×7目标格上，110个实材料参数均未知。每个频率装配：
L(chi) = L0 - k² Pᵀ W diag(chi) P,
b(chi) = k² Pᵀ W diag(chi) e_inc,
L u_sc = b, d = d_inc + R u_sc。
因此db不能丢失。P为真实P3插值，W为实际材料积分权重。所有源与频率共享epsilon_infinity/loss。

数据目标为0.5 sum|d-y|²/sigma²。对每个源一个加权伴随lambda=L^(-H) R^H(d-y)/sigma²，材料实梯度为Re sum_s k² W conj(P lambda_s) (e_inc_s+P u_s)。loss方向乘i/f后取实部。主基线加入相同尺度TV、材料box和墙properGaussian正则，采用逐频continuation和两个公开均匀起点。最终仅按完整三频测量目标+相同先验选支。

## 原生固定GS电流表示 B1

j(q,s)=chi(q)[e_inc(q,s)+GD0 j(q,s)]，
GD0 j = P L0^(-1) k² Pᵀ W j，
GS0 j = R L0^(-1) k² Pᵀ W j。
y_sc=GS0 j。L0只依赖公开参考介质/网格/PML/频率，迭代中固定；材料Jacobian不固定。有限元积分点上的j为物理contrast current，不是R3低秩forcing chart。

## 精确吸收 B2

材料A用LA=L0-VA及bA求背景。uA=LA^(-1)bA，e_eff=e_inc+P uA。剩余电流jR的场为uA+LA^(-1) k² P_Rᵀ W_R jR。原始观测同时保留背景响应、有效照明、目标返回墙再到接收器的路径；也可重构jA与jR后由原始GS0测量。

## 内部谱 B3 与可逆动作 B4

内部GD的数值模态用实际物理Green actions构建。投影后的M=I-chi GD分成强H/弱L块；精确Schur保留M_HL M_LL^(-1) M_LH。投影内恒等不代表完整空间准确，必须另与完整FEM比较。降低rank不自动代表省FE/RHS/CPU。

吸收仅暂停对应材料更新；保留材料状态与恢复路径。出现墙梯度/全波目标冲突时REOPEN。局部不确定性和风险只作为经验筛选，未建立完整posterior或严格安全证书。

## B3实际做的是全电流空间的模态Schur预条件

把电流改写为w=√Wj后，内部算子C=√W GD W^(-1/2)。从这个真实Green action的右奇异向量取得公开Q（32列）。每次材料变化都重新计算Mq=I−Q*diag(chi)CQ；强8/弱24的Schur补偿只用于求Mq的逆Z。预条件器为D=I+Q(ZZ*−I)Q*。Q外的全部电流方向仍然保留；没有删去完整空间的弱电流，也没有把材料锁成弱坐标。若弱块或核心不稳定，记录并退回普通方向。

实际方向是物理加权梯度的负SPD作用，吸收区域用两侧相同活动投影处理；在当前固定材料下，用完整GS/GD计算精确二次线步。这个Schur预条件与普通同空间方法相同，B5使用完全同样权限。B3的可检验增量来自是否能减少达到同一目标质量的实际工作，而不是Schur恒等式本身。

## 经验目标风险筛选的身份

B4先看材料变化、条件状态LS曲率、完整数据目标梯度和冻结场的GS影响。仅有候选时，额外付费构造完整材料Jacobian，消去剩余墙nuisance，并估计δx≈−Sx^(-1)BᵀQ_R A_A δw_A。这里的材料标准差只来自条件曲率代理；活动box、全非线性及真实posterior并未认证。所有记录保留LOCAL_UNCONSTRAINED_GAUSSIAN_TANGENT_WITH_CONDITIONAL_SD_PROXY，禁止由此声称目标全局安全。

## 公共参考复用

九个parent共享公开网格、源/接收器和频率。因此空气参考L0、P/W、GS0、GD模态可精确复用，材料、当前总场、电流、分支选择和测量仍逐case重新初始化。batch将一次冷setup和每case的warm开销分开；普通B5可用相同缓存。独立部署仍从空缓存构造，历史bank不提供免费物理信息。
