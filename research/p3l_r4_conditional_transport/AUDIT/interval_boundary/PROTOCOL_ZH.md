# 较大步长四臂适用边界：已预注册，未启动

`interval_chain_boundary_v1` 是看结果之前登记的研究范围扩展。原 P1 目标、噪声、full66 probit 先验及 likelihood-only pCN 判决不变。固定 DEV56–59、beta=.005/.01/.02、种子20261025181–183，每单元256步，共36单元、每臂9216次转移。只有4个原独立父对象，不能把半径/种子当成36个独立对象。

四个独立进程臂为 Full、Richardson、global GMRES、固定精度 global GMRES；每单元从同一原在线锚点重新开始，同种子/父对象的创新及 U 跨臂和半径重置相同。四臂均享有相同的原在线锚点缓存、接收器数据权限和 stored/canonical endpoint-once 数组合同。

每个近似臂的4个原 R0 各冷付576逆列 RHS，合计2304RHS；原锚点不变，故同父对象9个设置复用已经付费的不可变 R0。各臂独立付这笔费，不共享另一个臂的免费逆。每单元新建端点/basis 缓存，不跨材料使用 Krylov basis，不重基到接受近似点。每次真实小型复 LS 先保守计1RHS，已付 R0 的 BLAS 应用不冒充新 LU 求解。

组合硬上限300000复RHS。Full子额73728，另三臂各75424；每臂整个36单元受控上限1200秒、3GiB。顺序固定 beta、seed、parent 递增。余额不跨臂借用、不按表现重排，也不因为单元停止而扩额重跑。某臂到硬额/时间上限则明确保留 PARTIAL，余下单元 NOT_RUN；其他臂只保留自身预分配权利。异常错误停止后续启动，全部既有证据保留。预算不足不能作未知的近似判决。

两端点包络、最多4次修正和严格 native Full8源回退均沿用原协议。逆界失效或 MH 仍未解析必须 Full；固定精度臂还需两端 white-error-upper≤.05 且 Phiwidth≤.05。若付不起回退则停止，不能放松门限。普通 float64 包络不是定向区间证书。

每步立即保留候选 x/q/chi、同U、接受/拒绝、当前与候选 proof；已计算的末端场在硬停止时也尝试封存，不把 UNKNOWN 冒充接受。逐步记录真实解码材料 proposal RMS、接受后的移动 RMS、距原锚点 RMS。最终逐单元报告实际移动、接受、Full回退、冷/复用设置、全部RHS与整段外部时间。每臂时间包括冷取得、写数组、proof、封存；协调器总时间额外包括子进程启动、postclose对照、报告与所有状态/场/proof发布，不把已经计费的子作业时间重复加账。

只有显式匹配父释放并给出已关闭的 VI 控制器收据才可启动。`configs/INTERVAL_CHAIN_BOUNDARY_RELEASE_TEMPLATE.json` 当前 released:false。父线程未来完成调度后运行 `tools/run_r4_interval_chain_boundary.py pair --release <实际释放文件>`。程序不自动等待 VI 后开跑；本次只执行准备、语法/零物理门禁检查。

postclose 对所有完整/部分已判决前缀检查相同创新/U、拒绝候选、flags/x路径、所有势 proof 与固定精度条件；缺失单元不计通过。实际小LS与全网格RHS分列，原在线上下文拥有费用另披露。程序不要求 beta=.002 之外 Richardson 仍有收益，不排除普通 Krylov 解释收益，不创造后验参考、混合或 MCSE 资格。旧 as-run 源码快照不改；fast 只添加默认空的 case/步长/步数、原R0池与材料移动 hooks，VI源码未修改。
