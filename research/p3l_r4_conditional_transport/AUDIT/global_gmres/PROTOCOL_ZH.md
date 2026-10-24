# 共同R0权限的普通global GMRES强控制：只准备

新增src/p3r4/decision_refresh_global_gmres.py、tools/run_r4_interval_chain_global_gmres.py。fast runner仅增加缺省为空的opt-in类/配置/resolver/协议扩展hooks；原fast默认路径与全部封账pair结果不变。没有启动新的物理或训练。

旧SciPy控制每次LU预条件callback付费，不能代表拥有已付R0矩阵缓存时普通Krylov的最佳控制。本修订让Richardson与global GMRES各自冷付576逆列RHS，再享有相同R0/G/S八源BLAS权限。R0 matrix apply不是新LU/RHS solve，不再错收source预条件费用；亦不从这些普通数值操作推导A17或current-space理论优势。

GMRES把N×8复电流矩阵看成C^(N*8)中的一个向量，内积为Frobenius共轭内积。初值与Richardson相同jinit=j0+R0diag(Δχ)E0。设ρ0=diag(χ)einc−Lχjinit，z0=R0ρ0；Arnoldi算子K(V)=R0[V−diag(χ)GV]。两遍modified Gram-Schmidt生成至多4个方向，最小化||βe1−Hbar y||，jhat=jinit+∑y_iV_i。它是普通left-preconditioned global GMRES，不是每源8套独立GMRES或block方法，也没有新物理查询。

只缓存当前与候选两个物理端点的basis/最多4次解，接受后保留同端点已有精度，拒绝候选的basis可释放；从不跨不同χ偷用basis或做状态自适应proposal。相同χ的已付小解可以重放而不重复收费，新LS每次必须经book.reserve_rhs(1,'samecache_global_GMRES_reduced_LS')保守计1个complex RHS，再实际执行≤5×4的一右端复least-squares；这不是576维FullLU。记录H/系数real-imag、真实小shape、纯数值LS时间、含ledger的callback时间、所有G/S/R0列动作、内积/copy/basis峰值及缓存hit。

提案仍full66 pCN beta=.002，所有方法独立进程、共同初值/创新/U。前瞻parents DEV48…55，每arm每parent1024步，seed20261025161。四arm为：

- optimizedFull：共同stored codec、Target native endpoint只引用一次。
- richardson：已付R0、原decision-drivenmax4修正。
- global_gmres：相同R0、max4 globalArnoldi、同U两端势区间决策。
- global_gmres_accuracy：当前和候选均需white_error_upper≤.05且Phiwidth≤.05才作判断；达到精度后MH若仍不确定，继续剩余迭代或付费Full。

所有arm均使用同一stored/canonical endpoint-once合同、原始R0和完整Full fallback，不删拒绝/失败。中间误差、Φ包络、LS proof标量append-only；最新current/field数组和Gaussianizedx路径封存。没有把接受近似点冒充精确cache。float64上界仍非定向interval证书。

组合硬上限300000RHS，固定无相互读取的子预算为Full65536、Richardson78154、globalGMRES78155、accuracy78155，每arm1200秒/3GiB。预算包括每arm冷R0、全部8源Full和实际reducedLS保守RHS。任一子预算耗尽即STOP，不借未声明费用或多加迭代救结果；不会因为剩余组总额多而暗改分配。

父调度器release：

```json
{"released":true,"experiment":"interval_chain_global_gmres_v1","campaign_id":"P3LR4-20261025"}
```

入口`<python> tools/run_r4_interval_chain_global_gmres.py --release RELEASE_PATH`依次启动4个独立物理arm及第5个postclose compare进程，只有新子进程启用stored。显式单arm也支持，需R4_ARRAY_CODEC=stored；正式wholewall比较使用协调器的完整外部进程计时。没有默认付费启动。

全部4Books CLOSED COMPLETE、PID不同、相同输入及数组身份后，Gold Full逐项审计全部flags/x paths、最终两端Φ区间和每个中间proof；accuracy另检查每次最后决策的双端精度条件。任何不一致NO_GO并停止。零miss只给PASS_SAME_KERNEL_FOUR_ARM_NUMERICS，时间/RHS/Fullfallback/小LS代价原样报告，由父调度器判断普通Krylov是否已解释收益；代码不要求Richardson胜，也不以精度控制更贵就过度宣传decision方法。没有posterior reference、混合或MCSE资格结论。

GLOBAL_GMRES_MATH_CHECK.json：任意14×14复线性系统、8 RHS、1…4次globalGMRES与独立SciPy flattened-system结果差约1.5e−15；40步双端accuracy-first mock审计0判决错、0accuracy违规、basis最多2端点，拒绝未计费LS。全部物理RHS/queries及训练为0。真实物理控制尚NOT_RUN。
