# DRAFT：固定材料纤维下的边缘乘积反例

**状态：父线程审查草稿；不是新Cauchy–Schwarz定理或原创性声明。物理脚本已准备但实施worker没有执行；尚无新Full响应结果。** 只讨论原A固定S3均匀材料、截断prior，V=QD。两个v来自已经保存的GH21 material likelihood mode，不是GT。原P1属性区间严格为正，所以本例vR>0；若推广到vR<0，下述协方差符号反转，不能不加条件写成≤0。

给定V=v≠0，δR=vR/Q。设m=E[Q|V=v]、h=E[1/Q|V=v]；原prior fibre在此非退化，且Q有正下界。则

    Cov(Q,δR|V=v)=vR−m·vR h=vR(1−mh)≤0。

最后一步是标准Cauchy–Schwarz：E[Q]E[1/Q]≥(E[√Q/√Q])²=1；Q非常数时严格。因此，即使两个一维边缘**分别完全精确**，独立耦合后的材料均值仍变成vR·mh>vR。复材料增量的均值按同一实因子mh径向偏移。这个反例来自丢掉同材料gauge内的copula，并没有制造第二个物理材料世界；“两个精确边缘”也不是“一个正确联合posterior”。

若额外假设某局部材料后验族随σ→0集中于非零v*，且相应条件prior-fibre矩在v*连续（有界Q和有界属性保证所需矩可控），则真材料均值→v*、材料方差→0；分开取Q与δR边缘后再独立相乘的均值趋向vR*·E[Q|v*]E[1/Q|v*]，只要极限纤维仍非退化，偏移不消失。这里**假设**材料集中，没有证明当前Maxwell模型全域可辨识、global tail受控或Full66后验存在同样极限。不能把局部GH15/21的数值稳定当成这些全局前提的证明。

已有POSTHOC_GL128 mode条件矩给出以下代数数值，均来自旧节点，**不是待运行物理试验的响应结果**：

- REGULAR: vR=0.1184110978，mode条件Cov(Q,δR)=-0.0026381958721，所以精确边缘独立乘积的实材料均值偏移=0.0026381958721。
- PURE_OVERLAP: vR=1.20196737309，mode条件Cov(Q,δR)=-0.0180278528714，所以精确边缘独立乘积的实材料均值偏移=0.0180278528714。

实际拟运行control先从GL128 transition fibre（权重φ(z)/(Z0 Q² rectangle_area)）与pure Q=1质量中取32个配对(Q,D=v/Q)，然后只作一次固定seed置乱Q相对于D。两组经验边缘完全保留，但单个有限置乱的平均乘积不要求恰等于经验边缘均值乘积；应单列两者差。不同置乱输出不是iid posterior draws，不能据此做Γ校准、覆盖率或posterior真实性判断。若纯相被抽到，其z仍取原截断Gaussian pure-tail条件分布，不能把pure源质量丢掉。

脚本 `copula_counterexample.py` 只对两个joint reference及最多64个factorized材料调用原Full forward，禁止adjoint，最多66×8=528 RHS；只对字节相同的decoded材料复用cache，不四舍五入合并不同材料。Full response差使用(Ffactor−Fjoint)/σ，并同时保存complex L2/RMS及sqrt2实白化范数；它是控制材料改变造成的物理响应差，不是观测校准统计。实际费用必须由root执行后的Book receipt填写，当前A总费用仍为157696 RHS，不能预记528已发生。
