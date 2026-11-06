# Scalar contraction versus Green feedback derivative

This is a zero-forward, zero-adjoint decomposition of the same ten saved endpoints. Existing START arrays supply w, a and the original Phi gradient. No new truth read or candidate selection occurred.

The scalar term dominates in norm at all ten states: ||g_beta||/||g_K|| ranges 2.288 to 5.390. At every mf3 endpoint g_beta opposes the constitutive gradient, whereas g_K aligns with it. Thus the strong early cancellation previously observed in g_metric is mainly scalar CIE contraction, not Green-specific cancellation. This rejects a Green-only shielding interpretation of that evidence.

|State|beta/K norm|cos(const,beta)|cos(const,K)|fixed-T projected descent/Phi descent cos|
|---|---:|---:|---:|---:|
|r10_twod_0001_s0|3.787|-0.917|0.596|-0.353|
|r10_twod_0001_s3|4.266|-0.283|-0.724|0.059|
|r10_twod_0004_s0|5.054|-0.911|0.328|-0.408|
|r10_twod_0004_s3|4.725|-0.425|-0.095|0.015|
|r10_twod_0022_s0|4.316|-0.972|0.849|0.160|
|r10_twod_0022_s3|2.497|-0.146|0.379|0.102|
|r10_conditional_0000_s0|5.390|-0.951|0.789|-0.342|
|r10_conditional_0000_s3|2.288|-0.205|-0.722|0.126|
|r10_conditional_0003_s0|3.328|-0.824|0.691|-0.430|
|r10_conditional_0003_s3|2.868|-0.314|0.791|-0.023|

The separate g_beta and g_K sum to the prior g_metric to floating-point roundoff. Each term has a saved norm, raw cosine with the original Phi gradient, work for a unit negative-gradient step, and projected-gradient work. These are local directional diagnostics, not finite-step Phi guarantees.

Fixed-T means z=m/(1+beta m), T=I-D_z K, and w=T^-1(J-D_z(Q+beta J)). Its mapped-m gradient is -rho beta² Re sum conj(a)(Q+beta J)/(1+beta m)+g_prior with a=A^-*(w Omega). It is neither freeze-A nor simply g_const+g_beta. The fixed-T projected descent points against physical descent at four of five mf3 endpoints, including the annulus and shell. Finite contracted-readout candidates remain a separate authorized control.

Cost: zero new physical FW or adjoint; final calculation used ten low-rank Green actions. One preliminary coding attempt used an incorrect simplified fixed-T gradient and was corrected before handoff; total including it is twenty low-rank actions / 240 auxiliary complex RHS. g_beta/g_K values were unaffected.
