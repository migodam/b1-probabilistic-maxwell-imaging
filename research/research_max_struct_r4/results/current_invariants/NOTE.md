# Saved-current collective invariants

No GT, inverse run, or new physical state. Parent-supplied reciprocal/real-material identities are evaluated only on saved auxiliary currents. No success threshold or theorem verification is introduced.

Retained receiver rank at relative1e-3: 17.
E.T rowspace raw relative distance=0.000232756841; orthonormal-row distance=0.00028598901; largest principal-angle sine=0.000694920262.
E.conj().T rowspace raw relative distance=0.000232756841; orthonormal-row distance=0.00028598901; largest principal-angle sine=0.000694920262.
Gamma leakage outside P Gamma P, relative Frobenius=1.12973662e-07.

| Stored current | Free family | Reciprocity defect / ||E.T J|| | Optical defect / sum sides | Minimum ||deltaJ||/||J|| | Remaining reciprocity constraint fraction | Data change / ||Y|| |
|---|---|---:|---:|---:|---:|---:|
|R2_native stage0|native_mf3|0.000815031|0.000244839|32034.1|2.21306e-05|0.0616933|
|R2_native stage1|native_mf6|0.000814936|0.000244676|3569.6|1.87967e-07|0.0402807|
|R2_native stage2|native_mf10|0.000814248|0.000244341|2639.95|1.54919e-07|0.0386269|
|R2_global stage0|native_mf3|0.000814908|0.000244899|31021.5|2.13888e-05|0.0616933|
|R2_global stage1|native_mf6|0.000814939|0.000244751|3493.67|1.90823e-07|0.0402807|
|R2_global stage2|native_mf10|0.000814089|0.000244322|2716.63|1.49865e-07|0.0386269|
|R3_noisy_strong_release stage0|native_mf3|0.000815031|0.000244839|32034.1|2.21306e-05|0.0616933|
|R3_noisy_strong_release stage1|native_mf6|0.000814936|0.000244676|3569.6|1.87967e-07|0.0402807|
|R3_noisy_strong_release stage2|released_DPCS_mf10|0.0100317|0.00322607|0.0013117|2.64919e-15|0.00663076|

Each correction uses an orthonormal physical-current basis and the66-complex-constraint Gram pseudoinverse. Its reported minimum applies to the retained numerical Gram range; discarded constraint components and actual post-correction defects are preserved. Basis orthonormalization changes only this diagnostic norm representation, not the historical CIE optimization.

All four constraint Grams retained numerical rank66. Their largest/smallest singular-value ratios are approximately1.069e12,4.785e9,3.848e9 for native mf3/6/10 and207.0 for the released family. Thus the remaining fractions in the table are numerical solve residuals, not evidence of an exact incompatible constraint component. The actual post-correction reciprocity ratios are about1.9e-8 at mf3 and1.4–1.6e-10 at mf6/10. Original native/global optical ratios near2.44e-4 increase to about.00380/.00170/.00164 after the corresponding reciprocity-only correction; it does not enforce the optical constraint.

Native APCS corrections preserve the locked DPCS coordinates to floating arithmetic, while weak/unretained receiver channels can change. R3 final released-stage corrections use its actual free[Vd,APCS] family and therefore can alter locked data components. These controls are not interchangeable.

Small rowspace/Gamma leakage means these quantities are largely visible in the retained measurement subspace for this acquisition; the measured minimum correction describes the cost of using the remaining freedom. Neither fact alone establishes feasibility of the full real-material state equation or an improved reconstruction. Parent owns the algorithm decision.

All detailed norms, optical changes after reciprocity correction, numerical ranks, SVDs, Gram solves and timings are in RECEIPT.json; matrices and corrections are in ARRAYS.npz. New physical LU/RHS and inverse runs are all zero. No source/current coefficient tuning or material scoring was performed.

Measured diagnostic wall=.741666583s and CPU=.719651s. Costs include one receiver SVD, four row-space SVDs, three native APCS SVDs plus three diagnostic orthonormalization SVDs, one released-family SVD, four66×66 constraint-Gram SVDs and nine constraint RHS solves. Shared families are reused between streams; R3's first two stored stages reproduce the native saved currents and are not independent evidence.
