# R4 working formula collection

Only derived computational formulas are listed. Conventions and complete assumptions/proofs are in [MATHEMATICAL_THEORY.md](MATHEMATICAL_THEORY.md); historical provenance is in [FORMULA_AND_CODE_REUSE.md](FORMULA_AND_CODE_REUSE.md). Fixed geometry, independent real scalar materials, and invertible finite VIE states are essential throughout.

| ID | Derived formula | Computational consequence and actual status |
|---|---|---|
| R4-F01 | $F(m+h)-F(m)=SM_{m+h}^{-1}D_hQ_m$ | Exact finite nonlinear material response; explains why a finite update changes receiving resolvent as well as material. Used in theory and the observed finite-step diagnosis. |
| R4-F02 | $E_m(h)=SM_{m+h}^{-1}D_hGM_m^{-1}D_hQ_m$ | All higher-order error, with multiplication order fixed. Actual finite-step remainder measured; no derivative precision experiment. |
| R4-F03 | $\beta_F^2=\sum_i h_i^2\|(GM_m^{-1})_{i,:}\|^2$ and $\|E_m(h)\|\le\|SM_m^{-1}\|\beta_F\|D_hQ_m\|/(1-\beta_F)$ when $\beta_F<1$ | A genuine sufficient ellipsoidal step certificate. Dense preparation costs $N$ RHS plus receiving factors. Tested once as an algorithm qualification experiment; conservative, not the global bridge gate. |
| R4-F04 | $g_i=\Re\sum_s\bar\Lambda_{is}Q_{is}/\|Y\|^2$, $M^*\Lambda=S^*(F-Y)$ | Source-adjoint data gradient: $s$ RHS after the forward. Implemented for every data L-BFGS correction. |
| R4-F05 | $Hh=\Re\sum_s(\overline{\dot\Lambda}Q+\bar\Lambda GU)/\|Y\|^2$, $U=M^{-1}D_hQ$, $\dot\Lambda=M^{-*}(G^*D_h\Lambda+S^*SU)$ | Exact residual curvature at $2s$ RHS per action, sharing the LU. Implemented in the three-endpoint Hessian diagnosis. |
| R4-F06 | $F(m+tp)=F(m)+tS(I-tM^{-1}D_pG)^{-1}M^{-1}D_pQ$ | Exact rational ray queries with one dense kernel/Schur preparation. V1 and V2/V3 uniform ray use it; all full-rank setup and triangular work charged. Finite scanned minima are not all roots. |
| R4-F07 | $F(c\mathbf1)=c(SU)(I-cT)^{-1}(U^*E_0)$ for $G=UTU^*$ | Special uniform-material family can be cached across targets. Executed implementation charges a cold reference-$m=1$ setup per target; direct shared $G$-Schur cache is a theoretical optimization, not a measured speedup. |
| R4-F08 | $DT(m_*)=I-J(m_*)^\dagger J(m_*)$ at an interior constant-rank exact root | Observable normal directions contract locally; null directions remain neutral. No entrance-to-basin or GT guarantee. Mathematical analysis only. |
| R4-F09 | $F(n)-F(m)=B_{n,m}(n-m)$ with $B_{n,m}h=SM_n^{-1}D_hQ_m$ | A uniform two-point singular-value lower bound would imply global inverse stability. A current-point Jacobian does not supply it. Reuses R10 bilateral structure; not certified in the 576-cell experiment. |
| R4-F10 | $\psi=\operatorname{tr}[R(Q^*Q)^{-1}R^*]$ | R3 source-field weighted objective, implemented as a comparator. It changes nonzero-residual geometry, supplies no new measurements, and is never used to certify a data solution. |

The block-Krylov alternative in [theory/RESOLVENT_GEOMETRY.md](theory/RESOLVENT_GEOMETRY.md) gives a $2s$-dimensional common rational proposal from at most $3s$ physical RHS and a residual-defect identity. It was derived but **not run**. No speed or GT claim is attached to it.

R2 Q6/Q20 were deliberately not invoked: this sprint did not derive a useful full-observation block-selection criterion for a wrong complement. That negative value gate prevents repeating the earlier unproductive root-count experiment.
