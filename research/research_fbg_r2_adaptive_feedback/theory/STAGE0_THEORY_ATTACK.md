# Stage 0 algebra and claim-boundary audit

**Scope.** Finite-dimensional complex contrast-current equations, with a real diagonal material vector and fixed incident fields. The companion `experiments/stage0_algebra.py` checks one deterministic small matrix instance. This note is a mathematical attack on a candidate feedback design, not a Maxwell experiment, material-recovery result, or novelty ruling. The actual scalar-TM implementation uses `j=(I-XG)^{-1}XE`, `Q=E+Gj`, and `y=Sj` (`research_z_alg_r3_digits/src/physics.py`, `forward`). The R3 small-Gram proposal is `psi_0=tr[R(Q*Q)^{-1}R*]` (`research_z_alg_r3/MATHEMATICAL_THEORY.md`, R3-05); R4 uses it with a separate acquisition-dependent scale and a material regularizer (`research_z_alg_r4_basin/src/basin_solver.py`, `correct`). Here `*` denotes conjugate transpose.

## 1. Exact gauge, with the source convention fixed

Write `X=diag(m)`, `L=I-XG`, `j=L^{-1}XE`, `Q=E+Gj`, and `R=Sj-Y`. The *incident* field `E` is fixed. For a fixed matrix `B` with `A=I+XB` invertible, define

```
Xtilde = A^{-1} X = X (I+BX)^{-1},
Gtilde = G+B,
Qtilde = (I+BX) Q.
```

`Qtilde` is the **new solved total field**, not a replacement incident source. Indeed,

```
Ltilde = I-Xtilde Gtilde = A^{-1} L,
Ltilde^{-1} Xtilde E = L^{-1} X E = j,
Qtilde = E+(G+B)j = Q+B X Q,
Xtilde Qtilde = j.
```

Thus `Sj`, `R`, the original data loss, and its physical material solution set are identical at the same `m`. This is an exact change of representation on the joint invertible domain, not a new physical acquisition. A mistaken substitution `E→Qtilde` on the transformed incident right-hand side breaks the identity. The equality `A^{-1}X=X(I+BX)^{-1}` does **not** permit exchanging `X` and `B`: `(I+BX)^{-1}X` is generally different. A dense `B` makes `Xtilde` dense even though `X` is diagonal; `Xtilde` is a reparameterized nonlocal operator, not a new local material image.

For fixed `B`, the derivative of the transformed material operator is

```
dXtilde = (I+XB)^{-1} (dX) (I+BX)^{-1}.
```

Omitting either factor or pretending the derivative remains diagonal is generally wrong. The map requires both `I-XG` and `I+XB` invertible. `I+XB` and `I+BX` have the same determinant, but their singular-value conditioning need not match.

## 2. Feedback Gram is a changed proposal metric

Let `F=Qtilde=(I+BX)Q`, `H=F*F`, and assume `F` has full column rank. The proposed feedback objective is

```
psi_B(m) = tr[R H^{-1} R*] = ||R H^{-1/2}||_F^2.
```

At `B=0` this is the R3/R4 unscaled small-Gram proposal. At any admissible fixed `B`, `psi_B=0` exactly when `R=0`, but away from zero it changes the material-space landscape. The physical data likelihood and real material box do not change. Consequently a lower `psi_B` at the same material does not establish a better fit or recovery. If `M=I+BX`, then, for the same `R`,

```
psi_0 / ||M||_2^2 <= psi_B <= ||M^{-1}||_2^2 psi_0.
```

The bounds show how transformed-field conditioning can create a large apparent objective change without a physical residual change. They are pointwise and do not prove useful basin changes.

The displayed `psi_B` is a **defined transformed-field metric**. R3 derived `psi_0` by eliminating an auxiliary all-receiver constitutive residual under its local diagonal material equation. An arbitrary dense `B` makes `Xtilde` nonlocal, so the same simple diagonal-variable-projection derivation cannot simply be copied with `Q→F`; any claim that `psi_B` is the eliminated residual of a new constitutive model needs a separate derivation with its norm and constraints specified.

Two unrestricted choices expose the issue:

* `B=-G` gives `Gtilde=0`, `F=E` for **every** material, and `dF/dm=0`. The physical current and nonzero physical data gradient remain. Thus minimizing a field-variation or projector-variation score alone chooses an algebraic cancellation with no evidence of physical-search benefit. This choice generally makes `Xtilde=(I-XG)^{-1}X` dense and might be excluded by a declared structural class.
* `B=tI`, `t→+∞`, with positive nonzero material and full-column-rank `j`, gives `F=Q+t j` and `psi_B=O(t^{-2})→0` even when `R` is fixed and nonzero. `I+tX` stays invertible. The *raw* Gram objective is therefore scale-degenerate as an unrestricted outer objective. An acquisition reference scale that depends on `B` can alter this asymptotic statement; such a scale must be declared and held consistent before making any comparative claim.

## 3. Exact fixed-`B` material gradient

For a real material perturbation `h`, set `X_h=diag(h)`. The physical current tangent `T=dz/dm[h]` (where `z=j`) solves

```
L T = X_h Q,       dQ = G T,       dR = S T.
```

Differentiate `F=(I+BX)Q` without dropping the direct term:

```
dF = B X_h Q + (I+BX) G T
   = B (I-XG)T + (I+BX)G T
   = (G+B)T.
```

The second line uses the *same physical tangent equation*. It does not hold if one replaces `T` by a reduced or approximate tangent without a corresponding residual term. If `B=B(m)`, add `(dB)j` to `dF`; the following fixed-`B` gradient cannot be used across such an update unmodified.

Set `K=H^{-1}`, `W=R K` and `V=F K R*R K` (so `V=F W*W`). Matrix differentiation gives

```
d psi_B = 2 Re <W, S T>_F - 2 Re <V, (G+B)T>_F.
```

With

```
P = L^{-*} [S*W - (G+B)*V],
grad_i psi_B = 2 Re sum_source conj(P_{i,source}) Q_{i,source},
```

one retained physical forward factorization and one adjoint application suffice for this algebraic gradient. It is the derivative of the **proposal metric**, including the Gram derivative, not the derivative of the original data loss. The formula assumes a fixed source set, real material coordinates, rank-full `F`, and differentiability inside the feasible box. It must be realified or otherwise rederived for independent real/imaginary contrast parameters.

The full-column-rank source-span projector `P_F=F(F*F)^{-1}F*` has the distinct derivative

```
dP_F = (I-P_F)(dF)Fdagger + Fdagger* (dF)* (I-P_F),
Fdagger=(F*F)^{-1}F*.
```

For fixed `B`, insert `dF=(G+B)T`. This projector can move even though physical current is gauge-invariant. Its derivative vanishes for `B=-G`; a small projector derivative by itself is not evidence of an improved material step.

## 4. Source basis and production normalization

For an invertible mixing `C` of the **same** incident illuminations, transform `E→EC` and `Y→YC`. Then `Q→QC`, `R→RC`, and `F→FC`. Exact algebra gives

```
tr[(RC)((FC)*(FC))^{-1}(RC)*] = psi_B,
P_{FC}=P_F.
```

The raw Gram proposal and source-span projector are basis-invariant on the full-rank domain. The raw Frobenius physical data loss `||RC||_F^2` is generally changed by a nonunitary `C`; it is invariant under unitary source mixing. If the source-noise covariance matters, it must be transported with the data rather than silently treating the mixed columns as independent equal-noise measurements.

R4's production objective is not just raw `psi_0`: its factor is `0.5 * reference_field_mean / ||Y||_F^2`, followed by a material regularizer. Under nonunitary `C`, a freshly recomputed reference-field mean and `||YC||_F^2` generally change by different factors; the normalized objective and its balance against regularization change. Holding the **original acquisition's fixed normalization constant** during a source-basis algebra check preserves the raw invariance statement. Any feedback method must freeze its normalization convention and physical ranking criterion before comparing methods.

## 5. Pole and nonnormal safety

Algebraic invertibility is insufficient for a stable inner solve. Monitor `sigma_min(I+BX)`, `sigma_min(I+XB)`, `sigma_min(I-XG)`, and `sigma_min(F)` or the source Gram condition, along the material trajectory. A norm bound `||BX||_2<1` gives `||(I+BX)^{-1}||_2 <= 1/(1-||BX||_2)`; it is sufficient and conservative. Spectral radius `<1` alone is unsafe for nonnormal matrices: `BX=[[0,20],[0,0]]` has spectral radius zero while `cond(I+BX)` is about 402. A diagonal near pole can make a transformed material coefficient of order `10^6` while the physical current remains unchanged. Such cases must not be accepted by a small transformed objective or a finite-value check alone.

## 6. Exact comparison with A17's current Galerkin Schur identity

The actual A17 closure sources are `Gaussian/A17/CLOSURE_R1/manuscript/CORE_THEORY.tex` (Proposition 1 and Theorem 1) and `CONDITIONAL_EXCHANGE_THEORY.tex` (common-core exchange and step gain). A17 uses a **current-space Galerkin approximation** `R_0=V(V*LV)^{-1}V*` for `L=I-XG_D` and a material tangent injection `B_p`; this `B_p` is unrelated to the feedback-design matrix `B` above. For a candidate current bank `C`, A17 defines

```
Cbar=(I-R_0 L)C,
Y_C=S Cbar,  N_C=C*(B_p-LR_0B_p),  Gamma_C=C*L Cbar,
J_C-J_0=Y_C Gamma_C^{-1} N_C.
```

For old and new banks `A,C` around the same retained core `V`, the exact Galerkin endpoint difference is `J_C-J_A=Y_C Gamma_C^{-1}N_C - Y_A Gamma_A^{-1}N_A`. This can alter both the material normal curvature and the residual-dependent linear term. The manuscript's conditional gain is evaluated against an unchanged full quadratic and explicitly distinguishes a numerical accepted move from a stopping certificate and from material-truth improvement.

The feedback gauge instead preserves the **full** current and data at every admissible `B`, changing only a proposal metric. `R_0L` is a current-space oblique projector; `P_F` above is a Euclidean projector onto a transformed illumination span. Their derivatives, domains, and roles are different. A17's DDA implementation also uses contrast-dependent polarizability and vector local fields; the scalar-TM toy algebra here is no substitute for that model. Neither A17's finite Galerkin identity nor this exact gauge implies improved physical material recovery, a continuous-PDE result, or a global basin guarantee.

The separate Stage 0 tangent-span diagnostic evaluates A17's **finite** current exchange `K_C-K_A=Cbar Gamma_C^{-1}N_C-Abar Gamma_A^{-1}N_A` from the shared core and compares it with two direct Galerkin endpoint solves. It also checks `P(K_A+Delta K)-P(K_A)` against `P(K_C)-P(K_A)` for two material columns. That is a projector onto **current-tangent ranges** and is unrelated to the requested illumination-projector derivative.

The actual A17-to-illumination bridge fixes the **exact** total-field anchor `Q` and one material direction `h`, uses the scalar contrast-source injection `D_h=diag(h)Q`, and forms two Galerkin current tangents `T_A,T_C` around the same retained `V`. Define `P_Q=Q Qdagger` and the linear derivative map

```
mathcal D_Q(T) = (I-P_Q) G T Qdagger
               + [(I-P_Q) G T Qdagger]*.
```

At this common anchor, the exact Schur exchange `Delta T=Cbar Gamma_C^{-1}N_C-Abar Gamma_A^{-1}N_A` equals `T_C-T_A`; therefore `mathcal D_Q(T_C)-mathcal D_Q(T_A)=mathcal D_Q(Delta T)`. The code compares the left side from independent direct Galerkin solves with the right side from Schur factors. In this one instance, the enriched derivative still differs from the **full** derivative by about 75% relative norm, a recorded Galerkin tail. The identity only propagates the exchange through a linear derivative map; it does not reduce that tail or validate the actual vector DDA derivative injection in A17.

For the proposed feedback fit, the fixed-anchor ridge problem is `min_B ||(G+B)T||_F^2+lambda||B||_F^2` on source-whitened physical current tangents `T`. The diagonal normal equation is `b_i (||T_i||²+lambda)+sum_j (GT)_{ij} conj(T_{ij})=0`; each disjoint block obeys `B_k(T_kT_k*+lambda I)+(GT)_k T_k*=0`. Under an invertible source change `C`, whitening changes the tangent columns by a right unitary matrix, so the same exact ridge minimizer is expected when all sources and directions are mixed consistently. This invariance applies to the *fit* with its fixed ridge convention, not to R4's recomputed production normalization or equal-noise physical data loss.

## 7. Stage 0 disposition

The companion JSON records gauge/current equality, ordering, full fixed-`B` directional gradient, a direct Torch CPU check of the production gradient helper against the independent NumPy calculation, source-mixing caveat, two unrestricted-design counterexamples, conditioning witnesses, A17 shared-core Schur/step identities, the distinct tangent-span projector difference, the common-anchor illumination-projector derivative bridge and its full-model tail, and the production diagonal/block ridge equations with source-mixing invariance. These are necessary algebraic checks. A substantive next claim requires a predeclared admissible `B` class and normalization, a physical pilot with original-data objective and matched compute accounting, and an out-of-sample comparison; algebra alone cannot choose a scientific direction.
