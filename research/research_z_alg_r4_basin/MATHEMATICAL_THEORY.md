# Finite Maxwell steps, correction neighborhoods, and what they do not imply

This is the parent researcher's derivation. All materials are independent real scalar contrasts, geometry and complex operators are fixed, and the statements concern a finite discrete VIE. They do not automatically extend to nonlinear decoders, geometry parameters, tensor media, or the continuous vector Maxwell PDE.

## T1. Exact finite-step response and remainder — standard resolvent algebra

Assume both $M=I-D_mG$ and $M_h=I-D_{m+h}G$ are invertible. Let $J=M^{-1}D_mE_0$, $Q=E_0+GJ$, and $F=SJ$. Subtracting the two current equations gives

$$M_h[J(m+h)-J(m)]=D_h(E_0+GJ)=D_hQ.$$

Consequently

$$F(m+h)-F(m)=SM_h^{-1}D_hQ,\qquad DF(m)[h]=SM^{-1}D_hQ.$$

Using $M_h^{-1}-M^{-1}=M_h^{-1}D_hGM^{-1}$ gives the exact higher-order term

$$\boxed{E_m(h)=SM_h^{-1}D_hGM^{-1}D_hQ.}$$

The order of the factors matters. A small Euclidean material step alone does not control this term uniformly near a singular resolvent. A small observed tangent can also arise from receiver cancellation, so it does not bound the internal perturbation.

## T2. A certified sufficient step condition and its practical cost

Set $B_h=D_hGM^{-1}$ and $\beta=\|B_h\|_2$. Then $M_h=(I-B_h)M$. If $\beta<1$, $M_h$ is invertible and

$$\|E_m(h)\|_F\leq\|SM^{-1}\|_2\frac{\beta}{1-\beta}\|D_hQ\|_F=:e(h).$$

Writing $r=F(m)-Y$ and $d=DF(m)[h]$, the computable inequality

$$\|r+d\|_F+e(h)<\|r\|_F$$

is sufficient for strict original-data descent. A regularizer requires its actual change to be included separately.

An explicit upper bound is

$$\beta\leq\beta_F=\|D_hGM^{-1}\|_F=\left(\sum_i h_i^2\|(GM^{-1})_{i,:}\|_2^2\right)^{1/2}.$$

For a fixed current state this defines a physics-dependent ellipsoid in material steps. Constructing $GM^{-1}$ costs a dense inverse action, namely $N$ complex RHS, even if an LU already exists. It is not a free per-iteration diagnostic. A power estimate or an illuminated-direction secant is not a certified upper bound. The ratio $\|E_m(h)\|/\|d\|$ has no uniform output-relative bound from $\beta$ alone, because $S$ can cancel $d$.

Algorithmic consequence: the sufficient bound can qualify a correction step, but its conservatism and cost must be measured. It must not force arbitrarily small steps and then be called a robust global solver.

## T3. Exact rational response along a complete feasible material ray

For a fixed real direction $p$, define $K=M^{-1}D_pG$ and $V=M^{-1}D_pQ$. For any real $t$ with $m+tp$ inside the material box and $I-tK$ invertible,

$$\boxed{F(m+tp)=F(m)+tS(I-tK)^{-1}V.}$$

An eigenvalue or complex-Schur representation can reuse this directional transfer function over many values of $t$. For a small active set, the exact R2 low-rank kernel is preferable. For a dense direction, constructing $K$ requires $N$ RHS and its dense factorization is a real cost; this is not a low-dimensional material model.

Separated one-dimensional loss minima can be proposed even when an ordinary backtracking line search would stop at the first descent region. The formula is standard resolvent/R10 amplitude algebra. Whether correcting several separated minima improves imaging is an experimental algorithm claim, not a consequence of the identity.

## T4. Local identifiability does not imply accessibility — rational VIE counterexample

Consider the one-cell algebraic VIE with $E_0=S=1$, $G=1+0.1i$, legal material interval $[0.5,3]$, and true material $a=2$. Its response is $f(x)=x/(1-Gx)$. No real material is a pole, and $f'(a)=1/(1-Ga)^2\ne0$, so the two-real-component map is locally identifiable at the true material. Indeed the Möbius map is globally injective away from its pole.

The squared data error is a positive constant times

$$\ell(x)=\frac{(x-2)^2}{1-2x+1.01x^2}.$$

Its derivative has the sign of

$$ (x-2)(-1+1.02x). $$

For $0.5\leq x<1/1.02$, the derivative is positive. Thus the lower box endpoint is a strict one-sided constrained local minimum with nonzero residual, while $x=2$ is the unique zero-residual solution. Starting at $x=0.5$, a projected descent or bounded local least-squares method can stop at the wrong endpoint. A global rational scan reaches $x=2$.

This proves the logical separation inside the same rational current/resolvent algebra. This chosen one-cell $G$ is **not asserted to be the self-interaction of the disk-pulse Hankel geometry used in the experiments**. The actual 576-cell evidence must be assessed independently; this counterexample is not a new Maxwell PDE nonconvergence theorem.

## T5. What transverse contraction certifies — conditional local statement

Let $f(m)=\operatorname{pack}_{\mathbb R}F(m)$ and consider an interior, constant-rank, unregularized update

$$T_y(m)=m-J_{\mathbb R}(m)^\dagger[f(m)-y].$$

At an exact solution $m_*$, differentiation gives

$$DT_y(m_*)=I-J_{\mathbb R}(m_*)^\dagger J_{\mathbb R}(m_*).$$

It is zero on the observable normal space and identity on the Jacobian nullspace. If the Jacobian has full column rank, it is zero in every material direction. Under the usual local smoothness and neighborhood assumptions this explains local correction, not entrance into that neighborhood from an arbitrary point. For identity-metric damping $\lambda>0$, the normal eigenvalues become $\lambda/(\sigma_i^2+\lambda)$; null directions retain eigenvalue one.

Box projection is generally nonsmooth where the active set changes. A smooth fixed-active-set analysis must use the free variables. A regularized update is not generally fixed at the exact material solution unless the regularizer gradient also vanishes there. These qualifications prevent a local contraction statement from becoming a GT-recovery guarantee.

## T6. Exact residual Hessian action and a measurable alternative to the GN model

Let $\ell(m)=\|F(m)-Y\|_F^2/(2c)$, where $c=\|Y\|_F^2$ is fixed, and retain a matching physical LU. Solve

$$M^*\Lambda=S^*R,\quad R=F-Y.$$

The real material gradient is $g_i=\Re\sum_s\overline{\Lambda_{is}}Q_{is}/c$. For a real direction $h$, differentiating the two physical equations gives

$$U=M^{-1}D_hQ,\quad \dot Q=GU,$$

$$\dot\Lambda=M^{-*}(G^*D_h\Lambda+S^*SU),$$

$$\boxed{[\nabla^2\ell(m)h]_i=\Re\sum_s(\overline{\dot\Lambda_{is}}Q_{is}+\overline{\Lambda_{is}}\dot Q_{is})/c.}$$

Each action needs $s$ tangent and $s$ adjoint RHS using the same LU; the base adjoint is another $s$ RHS once. This is a standard second-order adjoint differentiation, derived here to expose the exact residual-curvature term that $J_{\mathbb R}^TJ_{\mathbb R}$ omits. No finite differences are required. The common TV/smooth regularizer adds

$$D^T\operatorname{diag}\!\left[\frac{2\lambda_s+\lambda_{TV}\epsilon^2/(d_i^2+\epsilon^2)^{3/2}}{n_{edge}}\right]Dh,\quad d=Dm.$$

On a fixed free-variable face of the box, a negative Rayleigh quotient demonstrates a local descent mechanism invisible to the positive-semidefinite GN Hessian. A numerical smallest-eigenvalue estimate is not a positive-definiteness certificate, and failure to find a negative direction does not prove a local minimum. Box-active directions need feasible-cone treatment; a free-face test only diagnoses that face. This calculation is being used as a mechanism experiment rather than a derivative precision test.

## T7. A genuinely fixed uniform-material transfer curve, and its limitation

For the public family $m=c\mathbf1$, the diagonal material is $cI$, so it commutes with $G$. If $G=UTU^*$ is a complex Schur decomposition, then

$$\boxed{F(c\mathbf1)=c(SU)(I-cT)^{-1}(U^*E_0).}$$

This is a special, exact fixed-operator family. The Schur decomposition, the port factors, and responses on a fixed amplitude grid are independent of the measured target $Y$. They can be reused across targets. Arbitrary spatially varying diagonal material does not commute with $G$, so this property cannot be transferred to all 576 independent pixel updates.

For a fixed response grid $F(c_j\mathbf1)$, each new target only changes the cross term in

$$\|F(c_j\mathbf1)-Y\|_F^2=\|F(c_j\mathbf1)\|_F^2-2\Re\operatorname{tr}(Y^*F(c_j\mathbf1))+\|Y\|_F^2.$$

This offers a reusable, data-conditioned source of public initial materials. A uniform seed is not a digit template, a material decoder, or a solution of the observed image. Several detected minima can correspond to different subsequent correction behavior, but no theorem makes a raw curve minimum a good GT candidate.

The executed V2/V3 implementation uses the equivalent paid reference-$m=1$ ray kernel and charges a cold full-rank Schur setup **for each target**. The direct $G$-Schur cross-target cache above has not been implemented or timed in this sprint; it must not be counted as an achieved amortized speedup. Its component identities are standard matrix-function algebra.

## T8. Two-point observability is the condition a local Jacobian cannot supply

For two legal nonsingular materials $m,n$, T1 gives the exact finite secant

$$F(n)-F(m)=SM(n)^{-1}D_{n-m}Q(m).$$

Define the real matrix $B_{n,m}$ by applying the operator $h\mapsto SM(n)^{-1}D_hQ(m)$ to the independent real material coordinates and packing real/imaginary data. This is an actual field-product matrix: its $(r,s),i$ complex entry is $[SM(n)^{-1}]_{ri}Q(m)_{is}$. It uses the receiving field at one endpoint and the illuminating field at the other, as in R10's bilateral response identity.

If a set of materials has the **uniform two-point** condition $\sigma_{\min}(B_{n,m})\ge\mu>0$ for every pair in that set, then

$$\|F(n)-F(m)\|_F\ge\mu\|n-m\|_2.$$

It follows immediately that the measurement map is injective and uniformly inverse-stable on that set. A lower bound only on the diagonal $B_{m,m}=J_{\mathbb R}(m)$ does not provide this two-point condition; field changes at separated endpoints matter. Checking finitely many pairs cannot certify the full box condition. The scalar counterexample in T4 already shows that even injectivity and nonzero local derivative do not ensure accessibility by local descent.

This is a sufficient conditional statement obtained from a standard exact secant identity, not an experimentally established global inverse bound for the current geometry. R4's singular tails and unsuccessful reconstructions do not establish the premise. An observed ratio $\|n-m\|/\|F(n)-F(m)\|$ only gives a lower bound on the inverse-stability constant needed for that observed pair.

## Research hypotheses, not theorems

1. Poor finite-step agreement and observable-frame rotation can identify when a local correction model has become uninformative. Thresholds are empirical unless separately certified.
2. Correcting several separated rational-direction minima before original-data ranking can retain entries that an early raw-score filter discards.
3. The resulting candidate set may improve finite-budget accessibility; it need not contain GT, represent all branches, or carry probability weights.

These hypotheses will be retained or rejected from actual experiments. No new measurements are introduced.
