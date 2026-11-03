# R4 ongoing derivations and route qualification

## 1. Exact port lifting collapses after full-current elimination

Let A(X)=I-XG invertible, j*=A^-1 X Einc, F=Sj*, and r=Aj-XEinc.
Then SA^-1r=Sj-F. This is the session's exact identity, not a new result.
For fixed X consider L(j,X)=1/2||Sj-y||²+lambda/2||Sj-F||², lambda>0.
Write P as the orthogonal projector onto range(S), z=Sj. Since F is in range(S),
the optimum is z=(Py+lambda F)/(1+lambda), and

min_j L = lambda/[2(1+lambda)]||F-Py||² + 1/2||(I-P)y||².

Thus the profiled material landscape is a positive rescaling of the original
data objective plus a material-independent constant. Under an invertible fixed
whitening operator the same proof uses range(WS). If a prior is appended without
the same rescaling, the effective prior changes; that is not a pure structure gain.
Restricting current capacity, leaving currents inexact, or adding eta||r||² can
change the landscape, but each must be tested as that different mechanism.
Status: PROVED_UNDER_ASSUMPTIONS. This qualifies the proposed pure port-weighting
route and avoids spending experiments on an exact algebraic restatement.

## 2. Receiver-visible scattering feedback and cooperative curvature

At a physical anchor let R=A^-1, C=GR, O=SR, E=Einc+Gj, residual e=F-y.
For D=diag(d), the full finite response is deltaF=O(I-DC)^-1 D E.
The first two coefficients are ODE and ODCD E. The least-squares Hessian is

d^T Hdata d = ||ODE||² + 2 Re <e,ODCD E>.

After division by ||y||² and adding the exact inherited prior Hessian, this is
the original physical Hessian. It is already present in R4-BASIN HVP code.
Here the residual-feedback matrix has entries
T_ij=C_ij sum_rs conj(e_rs) O_ri E_js,
Hfeedback=Re(T+T^T)/||y||². The GN matrix has entries
Re[(O*O)_ij sum_s conj(E_is)E_js]/||y||².
This supplies a direction criterion involving source, receiver, phase and current
data mismatch; a large |C_ij| or interaction energy alone omits that alignment.

The first diagnostic asks whether the already-near-stationary R3 failure has
negative cooperative curvature on its free face. Positive free-face curvature
would reject that local explanation, without certifying the whole box problem.
The dense contraction is an implementation of known derivatives; novelty and
algorithm benefit remain unestablished.

## 3. Illuminated localization of a completed T needs no full T inverse

The DCTMC control shows data-compatible nonlocal T is not a local material.
For any completed T with I+TG invertible let V=(I+TG)^-1 T, J=T Einc,
Q=Einc+GJ. Then VQ=J. Rather than approximating all entries of V by diag(m),
one may minimize the illuminated discrepancy ||(diag(m)-V)Q||_Omega².
Its real box-constrained row LS is m_i=clip(Re sum_s w_s conj(Q_is)J_is /
sum_s w_s |Q_is|²); the same public spatial prior uses the existing material
inner solver. This is the classical CSI material projection, not a new readout.

For a physical T0, receiver S=Ua Sa Va*, incident Einc=Ub Sb Vb*, and a fixed
retained mask M_ij on products sa_i sb_j, the data-completion correction is
DeltaT=Va [M odot (Ua*(Y-F0)Vb)/(sa_i sb_j)] Ub*.
Consequently DeltaJ=DeltaT Einc = Va [M odot Ua*(Y-F0)Vb / sa_i] Vb*.
The incident singular values cancel exactly in this action; their sole remaining
role is the frozen mask. Thus illuminated localization can be executed with the
physical current, a truncated data projection and one Green action, without
forming/inverting T or V. If the mask is source-independent this is the familiar
receiver pseudoinverse current projection. It is a simplification/known-method
connection, not automatically a novel algorithm or an imaging gain.

The decisive comparison keeps the same T-data mask and 30 outer steps but changes
locality projection from operator-Frobenius diagonal extraction to illumination-
weighted shared-real constitutive fitting. Both select only by original Phi and
use the same original corrector. This tests why the nonlocal completion fails.

## 4. Collective current invariants and observational redundancy

For real diagonal eta=X^-1 and reciprocal G, (diag(eta)-G)J=Einc implies
Einc^T J=J^T Einc. With Gamma=(G-G*)/(2i), it also implies
Im_H(Einc*J)=J*Gamma J. These are classical reciprocity and discrete optical-
theorem relations, not new discoveries. Known nonuniform cell volumes require
their mass weights; this study uses the existing equal-volume scalar grid.

If Einc^T=L S, the reciprocity relation is fixed by SJ alone. If Einc*=K S and
Gamma=S* Z S, the optical theorem is also an equation in SJ alone. With exact
data locking these constraints then cannot remove hidden-current ambiguity;
incompatible noisy data instead make the constraints inconsistent. Approximate
factorization can make enforcing them expensive in weakly visible currents.
The diagnostic measures these defects in existing saved CIE states before
adding a new constraint or claiming extra physical information.

## 5. A: local defect transport, exact block absorption, and its one-exchange reduction

This is a finite equal-cell scalar-TM statement. Let X=diag(m)>0, Q=Einc+GJ,
r=J-XQ, and partition the actual Green matrix into disjoint four-cell spatial
blocks Gin and the remaining interactions. The state penalty uses

K_beta=I+X(beta I-Gin), R_beta=beta K_beta^-1 r, beta>0.

Gin=0 gives ordinary CIE. Its zeros are the original Maxwell state equation
whenever K_beta is invertible; this does not make its off-shell objective or
optimizer equivalent to ordinary CIE. The same material prior is retained.
With B=beta I-Gin, differentiation at fixed J,Q gives
dR=-K_beta^-1 diag(dm)(B R+beta Q). For a fixed Hermitian source metric Omega,
set Z=K_beta^-* (R Omega); the exact real material gradient of
rho/2 tr(R Omega R*)+Rprior is

g_i=-rho Re sum_s conj(Z_is)(B R+beta Q)_is + (gprior)_i.

The derivative of K_beta must be included. Freezing that derivative changes the
method. Our implementation includes it. A sufficient pole exclusion is
Im_H Gin positive definite in each block: writing K_beta=X H with
H=X^-1+beta I-Gin, H v=0 would imply v* Im_H(H)v<0, a contradiction.
Semidefinite radiation alone does not give this strict argument. This is a
finite-matrix sufficient condition, not a continuous Maxwell uniqueness proof.

Separate the diagonal self term: D=I+X(beta I-diag G),
L=D^-1 X Goff, Goff=Gin-diag G, u=beta D^-1 r. Then

R_exact=(I-L)^-1 u,       R_one=(I+L)u.

For ||L||<1, the Neumann identity gives
R_exact-R_one=(I-L)^-1 L^2 u, so its norm is at most
||L||^2/(1-||L||) times ||u||. Also I+L is invertible, hence the simplified
penalty has exactly the same state zeros. These statements hold in the norm
used for the bound; an infinity-norm bound is not silently an H(curl) or
Euclidean estimate. On the development public grid the analytic whole-box
infinity-norm bounds are 0.16810 (beta=.2) and 0.03784 (beta=1), obtained from
max_{m in [.02,20]} m/|1+(beta-Gii)m| times the off-diagonal block row sum.
The maximum is checked at the endpoints and the legal critical point
m=-1/Re(beta-Gii); no material truth is used.

The exact Euclidean identity
||R_one||^2=||u||^2+2 Re<u,Lu>+||Lu||^2
identifies what the penalty changes: complex, signed coupling of neighboring
constitutive defects, together with the quadratic term. The diagonal control
removes both terms and did not escape the development failure. It does not
separately establish the benefit of the signed cross term or of the actual
Green phase versus another norm-matched coupling. The one-exchange method retained the
physical fit, so summing infinitely many local loops is not required by that
evidence. It is not a truncated forward solver: final fields and selection use
the unchanged full Maxwell resolvent, including all scattering orders.

For b_i=beta-Gii, d_i=1+b_i m_i, a_i=m_i/d_i,
u_i=beta(J_i-m_i Q_i)/d_i, the simplified derivative is
du_i/dm_i=-beta(Q_i+b_i J_i)/d_i^2, da_i/dm_i=1/d_i^2.
With W=R_one Omega and Z=W+Goff*(conj(a) W),
g_i=rho Re sum_s [conj(Z_is)du_is/dm_i+
conj(W_is)(Goff u)_is/d_i^2]+(gprior)_i.
This costs fixed-size block products and scalar divisions, no block inverse.
At fixed four-cell block size the added work is O(Ns); it does not reduce the
cost of the global Fourier current updates or full physical corrector.

Status: DERIVED, with EMPIRICAL development evidence. Diagonal/self
renormalization, modified Green CIE/NIE and T-matrix changes are known; see
literature/BLOCK_SELF_CIE_PRIOR_ART.md. Novelty, if any, must lie in the specific
off-shell transport and its demonstrated behavior, not the generic resolvent
identity. The current implementation remains a CIE modification, not a
replacement theory of inverse scattering. Neither the bound nor a successful
case proves basin enlargement for all materials.
The two-beta current/material sweeps do not optimize one common auxiliary
objective, so a single-objective descent or no-cycle theorem is not claimed.

## 6. B: finite state regularization yields a resolvent-dependent profile metric

Let e=F-Y, q=A(j-j*), C=SA^-1. The tested lifted objective, without the prior
and common normalization, is

L(q,m)=1/2||e+Cq||^2+lambda/2||Cq||^2+eta/2||q||^2, eta>0.

Its unique q is -C*[eta I+(1+lambda)CC*]^-1 e and its profile is
1/2 <e,M e>, where

M=(eta I+lambda CC*)[eta I+(1+lambda)CC*]^-1.

For lambda=0 this is the familiar wavefield/state-relaxation profile. As
eta tends to zero on the range of C, M tends to lambda/(1+lambda), recovering
the collapse in section 1 for lambda>0. Finite eta changes the material
landscape through C(m); this is a real modification but overlaps known
variable-projection/WRI/CC-CSI structures. Holding M fixed while changing m
would omit a derivative. The implementation uses the exact envelope gradient.

With v=Cq, z=A^-* S* v, Qhat=Einc+G(j*+A^-1q), the material gradient is
-lambda Re sum_s conj(z_is)Ephysical_is
-eta Re sum_s conj(q_is)Qhat_is, plus the original prior, all data terms divided
by ||Y||^2. Lambda 0 and 1 use the same eta frozen from the public host and
identical iteration budgets, then the same original-Maxwell corrector.
The lambda=1 lower final residual accompanied worse material on the disclosed
development object. Thus a changed resolvent metric is established, an imaging
advantage is not. Status: DERIVED / EMPIRICAL, with no originality claim.

## 7. Why nonlocal residual transport can change the material readout

This is a conditional weighted-least-squares consequence, not a new abstract
projection theorem. It sharpens the already known R3/D1 diagonal-weight
cancellation. Use realified coordinates, absorb the fixed source metric into
the fields, and linearize at an exact constitutive state with an interior,
unregularized material update. Write the state defect as b-C dm, where C has
full column rank. Its columns contain each cell's illumination vector and
have disjoint spatial support. Let M=P^T P>0 be the frozen residual metric.
At a zero defect, derivatives of the material-dependent P do not enter this
first-order readout. The two material maps are

D_M=(C^T M C)^-1 C^T M,       D_0=C^dagger,
Pi=C(C^T C)^-1 C^T.

Directly using D_M C=I gives the exact difference

D_M-D_0=(C^T M C)^-1 C^T M(I-Pi).

Therefore D_M=D_0 for every right-hand defect b if and only if
C^T M(I-Pi)=0. By symmetry of M this is equivalent to invariance of range(C)
under M, or [M,Pi]=0. A positive scalar weight shared across all sources of a
cell has this invariance and cancels from the unregularized readout. Spatial
off-diagonal transport generally does not: it changes how the part of b that
no local real material perturbation can explain feeds into the material step.
It need not improve that step; the sign and phase of the transported defect
can matter. A difference for the actual reachable b additionally requires
C^T M(I-Pi)b!=0; non-invariance only guarantees that some right-hand defect
has a different readout. This predicts a possible difference between +L and -L despite their
identical amplitude and quadratic-term contributions.

The condition also supplies a decisive limitation: if C spans the entire
residual space, then Pi=I and every invertible weighting cancels. For instance
a one-real-source, one-real-equation-per-cell toy model with nonzero field
cannot demonstrate this first-order mechanism. Even one complex illumination
with a shared real material can leave a constitutive residual complement;
multiple complex illuminations generally enlarge it.
The full Maxwell experiment retains that coupling; a scalar one-equation
analogue would discard it.

An active box, nonzero prior, nonzero state defect, incomplete inner solve, or
different beta in the two sweeps invalidates a direct use of this formula as
the derivative of the complete implemented iteration. These conditions are
present in parts of the experiments, so this is a structural explanation and
failure boundary, not a convergence or causality proof for the observed gain.
Status: DERIVED, under the stated assumptions. No machine-precision identity rerun was
needed, and no additional observation or material prior is supplied by it.

## 8. A structural refinement after failed fresh transfer: average the partition

The six-object confirmation did not show a success-count advantage. One
remaining, concrete implementation hypothesis is dependence on the arbitrary
origin of disjoint2x2 blocks. This motivates one new development variant, not a
post-hoc modification of the frozen confirmation method.

Let Goff^(q), q=1,...,4, denote the four parity translations of the2x2 partition,
with boundary singleton/partial blocks and no periodic wrap. Set
Gbar=(1/4) sum_q Goff^(q). Its horizontal/vertical neighbor entries are Gij/2;
diagonal neighbors are Gij/4; all other entries vanish. These weights are the
exact partition-inclusion probabilities, not tuned parameters. With the same
material-dependent a and u as before,

Rbar=[I+diag(a)Gbar]u = (1/4) sum_q [I+diag(a)Goff^(q)]u.

This averages the residual transport, not the four squared objectives. The
diagonal self transform, material prior and full physical forward stay fixed.
The convex norm bound gives ||diag(a)Gbar|| <= (1/4)sum_q||diag(a)Goff^(q)||,
so any common q<1 suffices for original state-zero equivalence. The actual
finite-grid infinity bound is computed analytically from row sums. The same
gradient derivation applies to Gbar and its exact adjoint, implemented as a
sparse8-neighbor action with O(Ns) work.

This removes the interior parity-origin choice; it does not remove the finite
spatial neighborhood or prove a better search direction. The initial minimal
test uses the disclosed R2 development object and the lexicographically first
public-data-ineligible confirmation ID, obs_42cb233e4f224bb3. It never uses a
GT direction. The latter is now POST_CONFIRMATION_DEVELOPMENT, excluded from
new confirmation statistics. If this fails, it rejects this particular
partition-origin repair, not every spatial transport or nonlocal algorithm.

Observed outcome: both prescribed development cases failed. On the old positive,
the physical relative residual became .107660 and material NRMSE .389719. On
42cb the NRMSE decreased from .761368 to .641395 while the residual increased
to .143491. This is a mixed material/data change, not a recovered image or a
GT-selectable online improvement. See `results/DEVELOPMENT_POST_CONFIRMATION_ADDENDA_V1.md`.
