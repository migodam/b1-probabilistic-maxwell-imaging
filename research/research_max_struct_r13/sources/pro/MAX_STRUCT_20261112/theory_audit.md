# MAX-STRUCT independent theory audit: finite constitutive branch events and representation-induced material penalties

Date: 2026-11-12 (Singapore). Repository snapshot reviewed: `ec3a37da6f3ee90beb61e753067d4f639f4118f5`.

## Executive judgment

The most useful additional theoretical result developed in this audit is a **finite-distance branch-switch law for the shared-real-material constitutive projection**, together with an exact account of how DR normal memory traverses that distance at an AP stagnation pair. Its key quantity is not a Green singular value or a local tangent angle: it is the competition between the energy of a shared constitutive field and the incompatible normal component accumulated by reflection.

This does **not** prove global material recovery. It predicts when a particular nonphysical, mutually closest field/current pair can no longer retain its old material shadow. It therefore addresses a specific search obstruction, while leaving data-feasible wrong materials, statistical ambiguity, and model mismatch separate.

A second useful result is an exact finite-material decomposition of a restricted auxiliary-current profile objective into its full-current counterpart plus a nonnegative representation penalty. That penalty has a computable material derivative and, at a compatible root, a positive-semidefinite Hessian contribution. This is sharper than reporting raw current projection error. Its underlying least-squares algebra is standard; a paper must earn novelty through Maxwell-specific predictive consequences and an efficient intervention.

The repository already contains the exact resolvent derivative, source/receiver product fields, shared-material least squares, dark-current material lifts, finite secant bounds, DPCS/Fourier projection-loss decomposition, physical-complement warm starts, and much of the relevant local geometry. Repeating them would not be a new contribution. The normal-drift mechanism itself is also established in difference-map/DR literature; the bounded-material event formula, its perturbation limits, and useful Maxwell predictions are the potential new content.

---

## 1. Model, realification, and scope

Use the repository convention

\[
J=D_m Q,\qquad Q=E+GJ,\qquad Y\simeq SJ,
\]

with real isotropic material \(m_i\in[\ell_i,u_i]\), \(0\le\ell_i<u_i<\infty\), and a fixed positive field/current scale \(s\). R12 uses \(X=J/s\) and metric

\[
\|X\|_F^2+\|Q\|_F^2.
\]

At cell \(i\), collect all illuminations into a complex row; in vector Maxwell also collect the three field components. Realify that row so that it is a real vector of length \(d=2L\) for scalar TM or \(d=6L\) for isotropic vector Maxwell. The constitutive constraint is

\[
\mathcal C_i=\{Z_i=u(\theta_i)v_i^T:
\theta_i\in[\theta_{i,-},\theta_{i,+}],\ v_i\in\mathbb R^d\},
\]

where the two rows of \(Z_i\) are \((X_i,Q_i)\),

\[
u(\theta)=(\sin\theta,\cos\theta)^T,
\quad\theta_i=\arctan(m_i/s),
\quad\theta_{i,\pm}=\arctan((\ell_i,u_i)/s).
\]

Thus each cell is a union of rank-one linear subspaces with a bounded material orientation. The interval width is less than \(\pi/2\) under the stated nonnegative finite material bounds.

The propagation/data set

\[
\mathcal A_\delta=
\{(X,Q):Q=E+sGX,\ \|W(sSX-Y)\|\le\delta\}
\]

is closed and convex when nonempty. A global noise ball is permitted; the analysis does not require an affine data equality. The material set is a product \(\mathcal C=\prod_i\mathcal C_i\), generally nonconvex.

These real-angle formulas do not automatically extend to unknown complex lossy contrast or anisotropic tensor contrast. The isotropic vector case does extend by enlarging the per-cell illumination/component vector. Lossy and tensor cases require their own local constraint geometry.

---

## 2. Exact normal-memory drift at a nonfeasible AP pair

### Proposition 1: drift along a proximal normal ray

Let \(\mathcal A\) be a nonempty closed convex set. Let \(c,a\) satisfy

\[
c\in P_{\mathcal C}a,\qquad a=P_{\mathcal A}c,
\qquad d=a-c\ne0.
\]

Define standard DR

\[
Tz=z+P_{\mathcal A}(2P_{\mathcal C}z-z)-P_{\mathcal C}z.
\]

For every \(t\ge0\) for which \(P_{\mathcal C}(c+td)=c\) with the chosen projection branch,

\[
\boxed{T(c+td)=c+(t+1)d.}
\]

**Proof.** Since \(a=P_{\mathcal A}c\), the convex projection characterization gives \(c-a=-d\in N_{\mathcal A}(a)\). Hence

\[
P_{\mathcal A}(c-td)
=P_{\mathcal A}(a-(t+1)d)=a,
\]

because positive multiples of a convex normal remain in the same normal cone. Substitution in DR proves the formula. \(\square\)

AP remains at \(a\mapsto c\mapsto a\). DR retains the mismatch in its state instead of resetting to one of those two shadows. This general normal drift is prior art, including Elser, Rankenburg & Thibault (2007) and Elser (2017). Its use here is a foundation for the material-specific finite event law below.

### Important limitations

1. The proposition does not say that every stalled physical optimizer is an AP best pair.
2. It does not say that every real DR iterate lies on this ray.
3. A nonzero gap may reflect model inconsistency, rather than an accessible good material.
4. At a wrong but feasible physical material, \(a=c\), so \(d=0\): there is no normal drift and no automatic escape.
5. Near a noise-ball interior point the data inequality is inactive. Its tangent geometry cannot stand in for statistical material identifiability.

---

## 3. Exact first material-branch event

Assume \(c=P_{\mathcal C}a\). For one cell write

\[
c_i=u_i v_i^T,\qquad d_i=u_{i,\perp}w_i^T,
\qquad u_{i,\perp}=(\cos\theta_i,-\sin\theta_i)^T.
\]

The second expression follows from first-order stationarity with respect to the unconstrained amplitude vector \(v_i\): \(u_i^Td_i=0\). This amplitude-normal condition is valid at material bounds as well as in the interior.

The computable quantities are

\[
a_i=\|v_i\|^2,\qquad b_i=\|w_i\|^2,
\qquad \gamma_i=\langle v_i,w_i\rangle.
\]

Using complex rows without explicit realification,

\[
v_i=\sin\theta_i\,X_{c,i}+\cos\theta_i\,Q_{c,i},
\]

\[
w_i=\cos\theta_i\,d_{X,i}-\sin\theta_i\,d_{Q,i},
\qquad \gamma_i=\operatorname{Re}\langle v_i,w_i\rangle.
\]

A stationary interior material has \(\gamma_i=0\). At the lower bound \(\gamma_i\le0\); at the upper bound \(\gamma_i\ge0\).

For a competing material orientation \(\theta_i+\varphi\), the squared projected energy at \(Z_i(t)=c_i+td_i\) is

\[
\mathcal E_i(\varphi,t)
=a_i\cos^2\varphi
+2t\gamma_i\sin\varphi\cos\varphi
+t^2b_i\sin^2\varphi.
\]

The old shadow has energy \(a_i\). Therefore it remains globally closest exactly when

\[
(t^2b_i-a_i)\sin^2\varphi
+2t\gamma_i\sin\varphi\cos\varphi\le0
\]

for every allowed \(\varphi\).

### Proposition 2a: interior material

If \(\theta_i\) is interior, \(\gamma_i=0\), and \(b_i>0\), the first event is

\[
\boxed{t_{*,i}=\sqrt{a_i/b_i}.}
\]

For \(0\le t<t_{*,i}\), the material shadow and its amplitude remain exactly \(c_i\). At \(t=t_{*,i}\), the real \(2\times2\) Gram becomes isotropic and all admissible material orientations tie. For \(t>t_{*,i}\), the maximizing orientation is whichever material bound has the larger angular distance from \(\theta_i\); if these distances tie there are two maximizing bounds.

**Proof.** The energy difference becomes \((t^2b_i-a_i)\sin^2\varphi\). Its sign gives the first assertion. Once positive, its maximum is attained at the largest \(|\varphi|\), because the allowed angular separation is less than \(\pi/2\). \(\square\)

This is an exact singular-value crossing of a rank-one constitutive cone. It is not a novel SVD theorem. The interpretation is that incompatible, cross-illumination normal memory becomes as large as the shared constitutive component.

### Proposition 2b: bound-active material

Suppose \(\theta_i=\theta_{i,-}\) or \(\theta_{i,+}\), and let

\[
\varphi_{i,f}=\theta_{i,\mathrm{other\ bound}}-\theta_i.
\]

If \(b_i>0\), the first competing branch is the opposite bound, and

\[
\boxed{
 t_{*,i}
 =\frac{-\gamma_i\cot\varphi_{i,f}
 +\sqrt{\gamma_i^2\cot^2\varphi_{i,f}+a_i b_i}}{b_i}.
}
\]

Equivalently this is the positive root of

\[
b_i t^2+2\gamma_i\cot\varphi_{i,f}\,t-a_i=0.
\]

**Proof.** Divide the energy difference by \(\sin^2\varphi>0\). The condition is

\[
b_i t^2-a_i+2t\gamma_i\cot\varphi\le0.
\]

At the lower bound, \(\varphi>0\), \(\gamma_i\le0\), and the left side increases with \(\varphi\), hence is maximal at the opposite bound. At the upper bound the analogous statement holds with both angular and normal signs reversed. The positive quadratic root follows. \(\square\)

Unlike the interior case, when \(\gamma_i\ne0\) this is generally a **global tie between two boundary branches**, not a loss of local stability of the original endpoint. The local material derivative still points outward at the original bound. Calling every event a local Hessian instability would be incorrect.

For exact projection data \(c=P_{\mathcal C}a\), all \(t_{*,i}\ge1\). A materially smaller value indicates a sign, scaling, projection, or numerical issue. It is not a new recovery gate.

If \(b_i=0\), assign \(t_{*,i}=\infty\). If \(d_i\ne0\), then \(a_i>0\): with two distinct allowed slopes the union of constitutive lines spans the two-row space, so a nonzero input cannot project to zero on every allowed line.

### Corollary: finite escape from a constant nonphysical shadow

At an exact nonfeasible AP best pair with the above nontrivial positive real material intervals,

\[
t_* = \min_{i:b_i>0}t_{*,i}<\infty.
\]

Thus DR cannot preserve the old constitutive shadow indefinitely. It must change a material branch by the first integer ray parameter strictly exceeding \(t_*\), with tie behavior at an integer \(t_*\) determined by the implementation.

This is a finite-distance statement, not only an asymptotic local convergence rate. It still provides no guarantee about the quality of the branch reached after the switch.

---

## 4. Near-AP pairs: a computable perturbation bound

Real saved states are usually approximate, not exact, AP best pairs. Let

\[
a\in\mathcal A,\quad c=P_{\mathcal C}a,\quad d=a-c,
\quad e=P_{\mathcal A}c-a,\quad\varepsilon=\|e\|.
\]

The constitutive ray event formula remains exact because it needs only \(c=P_{\mathcal C}a\). What becomes approximate is the assertion that the real DR trajectory follows that ray.

### Proposition 3: one-step ray defect

For every \(t\ge0\) before the constitutive ray event,

\[
\boxed{
\|T(c+td)-[c+(t+1)d]\|\le(t+1)\varepsilon.
}
\]

**Proof.** Put \(a'=a+e=P_{\mathcal A}c\). Then \(-(d+e)\in N_{\mathcal A}(a')\), so

\[
P_{\mathcal A}\big(a'-(t+1)(d+e)\big)=a'.
\]

The argument of this projection is \(c-td-te\). Nonexpansiveness of the convex projection yields

\[
\|P_{\mathcal A}(c-td)-a'\|\le t\|e\|.
\]

Add \(\|a'-a\|\) and use the DR formula. \(\square\)

If the true trajectory remains in a region where \(T\) is Lipschitz with constant \(L_k\), and \(E_k\) denotes distance from the reference ray at \(t_k=t_0+k\), then

\[
E_{k+1}\le L_k E_k+(t_k+1)\varepsilon.
\]

For a product of strict fixed bound branches, \(P_{\mathcal C}\) is locally an orthogonal linear projector; with convex \(\mathcal A\), DR is nonexpansive on that common branch. As long as actual and reference states retain these branches,

\[
E_K\le E_0+\frac{\varepsilon K(2t_0+K+1)}2.
\]

This is a conditional bound: it does not authorize ignoring an intervening branch change.

For an interior cell on its normal ray, the differential of the nearest rank-one projection has norm

\[
\|DP_{\mathcal C_i}\|
=\frac1{1-t\sqrt{b_i/a_i}}
\]

while \(t<t_{*,i}\) and the material angle remains interior. The expression follows by reducing the state to the two singular directions \(\sigma=\sqrt{a_i}\), \(\tau=t\sqrt{b_i}\); the two coupled tangent eigenvalues are \(\sigma/(\sigma-\tau)\) and \(\sigma/(\sigma+\tau)\). It diverges near the event. Hence an arbitrarily small AP residual need not produce an accurate long-range event-time forecast if a trajectory approaches a projection degeneracy.

**Practical implication:** report AP defect, gap, event margin, active material status, and normal-ray deviation together. A threshold or Gram eigengap alone is not a success predictor.

---

## 5. Algorithm candidate: event-driven reflection, with two distinct semantics

### 5.1 Exact fast-forwarding at an exact AP pair

If the AP pair is exact and the old projection branch is unchanged, all intermediate DR steps add the same \(d\). For an integer ray offset \(n_0<t_*\), starting at \(z=c+n_0d\), fast-forward directly to

\[
z_{\mathrm{event}}=c+\lceil t_*\rceil d
\]

through all preceding unambiguous steps. This skips exactly `ceil(t*) - n0` DR updates and ends at the first integer state that may have changed projection. If `t*` is an integer, this target is the tie state, not a state guaranteed to have switched; call the actual constitutive projector there and preserve its tie rule. For a general noninteger initial ray parameter `t0`, the correct target is `t0 + ceil(t* - t0)`, not necessarily `ceil(t*)`. A branch change is guaranteed only strictly beyond the threshold. This reproduces the future state of ordinary DR and skips redundant projection pairs. It does not change which branch DR ultimately chooses. The step indexing counts the initial projection at `z0` as step 0; with `t0=1` and noninteger `t*`, the first possible changed projection is at step `ceil(t*) - 1`.

Computing every threshold costs \(O(NL)\), or \(O(3NL)\) for vector fields, and one pass over existing state arrays. It requires no material-dependent Maxwell LU, no Jacobian, and no additional Green rank. At an exact pair this can be a provable efficiency improvement, not a recovery improvement.

### 5.2 A finite jump at an approximate AP pair

If \(e\ne0\), setting

\[
z_{\mathrm{new}}=c+(1+\eta)t_*d
\]

is a **new finite-memory proposal**. It is not exact fast-forwarding and may change the basin. The three-case Maxwell probe used \(\eta=0.001\); that is an implementation choice, not a theorem.

A minimal implementation can use an already available AP-like near-pair, compute the event, make one proposal just across it, and continue ordinary DR. It should retain the prior physical candidate so a new branch never erases the only available incumbent. All physical audits, initialization and projection setup must count toward actual cost.

Do not automatically apply a per-cell threshold independently to each cell: a spatially varying scaling of \(d\) generally destroys the global propagation/data normal-ray identity. The global first event uses one scalar \(t_*\).

Do not select a branch by requiring immediate physical-objective decrease. The R11 endpoint paths already show that such a condition can forbid the very finite excursion required to leave a bad physical basin. Candidate selection after search must remain separate from the event geometry.

### What must be tested next

- Does the first-event law predict large actual material jumps, beyond direct projection along a prescribed ray?
- Does an event proposal outperform equal-cost ordinary DR, reset DR, AP, and a magnitude-matched random/non-geometric normal proposal?
- Is any benefit explained by the chosen shared-material direction, rather than merely a large perturbation?
- Can the near-pair be identified cheaply during ordinary operation rather than by paying 3000 preparatory AP steps?
- Does the event mechanism survive vector Maxwell, changes of scale, and unseen geometries?
- For already data-feasible wrong materials, the expected result is **no escape mechanism from zero gap**. This negative control is necessary.

---

## 6. Completed evidence and what it does not establish

### 6.1 Analytic constitutive mechanism model

`normal_ray_demo.py` constructs one real cell with two shared channels, a positive bounded material slope, and an affine observation set. AP is exactly trapped at a nonfeasible best pair of gap 1. Two cases were run:

| Constitutive status | Predicted threshold | First evaluated ray parameter with switch | AP final gap | DR final feasibility gap |
|---|---:|---:|---:|---:|
| Interior material | 3.2 | 4 | 1 | below \(5\times10^{-16}\) |
| Lower-bound material | 3.4802491312 | 4 | 1 | 0 in floating point |

The interior shadow first jumps to the lower bound; the active lower-bound shadow first jumps to the upper bound. Both ordinary trajectories later reach an intersection. Skipping the three unchanged-shadow projection pairs reproduces the same trajectory up to floating-point error.

This is a mechanism counterexample to the claim that identical projectors imply identical AP/DR behavior. It is not a Maxwell imaging validation, and the affine constraint admits multiple materials, so feasibility does not identify a ground truth.

### 6.2 Real Maxwell states, run by the parent research thread

Three already exposed repository instances were rerun to 3000 AP steps and then analyzed. The stored outputs are under `work/normal_memory/`. No ground truth was used for the construction or intervention.

| Instance | AP pair defect / gap | Predicted first ray event | First event material |
|---|---:|---:|---|
| Old shell `A__r10_twod_0004` | about 0.000948 | 5.7317 | lower 0.02 to upper 20 |
| R12 negative case `B__r12_confirm24_0010` | 0.00090698 | 12.495138 | lower 0.02 to upper 20 |
| NN case `C__r8_test_0001` | 0.00598649 | 14.079063 | lower 0.02 to upper 20 |

Probing the prescribed normal ray immediately below and above the formula confirmed the predicted first material branch event in each instance. This checks the relevance of the exact local constituent structure in actual saved Maxwell geometry; it is not itself a new recovery result.

Actual ordinary DR from the approximate AP pair exhibits small material movement from its first step, as predicted by the near-pair caveat. For the shell a large jump occurs near ray parameter 6, consistent with the 5.73 forecast. The formula must not be advertised as an exact DR iteration forecast on all three approximate states.

The original 500-step probe is an intermediate observation, not the final outcome. It was followed by **3000-step continuations** and an integer-event control. The final selected candidates (minimum original Phi over the original incumbent and scheduled physical audits, no corrector, no truth used for selection) are:

| Instance | Ordinary DR: residual / delta | Fractional event: residual / delta | Integer event: residual / delta |
|---|---:|---:|---:|
| Old shell `A__r10_twod_0004` | 0.957701 | 0.973092 | 0.955139 |
| Negative case `B__r12_confirm24_0010` | 30.079656 | 30.079656 | 30.079656 |
| NN case `C__r8_test_0001` | 24.095912 | 0.929641 | 0.929749 |

The shell therefore reaches the data ball under **all three** continuations; the earlier statement that shell failed under both 500-step prefixes must not be presented as the final conclusion. B0010 fails under all three, and all three preserve the same original incumbent. NN0001 reaches the data ball with both fractional and integer event proposals, while ordinary DR remains outside under the declared candidate-selection rule.

The full ordinary traces also locate the first full lower-to-upper box jump at zero-based steps **5, 12 and 14**, respectively, including the predicted cell. These match `ceil(t*) - 1` for an initial state `z0 = c + d`. Small interior material movement still starts at step 1, so this is agreement for the first full box-flip event, not exact constancy of the complete approximate-pair trajectory.

These are exposed-state exploratory continuations, not independent confirmation. The AP3000 preparation, projector setup, physical audits and all branch runs belong in total research cost. Branch incremental CPU excludes shared preparation and cannot by itself establish a complete solver speedup. Material accuracy was subsequently scored offline, after the branch outputs were fixed; see the table below. The integer control reduces concern that the one positive case depends only on the chosen 0.1% overshoot, but it remains one positive instance; it does not establish robustness to general perturbations or unseen objects.

The completed offline amplitude-NRMSE audit is stored in `work/normal_memory/combined_offline_scores.json`:

| Instance | Original input | AP3000 material | Ordinary DR | Fractional event | Integer event |
|---|---:|---:|---:|---:|---:|
| A shell | 0.515137 | 0.521084 | 0.187079 | 0.177757 | 0.186115 |
| B0010 | 0.365795 | 0.368014 | 0.365795 | 0.365795 | 0.365795 |
| C NN0001 | 0.273914 | 0.334856 | 0.334856 | 0.133181 | 0.128761 |

A/B use the inherited amplitude normalization 17.9 and C uses 12, so the table is for within-instance comparisons. B's identical outputs are the retained original incumbent, not newly recovered branches. C confirms actual material improvement over the AP-damaged state in this exposed case. It does not show superiority to the already published DR from the original NN initialization, whose reported NRMSE is 0.077988. The three observed cases cannot establish an independent success rate or eliminate dependence on exposed-case choices.

A crucial theoretical limit follows: **exact fast-forwarding preserves the subsequent DR trajectory, so it cannot improve recovery**. The observed NN improvement belongs to the approximate-pair event proposal, which intentionally changes the state/trajectory. The finite event law and its perturbation bound help define that proposal; they do not prove its recovery benefit. Skipping only 5, 12 or 14 projection pairs cannot explain a large general speedup over a 3000-step run.

---

## 7. Exact finite-material profile penalty from restricted current representation

This section is a second, distinct contribution candidate. It should not be conflated with the DR event law.

At fixed material \(m\), let

\[
A_m=I-D_mG,
\]

and consider a fixed-material quadratic auxiliary objective

\[
\phi(m,J)=\tfrac12\|W(SJ-Y)\|^2+
\tfrac\rho2\|L_m(A_mJ-D_mE)\|^2.
\]

The matrices may represent one source, or a vectorized all-source state with the actual positive source metrics. \(L_m\) can be the scalar CIE transform or a frozen-at-this-material invertible G19 transform. Define

\[
C_m=\begin{bmatrix}WS\\\sqrt\rho L_m A_m\end{bmatrix},
\qquad b_m=\begin{bmatrix}WY\\\sqrt\rho L_mD_mE\end{bmatrix},
\qquad H_m^{\mathrm{cur}}=C_m^*C_m.
\]

Assume \(C_m\) has full column rank. Write \(J_*(m)\) for the unconstrained full-current minimizer. Restrict the current to an affine family \(\mathcal V=J_D+\operatorname{range}V\). The reduced profiles obey

\[
\boxed{
\phi_{\mathcal V}(m)
=\phi_{\mathrm{full}}(m)
+\tfrac12\operatorname{dist}_{H_m^{\mathrm{cur}}}^2(J_*(m),\mathcal V).
}
\]

**Proof.** The full normal equation is \(C_m^*(C_mJ_*-b_m)=0\). Pythagorean expansion gives

\[
\|C_mJ-b_m\|^2
=\|C_mJ_*-b_m\|^2+\|C_m(J-J_*)\|^2.
\]

Minimize the second term over the affine family. \(\square\)

This is a finite-\(m\), exact statement. The extra term is not merely an initialization error. It changes with material, so its gradient can favor one material over another even if the full-state objective would prefer the opposite ordering.

It is nonetheless standard least-squares/variable-projection geometry, not a new Maxwell identity. The research contribution would be proving which material directions it suppresses in the actual representation, diagnosing failures prospectively, and removing that suppression at lower cost than unrestricted current optimization.

### Legacy G19 qualification

R11/R12 legacy G19 retains a scalar CIE current block, nonlocal G19 material readout, two different betas, and PR history. It is not alternating minimization of the single common \(\phi\) above. Therefore the formula is exact for a specified fixed-current block/profile diagnostic or a new matched objective, **not an exact description of the complete legacy G19 trajectory or its Hessian**. This caveat is essential.

### A local positive-semidefinite artificial material stiffness

At an exact zero-residual compatible root (both data and constitutive residuals vanish) with current anchor represented, let \(B_Qh=D_hQ\) and hold the fixed-material transform used in the linearized residual. Define

\[
C=\begin{bmatrix}WS\\\sqrt\rho L A_m\end{bmatrix},
\qquad B=\begin{bmatrix}0\\\sqrt\rho L B_Q\end{bmatrix}.
\]

With orthogonal projectors onto the indicated column spaces, the two linearized reduced Hessians are

\[
H_{\mathrm{full}}=\operatorname{Re}\{B^*(I-P_C)B\},
\]

\[
H_V=\operatorname{Re}\{B^*(I-P_{CV})B\}.
\]

Hence

\[
\boxed{H_V-H_{\mathrm{full}}
=\operatorname{Re}\{B^*(P_C-P_{CV})B\}\succeq0.}
\]

This PSD Hessian assertion concerns the zero-residual linearized model. Away from such a root, the finite representation penalty remains nonnegative, but its exact material Hessian need not be PSD: residual curvature and derivatives of the moving metric matter. For nested Fourier/APCS/full spaces the linearized gaps split into separate PSD contributions. A high current representation error does not imply high material stiffness; this expression pulls the residual obstruction back through the shared material tangent. Conversely a small total projection error can conceal a severe restriction in one important material direction.

### Minimal directional enrichment

For a chosen material direction \(h\), compute the full linearized optimal current response

\[
x_h=C^\dagger Bh.
\]

Adding its component not already represented by \(V\) removes the profile gap exactly along that one direction. With a common spatial basis and many illuminations the addition is a block of up to \(L\) vectors, not always a single spatial vector. A truncated SVD gives an approximate low-rank completion with a residual bound in the \(C\)-metric.

This is not a global-basin theorem. A basis that includes all currently useful local directions can still fail at finite distance, and the added solve can cost as much as a forward/adjoint block. Use LSQR/QR or stable factorizations, rather than explicitly forming unnecessarily squared normal equations.

---

## 8. Elastic DPCS: useful control, not a new conceptual discovery

For receiver SVD \(S=U\Sigma V^*\), let \(V_D\) be the retained receiver subspace, \(J_D=V_D\Sigma_D^{-1}U_D^*Y\), and \(B_F\) the Fourier APCS family in \(\ker V_D^*\). Introduce retained-mode slack

\[
J=J_D+V_Db+B_Fa.
\]

Because the receiver subspaces are orthogonal,

\[
\boxed{
\|SJ-Y\|_F^2
=\|\Sigma_Db\|_F^2
+\|SB_Fa-Y_\perp\|_F^2.
}
\]

This makes a full noise-ball constraint explicit without fixing noisy retained coefficients exactly. It adds only \(r_D\times L\) complex coordinates. The allowable slack depends on how much of the error budget the complementary residual consumes; it is not a free independent radius for every retained mode.

R8 already identifies the conflict between noisy fixed DPCS and learned physical initialization and recommends addressing it. Therefore this is a warranted implementation/ablation, not a new research direction discovered here. Its real value is disentangling representation freedom from reflection dynamics under identical noise semantics.

At fixed \(m\), minimizing the state residual under the global data ball is a convex quadratic constrained least-squares problem. Its scalar multiplier produces a monotone residual curve. Existing R12 data-ball projection machinery and low-rank receiver structure can be reused. Releasing DPCS changes the feasible auxiliary family and, if replacing a penalty by a hard ball, also the objective contract; factorial comparisons must distinguish those changes.

---

## 9. Observability, angles and corrector limits

The physical tangent remains

\[
\delta J=(I-D_mG)^{-1}D_hQ,
\qquad DF(m)[h]=S(I-D_mG)^{-1}D_hQ.
\]

This is already in the atlas. At an exact equality-data intersection, eliminating fields maps the tangent intersection of \(\mathcal A\) and \(\mathcal C\) to \(\ker DF(m)\), assuming the physical resolvent is invertible. Nonzero principal-angle values require the actual state metric and material lifting; they are not raw singular values of \(DF\).

For a strict data-ball interior point, the inequality is locally inactive and every nearby physical material tangent belongs to the feasibility tangent intersection. Therefore a DR feasibility angle cannot quantify noise-level statistical identifiability there. A noise-weighted material Jacobian, finite candidate separation and independent measurement validation are separate objects.

At a generic nonzero-normal DR fixed state \(z_*=x_*+n_*\), the derivative of the nonconvex projection includes curvature of the normal-offset manifold. Replacing it by a tangent-space projector is generally wrong. The familiar DR \(\cos\theta\) subspace rate applies to a much narrower linear/zero-offset situation; it cannot explain all saved trajectories.

A discrepancy-preserving proximal physical repair,

\[
\min_m\tfrac12\|m-m_c\|_M^2
\quad\text{s.t.}\quad\|W(F(m)-Y)\|\le\delta,\ m\in\mathcal M,
\]

has the simple guarantee that an already feasible \(m_c\) remains unchanged. It separates repair from prior-based branch choice. This is established constrained regularization, not novel, and a local nonlinear solver does not ensure the globally closest feasible material. R11 also showed that simply stopping at the first noise-ball crossing did not fix its tested failures; those negative results must be retained.

---

## 10. Originality and paper risk

### Established foundations

- Exact resolvent/secant identities and distorted source/receiver products: inverse scattering and difference factorization; also Lippmann–Schwinger–Lanczos work.
- Shared-material rank-one or total-least-squares projection: R12 already derives the quadratic slope equation and notes Gram degeneracy.
- Normal drift from nonfeasible AP pairs: difference-map/DR literature, explicitly Elser and coauthors.
- Unconstrained singular-value crossing and rank-one projection curvature: standard matrix geometry.
- Profile penalties, Schur complements, variable projection and reduced-basis errors: standard optimization/model-reduction tools.

### Potential new, narrower contribution

A bounded shared-material first-event law with explicit active-bound formulas; its exact and approximate DR semantics; a finite state-dependent prediction tested on actual Maxwell states; and a demonstrably cheaper or better controlled event-driven search. The representation profile theorem can serve as a complementary attribution/control framework rather than another main algorithm.

### Strongest reviewer objections

1. “This is an application of rank-one projection and difference-map drift.” Correct unless independent Maxwell mechanism prediction and a useful algorithmic consequence are established.
2. “The event predicts a box flip, not recovery.” Correct. The paper must separate escape, feasibility, and material accuracy.
3. “You paid thousands of AP steps to find the event.” Correct for the current probe. This is mechanism evidence, not acceleration evidence.
4. “A large random perturbation would do the same.” Requires matched-norm, matched-cost non-geometric controls.
5. “The observed benefit is from full currents, not reflection.” Requires full/restricted representation factorial tests including elastic measured modes.
6. “The event law only holds at exact AP stagnation.” The perturbed bound is explicit and becomes weak near curvature singularities; that limit should be reported, not hidden.
7. “The difficult examples are already exposed.” Current evidence is development evidence. New confirmation must be generated and scored independently after the algorithm and forecast rules are fixed.

### Suggested allocation

Prioritize the finite constitutive-event mechanism and a small, precise algorithm test. Keep profile-penalty and elastic-DPCS work as the necessary representation control. Reduce investment in unsupported universal claims about Green rank, local condition number, and local DR angle. Do not commit to a high-impact recovery paper until the mechanism survives independent objects and the intervention improves useful outcomes at full cost.

---


## 11. Adversarial audit of the finite-event statement

The following conditions are part of the theorem, not optional implementation details.

### Dark cells and degenerate material intervals

- If a cell has two distinct allowed material slopes and its exact projected state is zero, then its input state must also be zero. Indeed, a zero projection onto every allowed line is orthogonal to two independent row directions and hence to the full two-row space. Thus an exactly dark projected cell has `a=b=gamma=0`; its material label is unidentifiable and must not be counted as a physical material branch event. Assign no event to this cell; do not evaluate `0/0` or interpret a default material label change as recovery.
- If `a>0,b=0`, there is no normal movement in that cell and its threshold is infinity.
- With a nontrivial allowed interval, `a=0,b>0` cannot arise from an exact global constitutive projection. If it appears numerically, first inspect amplitude stationarity, projection accuracy, and the state convention. An artificial small positive denominator changes the problem.
- A fixed known material, `lower=upper`, gives a linear subspace. A nonzero normal ray can remain projected onto that subspace forever. Therefore the finite-escape corollary requires at least one nonzero-gap cell with a nontrivial interval, and a version claiming finite escape for every nonzero gap requires nontrivial intervals at every cell where such a gap may occur. It does not automatically cover fixed background/host cells.
- A lower material bound of zero is allowed: `theta_lower=0` is regular and the derivation never divides by material. An infinite upper bound is outside the stated finite-material theorem. Adding the projective direction `theta=pi/2` would add limiting states with `Q=0,J≠0`, which are not finite-contrast constitutive states.

### Dimensional and metric assumptions

The real isotropic vector case is covered only because the same scalar material multiplies every field component and illumination at a cell. Concatenating three complex components changes the row length from `2L` to `6L`, not the two-row rank-one geometry. An unknown complex lossy scalar, a tensor material, magnetic/electric coupled unknowns, or source-dependent material does not obey this bounded real-angle theorem without a separate derivation.

The metric must remain fixed during the ray. A fixed cellwise rescaling and a common positive source metric for both rows can be absorbed by coordinate transformations. A metric that couples cells destroys the product projection; different nonproportional source metrics for current and field do not preserve the simple common two-row director; changing the scale during the jump changes the ray, its normal cone, and the threshold. A common scalar multiplication of the whole metric leaves the projection and threshold invariant, but changing the relative current/field scale generally does not. Therefore the threshold is an algorithm-state quantity, not a coordinate-free physical observable.

### Projection ties and numerical semantics

At an interior event all allowed slopes tie. At a bound-active event with nonzero gamma the two endpoint branches tie. Multiple cells may attain the same minimum threshold. One must preserve the actual projector's deterministic tie rule when claiming exact trajectory equivalence. For a symmetric interior material orientation, the two farthest endpoints can remain tied even just beyond the event.

Use an explicit initial-state indexing convention. For `z0=c+d`, noninteger `t*` predicts the first potentially changed projection at zero-based index `ceil(t*)-1`. At an integer `t*`, the tie may retain the old shadow, so the guaranteed change can occur one step later. Exact fast-forwarding ends at the first tie/crossing state and lets the actual projector decide; it must not silently force the far endpoint at an unresolved tie.

The approximate-pair bound assumes that `c=P_C(a)` is itself exact. An inexact constitutive projection requires additional error terms. Audit `u^T d`, interior gamma or active-bound sign, and the candidate objective values separately from `||P_A c-a||`. Near dark cells, nearly repeated thresholds, and projection ties, floating-point classifications can change the first event. Reporting only a rounded material value and a threshold is insufficient for a numerical theorem claim.

### What is and is not new beyond R12 and Elser

R12 already presents the shared-material quadratic slope equation, real 2-by-2 Gram degeneracy, possible projection branch switches, and the caveat that local smooth rates cannot explain full trajectories. Elser and coauthors already explain normal drift away from nonfeasible AP pairs. The new derivation here combines those ingredients to obtain an explicit bounded-material first-event formula and a near-pair error accounting. The active-bound formula is a short consequence of the existing projection formula, so its algebra alone is unlikely to carry a strong paper.

Moreover, the event theorem applies to any bounded shared-scalar bilinear constraint with a convex partner set. Its statement does not use a special property of the Maxwell Green operator. Maxwell-specific content must come from identifying which propagation/receiver fields create the normal state, proving useful consequences for those fields, and testing predictive or algorithmic value across physical regimes.

The strongest conceptual threat is that **exact skipping cannot improve recovery**, because it reproduces the original trajectory. The observed NN benefit must be attributed to the approximate-pair finite proposal. Its mechanism requires further study; it cannot borrow a recovery guarantee from exact fast-forward equivalence. A matched-length random or generic extrapolation control is necessary before claiming the derived direction matters.

## 12. Maxwell source/receiver content of the AP normal field

The event law is generic geometry, but the normal field that feeds it has an exact Maxwell characterization. This supplies a useful bridge to illumination and reception without pretending that it is new observability information.

Assume real scalar material, fixed scale `s`, exact AP best pair `(a,c)`, invertible `A_m=I-D_m G`, and a strictly feasible current for the data-ball constraint (Slater condition). Let

\[
a_Q=E+sG a_X,\qquad
r=W(sS a_X-Y),\qquad
\|r\|\le\delta.
\]

Use the half-squared ball inequality and let `lambda>=0` be its scalar KKT multiplier. The projection of `c` onto `A_delta` gives

\[
d_X+sG^*d_Q+\lambda sS^*W^*r=0.
\]

The constitutive amplitude-normal condition, valid even at material bounds, gives

\[
d_Q=-D_{m/s}d_X.
\]

Consequently

\[
\boxed{
A_m^*d_X=-\lambda sS^*W^*r.
}
\]

Define the residual-driven adjoint field

\[
p=A_m^{-*}S^*W^*r.
\]

Then

\[
\boxed{d_X=-\lambda s p,\qquad d_Q=\lambda D_m p.}
\]

These signs follow directly from differentiating the nearest-state Lagrangian with equality `Q-E-sGX=0`; the data residual is predicted minus measured data. If the opposite residual convention is used, the adjoint sign must be changed consistently.

### Consequence 1: a nonzero-gap AP trap requires active data pressure

If the data constraint is inactive, `lambda=0`. Invertibility of the Maxwell resolvent then implies `d_X=d_Q=0`. Therefore, under these assumptions, an exact nonfeasible AP stagnation pair must have a strictly positive data-ball multiplier and lie on the ball boundary. A propagation-only incompatibility cannot support such a pair for a regular material; the unique physical field at that material removes it.

This statement fails without the invertibility/qualification assumptions. At zero noise radius, the norm-squared inequality gradient vanishes on equality data and a finite scalar multiplier description is unsuitable; use an unrestricted data-equality multiplier instead. At an unattainable or minimally attainable ball radius, Slater may fail. At a singular resolvent, a nonzero normal field may lie in `ker A_m^*` even with zero data multiplier. None of these cases should be hidden under the same formula.

### Consequence 2: the threshold depends on source and receiver fields jointly

At cell `i`, write `q_c=Q_c`. The normal coordinates become

\[
v_i=\sqrt{1+(m_i/s)^2}\,q_{c,i},
\qquad
w_i=-\lambda s\sqrt{1+(m_i/s)^2}\,p_i.
\]

Thus an interior material satisfies

\[
\operatorname{Re}\langle q_{c,i},p_i\rangle=0,
\qquad
\boxed{t_{*,i}=\frac{\|q_{c,i}\|}{\lambda s\|p_i\|}.}
\]

For a bound-active material, put

\[
Q_i=\|q_{c,i}\|^2,\quad P_i=\|p_i\|^2,
\quad z_i=\operatorname{Re}\langle q_{c,i},p_i\rangle.
\]

The threshold can be written

\[
\boxed{
 t_{*,i}=
\frac{z_i\cot\varphi_{i,f}
+\sqrt{z_i^2\cot^2\varphi_{i,f}+Q_iP_i}}
{\lambda sP_i}.
}
\]

At the lower bound `z_i>=0`, while at the upper bound `z_i<=0`. Both the source-dependent shared field and the receiver-adjoint residual field enter. A Green singular value alone omits the actual source excitation, receiver coverage, noise metric, active data pressure and material resolvent. This is a structural reason why static Green energy cannot forecast the event.

`q_c` is an auxiliary shared-material field, not automatically the physical total field `Q(m)`. Likewise `p` is driven by the auxiliary data residual. The formula is an event predictor for an optimization state, not a Fisher-information or GT-recovery formula.

### Consequence 3: AP stagnation is a noise-constrained extended-source stationary point

Let

\[
T_m=D_{s^2+m^2},\qquad J_m=A_m^{-1}D_mE,
\]

and `a_J=s a_X`. The constitutive defect at the propagation/data point is

\[
a_J-D_m a_Q=-\lambda T_m p.
\]

Therefore

\[
a_J=J_m-\lambda A_m^{-1}T_m p.
\]

With physical residual `r_phys=W(SJ_m-Y)`, define the positive-semidefinite receiver-space operator

\[
H_m^{\mathrm{rec}}=WSA_m^{-1}T_mA_m^{-*}S^*W^*.
\]

Eliminating the extended state yields

\[
\boxed{r=(I+\lambda H_m^{\mathrm{rec}})^{-1}r_{\rm phys},}
\qquad
\boxed{\|d\|^2=\lambda^2\langle r,H_m^{\mathrm{rec}}r\rangle.}
\]

For a common receiver map this is one receiver-space operator acting on all illumination columns, with a single global multiplier for the Frobenius data ball. More general noise covariance uses the vectorized observation operator. The scalar multiplier is determined by the ball boundary, subject to any unreachable receiver residual floor.

Equivalently, for fixed material the joint closest-state objective is

\[
\psi_\delta(m)=
\min_{\|W(SJ-Y)\|\le\delta}
\frac12\left\|T_m^{-1/2}
(A_mJ-D_mE)\right\|^2.
\]

The fixed-material problem is convex and, with invertible `A_m`, strictly convex in current. At an AP best pair its current is the unique minimizer. Where the envelope is differentiable, its material derivative is

\[
\partial_{m_i}\psi_\delta
=\lambda\operatorname{Re}\langle q_{c,i},p_i\rangle,
\]

with the corresponding material box KKT signs. This can be checked directly: differentiating the factor `1/(s^2+m_i^2)` creates a term that cancels the part from `a_Q=q_c+lambda D_m p`.

This is a concrete Maxwell realization of established extended-source/WRI and variable-projection ideas. It explains why AP can be stationary in a profiled compatibility objective while the direct physical objective has a different gradient: the two use different state and residual fields. It does not turn the profile into a convex material problem. Every physical material already inside the data ball gives `psi_delta=0`, so this profile cannot choose GT among such materials.

### Assessment of added value

The strongest additional result is the exact restriction that a regular nonzero-gap AP trap must be supported by an active receiver-data multiplier, together with the source/receiver field expression for the branch threshold. This connects the otherwise generic event law to Maxwell propagation and helps design controls: remove or change receiver constraints, illumination or noise radius and observe how the *actual* dual field and event change. It does not justify a monotonic claim in noise radius, because the auxiliary field, multiplier and AP pair all change together.

The resolvent adjoint, receiver-space filter and extended-source profile are familiar structures. They should be presented as a rigorous bridge and a mechanism hypothesis, not as a new inverse-scattering identity. The decisive open question remains whether these state-dependent quantities predict useful changes in material recovery on independent instances, beyond predicting a generic projection branch flip.


## References for parent verification

1. Elser, Rankenburg & Thibault, *Searching with iterated maps*, PNAS (2007), DOI: https://doi.org/10.1073/pnas.0606359104 . Generic normal drift / difference-map search, not Maxwell event novelty.
2. Elser, *Matrix product constraints by projection methods*, arXiv: https://arxiv.org/abs/1601.01003 . Matrix product projection and proximal-pair escape discussion.
3. Hesse & Luke, *Nonconvex Notions of Regularity and Convergence of Fundamental Algorithms for Feasibility Problems*, SIAM J. Optimization (2013), https://doi.org/10.1137/120902653 . Local regularity scope.
4. Lewis & Malick, *Alternating Projections on Manifolds* (2008), https://doi.org/10.1287/moor.1070.0291 . Local manifold angles.
5. Luke & Martins, *Convergence Analysis of the Relaxed Douglas–Rachford Algorithm*, https://doi.org/10.1137/18M1229638 . Inconsistent/nonconvex regularity.
6. van Leeuwen & Herrmann, *Mitigating local minima in full-waveform inversion by expanding the search space*, GJI 195 (2013), https://doi.org/10.1093/gji/ggt258 . Extended wavefield search and penalties.
7. Zhong & Chen, *Twofold subspace-based optimization method for solving inverse scattering problems*, Inverse Problems 25 (2009) 085003, https://doi.org/10.1088/0266-5611/25/8/085003 . Fixed current subspaces and staged refinement.
8. Druskin, Moskow & Zaslavsky, *Lippmann–Schwinger–Lanczos algorithm for inverse scattering problems*, https://arxiv.org/abs/2101.12317 ; https://doi.org/10.1088/1361-6420/abfca4 . Finite difference product fields and data-driven reduced internal fields; no novelty claim for the secant identity alone.

Repository originals additionally inspected include `research_max_struct_r11/THEORY_DERIVATIONS.md`, R12 theory and `src/field_current.py`, supplement source `run_paired.py` and audit, R8 `cie_fft.py` and neural warm-start reports, R6 Fourier/APCS implementation, and formula atlas sections on shared material, material lifts, bilateral modes and variable projection.
