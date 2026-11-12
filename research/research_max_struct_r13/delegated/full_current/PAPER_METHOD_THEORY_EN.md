# Bounded material events and paid physical candidate selection

## Physical model and lifted geometry

We use the time convention \(\exp(-i\omega t)\) and reconstruct a real isotropic susceptibility \(m=\chi=\varepsilon_r-1\) within a prescribed box \([\ell,u]\), with \(0\leq\ell<u<\infty\). All sources, and all three field components in a vector cell, share one material value. The discrete equations are

\[
J=D_mQ,\qquad Q=E+sGX\;(J=sX),\qquad F(m)=SJ.
\]

Here \(J\) is a normalized contrast source, not an SI current measured in amperes. The scalar-TM disk-pulse VIE and vector CM+RR dyadic DDA retain their respective Green functions and self terms. Shared algebra does not identify their continuum errors. Source-dependent measured reception uses its actual block operator and adjoint; absent channels are not zero-filled.

For the scalar synthetic acquisition, the background wavelength is one, \(k=2\pi\), and the inverse region is a \(0.96\)-wavelength square discretized into \(24^2\) independent real materials. Six or twelve cylindrical line sources lie at radius \(1.4\), with incident field proportional to \((i/4)H_0^{(1)}\); these are not plane waves. Reception uses the declared full or partial aperture at radius \(1.1\). Classical plane-wave, vector-dipole and measured acquisitions remain separate physical contracts. The scalar synthetic host is \(m=0.1\), rather than vacuum support.

Fix \(s>0\) and equip \(z=(X,Q)\) with

\[
\|z\|^2=\|X\|_F^2+\|Q\|_F^2.
\]

The two constraint sets are

\[
\mathcal C=\{(X,Q):sX=D_mQ,\ m_i\in[\ell,u]\},
\quad
\mathcal A=\{(X,Q):Q=E+sGX,\ \|W(sSX-Y)\|_F\leq\delta\},
\]

where \(W\) is fixed and declared; the raw synthetic ball uses \(W=I\). The nonempty set \(\mathcal A\) is closed and convex, whereas \(\mathcal C\) is generally nonconvex. The common projection metric uses \(s=6\) in the stated experiments. Changing this scale changes the event geometry.

At cell \(i\), concatenate sources and vector components and realify complex entries, obtaining \(Z_i\in\mathbb R^{2\times d}\). With

\[
\theta=\arctan(m/s),\quad
u(\theta)=(\sin\theta,\cos\theta)^T,
\]

the material set is \(Z_i=u(\theta_i)v_i^T\). At fixed angle, the optimal amplitude is \(v_i^T=u(\theta_i)^TZ_i\). Thus material projection maximizes \(u(\theta)^TK_iu(\theta)\), where \(K_i=Z_iZ_i^T\), over the legal angle interval. Endpoints and stationary angles suffice; no angular grid search is required. This representation applies to bounded shared real materials, not arbitrary loss or anisotropic tensors.

## Finite events on a constitutive normal ray

Let \(c=P_{\mathcal C}(a)\), \(d=a-c\), and suppress a cell index. Amplitude stationarity gives \(u^Td=0\), so

\[
c=uv^T,\qquad d=u_\perp w^T,\qquad
a_0=\|v\|^2,\ b_0=\|w\|^2,\ \gamma=\langle v,w\rangle.
\]

For the ray \(c+td\), compare the current angle with a relative angle \(\phi\). Substitution of \(u(\theta+\phi)=u\cos\phi+u_\perp\sin\phi\) gives the exact projection-energy difference

\[
\Delta(\phi,t)=(t^2b_0-a_0)\sin^2\phi
                  +2t\gamma\sin\phi\cos\phi.
\tag{1}
\]

For an interior material, angle stationarity implies \(\gamma=0\). Since \(c\) is optimal at \(t=1\), \(a_0\geq b_0\). If \(b_0>0\), the first competition threshold is

\[
t^*=\sqrt{a_0/b_0}.
\tag{2}
\]

Before the threshold, every nonzero legal angle change has smaller energy. After it, the maximum is attained at a legal endpoint farthest from the current angle. Equal-distance endpoints and equality at the threshold may be tied; uniqueness is not asserted there.

At the lower endpoint, one-sided optimality gives \(\gamma\leq0\). For \(0<\phi\leq\phi_{\rm far}\),

\[
\frac{\Delta(\phi,t)}{\sin^2\phi}
=b_0t^2-a_0+2t\gamma\cot\phi.
\]

Because the legal interval has width less than \(\pi/2\), this expression is largest at the opposite endpoint. The upper-endpoint argument is analogous. The first boundary competition therefore solves

\[
b_0t^2+2\kappa t-a_0=0,
\qquad \kappa=\gamma\cot\phi_{\rm far}.
\tag{3}
\]

Its positive root is \((\sqrt{\kappa^2+a_0b_0}-\kappa)/b_0\); when cancellation is possible it is evaluated as \(a_0/(\sqrt{\kappa^2+a_0b_0}+\kappa)\). A zero \(b_0\) gives an infinite threshold, without an artificial denominator floor. Dark states, persistent ties and threshold equalities are recorded as degeneracies, not identifiable material events.

The first global event is \(\min_i t_i^*\). Its arithmetic requires field-energy contractions and a minimum, rather than a material Jacobian or Hessian. Obtaining the underlying projected state still costs propagation/data projections. This finite-dimensional law belongs to bounded bilinear projection geometry; it is not Maxwell-exclusive and does not identify a correct material branch.

Excluding dark cells, persistent symmetric ties and threshold equalities, the material readout is constant between consecutive distinct events on this ray, while the projected amplitudes are affine in t within a fixed material branch. This does not identify an interval of equivalent future reflection trajectories.

## Exact drift, approximate proposals and retained memory

Use the Douglas–Rachford map

\[
T=I+P_{\mathcal A}(2P_{\mathcal C}-I)-P_{\mathcal C}.
\]

Suppose \(a=P_{\mathcal A}(c)\), \(c=P_{\mathcal C}(a)\), and \(d=a-c\). While \(P_{\mathcal C}(c+td)=c\), convex normal-cone geometry gives \(P_{\mathcal A}(c-td)=a\), hence

\[
T(c+td)=c+(t+1)d.
\tag{4}
\]

Skipping repeated unchanged-shadow projections on this exact ray deletes redundant operations; it does not create additional recovery capability. Normal drift and discrepancy accumulation have established projection-method precedents, including [Elser's matrix-product formulation](https://arxiv.org/html/1601.01003v1).

For an approximate pair, let \(a\in\mathcal A\), \(c=P_{\mathcal C}(a)\), and \(e=P_{\mathcal A}(c)-a\). Put \(a'=a+e\). The normal-cone property yields

\[
P_{\mathcal A}\bigl(a'+(t+1)(c-a')\bigr)=a',\qquad t\geq0.
\]

Its argument equals \(c-td-te\). Nonexpansiveness of \(P_{\mathcal A}\) and the triangle inequality consequently bound the one-step deviation from (4) by \((t+1)\|e\|\), provided the same constitutive shadow is retained. This is not an unconditional multi-step trajectory or event-time bound. Fixed-branch nonexpansive arguments require the relevant iterates to remain in that branch; interior projection sensitivity can deteriorate as the Gram eigenvalue gap closes.

Finite AP preparation and a DR-state construction are therefore treated as distinct proposal mechanisms. The former applies a declared finite number of AP steps; the latter forms one \(\mathcal A\)-shadow from the current DR constitutive shadow. Neither asserts exact future fast-forwarding. Integer targets use the first integer strictly beyond the tie; fractional targets use a predeclared overshoot. Matched-norm random displacements, a same-direction \(t=2\) action, no-event shadow reset, physical-memory reset and Anderson controls isolate different effects. They do not inherit a convergence guarantee merely from their names. Established nonconvex projection convergence results require explicit regularity assumptions; see [Hesse and Luke](https://doi.org/10.1137/120902653).

If the initial physical state satisfies the legal constitutive constraint, \(z_0\in\mathcal C\), then AP and DR have the same first step:

\[
P_{\mathcal A}P_{\mathcal C}z_0=P_{\mathcal A}z_0=Tz_0.
\]

This equality does not extend automatically to subsequent iterates.

### An exact same-shadow counterexample

Identical material readouts do not generally define a closed reduced dynamical system. Consider the abstract real geometry

\[
s=1,\quad \mathcal C=\{(mv,v):0\leq m\leq1,\ v\in\mathbb R\},\quad
\mathcal A=\{(X,Q):Q=1-X,\ |X|\leq3\}.
\]

The coefficient \(G=-1\) here is abstract; it is not a Maxwell Hankel or dyadic self term. Take \(z_1=(-1,1)\), \(z_2=(-11/5,1)\). For \(z=(-a,1)\), projection energy is \((1-am)^2/(1+m^2)\). The endpoint energies are \((1,0)\) for \(a=1\), and \((1,18/25)\) for \(a=11/5\); an interior stationary point cannot improve them. Thus both have the strict identical full shadow \(P_{\mathcal C}z_j=(0,1)\).

The reflected states are \((a,1)\). Their projections onto the affine line are \((a/2,1-a/2)\), inside the segment. Consequently

\[
Tz_1=(-1/2,1/2),\qquad Tz_2=(-11/10,-1/10).
\]

The next material projections select respectively \(m=0\) and \(m=1\). Even an identical complete constitutive shadow cannot determine the next material readout. This counterexample concerns discarded state memory, not which branch is closer to truth, nor whether resetting memory is always harmful.

## A conditional state-ball branch bound

For a realified cell state \(Z\) and perturbation \(E\), assume \(\|E\|_F\leq r\). Expanding the Gram product gives

\[
\|(Z+E)(Z+E)^T-ZZ^T\|_2
\leq2\|Z\|_Fr+r^2=\rho.
\]

For any closed legal angle interval \(I\), define \(h_I(K)=\max_{\theta\in I}u(\theta)^TKu(\theta)\). Unit-length \(u\) implies \(|h_I(K+\Delta K)-h_I(K)|\leq\|\Delta K\|_2\). Split the legal angle interval into low/high halves and let \(M=h_H-h_L\). Then

\[
|M(Z+E)-M(Z)|\leq2\rho.
\tag{5}
\]

If \(M>2\rho\), every optimal projection in that state ball lies outside the low half; the opposite sign excludes the high half. Strictness also excludes a midpoint-only tie. If \(r\) is an observed deviation between actual DR and a reference ray, this is an a posterior branch certificate. Without a valid future radius bound it is not a future-event certificate. An inconclusive margin remains inconclusive.

## Maxwell-specific dual pressure and its approximate-pair correction

Let \(A_m=I-D_mG\) be invertible. Assume exact amplitude stationarity for \(P_{\mathcal C}\), compatible propagation, a positive-radius ball satisfying Slater's condition on the propagation graph, and a valid finite multiplier \(\lambda\) for the half-squared ball constraint. For \(a'=P_{\mathcal A}(c)=a+e\), define

\[
r'=W(sSa'_X-Y),\qquad
p'=A_m^{-*}S^*W^*r'.
\]

Propagation gives \(e_Q=sGe_X\), while constitutive amplitude normality gives \(d_Q=-D_{m/s}d_X\). Projection KKT acts on \(d+e\):

\[
d_X+e_X+sG^*(d_Q+e_Q)+\lambda sS^*W^*r'=0.
\]

Therefore

\[
\boxed{d_X+\lambda sp'
=-A_m^{-*}(e_X+sG^*e_Q)
=-A_m^{-*}(I+s^2G^*G)e_X.}
\tag{6}
\]

For an exact pair, \(e=0\), this reduces to \(d_X=-\lambda sp\), \(d_Q=\lambda D_mp\). At a nondegenerate interior event,

\[
t_i^*=\frac{\|Q_{c,i}\|}{\lambda s\|p_i\|}.
\]

These are auxiliary residual-adjoint fields, not truth fields. Invertibility and \(\lambda=0\) imply zero exact gap. Zero-radius, unreachable or non-Slater limits and singular Maxwell matrices require separate treatment.

For approximate pairs, scalar \(\|e\|/\|d\|\) omits the directional resolvent action in (6). A retained material LU can apply \(A_m^{-*}\) to the receiver and correction RHS in one batch. The resulting directional gain is neither an operator-norm certificate nor a recovery guarantee. Receiver-Gram and adjoint representations also arise in established extended-field elimination; see [van Leeuwen and Yang](https://arxiv.org/html/2403.15292v2).

The Q-component correction is \(d_Q-\lambda D_mp'=-D_{m/s}(d_X+\lambda sp')\). Diagnostics use the actual finite multiplier returned by the data-ball projection, not a coefficient fitted to the gap. A bound on \(\|A_m^{-*}[I,sG^*]\|\), if independently justified, would control this correction by \(\|e\|\); measuring its action on one RHS does not supply that bound.

## Fixed-current profiles and baseline distinctions

Fix one explicitly specified single-\(\beta\) current block \(\tfrac12\|C_mJ-b_m\|^2\), with full-column-rank \(C_m\). If \(J^*\) is its unrestricted minimizer, the normal equations eliminate the cross term, giving, for an affine current family \(\mathcal V\),

\[
\varphi_{\mathcal V}(m)=\varphi_{\rm full}(m)
+\tfrac12\operatorname{dist}_{C_m^*C_m}^2(J^*(m),\mathcal V).
\tag{7}
\]

A common material prior adds equally to both sides. The representation penalty depends on m and differs from an initial Euclidean current-projection loss. Equation (7) concerns that fixed block, not a common Lyapunov function for an entire dual-\(\beta\) G19 trajectory.

The matched G19 controls retain its PR current loop, material denominator/readout, prior and beta schedule, changing only \(JD+PF_{\rm mf}\alpha\), full APCS \(JD+PF_{\rm full}\alpha\), or full current \(JD+F_{\rm full}\alpha\). On \(24^2\), the original final mf12 already spans full APCS. Early truncation, final capacity and fixed DPCS must be separated. Independent real current capacity is \(2n_s\operatorname{rank}(PF_{\rm mf})\), not the redundant coefficient count; material capacity remains N.

CIE/Dual-CIE retain the declared real-material/common-prior adaptations of [the contraction-integral formulation](https://doi.org/10.3390/jimaging5020027). A global PR coefficient in Dual-CIE is distinct from **physical GlobalDual**, which optimizes the material-only total-field source-Gram objective at \(\beta=0\); \(\beta=.2\) denotes a different feedback variant. Full-J TRF uses the inherited physical residual/Jacobian with \(x_{\rm scale}=\mathrm{jac}\), whereas physical L-BFGS-B uses the inherited analytic gradient. For scalar prior edges \(d\), \(q=\sqrt{d^2+\epsilon^2}\), the prior residual \([\sqrt{2\alpha/M}d;\sqrt{2\eta/M}d/\sqrt{q+\epsilon}]\) has half-squared norm equal to \(\alpha\operatorname{mean}(d^2)+\eta\operatorname{mean}(q-\epsilon)\). This verifies the residual objective, not an exact nonlinear Hessian.

For the scalar native schedule, masks mf3/6/10/12 permit 500/500/500/1500 updates. Single-beta CIE uses current/material beta1/.2/.2/.2; Dual-CIE retains current beta1 while changing material beta. CIE/Dual/A+ initialize zero APCS coefficients above JD. The G19 controls project the same initial physical current into their respective admissible families. Material variables are not reduced by the current basis. Setup and warm-projection costs are charged.

More explicitly, physical GlobalDual uses

\[
F_{\rm GD}(m)=\frac{c_0}{2}\operatorname{tr}\!\left[
\mathcal E(Q^*Q)^{-1}\mathcal E^*\right]+R(m),\qquad
\mathcal E=F(m)-Y,
\]

where \(c_0=\|Q_{\rm ref}\|_F^2/(n_s\|Y\|_F^2)\) and the reference is fixed. Under the declared uniform0.1 initialization, it matches the historical uniform reference. Its analytic derivative includes the changing Gram; it is not ordinary LS with a frozen weight. TRF and physical L-BFGS-B instead optimize the original physical objective through different numerical models. Finite CPU caps can truncate any native schedule, so configured iterations are not completed iterations.

The complete-current NFFT-equations control uses ordinary LSIE, \(JD+F\alpha\), data weights \(1/\|Y_s\|^2\), state weights \(1/\|JD_s\|^2\), the declared prior/box and joint L-BFGS-B. It is explicitly an adaptation of [Wei–Chen's equations](https://usncursi.org/archive/aps-ursi/2018/papers/1927.pdf), **not an original NFFT-SOM author-algorithm reproduction**.

Its joint objective is \(\sum_s[\|SJ_s-Y_s\|^2/\|Y_s\|^2+\|J_s-m(E_s+GJ_s)\|^2/\|JD_s\|^2]+R(m)\), without half factors. Zero initial correction means J=JD, not zero current. Complete unitary Fourier and spatial-correction coordinates share the global Euclidean metric; no current/material rescaling is introduced. This normalization and optimizer must not be silently equated with the CIE current block or the underspecified original conference algorithm.

These GlobalDual, NFFT-equations and TRF control specifications concern the common-port real scalar implementation; vector and measured models retain their separately qualified interfaces.

## Paid audits, correction and provenance

Online methods use only public inputs and frozen configurations. Raw materials undergo the original physical audit. The primary estimator minimizes paid, legal, in-budget \(\Phi=D+R\); the secondary minimizes R among discrepancy-feasible candidates, with absence reported and primary fallback distinguished. Raw, corrected and primary are retained separately. The common corrector optimizes the same physical objective; lower Phi does not imply lower material error.

Synthetic circular-complex noise uses \(\delta=\sigma\sqrt{\chi^2_{2d_{\rm obs}}(.99)/2}\), where \(d_{\rm obs}\) is the number of effective complex observation channels. Cross-grid error is outside this noise model. Measured/classical diagnostic balls are not relabeled CN probability certificates. Scalar and vector physical objectives retain their different inherited normalizations.

Proposal audits follow declared absolute process-CPU buckets after completed updates, preserving and charging native audits. Physical optimizer trials already paid for forward evaluation remain eligible. Atomic operations can overshoot; spent work is counted and late candidates excluded. Prefix selection cannot use future candidates or truth. CPU, wall time, LU/RHS, projection factors, FFT/Gram work, setup and failures remain separate; hosts are not pooled into a fictitious speedup. Source shipping/config references are distinguished, missing costs remain unknown, and truth scoring follows completed-queue prediction sealing.

Archives contain the material candidate pool, audit metadata and selected physical field snapshots. They do not contain every CIE/G19 auxiliary state or PR restart vector. Separate profile/short-state diagnostics do not fill that gap or establish exact optimizer continuation. These coverage limits accompany any subsequent interpretation.
