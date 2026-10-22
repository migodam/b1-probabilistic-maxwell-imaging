# Read-only probability and source-adjoint review

Status: **CONDITIONAL_PASS_STATIC_REVIEW; PHYSICAL_FINITE_DIFFERENCES_NOT_RUN**.
The source-adjoint formula and ideal probit prior transport are mathematically
consistent. Two concrete helper/cache API bugs and several integration
conditions require attention before a general reusable reference wrapper.
The parent Codex retains final mathematical judgment and owns paid checks.

Scope: `src/p3l/core.py`, relevant `src/p3l/kernels.py` interfaces, and the
already existing A30 adapter's `solve_anchor`/`receiver_transfer` contracts.
No core, ledger, protocol, or physics file was edited. No `Target`, `Book`,
backend, state solve, Jacobian or physical gradient was instantiated/called.
Only NumPy/SciPy probability helpers were evaluated. Maxwell/RHS cost: **0**.

## 1. Likelihood and source adjoint

Use complex Frobenius inner product

\[
\langle A,B\rangle_F=\operatorname{tr}(A^*B),\qquad
r=Y(q)-y,\qquad \Phi(q)=\|r\|_F^2/\sigma^2.
\]

The current-form state satisfies

\[
L j=X e,\qquad L=I-XG,\quad X=\operatorname{diag}(\chi),
\quad E=e+Gj,\quad Y=Sj.
\]

Differentiating in a real latent direction \(dq\) gives

\[
L\,dj=\operatorname{diag}(d\chi)E,
\qquad dY=SL^{-1}\operatorname{diag}(d\chi)E.
\]

Define the eight-column source adjoint by

\[
L^* a=S^*r.
\]

Then

\[
d\Phi=\frac{2}{\sigma^2}\operatorname{Re}
\sum_{i,s}\overline{a_{is}}E_{is}\,d\chi_i.
\]

Let \(D_{ik}=\partial\chi_i/\partial q_k\), with real \(q_k\), and
\(h_i=\sum_s\overline{a_{is}}E_{is}\). Therefore

\[
\boxed{\nabla_q\Phi=\frac{2}{\sigma^2}\operatorname{Re}(D^T h).}
\]

This agrees exactly with `core.py:169–176`. In particular **`D.T` without
complex conjugation is correct**: the conjugation has already entered through
\(\overline a\). Replacing this by `D.conj().T @ h` would be wrong.
There is no missing material-volume or host factor: these are already in the
frozen state law and the decoder derivative. This statement depends on the
adapter actually using the stated current equation, as its source indicates.

`Target.phi_gradient` returns the gradient of positive likelihood potential
only. It is neither a posterior gradient nor a log-target gradient. In original
coordinates, inside the property box,

\[
\nabla_q\log\pi(q\mid y)=
-\nabla_q\Phi(q)-\begin{bmatrix}z\\0_2\end{bmatrix}.
\]

The uniform-property prior has zero interior score; it does not contribute
an invented Gaussian property precision in original coordinates.

## 2. Whitened Jacobian consistency and packing

`core.py:156–167` uses a real packing operator \(P\) that orders each source's
32 receiver real parts followed by its 32 imaginary parts. It is a real
isometry from the complex Frobenius norm, independent of this ordering:

\[
r_w=\sqrt2\,P(r)/\sigma,\qquad
J_w=\sqrt2\,P(D_qY)/\sigma,\qquad
\Phi=\tfrac12\|r_w\|^2.
\]

Consequently the independently useful identity is

\[
\boxed{\nabla_q\Phi=J_w^T r_w.}
\]

The intermediate einsum has indices receiver/source/pixel; flattening its
first two axes and then reshaping each derivative column back to `(32,8)`
before applying `white` preserves the declared source-major packing. No
packing discrepancy was found by source inspection. The parent should still
pay a small set of substantive finite differences and adjoint/Jacobian
comparisons before relying on either implementation for reference production.

## 3. Exact ideal probit transport

Let \(x=(z,\eta)\in\mathbb R^{66}\), with

\[
T(x)=(z,u),\quad u_j=\Phi_N(\eta_j),\qquad
DT(x)=\operatorname{diag}(1_{64},\varphi_N(\eta_1),\varphi_N(\eta_2)).
\]

Since the original property density is uniform on \((0,1)^2\), change of
variables gives

\[
p_{0,x}(x)=(2\pi)^{-32}e^{-\|z\|^2/2}
\prod_{j=1}^2\varphi_N(\eta_j)
=(2\pi)^{-33}e^{-\|x\|^2/2}.
\]

Thus the ideal transformed posterior, relative to \(dx\), has

\[
\boxed{\log\pi_x(x\mid y)=-\Phi(T(x))-	frac12\|x\|^2+\mathrm{const}.}
\]

The gradient is

\[
\boxed{\nabla_x\log\pi_x(x\mid y)
=-DT(x)^T\nabla_q\Phi(T(x))-x.}
\]

`core.py:192` supplies the diagonal of \(DT\) correctly. Multiplying the
likelihood gradient componentwise by this diagonal and subtracting **all 66
coordinates** is the proper transformed log score. Once the full standard
Gaussian prior is included, adding another log-Jacobian term would double
count the property prior transport.

`Target.objective` at `core.py:158` is the **q-space** potential
\(\Phi+\|z\|^2/2\). It is correct for its original coordinates but is not the
transformed x-space target. A wrapper must explicitly use
`-target.phi(point,y) - .5*x@x` with `point.q = from_probit(x)`.
The 66-dimensional MALA interface in `kernels.py` allows an unbounded
transformed target and expects its **full log-target gradient**. Original
property reflection/slice kernels must not be applied to \(\eta\) coordinates.

The two MAP notions differ: MAP relative to \(dq\) and MAP relative to \(dx\)
need not map into one another. Probit is a probability-preserving inference
coordinate choice, not a claim that transformed mode equals the original MAP.

## 4. Initialization clipping and numerical tails

`to_probit` currently clips property values to `[1e-4,1-1e-4]`, yielding
\(\eta\) approximately in `[-3.71901648545568,3.71901648545571]`.
This is acceptable **only as a declared finite initialization rule**, including
for a q-space MAP exactly on a property boundary. It changes the initial
state; it must not clip proposals, constrain production support, transform
every current state, or silently define the target prior.

No eta clipping or target truncation appears in `from_probit` itself. Thus
initialization clipping does not, by itself, truncate the production target.
The existing `to_probit` is nevertheless **not a general inverse**:

```text
eta input:                            [ 5.0, -5.0 ]
to_probit(from_probit(x))[64:]:        [ 3.7190164854557084, -3.7190164854556804 ]
```

For an exact mathematical inverse on interior property states, use un-clipped
\(\Phi_N^{-1}(u)\); boundary values correspond to infinite eta, so they require
the separate finite initialization policy. Renaming/documenting the helper
as initialization-only would make this distinction reviewable.

Float64 `ndtr` itself has unavoidable tail limitations. The local runtime gave:

| eta | computed u |
|---:|---:|
| 8.0 | 0.9999999999999993 |
| 8.2 | 0.9999999999999999 |
| 8.3 | 1.0 |
| 9.0 | 1.0 |
| -37.0 | 5.725571222523923e-300 |
| -38.0 | 0.0 |

Therefore the implemented map ceases to be numerically one-to-one in extreme
tails, even without clipping. Keep the chain state in x-space; do not recover
eta from rounded u. Evaluate its Gaussian log prior directly in x-space,
rather than forming products of tiny Jacobian densities or `log(ndtr(...))`.
The analytic `probit_derivative` is not the derivative of the rounded stepwise
float64 CDF in these extreme tails, so a numerical smoothness claim there
would be false. A fixed exact-MH correction using the implemented target still
corrects proposal drift imperfections, provided forward/reverse proposal
densities use the drifts actually used to generate proposals. Tail rounding
remains a numerical approximation to the ideal transported target, not a new
authorized hard support or a reason to clip production eta.

## 5. Concrete API/cache bugs

1. **Observation-dependent gradient cache lacks observation identity.**
   `Point.gradient` is cached once at `core.py:170–176`, but it depends on y.
   Calling `phi_gradient(point,y1)` and then `phi_gradient(point,y2)` returns
   the first gradient unchanged, while `phi(point,y2)` computes a new
   likelihood. Points reused across observations/noise replicas can therefore
   generate inconsistent drift/target pairs. Scope the cached gradient to one
   observation, include observation identity, or enforce that a point cannot
   cross observations. Current single-observation ownership can avoid the
   failure but is not enforced by this API. No physical reproduction was run.

2. **Integer arrays silently truncate probit outputs.** Both `to_probit` and
   `from_probit` use `np.asarray(...).copy()` and retain integer dtype. The
   actual pure-helper check gives `from_probit(zeros(66,dtype=int))[64:] ==
   [0,0]`, whereas a floating input gives `[.5,.5]`. Both conversion helpers
   should normalize to floating dtype and validate the 66-coordinate shape.
   This is independent of any physics budget.

3. **Validate observation shape where the target accepts it.** `phi` and
   `phi_gradient` do not enforce `(32,8)` finite complex observations; NumPy
   can silently broadcast a `(32,1)` observation. The saved-observation
   metadata check currently validates keys/target/sigma, not array shape.
   If the input pipeline guarantees the frozen shape this is contained, but
   the current target interface does not enforce it.

## 6. Cost accounting and cache reuse

The source adjoint solves eight RHS against an **existing** state factorization.
The A30 adapter reserves `columns.shape[1]` before this solve, under role
`likelihood_gradient_adjoint`. It does not construct a fresh state or LU.
Therefore a fresh forward plus one source gradient costs **8 + 8 = 16 RHS**;
a gradient on a cached same-material state costs **8 incremental RHS**. An
already cached gradient costs zero RHS, subject to the observation identity
condition above.

The receiver transfer instead solves 32 adjoint RHS. A fresh forward plus
full receiver Jacobian costs **8 + 32 = 40 RHS**; a new full Jacobian on an
existing state costs **32 incremental RHS**. Transfer construction/cache
ownership is explicit in A30. Matrix multiplications and residual audits are
wall-time work, not new linear-system RHS.

One avoidable extra cost remains: when a point already has its full whitened
Jacobian, `phi_gradient` still performs eight source-adjoint RHS. The identity
\(\nabla_q\Phi=J_w^Tr_w\) can reuse that Jacobian at zero additional RHS.
This is an optimization opportunity, not undercounting or an algebraic bug.
It should not trigger building a 32-RHS Jacobian merely to save an 8-RHS
source gradient.

`physical(...)` records elapsed time; it is not itself a reservation. The
actual `reserve_rhs` call is inside the unchanged adapter, before the solve.
Do not add a second eight-RHS reservation at the wrapper or count the same
elapsed interval twice when summing campaign costs. Production MALA should
retain a rejected current state's gradient and state factorization so that
every rejection does not incur an unnecessary recomputation.

## Parent verification items

Before large reference work, the most decisive paid checks are a few real
directional finite differences for the likelihood, agreement of the source
gradient with \(J_w^Tr_w\), and chain-rule checks in moderate eta where the
float64 CDF remains resolved. Include imaginary-property perturbations to
detect misplaced conjugation and verify the factor of two from circular
complex noise. Check the ledger increments for these same operations.
These checks were deliberately **NOT RUN** by this read-only reviewer.
