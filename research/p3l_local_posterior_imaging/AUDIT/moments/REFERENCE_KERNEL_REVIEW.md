# Full-reference runner: independent read-only kernel review

Verdict: **NO_FATAL_DEFAULT_MALA_KERNEL_DEFECT_FOUND_BY_STATIC_REVIEW**.
The default probit MALA runner targets the intended full posterior, retains
rejected states and their current gradients, and freezes production tuning.
This verdict is not a mixing, finite-sample reference-quality, or physical
solver certificate. The parent Codex owns the final mathematical decision.

Reviewed: current `src/p3l/reference.py`, `src/p3l/core.py`, and the PCN,
property-reflection, elliptical-slice and MALA definitions in
`src/p3l/kernels.py`. No runner/core/ledger/physics edits or physical calls were
made. Maxwell solves/RHS during this review: **0**. Earlier physical contract
figures supplied by the parent (source-adjoint/J relative difference
`4.08e-16`, finite-difference discrepancy `1.42e-10`, five paid states) are
context, not independently rerun evidence from this reviewer.

## 1. Full probit target and gradient: consistent

For \(x=(z,\eta)\in\mathbb R^{66}\), \(u_j=\Phi_N(\eta_j)\),
the original Gaussian64 × Uniform2 prior pushes back to Gaussian66:

\[
\pi_x(x\mid y)\propto
\exp\{-\Phi(z,\Phi_N(\eta))-\tfrac12\|x\|^2\}.
\]

`Posterior.__call__` (`reference.py:22–26`) uses exactly this density,
relative to dx, when `probit=True`. It does not reuse the q-space
`Target.objective`, omit the eta prior, or multiply the Jacobian twice.
There is no local-ball restriction, likelihood cutoff, or production clipping.

The supplied gradient (`reference.py:27–30`) is

\[
g(x)=\nabla_x\log\pi_x(x\mid y)
=-\operatorname{diag}(1_{64},\varphi_N(\eta))
\nabla_q\Phi(T(x))-x,
\]

which is the correct chain rule. The core now distinguishes cached likelihood
gradients by observation array and uses a copied observation, fixing the prior
cross-observation cache issue. The conversion helpers now produce floating
arrays, and `to_probit(..., initialization=True)` explicitly names the inward
initialization rule. Only initialization uses it in this runner.

Finite CDF tails still round to u=0 or u=1, as discussed in the earlier core
review. The runner correctly retains x as the chain state and does not recover
eta from rounded u. This remains a floating-point limitation, not an imposed
eta support restriction.

## 2. MALA detailed balance with frozen tuning

For fixed positive h and SPD C, the implemented proposal is

\[
x'=x+\tfrac h2 Cg(x)+\sqrt h L\xi,
\quad LL^T=C,\quad\xi\sim N(0,I_{66}).
\]

Let \(d_f=\tfrac h2 Cg(x)\), \(d_r=\tfrac h2 Cg(x')\),
and \(\Delta=x'-x\). Then

\[
\log Q(x\mid x')-\log Q(x'\mid x)
=\frac{\|L^{-1}(\Delta-d_f)\|^2
-\|L^{-1}(-\Delta-d_r)\|^2}{2h}.
\]

This is the sign and scaling in `PreconditionedMALA.step`. Adding
`logpi(x')-logpi(x)` gives the correct MH ratio. Since C is the same at both
ends of a production transition, its determinant cancels. No missing
state-dependent metric determinant or derivative is needed in this **fixed**
proposal. It would become wrong to refresh C between a proposal and its
reverse-density calculation; this runner does not do that.

`metric` (`reference.py:33–36`) constructs
\(C=(J_x^TJ_x+I_{66})^{-1}\), symmetrized before Cholesky. In ideal arithmetic
this is SPD. The identity prior term is appropriate in probit coordinates.
The metric is an efficiency choice; it does not replace the full target.

One initialization approximation should be named precisely: `center` can move
boundary attributes inward, while the first J is evaluated at the original
anchor, then multiplied by the derivative at the moved center. Thus the first
matrix need not be the exact GN metric of `center`. It remains a legal SPD
preconditioner, so this mismatch does **not** invalidate MH.

## 3. Warmup and production boundary: frozen correctly

The warmup metric updates occur only at i=256,512,768 and only when
`i < warmup` (`reference.py:62–64`). Step-size updates occur only after
warmup transitions with `i < warmup` (`reference.py:78–83`).

At i=warmup the MALA kernel is constructed once with the last warmup h and C.
For every later iteration, `i <= warmup` is false, so the same kernel instance
is reused. No production state-dependent tuning is present. For warmup=0,
the initial kernel is constructed and remains fixed, which is also consistent.

Every saved array includes warmup states, and the initial state before the
first transition is not saved. Therefore production samples are precisely
`q[warmup:]` / `x[warmup:]`; discard the first `warmup` **rows**, not one row
per accepted move or one row per physical target evaluation.

Freezing a metric learned from the warmup does not make the first production
state stationary. Warmup draws are adaptive and must not be relabeled exact
posterior draws. Mixing diagnostics, cross-chain agreement and reference
qualification remain necessary, particularly under a fixed RHS cap.

## 4. Rejections, current gradient and physical-point pinning

MALA initializes `ChainState` with the current full x-space gradient. The
kernel stores a proposed gradient on acceptance and reconstructs the rejected
current state with its original/current gradient on rejection. It also retains
that gradient on the -infinity-target rejection branch. Thus rejection does
not replace a current gradient with a proposed gradient or erase it, and
repeated rejections do not automatically pay another current adjoint solve.

The reference loop (`reference.py:76–86`) records **every completed transition**,
including repeated q, x, potential, log density and endpoint key when rejected.
It neither deduplicates states nor likelihood-softmaxes them. Rejected proposals
are not inserted into the posterior state sequence.

`Posterior.pin(current.q)` protects the current x byte key. During proposal
evaluation, eviction excludes both that old current key and the newly evaluated
proposal key. The proposal therefore survives through its gradient calculation;
the old current point survives a rejection. After acceptance the new point is
pinned. Checkpointing reuses the pinned current point. The five-point limit is
sufficient for this MALA current/proposal lifecycle. The dictionary is FIFO
with protected entries, rather than strict LRU, but that does not change the
target or acceptance ratio.

The core Point cache stores physical state quantities at q=T(x), while the
reference state/log density and cache key use x. This distinction is correct:
even if two far-tail x values round to the same physical u, they retain
different Gaussian prior penalties. Their different x keys can cause redundant
physical calls but cannot collapse their probability weights.

The proof assumes the observation is fixed during a reference run and target
evaluations are deterministic to the stated solver tolerance. The runner
maintains that ownership pattern. The broader core still allows mutating arrays;
callers must not mutate an active observation or stored state behind the runner.

## 5. Other production kernels: invariant, with an ESS reporting caveat

PCN proposes a prior-reversible field update and an independent fixed-scale
coordinatewise reflected property increment. Its full-target MH ratio includes
the Gaussian proposal/prior correction, leaving the expected likelihood ratio.
The final beta/property scales freeze after warmup, despite reconstructing the
PCN object each step. No probability issue was found.

The ESS branch uses an uncapped field elliptical slice at fixed properties,
followed by symmetric property-reflection MH. Both conditional kernels preserve
the same full target; their ordered composition therefore preserves it. The
whole sweep is **not generally reversible** unless the component kernels
commute. Invariance is sufficient, and the runner's status claims invariance,
not whole-sweep detailed balance. Internal rejected ESS bracket angles are
proposal evaluations, not extra Markov states, so excluding them from q rows
is correct.

**Concrete nonfatal reporting issue:** `tr` is the property transition. Hence
the saved `accepted` flag and MH/proposal log ratios describe the property
substep only. `accepted=False` can coincide with a changed full q because the
preceding field ESS move was accepted. The `repeated_state` value, when inspected
in memory, also refers to the property substep. Do not interpret the ESS saved
acceptance fraction as a whole-sweep movement fraction or use it to reconstruct
rejection repeats. The saved q rows themselves are correct.

`field_rejected_proposals` is added to the in-memory metadata but is not saved
in CHAIN.npz or the final summary. Total target-evaluation counts do include
field trials plus the property call. Thus physical accounting remains in the
ledger, but a detailed ESS bracket-efficiency report cannot be recovered from
this saved metadata alone.

## 6. Scope, qualification and cost limitations

- `total_target_evaluations` in the transition summary counts kernel target
  invocations and excludes initialization and geometry setup. A target-cache
  hit still counts as an invocation but has no new RHS. It must never be used
  as a substitute for actual ledger RHS or total runtime.
- Default MALA initialization uses a full GN-preconditioned random perturbation,
  with four scale choices. This is more dispersed than one repeated point but
  is not a certificate that all relevant full-posterior basins are covered.
  The PCN/ESS branch starts all properties at their anchor values and only
  perturbs z by .03; it must not be described as proven overdispersed across
  bounded properties. This affects reference quality, not invariance.
- Each warmup Jacobian point can remain referenced by the A30 receiver-transfer
  cache even after eviction from the five-point Posterior cache. A single
  run has only the specified few geometry points; reusing one Target/backend
  across many runs could accumulate those dense state caches. The five-point
  limit does not bound the backend's complete memory footprint.
- CHAIN.npz checkpoints are partial traces, not evidence of completion. If a
  budget stop occurs, analyze only completed production rows, retain the
  stopping/failure status, and do not silently claim requested production
  length or reference qualification. A production stop selected by current
  diagnostics can also change finite-sample estimator behavior.
- No general resume path was reviewed or certified. Saving RNG/current-state
  files does not by itself authorize replay of paid work or prove that a
  restarted job uses the identical frozen kernel.

## Final review judgment

The default MALA implementation has the expected target, proposal correction,
production freeze, state repetitions, current-gradient retention and point
pinning. No fatal sampling issue was identified in those mechanisms. The
reference label must still be conditioned on actual material-function ESS,
MCSE, split-Rhat, cross-mechanism comparisons, complete costs and any budget
termination. ESS output acceptance metadata needs the narrower interpretation
above. The parent should not infer that this source review or the earlier
physical derivative checks certify full-posterior mixing for sixteen objects.
