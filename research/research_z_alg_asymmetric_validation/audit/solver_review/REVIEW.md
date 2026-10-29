# Independent solver implementation review

Review scope: `EXPERIMENT_CONTRACT.json`, `src/accounting.py`, `src/physics.py`, `src/charts.py`, and `src/probe.py`. No real `Model` or `Book` was constructed by this reviewer. No Hankel kernel, campaign geometry, physical forward solve, or cost-ledger mutation was executed. The only numerical work uses deterministic random dense matrices, an in-memory fake ledger, and unbound implementation methods. The parent retains scientific interpretation and final acceptance.

## Core solver assessment

**PASS on static algebra and fake-matrix implementation checks. All prephysics implementation notes below are closed by static re-review.** This is not a physical-model validation receipt.

- For the declared outgoing `exp(-i omega t)` convention, the point incident Green is `i H0(kr)/4`. Integrating `k^2 i H0(k|r-r'|)/4` over a circular pulse gives the external coefficient `i*pi*k*a*J1(k*a)*H0(k*r)/2`. The centered self term is `i*pi*k*a*H1(k*a)/2 - 1`; the subtraction is the lower integration endpoint, not an omitted self interaction. These match `Model.external` and the diagonal assignment. The actual high precision/physical Hankel comparison remains the parent's pending work.
- The current equation is `L j = diag(m) E0`, with `L = I - diag(m) G` and `E = E0 + G j`. Its material derivative is `delta j = L^-1 diag(E_s) delta m`. The code forms `T = S L^-1` through `L.T` solves and then `J[s,r,i] = T[r,i] E[i,s]`. `trans=1` is correct; a conjugate transpose would be wrong.
- Data packing is source major, then receiver, then real/imaginary component. A two-source/two-receiver labeling check passes exactly. The material Jacobian uses the corresponding row order.
- `CurrentChart` retains the complete receiver kernel `N`; it parameterizes `j_s = j0_s + scale*N*c_s`, with one shared real material vector. The constitutive residual `j - m E` uses all complex source equations and no division by a potentially vanishing field. Its current blocks are the realification of `scale*(N-diag(m)GN)`, and its material columns are `-E`. Static inspection and fake central differences agree.
- Current coordinates contain explicit real and imaginary current components followed by real dependent material coordinates. Free coordinates fill the remaining real material entries. No independent complex material per source is introduced.
- The current implicit material tangent agrees with the direct material graph tangent, including signs and scaling, on a fake regular graph.
- Each new `Model.forward` clears the previous instance-level LU and charges `ns` complex Maxwell RHS and one LU before execution. Independent validation uses a fresh `numpy.linalg.solve` and leaves no persistent factor. Material Jacobians charge `nr` transpose RHS against the live LU; cached Jacobians correctly trigger no new solve. Failed operations retain their pre-action RHS charge.

## Numerical fake-matrix evidence

`ALGEBRAIC_QA.json` reports material-J central-difference relative errors of `5.23e-12`, `4.08e-11`, and `4.78e-10` at steps `1e-4`, `1e-5`, and `1e-6`. The explicit inverse check differs by `1.63e-16`; a deliberately wrong Hermitian construction differs by `0.209`. The fake ledger records two transpose RHS for two receivers.

`CURRENT_ALGEBRAIC_QA.json` reports current-J relative errors of `1.33e-12`, `1.48e-11`, and `1.33e-10`. The center constitutive residual is `1.59e-15`; receiver compatibility is `2.18e-15`; the current/direct material-tangent relative gap is `1.64e-15`. The fake current system has 24 real variables and 24 real constitutive equations, with two free material coordinates. Both gates pass.

Reproduction: run `audit/solver_review/algebraic_qa.py` with the existing Gaussian `.venv_nn` Python. The system Python's missing NumPy was a tooling failure only; it executed no model work.

## Findings communicated to parent

1. **CPU continuation guard — FIXED:** `Book` now rejects any previous `RUN_START` using the same run name before a new append/action. Both interrupted and completed runs are single-use, preventing silent CPU resets. Pre-action RHS payment remains intact.
2. **Ordinary algebra overhead — FIXED:** `solve_real` records the actual number of RHS columns. A failed direct solve now records a separate `FALLBACK_SVD` decomposition and RHS. `CURRENT_CHART_SETUP` records one algebraic complex factorization, `ns` algebraic complex RHS, and receiver matrix dimension. These remain distinct from Maxwell resolvent RHS.
3. **Probe gates — FIXED:** probe PASS now requires current center residual at most `2e-13` and a full-rank square current Jacobian at relative `1e-12`. Material finite differences now cover every coordinate at every declared step. The current directional checks are independent derivative diagnostics, not a step-convergence estimate along one common direction.
4. **Endpoint acceptance remains a caller obligation:** `CurrentChart.correct(success=True)` means that the constitutive residual converged. It does not itself enforce full box legality, node-field threshold, resolvent floor, response qualification, or independent physical validation. Callers must retain the frozen endpoint gates and reject illegal roots. `MaterialChart.correct` likewise returns convergence rather than full endpoint qualification. Current `quality` implements response/state/resolvent/box checks; node/rank gates require explicit caller handling.

## Limits

The review has not executed the physical probe, independent endpoint solves, mpmath comparisons, interval enclosures, chart-switch coverage, atlas continuation, posterior quadrature, or efficiency comparisons. The single-live-LU behavior is verified per `Model` instance; the top-level physical-worker code must avoid keeping two instances' factorizations live simultaneously. The accounting checks are static: actual campaign totals must be reconciled against the ledger after parent execution.

Refresh receipt: corrected accounting/charts/probe source read once after parent changes. No new numerical tests were added or run; no physical Model/Book was instantiated.
