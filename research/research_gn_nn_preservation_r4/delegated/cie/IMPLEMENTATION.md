# CPU CIE reference adapter

Status: **IMPLEMENTED / NUMERICAL EXECUTION NOT_RUN**. Syntax, import and CLI help passed; the prepared default plan contains exactly the frozen 13 mothers × clean/noisy = 26 public observations. Neither VAL nor TEST numerical CIE optimization has run. Runtime qualification and any scientific baseline conclusion remain with the parent. No new TEST truth was opened.

The old CPU source was absent from the worktree but present at commit `4baa5b0f88db1a49e9880c8e159e8ba336cc7da1`. `vendor/cie_tsom_baseline.py` is copied verbatim with `git show`; no frozen source was written. `vendor/PROVENANCE.json` records the exact source path and the existing baseline/source/E3 audits read for this adapter. The directly reused function is `run_cie`.

This is the already audited 2D real-material scalar-TM **explicit dense Fourier** adaptation. It is not an FFT complexity result, a new physics backend, or a literal 3D complex-material paper reproduction. Native projected Fourier coefficient coordinates, sourcewise PR+, non-descent reset, exact quadratic line step, zero-padding continuation and real modified-contrast readout are unchanged.

The schedule is frozen: mf3/beta1 → mf6/beta.2 → mf10/beta.2; maxiter500 per stage, miniter50, coefficient tolerance1e-8, DPCS rtol1e-3. Dual CIE uses beta_current1 in the two later stages. The source retains source-normalized current/data terms. Stage histories record actual mf, coefficient counts, APCS ranks, iteration stops and WALL_BUDGET/MAXITER; incomplete current schedules are not relabelled converged or weakened to fit the clock.

Each case starts cold at the public uniform m=.1, with all 576 real material variables free. It does not consume a GN anchor. All public operator/observation loading and Model construction are charged to the case's 55-second envelope; Python import/process startup is excluded from the per-case timer and included in queue elapsed. One retained frozen R4 `Model` supplies G/E0/S, full nonlinear forward and analytic Jacobian/counts. Prospective LU/RHS checks cap the whole current-plus-corrector pipeline at 600 LU and 36000 actual complex RHS. The current phase receives 55% of the cold envelope; the remaining time is available to full-material original-data bounded analytic TRF. Atomic operations can cross the wall polling boundary; actual overrun and completed operations are retained.

The spatial prior is explicitly **zero**, an adaptation to match NN data-only correction. The original-data objective is .5||F(m)−Y||²/||Y||², with box [.02,20] and clean/noisy physical gates .001/.003. TRF uses the same exact analytic physical J, x_scale=jac, exact trust-region subproblem and tight1e-12 ftol/xtol/gtol; its wrapper uses the remaining55s budget rather than the older R2 helper's hard-coded45s phase endpoint. It stops on first verified physical feasibility. The final selector uses the lowest original-data objective among physically feasible finite candidates and otherwise reports FAIL; no GT, support, label or image metric is used. `raw_1` transports the current-phase physical incumbent as a disclosed duplicate candidate slot. Source and full-physics costs are separately retained.

The resulting JSON/NPZ use the existing evaluation candidate/top1 schema. For this cold baseline the saved `anchor` is public uniform .1, so its before→after energy reference differs from the NN's common GN anchor. Compare selected absolute endpoint errors or explicitly retain this different initial-reference scope. The shared NN pipeline is GN45s plus post10s; nominal55s equality does not imply matched actual LU/RHS, setup or exploration.

Default CLI invocation only prints NOT_RUN:

```sh
python research/research_gn_nn_preservation_r4/experiments/run_cie_reference.py
```

After the parent freezes TEST methods and coordinates CPU resources, the concrete subset command is:

```sh
python research/research_gn_nn_preservation_r4/experiments/run_cie_reference.py --split test --scope figure13 --methods cie_dual --workers 1 --execute --test-methods-frozen
```

`--workers 2` permits at most two BLAS1 processes when the parent authorizes resources. No remote i9 queue is used. `--methods cie_single cie_dual` adds the separate single arm only when requested. `--scope all108` is reserved for the parent after core NN completion and only if time allows. TEST observations cannot be chosen ad hoc with `--obs-id`.

A cost-only public VAL observation can be executed separately after resource coordination:

```sh
python research/research_gn_nn_preservation_r4/experiments/run_cie_reference.py --split val --obs-id PUBLIC_VAL_OBSERVATION_ID --methods cie_dual --workers 1 --execute
```

The CLI never reads a truth path. Existing result files are preserved. Each new finite queue writes its complete plan, actual per-case receipts and queue summary. Public numerical failures retain HOLD/FAILED/BUDGET_STOP; no schedule or acquisition changes are used to rescue them. An interrupted source exception may prevent its internal nonphysical counters from returning; those cases are explicitly HOLD_INCOMPLETE_CURRENT_INTERNAL_COUNTERS while actual frozen-Model physical counts remain recorded.
