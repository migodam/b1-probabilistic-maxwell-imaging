# A0 adjoint polish: code ready, execution NOT_RUN

Frozen root `A0_ADJOINT_POLISH_SUPPLEMENT.json`; runner `tools/a0_adjoint_polish.py --root ROOT [--output PATH] [--ledger PATH]`, default output `results/a0/adjoint_polish`. Source compiled in memory only; no imports, tests or optimization executed by worker. Parent supplies `physics.weighted_material_gradient` before execution.

Fixed parents0,3,6,9 × two repeats, known intrinsic contrast3/radius.12/loss.02. Select starts across all three saved main arms by recomputed online hierarchical objective. No GT/holdout/error-based choices. Support remains that chosen public support; the saved-main no-op competes with every completed paid online trial under the unchanged proper objective.

L-BFGS-B runs in public normalized coordinates `u=(theta−LOWER)/(UPPER−LOWER)`, bounds[0,1]^6; returned physical gradient is multiplied by widths[4.5,.25,2,1,.02,.02]. This is numerical preconditioning, not a changed prior, objective or scientific gate. Known-m primary theta remains eps,ell,px,py,gain_offset,phase. maxfun/maxiter/maxls20, with a separate hard guard of20 unique joint-theta combined evaluations per condition; optimizer overrun stops before further PDE. Exact combined-cache repeats are unpaid, and current FE material state is reused for gamma-only changes.

Each new point uses forward2 RHS and parent weighted-adjoint2 RHS. Function contract: `weighted_material_gradient(ev,world,raw_residual,gamma) -> (gradient4, metadata)`, with raw complex residual `gamma*F−y` **unwhitened**. Gain and phase use analytic gradients with unchanged sigma.005 prior.

Only parent0/repeat0 first evaluation upgrades the same state to native full shared J (16 extra RHS). Compare weighted physical4 to `Re((gamma*J4)^H residual)/.01²`; relative scaled max threshold1e−9 using public physical scales[4.5,.25,2,1]. Persist both before comparing. Audit mismatch raises immediately, saves all errors, and prevents offline access; no adjustment/jitter to rescue it.

Headroom must be≥700RHS/180FE in A0, positive remaining CPU, no borrowing. Upper bound676RHS170FE:640/160 fitting,16/0 audit,16/8 independent outer holdout,4/2 fixed finer parent0 both repeats. Every action reserves through supplied Book/FEM; no separate manual physical charge.

Per-case artifacts: `INITIALIZATION.json`, all `ATTEMPTS.json`, `EVAL_XX.npz` with theta/F/prediction/raw residual/physical+normalized gradients/costs/metadata, `NATIVE_FULL_J_AUDIT.npz` and audit JSON only at fixed audit condition, current best/no-op, `ONLINE_RESULT.json`. All8 online selections become terminal and durable in root `ONLINE_COMPLETE.json` before first offline label or holdout loading. Fatal online errors leave remaining fixed conditions explicitly NOT_RUN and no offline scoring.

Independent outer holdout is **previously exposed development re-evaluation**, not a new blind test. Same frozen64-real noise envelope, pose/material metrics, fixed parent0 finer n384/grid3; primary h.17/n256/grid2. Offline scores separate and all nonconvergence retained. Main A0 files are never modified and no new Go gate appears. Starts depend on all three paid main arms plus common paid proposal bank, so this is ordinary optimizer-error diagnosis, **not a matched-cost speed win**.
