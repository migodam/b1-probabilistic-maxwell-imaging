# Pure finite probability API

Implementation only, tests NOT_RUN by worker. No physics, forward calls, network, or installation.

`build_grid()` defines **FINITE_PRIOR_ENUM4_V1**, an explicit discrete prior. It is **not** a continuous Uniform integration certificate. There are 14 legal 2×2 masks with normalized Bernoulli(.65) mass, CC3 material nodes and positive masses (1/6,4/6,1/6), and CC5 pose nodes/masses (1/30,4/15,2/5,4/15,1/30). Empty material is physically collapsed into one wall state, with unused eps=3.75, ell=.125 placeholders. Result: 118 wall groups × 25 poses = 2950 states. Pose is the outer loop; each consecutive 118-row block shares its pose.

Returned arrays: `support_id` (N,), `support` (N,4), `epsinf` (N,), `ell` (N,), `pose` (N,2), `prior_weight` (N,), `wall_group_id` (N,), `legal_supports` (14,2,2). Support columns are row-major 2×2 regions. Also returns labels, support-prior masses, node metadata, counts, and certificate flags. Wall group encodes occupancy plus material, independent of pose.

`build_coarse_grid()` yields 477 rows (53 walls × 9 poses) from specified two material nodes and 3×3 target bank. Label **PROPOSAL_BANK_ONLY**, `posterior_claims_allowed=False`. Its normalized bank weights do not establish posterior claims.

`update_logweights(prior, uniqueLLlist)` returns normalized log posterior. Prefer `{observation_id: loglikelihood_array}` or `(id,array)` pairs; duplicate IDs fail. Bare arrays are supported only with caller-guaranteed unique observations. Each likelihood is added once. NaN/+∞ fail; impossible states may use −∞. This routine is independent of solver and factor registry.

`summarize(weights, grid, quantiles=(.025,.5,.975))` returns occupancy probabilities, each occupied-region conditional joint/marginal material distribution (undefined when occupancy mass is zero), physical wall means/covariances, joint 10×10 covariance, pose moments, 14 support-mode probabilities, event p_x≥0, and coordinate mass/CDF/quantiles. Physical coordinates are `[z0*(eps−1),...,z3*(eps−1), z0*ell,...,z3*ell, px,py]`. Empty placeholders never enter material conditionals or physical wall variables. Distributions use arrays `values,mass,cdf,quantile_levels,quantiles`.

`conditional_target_only(oldw,newLL,wall_group)` returns `{label,weights,exact_new_target_posterior,interpretation}`. Label **STALE_WALL_MARGINAL**: old wall masses remain exactly fixed while each pose conditional updates by the likelihood. This is an approximation to the new joint posterior and is not labeled the exact new target posterior. An occupied old wall group with zero new likelihood everywhere fails because its new conditional is undefined.

`posterior_metrics(approx,ref,grid)` returns each coordinate CDF supremum error, p_x≥0 event error, maximum support-mode mass error among reference modes ≥.05, occupancy errors, target-pose Euclidean squared-loss actions, reference Bayes risk, approximate action risk evaluated under reference, regret, and approximate Bayes risk. There is no global posterior certificate. Primary action and all primary risks/regret use only the two target-pose coordinates. Wall-material uncertainty remains in the separate covariance and discrepancy diagnostics; it does not enter the action loss.

`importance_sample(grid,log_likelihood,n_particles=1000,seed=0)` draws the finite prior and returns self-normalized importance particle/state weights. `tempered_smc(...,temperatures=(0,.25,.5,.75,1),mh_steps=2,resample_ess_fraction=.5)` applies incremental likelihood powers, multinomial resampling when ESS falls below threshold, and exact-invariant prior-proposal independence MH at each temperature. MH log acceptance is beta*(LL_proposed−LL_current); prior proposal/target terms cancel. Resampling is representation change and never new data. Both return reproducible seed, particle indices/weights, full grid weights, and `global_posterior_certificate=False`. Sampling currently requires finite likelihood logs. SMC returns stage ESS/resampling/acceptance history and `observation_likelihood_insertions=1`.

Tests use standard unittest, small analytic four-state arrays for updates/moments/metrics, finite-prior mass/count checks, seeded sampling against an analytic four-state posterior, and duplicate-observation/invalid-schedule rejection. No forward model is imported.

## Paid bank runner (implemented, NOT_RUN by worker)

`tools/finite_bank.py --output PATH [--resume] [--ledger PATH]` is designed exclusively for the parent-serialized E2E3 wrapper. It calls qualified `Evaluator(h=.17,n=256,order=3,grid_refinement=2)` at f=1 with four declared TX and 16 RX. Planned cost 2950 LU / 11800 RHS, no derivatives or truth solves. All physical reservations are delegated to the supplied shared Book via FEM. Known target defaults radius=.12, contrast=3, loss=.02; gamma=1. Finite pilot supplement, when present, is embedded in the resume config.

Outputs under PATH:

- `online/FINITE_GRID.npz`: complete public candidate coordinates/weights and discrete-prior labels.
- `online/FINITE_BANK.npz`: `data(completed,4,16), residuals(completed,), completed,total,tx,rx,frequency,prior_id`. Aggregate checkpoint every 25 rows and on caught failure.
- `online/paid_rows/0000.npz` etc: each paid completed row is atomically persisted, including `row_index,data,residual,receipt_json`. These define the exact recoverable prefix. Missing uncertain/failed paid rows are never retried automatically; reservations are never refunded. A partial durable final block recovers without duplicate solver calls.
- `online/COST_RECEIPT.json`, root `CHECKPOINT.json`, and completion-only `COMPLETE.json` preserve counts/status. Existing output requires `--resume`; changed config or holes fail closed.
- `online/representative_fields/ROW_0000.npz`, `ROW_1475.npz`, `ROW_2949.npz`: only these three actual candidate states may save full material fields/mesh. Interrupted missing field exports are not regenerated with expensive solves.
- `online/OBSERVATIONS.npz`: raw `data(12,2,4,16)`, `base(12,2,2,16)`, `outer(12,2,2,16)`, IDs and known noise/target/gamma. `sequential_factor_id(12,2,3)` labels disjoint t1 TX0/RX8, t2 TX1/RX8, t3 both base TX with remaining 15 RX. Index arrays `sequence_t{1,2,3}_{tx,rx}` make channel selection explicit. Full-base and sequential factor IDs are alternative representations and cannot both be inserted.
- `online/FACTOR_IDS.json`: raw block IDs, disjoint sequential IDs/channel membership, and duplication warning.
- `offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz`: `truth_index(12,),independent_noise(12,2,4,16),clean_means(12,4,16),parent_id,seeds,known_gamma,prior_id`. Parents are sampled after bank completion from its explicit finite prior with seed2026102504, means reused by index; all noise is independent real/imag N(0,.01²), seed2026102505.

No scientific acceptance gate is executed or inferred. Completion is only bank construction; the finite prior is not a continuous A0 quadrature certificate.

## Final measured bank cost

`tools/finalize_bank_cost.py --root ROOT [--bank PATH]` is a read-only, no-physics finalizer, to run only after complete bank and all `finite-bank-*` wrappers terminate. Use a wrapper label **outside** the `finite-bank-` prefix (for example `bank-cost-final-001`) because the prefix denotes physical bank jobs. Every matching start must have one authoritative completion and matching durable receipt; unknown/nonterminal bank job directories fail closed. The finalizer unions their ledger intervals, counts reservation fields once, requires actual11800RHS/2950FE, sums only measured completion CPU/wall, and retains all raw RSS units/receipts including failed launch001.

Writes only `online/JOB_COST_FINAL.json` (`complete=true`, explicit prior ID, physical planned/actual counts, measured CPU, summed walls, raw failed/successful receipts, maxima by RSS unit). It preserves raw in-child `COST_RECEIPT.json`, `CHECKPOINT.json` and `COMPLETE.json`: those cannot include their in-flight wrapper CPU. There is no manual extra CPU charge; the finalizer's own wrapper is separate E2E3 cost, while total campaign accounting remains authoritative in Book. E5 now **requires** this completed final receipt before creating its output and attributes its measured setup CPU rather than the incomplete in-child snapshot. Source compilation only; finalizer NOT_RUN by worker.
