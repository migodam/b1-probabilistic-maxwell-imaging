# Operator-driver implementation handoff

Status: implementation delivered; physical driver `main()` **NOT_RUN** by this worker. No `results/g0_operators` directory was created. This handoff is not a whole-G0 PASS or a scientific acceptance judgment.

## Interface and fixed scope

- Entry point: `tools/g0_operators.py`. It requires the existing `tools/run_job.py G0 <unique-job-label> <python> tools/g0_operators.py` wrapper environment and all six single-thread settings. The wrapper owns aggregate process CPU exactly once.
- Reads the current `AUTHORIZATION.json`, `CONTRACT.json`, `GATES.json`, `G0_CASES.json` and prior `results/g0_qualification/RESULTS.json`. The first seven case IDs/order, full-2x2 background with epsinf=3.75/ell=.125, 2257 public ports, and strict thresholds are checked. Existing qualified output arrays never select the ROM basis.
- One public fixed P3 mesh, one `g0_qualify.background_system` assembly, one `FixedSpaceCache` background sparse LU/2 RHS, and one complete `build_ports()` setup/2257 RHS. No separate `schur()` call or exterior factorization/2259 RHS.
- Exactly seven `cache.full()` baseline target solves, two RHS each. `weak_center` and `left_small` run first to supply the predeclared training fields; their saved full baselines are reused in the comparison instead of paying twice.
- Calls the existing `woodbury`, `schur_from_paid_ports`, `krylov`, and `goal_galerkin` interfaces. Their actual sparse LU/RHS costs go through `BudgetBook`; Krylov charges every background-preconditioner application. Dense factorization dimensions remain in the ledger. Dense estimator/probe solves additionally record dense RHS and transpose, without claiming new physical PDE RHS.
- ROM basis: fixed-order two-pass modified Gram-Schmidt QR of 2 background primal + 4 training full + 16 background receiver-adjoint columns. The adjoints use the same paid factor with `factor.solve(R.conj().T.toarray(), trans="H")`, after a 16-RHS reservation. Fixed 1e-12 relative algebraic deflation only, normalized columns, no sweep, no held-out field selection. The five held-out cases stay held out even when ROM accuracy is poor.
- Stability: only predeclared `low_loss_high_contrast` gets one additional 2257-dimensional `I-CV` dense LU, fixed-seed `onenormest(t=2,itmax=5)`, and the actual resolvent-RHS LU-solve backward residual. No SVD, jitter, pseudoinverse, inferred inf-sup lower bound, or resonance claim.
- Full/woodbury/paid-Schur/Krylov strict field and response relative identities above 1e-9, or equation/core backward residual above 1e-10, save artifacts and stop. ROM full-equation residual and differences are accuracy measurements; only a broken reduced solve triggers strict failure. There is no retry/fallback or basis change.

## Output contract

`results/g0_operators` is exclusively created with `exist_ok=False`; an old run is never replaced. Each case/method has `.npz` with response, scattered internal field, equation residual, and field/response differences to full, plus `.json` with solver CPU/wall, verification/I/O costs, residuals, errors, strict status and ROM training/held-out status. The directory also records the frozen plan, memory estimate, setup timings, background and target operators/RHS, paid adjoints, QR provenance/basis, one core diagnostic, accumulated results and status. Exceptions save `STOP.json` with context and available failed-kernel timings. Partial rows remain inspectable.

Known fixed physical cost before Krylov: **2289 RHS** = background 2 + public ports 2257 + target full 14 + adjoints 16; **8 sparse LU** = background 1 + target full 7. Dense algebra is separately dimensioned. Prior G0 costs remain cumulative and are not reset. A final ledger snapshot precedes the wrapper completion charge and explicitly says to use the wrapper receipt for aggregate CPU.

## Verification performed; cost gap preserved

1. Static Python parse/compile of the driver passed using system `python3`; no execution of driver `main()`.
2. An initial direct, pure-synthetic helper attempt under system `python3` failed before numerical work because `numpy` was unavailable. No dependency was installed.
3. One direct zero-PDE helper run completed with the existing A25 interpreter. Command, run from the campaign root:

   ```sh
   PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="src:vendor:tools:../../../research-w-continuation/.deps" "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A25/a25-adaptive-om-feedback/.venv_a25/bin/python" research/delegated/operator_driver/check_helpers.py
   ```

   Evidence: `research/delegated/operator_driver/ZERO_PDE_CHECKS.json`; synthetic retained-rank/orthogonality, inverse-norm estimator on a known diagonal core, strict-failure artifact retention, inaccurate-ROM retention, manifest and memory-guard checks passed. QR Gram error was 3.6025491547985136e-16; synthetic inverse-norm relative error was 2.220446049250313e-16. Synthetic dimensions are not the real 2257-port problem.

   **Accounting gap:** this helper run did not use `run_job.py` and the helper did not actively measure process CPU/wall. It wrote no `BudgetBook` records. Exact helper CPU is **NOT_RECORDED** and must not be reconstructed or entered as zero. The tool reported about 0.5 s elapsed for that shell call; that is not a wrapper CPU receipt. Both helper attempts performed zero FEM instances, mesh assemblies, sparse LU, or physical RHS; `main()` stayed uncalled. No further helper run followed the parent's instruction to centralize checks through the wrapper.

The final failure-timing metadata patch was not numerically rerun; parent verification remains pending.

## Remaining limits and review points

- **NOT_RUN:** real cache setup, seven-condition operator comparison, real core estimate and real ROM accuracy. No claim about speed, numerical stability, or generalization has been made.
- The original 1.5 GB guard is retained and strengthened with a named-array comparison estimate (about 1.462 GB for ~28420 free DOFs/2257 ports). It excludes FEM/Basis storage, sparse LU storage and opaque library workspace; it is not a process-RSS guarantee. Actual background sparse LU storage is separately checked against the frozen cap. Peak process RSS comes from the wrapper.
- Exact interface calls do not expose their dense factors; the one stability diagnostic deliberately pays and records an additional dense LU instead of pretending it is free. Method solver timings and verification/I/O timings are separate.
- Whole G0 remains partial: finite perturbation/Jacobian, independent information/posterior stress, and parent scientific review are independent requirements. The script never writes a whole-G0 PASS.
- Only `tools/g0_operators.py` and this delegated output directory were written by this worker.
