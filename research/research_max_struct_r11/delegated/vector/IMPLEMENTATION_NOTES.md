# R11 true-vector readout implementation

Status: IMPLEMENTED; one public saved-state derivative check PASS; new inversion NOT_RUN pending root-approved numerical contract. No GT generation/access. R10 files are read-only.

## Physical and acquisition contract

Reuse R10 VectorGrid and DenseResidentVectorModel unchanged: 12³ real isotropic materials m∈[0,8], 3N=5184 interleaved complex vector-current components per source, same true dyadic CM+RR DDA (width .6λ, k=2π, s=-1/3+i k³V/(6π)). All60 actual external electric dipoles remain RHS of one retained physical Torch LU; sameLU adjoint. Same 60 receiver positions with three components (180×60 data). No scalar-TM substitution, stacked2D, support oracle, or direct projection of scalar materials into vector modes.

Use identical existing R10 public Y and scalar sigma, iid circular complex CN, expected3% Frobenius noise without realization renormalization. Auxiliary source weights are ws_s=1/||Einc_s||². The inherited physical selection objective Phi=mean_s ||F_s-Y_s||²/||Y_s||² remains a source-normalized physical objective, not iid-sigma ML. All material/current/LU/optimizer histories reset between arms; only immutable background and receiver geometry/DPCS shared.

## New material readout

`src/vector_denominator.py` uses U∈C^(3N×r), K=BU*, with B=U(U*GU) for Galerkin or GU for explicit right-only control. D_m repeats each scalar material three times. A(m)=I+βD_m-D_mK. Small rank-r Torch LU and all60 RHS are kept on the original device; no full current or gradient crosses to CPU inside a readout.

Full metric uses the material-dependent A at each trial and the root-supplied derivative, including sum over all3 components and source metric. Lagged factor uses A(m_old) once per material readout. Both compare old and ordinary-CIE seed under their own explicit objective. Lagged decrease is not a guarantee about the full material-dependent metric or original Phi and is not claimed to be MM.

Eight projected-gradient steps are always budgeted. Projection first produces the feasible direction clip(m-g/||g||∞,0,8)-m. Lagged R0 then computes the exact quadratic alpha on that direction constrained to [0,1]. Full Armijo has at most8 trials per step; only scalar acceptance crosses to host. This is a bounded inexpensive engineering solver, not optimizer-independent mechanism evidence.

## Basis audit

R10 top32 is available with correct5184×32 U/GU and5,308,416 uncompressed array bytes. Existing numerical Gram residuals support that returned finite range; no full spectrum or all-strong completeness certificate exists. s32/s1=.856 is insufficient to decide field relevance.

Public R10 v04/v09 GlobalDual stage3 frozen fields give Galerkin defect propagation relative errors .9255/.9501 and defect projection energy6.83%/1.95%. This motivates root review of one capacity correction; no rank sweep was started. Audit includes same-source weighted J/Q/defect, right-only and Galerkin residuals. `top32_audit/RECEIPT.json` carries CPU/wall/operator cost, zero new physical solves and zero GT access.

## New derivative evidence

One centered directional difference per new variant on realN m=.7 and actual3N×60 public saved J/Q, beta1, nonuniform all-source ws. Full relative derivative error6.95e-9; lagged7.59e-10. Both small-LU infos zero. No old forward/identity tests were repeated. See `derivative_check/RECEIPT.json`.

## Suggested inversion budget awaiting root freeze

Developer cases v04 and v09 only; GlobalDual/full/lagged matched global ordinary-CIE current and same MF2/3/4 beta schedules. Preserve old-stage endpoints. Common physical corrector at least320 accepted/800 strict requests with80/160/240/320 checkpoints and matched physical-LU/RHS/cumulative wall curves. Exact upper budget still root decision. Historical80-step primary is insufficient to isolate representation/readout from corrector truncation.

## Failures and costs

- Initial shell read omitted explicit worktree and found no relative R10 path; no scientific operation performed.
- Default `python3` had no numpy; failed before loading public data or any operator. Reused existing `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`; no install.
- Initial guessed `.venv_nn` directly below Inv_SLAM did not exist; no scientific operation performed. Exact failed startup OS CPU not instrumented, retained as NOT_MEASURED; tool wall was below one second each.
- R10 operator catalog was absent inside this worktree. Read original workspace path given by project instructions; no backend was copied or changed.
- All newly measured audits/readouts preserve full receipts; costs are nested views unless explicitly whole-process.

## Root-approved development amendments

Root rejected rank expansion as a remedy for the order-zero longitudinal dyadic contribution. One physics-motivated basis was authorized instead: `U_src=orth([Einc,G Einc])`, relative thin-SVD threshold1e-10, actual rank120, U/GU19,906,560 bytes. This uses only background physics and the actual60 sources. Its seed singular values are not the Green spectrum. It captures early uniform-contrast propagation but does not cover unknown D_m multiplication.

The source-basis audit captures Q energy98.49%/92.46% on public v04/v09 stage3, but the late defect propagation relative errors remain.8418/.9361. Full initial-alpha now uses the root-approved Gauss-Newton curvature along the feasible projected direction, reusing the currently retained smallLU; the true full metric is still verified by Armijo. Eight steps and at most8 trials remain frozen. The former top32 fixed-state trial is preserved in `fixed_state_cuda`, its source snapshot in `fixed_state_transport/SOURCE.zip`. The amended source120/GN trial is separate in `source_gn_cuda`, with its own source snapshot and transport costs. No old derivative check was repeated.

Root then froze exposed DEV on v04/v09: sourceCIE/sourceDual/GlobalDual/full/lagged, same R10 public Y, geometry, material box[0,8], m0=.25 and R0. Original MF2/3/4 stages retain500 cap/min20. Common physical corrector320 accepted/800 requests saves80/160/240/320 accepted checkpoints, actual paid physical-Phi events and LU/RHS/CPU/wall. Only the first three ordinary arms were launched initially while root supplies a possible sixth z-lagged formula. This is developer comparison, not new confirmation, and the worker never reads GT or scores.
