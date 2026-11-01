# Scalar-TM CIE-I / Fourier-TSOM implementation audit

Implementation: `src/cie_tsom_baseline.py`. Status: **RUNNABLE ON PARENT REAL MAXWELL MODEL**. This is a mature-structure 2D adaptation, not a reproduction of the 2019 3D numerical benchmark and not an inference of imaging superiority.

## API and public inputs

```python
result = run_cie(model, Y, m0, config)
# result: m, candidates, history, cost, last_joint_m, status,
#         method, receiver_rank
```

`model` supplies complex `G/E/S`, `shape`, and `state(m)` returning physical `j/Q/Y`. `model.cost` is sampled before and after. The solver reads no material labels, support, GT image metrics, or hidden arrays. Full material `m` has one real unknown per cell.

Required config:

```python
config = dict(
    method='cie',                 # 'tsom': LSIE Fourier-APCS; 'csi': full current LSIE
    dpcs_rtol=1e-3,                # or a preregistered dpcs_rank
    stages=[dict(mf=1, beta=1., maxiter=300),
            dict(mf=2, beta=.2, maxiter=300)],
    material_box=(.02, 20.),
    regularizer=model.regularizer, # returns (value, real gradient[, other])
    physical_sigma=np.linalg.norm(Y),
)
```

The example beta values are the declared local calibration schedule, not transplanted published constants or GT-optimized settings. For dual beta set second-stage `beta_current=1.` while `beta=.2` remains the material beta. `physical_sigma=||Y||` matches the parent's normalized original data objective. A supplied `physical_whiten(residual)` overrides this scaling. The joint baseline uses positive fixed source weights; **general correlated joint-data whitening is not implemented**. The current experiment is clean; a future correlated-noise comparison needs the same whitening in the current block, not only final scoring.

`result['m']` is the best physically scored incumbent among initial and stage endpoints, while `last_joint_m` is the raw final continuation endpoint. `candidates` retains each real material and its original Maxwell score. A common full-material corrector is applied by the parent uniformly, outside this function.

## Actual optimization

The joint objective for one beta is

`1/2 sum_s wd_s ||S j_s-Y_s||² + rho/2 sum_s ws_s ||beta j_s-r_beta(m)(E_s+(G+beta I)j_s)||² + R(m)`.

For TSOM and CSI the state residual is `j-m(E+Gj)`; this is not obtained by putting beta=0 in the unscaled CIE objective. TSOM fixes the receiver DPCS and restricts current to projected low-frequency Fourier columns. CSI allows every complex current entry and has no fixed DPCS, providing an unrestricted ordinary-lifting control.

One sourcewise PR-CG step uses the exact quadratic line minimizer for fixed material. Polak–Ribière coefficients are clipped below at zero and a non-descent direction resets to the negative gradient. Native Fourier coefficient coordinates are retained: orthonormalizing projected columns while leaving a Euclidean coefficient CG metric would change the iteration. Roundoff-only singular actions of projected Fourier columns are removed by SVD **while reconstructing an operator in the same original coefficient coordinates**. No nonzero representable mode is deliberately discarded, and no coefficient ridge is added.

Stage transport exactly zero-pads the old Fourier coefficients into nested masks. The state estimate is inherited, and the same physical material generates both `r_betaJ(m)` and `r_betaR(m)`. Dual beta uses two state metrics and is not advertised as block descent on a single scalar function.

The native stopping statistic is `mean_s ||delta_gamma_s||² / ||gamma_s||²`, with threshold `1e-3`, a minimum of two iterations, and the configured cap. A zero initial denominator is handled by a numerical floor. This is the original type of coefficient stopping rule; `SCHEDULE_COMPLETE` means scheduled stages completed under those rules, not a physical or material convergence certificate. Every stage endpoint is re-solved by original Maxwell physics.

## Shared real material and common prior

The paper permits complex contrast. This implementation intentionally fits one **real** modified contrast per cell, shared over all incident sources. It minimizes the weighted real LS block directly, rather than doing complex LS and taking the real part afterward.

For positive beta the physical box `[lo,hi]` maps monotonically to `[beta*lo/(1+beta*lo), beta*hi/(1+beta*hi)]`. The real LS result is clipped in that legal r interval **before** inverse transport; it cannot cross `r=1`.

Without a common spatial regularizer the fixed-current material block has the exact bounded real r-LS solution. With a common smooth/TV material prior, that pointwise formula is **not** the minimizer of the full material block. We solve the actual state term plus `R(m)` by bounded L-BFGS-B, initialized by the better of the old point and legal r-LS seed. Its analytic residual derivative is

`d[beta J-r_beta(m)(Q+beta J)]/dm = -beta*(Q+beta J)/(1+beta*m)^2`.

This equals `D_beta' e-D_beta Q` and includes the off-manifold `DP e` term. A failed or worsening inner solve retains the better feasible input. Inner evaluations, status and iterations are recorded. This modification preserves the **same material prior**, whereas a ridge on r would change it.

## Computation and costs

This small-grid implementation constructs IDFT Fourier columns explicitly and reuses dense `G@basis` / `S@basis` products. The representation is the Fourier-TSOM subspace; **it does not claim the CG-FFT speed or memory scaling of the original large-grid implementation**. There is no new approximate Maxwell backend. FFT calls are zero; explicit IDFT column generation, receiver/APCS SVD, Green/receiver column actions, reduced products, inner material evaluations and original model costs are separately counted. Actual LU and complex RHS costs are in `cost['model']`; source and basis propagation are not relabeled as free.

Wall timing includes solver preparation, failed/inexact inner attempts, physical endpoint evaluation and selection. It excludes the parent's already-created geometry/operator assembly, which the parent must retain in its total ledger. Wall caps are checked between operations, so an atomic SVD, material block or endpoint Maxwell evaluation can cross a cap; that paid completion remains in recorded wall time. No iteration-count-only speed comparison is justified.

## Parent-model execution and time-only calibration

`audit/literature_baseline_execution.json` records eight initial real 4×4 parent-model executions: single/dual CIE, TSOM, CSI, with/without common prior, all finite and within the box. Those checks used a declared public uniform calibration material solely to verify execution, not the E1 test set.

The final time-only calibration is reproduced by `python -m research.research_fbg_r3_experiments.audit.literature_calibration_run`. It opens only `G`, `E0`, `S`, `grid_shape`, `Y` in `data/dev_calibration.npz`; **no material array is opened**. Original data residual is reported without GT error or support readout.

At the public spectral gap, `dpcs_rtol=1e-3` retains **12** receiver modes, leaving APCS rank **4**. Relative receiver singular values drop from about `.02694` to `9.22e-5`; no threshold was retuned using material recovery. An early worker message incorrectly stated rank16; that statement was immediately corrected and is not evidence. The stored receipt consistently records rank12. The tentative `.05` suggestion was never run.

For 4×4 both masks `mf=1` and `mf=2` already span the same four-dimensional APCS; the second stage changes beta and coefficient representation, not current capacity. Fourier coefficient counts are 4 and 16, corresponding to 48 and 192 real source coefficients for 6 sources, but the independent current capacity is 4 complex spatial directions ×6 sources =48 real directions at both stages. The full CSI capacity is 16×6 complex =192 real directions. Material capacity is always 16 real unknowns. This grid limitation must remain visible in interpretation; the 6×6 schedule can supply a meaningful mask-capacity progression.

Final receipt: `audit/literature_calibration_timing.json`.

| Run | Actual alternating steps | Wall seconds | Final original relative data residual |
|---|---:|---:|---:|
| CIE single beta, paper stop, 300/stage cap | 22 | .05675 | .32659 |
| CIE transported dual beta, paper stop, 300/stage cap | 23 | .06532 | .078185 |
| Ordinary LSIE Fourier-TSOM, paper stop | 39 | .03150 | .85498; stage-1 physical incumbent .85493 retained |
| Full-current ordinary CSI, paper stop | 13 | .01094 | .84159 |
| CIE single beta, strict stop, 500/stage cap | 305 | 1.18616 | .00042738 |
| CIE transported dual beta, strict stop, 500/stage cap | 305 | 1.15870 | .00026615 |

All schedules stopped by their declared coefficient criterion. **The strict configuration is `coefficient_tol=1e-8`, `miniter=50`, `maxiter=500/stage`, with the SAME beta and mask schedule.** It is frozen as the adequate-duration strong baseline, while the original paper stopping variant remains separately reported. No third beta stage or parameter search was introduced. The sharp improvement under longer iteration explicitly prevents calling CIE weak from the 22/23-step endpoints. No material GT was opened to make this decision.

Three timing runs each performed **50 actual alternating steps**, with stopping disabled and minimum iteration count 50. The first median **T0=.1989769163 s** was reported before the strict runs. Repeating the reproducer while adding strict rows gave **T0=.1855031247 s** in the final receipt; this timing variation is exposed, and the parent may retain the initially frozen T0 for its budget grid. These are setup-inclusive measurements on this process and grid, not universal speed claims. Every full schedule's complete paid cost is a required additional budget point; strict CIE costs approximately 1.19 s here, exceeding 4T0. No 576-cell baseline, paper image reproduction, success rate, or baseline superiority is established.

## Remaining limits

The original FFT-TSOM 2011 full text and 2016 CIE paper were not independently reacquired here; 2019 directly supplies the implemented structure, and the supplied report supplements earlier context. Complex materials, vector 3D grids, general correlated joint noise, LUCIE and CC-CSI are not implemented in this module. A strong-baseline comparison must name these limits and include the shared physical corrector, full cost and frozen schedule. R2's 100-step prototype remains a separate historical control.

## Post-closure576 coverage and final provenance clarification

This later section supplements the earlier4×4 calibration record. A real576 budget-point comparison has now been completed through `src/cie_tsom_cuda.py`; the earlier statement that calibration alone establishes no576 result remains a statement about that earlier phase. The48 closed receipts and the independent post-closure scorer are recorded in `results/e3_all_scores.json` / `.csv` and `audit/E3_BASELINE_PORT.md`. Scorer status is `POST_CLOSURE_SCORED`, classification is `DISCLOSED_HISTORICAL_MATURE_BASELINE_DIAGNOSTIC`, and its truth-access record explicitly occurs after parent closure. This provenance check reads scorer outputs, not raw GT arrays, and launches no numerical solver.

For these twelve specific historical clean/amplitude12 objects, the joint gate is amplitude-normalized NRMSE≤.15, reconstruction-Otsu IoU≥.7 and original physical relative data residual≤.001. The selected online physical incumbent, not GT-ranked best@K, gives:

| Arm | Joint pass/12 | Data pass/12 | Exact finite-set coverage |
|---|---:|---:|---|
| Transported dual-beta CIE + original-data L-BFGS corrector |10 |10 |All nine specific digits1–9 and two_islands pass; ring and curved_bar fail |
| Single-beta CIE + same original-data L-BFGS corrector |8 |9 |Digits1,2,3,4,6,8,9 and two_islands pass; digit5 passes material/data but fails IoU(.67826); digit7, ring and curved_bar fail |
| Fixed beta=.2 feedback + original-data corrector |8 |8 |Digits2,3,5,7,8,9, ring and two_islands pass; digits1,4,6 and curved_bar fail |
| Native Full-J TRF |0 |0 |All twelve fail the joint gate |

“All nine digits” means those nine fixed raster objects in this disclosed set, not arbitrary digits, fonts, morphologies, noise levels or amplitudes. The amplitude12 denominator is verified in the scorer's row metadata and `material_metrics(..., meta['amplitude'], ...)` implementation. A stale generic `definitions.nrmse` string in the current JSON mentions “8 or18”; it must not override the actual E3 amplitude12 contract.

The solver always exposes all576 real material cells at its final24×24 optimization stage. Native and fixed-feedback arms reach and select stage24 in all twelve receipts; their earlier8→12 stages are public material multiscale initialization, not known target support. CIE material updates and its corrector are full576 throughout. The Fourier restriction acts on auxiliary current coefficients, not a digit/material decoder. No GT shape, label, mask, support or outside-background freeze enters any online update or candidate ranking. Only the independent closed scorer uses stored GT support for evaluation. No default Codex provider/model, research policy, physical operator, precision, prior, acquisition or frozen post-pilot solver configuration was changed.

The source fidelity is to the directly read2019 principles: receiver SVD DPCS, projected Fourier APCS, source-normalized data/state residuals, one-step sourcewise PR current CG, material alternation and staged continuation. The fixed-material current residual/gradient/line-minimization mathematics follows that published formulation, with disclosed PR+ / non-descent-reset safeguards. The single arm uses one beta for both blocks. The dual arm restores the same physical m between blocks and applies beta_current=1 versus beta_material=.2 in later stages, following2019's different-beta principle; it is not the same objective or a claim of common monotone block descent. Restricted real material plus the shared spatial prior requires the disclosed CPU regularized inner readout, rather than literal complex contrast LS from the paper.

These are explicit dense-Fourier scalar-TM adaptations, **not** literal3D/complex-material reproductions and **not** FFT/CG-FFT complexity or speed claims. They reuse the original complex128 Maxwell operators for final validation. LUCIE and CC-CSI remain outside this four-arm comparison; the parent's separate native-TSOM scoring batch needs its own closed provenance and is not silently merged into these48 receipts.

The most important remaining limit is duration: current schedules are wall-truncated in9/12 single and11/12 dual trajectories; three dual cases never reach mf10. Even `SCHEDULE_COMPLETE` permits500-step MAXITER stages. The reported finite-set result is a45-second budget point with a shared original-data corrector, not a mature convergence/performance ceiling, and not evidence of an independent N1 gain. N1's E3 scientific gate remains NOT_PASSED.

Budget equality is nominal, not matched actual compute or search opportunity. All arms share two FIFO lanes and the same physical caps, but their stage plans, optimizer stopping, current/proposal costs and corrector time differ. Method order was not randomized or counterbalanced, and co-running tasks share GPU/CPU resources; ordinary CUDA streams promise no fixed SM quota. Setup is paid separately and shared once per operator cache, not assigned identically per arm. Thus throughput/resource contention, warm-up/order and stopping rules remain comparison limits. The raw receipts, not the scorer's currently missing aggregate stop/cap metadata, are authoritative about actual LU/RHS stops and wall overruns.
