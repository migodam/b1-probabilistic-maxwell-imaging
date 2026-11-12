# Qualified NFFT-equations / joint L-BFGS-B control

Implementation: `src/nfft_control.py`. Qualification everywhere: **NOT_ORIGINAL_NFFT_SOM**. This is an explicitly requested Eq(3)–(6) adaptation added after the primary-source audit, not an assertion that the missing mature author algorithm has been reproduced.

## Declared objective and coordinates

For the original public receiver SVD and frozen rank, retain the original determined current JD. Optimize complete current corrections and one shared real material per cell:

`J=JD+Z`, `Q=E+GJ`, `e=J-mQ`, `D=SJ-Y`,

`f(Z,m)=sum_s wd_s ||D_s||² + ws_s ||e_s||² + R(m)`,

`wd_s=1/||Y_s||²`, `ws_s=1/||JD_s||²`.

There are **no half factors**, no incident-field state normalization, no CIE contraction and no separate material readout. The common physical prior R and material box are passed by the caller. The spatial correction Z equals `F alpha` for complete unitary Fourier F; this is an orthogonal change of the global joint Euclidean coordinates. It adds no scaling or preconditioner and preserves the admissible full-current family. Initialize Z=0, so the first joint current is JD; the material initializer is exactly m0. A zero/unsupported DPCS rank or zero determined-current source energy is rejected rather than repaired by an arbitrary epsilon.

The analytic gradients are

`gZ=2[S* (D wd) + e ws - G* (m e ws)]`,

`gm=-2 Re sum_s conj(e_s) Q_s ws_s + grad R(m)`.

Packing real(Z), imag(Z), m gives the real joint gradient used by SciPy L-BFGS-B. Only material coordinates have bounds. No tuning based on material truth or GT residual is present.

Explicit adaptations relative to the [2018 conference equations](https://usncursi.org/archive/aps-ursi/2018/papers/1927.pdf): real shared material, the common project box/prior, the observed-data interpretation of the printed data denominator, and joint L-BFGS-B in place of unspecified original CG/update details. This control tests whether a strong known optimizer on those equations explains a comparison; it does not satisfy the separate mature-original-NFFT benchmark qualification.

## API and budget/selection behavior

```python
from nfft_control import run_nfft
result = run_nfft(model, Y, m0, config, check, audit_callback, audit_every=50)
```

Config requires exactly one `dpcs_rank` or `dpcs_rtol`. It accepts `material_box`, `regularizer`, `maxiter` (3000 default), `ftol` (1e-12), `gtol` (1e-8), `maxls` (20), optional `maxfun`, `receipt_path` and `cost_callback(receipt)`. The default L-BFGS memory rule is unmodified. No stages, beta continuation, line-search rescue, current warm projection, physical corrector or hyperparameter search is added.

The optimizer uses the declared joint equations objective. Candidate selection remains the original physical Phi, through paid audits of initial material, every 50 accepted joint updates, and a terminal material if it differs from the last audit. A trial joint objective evaluation is never used as a physical candidate without that forward audit. Periodic audits do not reset optimizer state. The caller owns its absolute process deadline and subsequent corrector allocation.

Budget checks surround receiver setup, joint evaluations, callbacks and physical audits. A physical audit that finishes over budget is charged and excluded. Caller BaseException propagates and the finally receipt still reaches `cost_callback` and `model.r13_nfft_cost`. An interrupted solver does not return an arbitrary latest iterate; the parent retains already delivered eligible physical candidates. CPU includes all adapter work and failed/interrupted work; the enclosing process timer must additionally include imports, final receipt/output serialization and exit.

Every joint evaluation counts one dense G action, one G adjoint, one receiver action, one receiver adjoint, their complex RHS and one prior evaluation. Receiver SVD and physical selection/prior work are counted separately. Physical model deltas are already nested in final model counters and must not be added again. No physical adjoint is performed by this control.

Integration: use `config['cost_callback']=auxiliary_receipt.update` and the returned status. The full returned result contains complex `final_joint_J`, array candidates and joint history; do not blindly merge it into the existing JSON-only receipt encoder. Save joint arrays separately if required. Returned `m` is the paid physical-Phi incumbent; `final_joint_m` and `final_joint_J` are distinct diagnostic outputs. Optimizer success and termination message remain explicit; reaching an iteration cap is not convergence.

The adapter currently supports common receiver maps and real scalar material. It rejects vector material and source-specific Fresnel ports; neither model is silently transplanted.

## Narrow wiring evidence

`validate_nfft_equations.py` uses only the already exposed public A r10_twod_0004 initializer/operators/observations. Its receiver rank13 was previously frozen from public data. Two new derivative directions check real/imag current and real material independently through central differences: relative discrepancies `1.86e-9` and `4.59e-10`. This verifies new derivative wiring, not an old frozen operator identity.

The one bounded wiring run accepted exactly 50 joint updates / 53 objective-gradient evaluations and produced two physical audits (initial and update50). Its optimizer correctly reports an iteration limit and success=false. It used two physical LUs / 24 physical source RHS and no physical adjoint. The 53 joint evaluations charged 636 RHS to each G/S forward/adjoint category. Local total validation CPU was about0.138 seconds excluding imports, with BLAS1. Separate derivative-check setup/evaluations are recorded in the validation JSON and are not concealed inside optimizer counts.

Artifacts: `nfft_equations_wiring/VALIDATION.json`, `EXECUTION.json`, `WIRING_STATES.npz`. No full inversion campaign, fresh material scoring, truth access, remote execution, dependency installation or legacy edit was performed.

A separate controlled caller interruption immediately after the initial paid audit (zero optimizer updates) propagated correctly and retained its finally receipt/eligible audit. Its extra one LU/12 source RHS and CPU are recorded separately in `INTERRUPTION_VALIDATION.json` / `INTERRUPTED_EXECUTION.json`; they are not hidden in the 50-update receipt.
