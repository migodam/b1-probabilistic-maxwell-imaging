# R13 full-current control: implementation and qualification

Status: implementation and one known-public-state wiring exercise completed. No full reconstruction experiment, truth access, remote execution, commit, or publication was performed by this worker. Scientific conclusions and acceptance remain with the parent Codex thread.

## Files and API

New implementation: `research/research_max_struct_r13/src/full_current.py`.

```python
from research.research_max_struct_r13.src.full_current import run_full_current, run_audited_cie

result = run_full_current(
    model, Y, m0, config,
    mode='restricted',     # or 'full_apcs', 'full_current'
    budget_check=check,    # may raise the caller's existing BaseException
    audit_callback=record_paid,  # (m, state, original_Phi, metadata)
    audit_every=50,
)
```

`config` is the inherited R8 `arm_config('a_plus', common, Y, model.E, continuation_stages)` result, augmented exactly as in R12 supplement `src/run_paired.py`:

```python
config.update(
    fft_geometry=dict(coords=public['coords'],
                      radius=float(public['cell_radius']), k=float(public['k'])),
    regularizer=lambda m: model.reg.evaluate(m)[:2],
    aplus_regularizer=model.reg,
    physical_sigma=float(np.linalg.norm(Y)),
    receipt_path=str(new_output_dir/'G19_EXECUTION.json'),
    initialize_physical_current=True,
    cost_callback=collect_receipt,
)
```

The caller must freeze `dpcs_rank` or `dpcs_rtol` using its existing public-data rule and account for that rule's setup. The adapter does not select ranks, stages, priors, CPU budgets, or corrector iterations. It calls original `prepare_full` when `matrix_denominator` is absent, retaining its setup artifacts and G19 `full_material_readout`. An already-prepared denominator is accepted and identified as reused setup; the caller must count its earlier preparation. The same object must not be shared concurrently across invocations because the inherited material readout rebinds its cost counter.

`run_audited_cie(model,Y,m0,config,...)` passes `g19=False` and preserves the supplied CIE/Dual-CIE/Global material block, sourcewise/global PR flag, and original physical-current initialization flag (absent/false remains absent/false). This permits the same physical audit schedule for mature inherited controls. An inherited A+ config needs its original explicit readout supplied by the caller when using this convenience entrypoint.

Callbacks receive only already-paid states. `metadata` includes phase, stage, iteration, adapter CPU time, original Phi, physical relative residual, and physical model costs. Parent callbacks should copy arrays if retaining them. Callback execution and I/O are part of the enclosing caller process's CPU budget. No corrector runs here.

## What changes and what is retained

Let `F_I` denote the original native orthonormal Fourier columns, `Vd` the public receiver DPCS basis, `P=I-Vd Vd*`, and `JD` its unchanged affine determined current. The three controls are:

| Mode | Current family at every stage | Independent complex spatial rank |
|---|---|---|
| restricted | `JD + P F_I A` | measured using complement-rank formula; original masks |
| full_apcs | `JD + P F_full A` | `N-rD` |
| full_current | `JD + F_full A` | `N` |

The two full modes replace every stage's mask with the complete grid mask from the first update. They retain the original number of stages, iteration budgets, beta values, reset at beta/stage boundaries, dual-beta convention, global PR+, original current objective and exact quadratic line search, material denominator and complete derivative, shared-real material box, prior and physical forward model. The unrestricted mode retains `JD` as an affine offset; its Fourier correction can change the DPCS coordinates and cancel that offset, so the admissible family is all currents. It imposes no additional noise ball, soft DPCS penalty, or support constraint.

The actual R12/R8 24×24 continuation `mf=3→6→10→12` already reaches all 576 native Fourier columns at `mf12`. Thus legacy final capacity is already full APCS with fixed DPCS. Fourier restriction is a transient continuation issue; persistent fixed-DPCS restriction is a separate issue. A paid CPU prefix can end before that final stage has produced its inherited stage-end physical audit.

The physical-current warm start is retained in all three modes using the original `lstsq` and its unchanged `rcond=1e-12`. Its loss, numerical rank and singular values are recorded. The full-current warm start can retain the initial physical current to roundoff; full APCS still replaces its DPCS with `JD`.

## Reuse rather than a new optimizer/backend

The adapter binds the actual `research/research_max_struct_r8/src/cie_fft.py::_run_cie` and `run_cie_fft` function code to private per-call global dictionaries. Only the basis factory, observed NumPy least-squares call, and `_run_cie` reference are rebound. No legacy source, module global, PR step, readout body, or optimizer is copied/edited. The restricted mode invokes the original `NativeFourierAPCS`; full APCS uses that same class with the full mask; unlocked current passes an empty DPCS projector to that class. The finite disk-pulse Toeplitz Green FFT actions, dense receiver map, and supplied physical MaxwellModel are exactly the inherited implementations.

The private NumPy proxy observes the original warm-start least-squares solve; it does not replace the result. Its projection diagnostics cost one additional dense G application and one S application to all sources, plus a DPCS projection. These are explicitly charged in `adapter_diagnostic_cost`. Restricted stage rank uses

`rank(P F_I) = |I|-rD+rank(F_complement* Vd)`

with tolerance `1e-12`; its small rank SVD/FFT setup are separately charged. Full-mask rank is exactly the stated `N-rD` or `N`. Native Fourier coefficient counts can exceed independent rank, and all reported actual real source-current DOF use `2*rank*number_of_sources`; material DOF remain N.

The original receipt includes current objective/gradient evaluations, PR steps, Green forward/adjoint actions, receiver products, material inner evaluations, small denominator LU/RHS, physical LU/RHS and Fourier actions. Added fields include current_spaces, warm projection loss, diagnostic costs, fresh denominator setup SVDs/G RHS, whole-adapter CPU/wall times, physical-audit count and whole-adapter physical-model delta. Setup retained in a parent process or imported dependencies must additionally be included in that parent's process-level accounting. Receipt timing is sampled before its final augmented JSON serialization; the enclosing process timer remains authoritative.

## Paid physical audit fairness

R8 originally audits initial and stage-end materials. The adapter retains those and additionally uses the existing `iteration_callback` every `audit_every` material updates (default 50). It evaluates the original physical material objective

`Phi(m) = 0.5 ||Y(m)-Y||_F² / ||Y||_F² + R(m)`.

It does not restart or chunk the PR loop, update current coefficients from a physical audit, or change the material readout. Its returned `m` is selected from all eligible periodic and inherited physical audits, rather than stage-end candidates alone. The inherited `result['candidates']` is left intact for provenance; complete selection candidates live in `physical_audit_candidates` and the winning record is `selection`. Duplicate periodic/stage-end materials use MaxwellModel's existing cache and do not incur another LU.

Budget checks surround setup, warm projection diagnostics, physical audits and callbacks, and the original loop retains its per-iteration budget checks. A physical forward finishing after the caller's deadline is charged and excluded by the post-forward check. A caller exception propagates; already delivered callbacks retain the eligible incumbent and the finally receipt records `INTERRUPTED_BY_CALLER`. Setup interruption emits a separate setup receipt. Parent-level callback eligibility rules remain authoritative if the supplied observer also records states internally.

In-flight SVDs, physical LU and inherited material inner solves remain non-preemptible. A deadline can therefore be overrun by the current operation; the caller must count that excess and exclude late candidates. Final corrector allocation and whether to use a last eligible audit stay with the parent runner.

## NFFT qualification and existing compatible sources

The supplied R13 review identifies Wei et al. (2016, EIT) and Wei–Chen (2018, inverse scattering) as full-Fourier corrections to determined current, including `J=JD+F alpha`. This adapter is explicitly a **full-Fourier/NFFT current-space transplant into G19**, not a reproduction of original NFFT-SOM's full algorithm, stopping, initialization, or material update. It must not be labeled a mature NFFT-SOM benchmark.

Source/file-name searches of the compatible MAX-STRUCT and FBG-R3 sources found the mature CIE/Fourier-TSOM code and the ordinary CSI mode, but no complete independent NFFT-SOM solver. Existing `research/research_max_struct_r3/src/cie_strong_release.py` adds DPCS columns only at the final stage and hooks the dense mature solver; it is a related historical control, not the full-space-from-first-step matched G19 control. Ordinary `method='csi'` changes the state-residual convention to uncontracted LSIE and is not a substitution for this experiment. No independent original-NFFT implementation was fabricated or run. The parent literature worker retains responsibility for verifying original source papers and deciding how to add a qualified mature NFFT-SOM comparison.

## Narrow wiring evidence

`validate_wiring.py` reads only the public R12 input `A__r10_twod_0004.npz` (576 cells, 12 sources, 48 receivers), uses the existing public noise-ball rank rule (`rD=13`), and performs exactly two updates in one explicitly shortened wiring stage per mode. It also exercises inherited Dual-CIE audit scheduling and caller interruption. This is neither a frozen-budget inversion nor recovery evidence.

| Wiring family | Independent spatial rank | Actual real source-current DOF | Relative initial-current projection loss |
|---|---:|---:|---:|
| restricted G19, mf3 | 36 | 864 | 0.9558910890 |
| full APCS G19 | 563 | 13,512 | 0.2327604403 |
| full current G19 | 576 | 13,824 | 3.58e-15 |
| inherited Dual-CIE, mf3 | 36 | 864 | not projected: inherited default initializes at JD |

All four exercises retained exactly two inherited PR steps, four paid audit callbacks (initial, each of two updates, cached stage-end), and three physical LUs/36 complex RHS. Callback Phi matched the original model's cached physical_score to `1e-14`. Full APCS projection loss agreed with its determined-component mismatch. A separate caller interruption after the first periodic audit retained its eligible prefix and emitted the expected receipt. Full numerical receipts are in `validation/WIRING_VALIDATION.json` and per-run `validation/*/EXECUTION.json`.

Interpreter: existing `Gaussian/.venv_nn/bin/python`; BLAS threads=1; no dependency installation. No repeated frozen operator identity/error-floor tests were run. Legacy R8 and R11 source diffs remained empty.

## Remaining limits

- Real scalar-TM uniform even-square 2D only; current controls do not qualify vector 3D, nonuniform NFFT, or measured-data transfer.
- Periodic audit selection changes the available incumbent set relative to historical stage-end-only reports. Main matched controls must all use the same declared audit policy and pay its physics cost.
- Inherits original dual-beta/non-MM limitations; lower auxiliary loss alone does not prove lower physical Phi or better material recovery.
- No recovery, speedup, basin, theorem, novelty, or universality conclusion follows from the wiring exercise.
