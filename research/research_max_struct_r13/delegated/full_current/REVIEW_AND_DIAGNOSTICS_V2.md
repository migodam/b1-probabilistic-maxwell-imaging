# Runner review and bounded public diagnostics v2

## Technical review and coordinated correction

Reviewed `src/run_experiment.py`, `src/material_events.py`, and `src/finite_queue.py` before the baseline campaign. The projection-only development pilot was already scoped separately; its runtime was not edited.

Coordinated with the parent and corrected the new runner's initialization policy: `initialize_physical_current=method.startswith('g19')`. G19's three current-space variants retain R12's physical-current warm start; mature CIE/Dual/Global and original R7 A+ retain their original false setting. The adapter receipt now explicitly records the flag. This prevents a newly warm-started CIE variant from silently being presented as the historical method.

Other reviewed limits, left unchanged during the projection pilot:

- `random_event` matches the fractional event displacement norm, not the integer event's `floor(t*)+1` displacement. Integer-vs-random differences therefore cannot be attributed solely to geometry. Near `t*=1` the mismatch can be large.
- `normal_double` is documented as a threshold-independent `t=2` control, but `event_start` first returns unchanged when no finite threshold exists. This conditionally disables that control; it needs an explicit future code-version decision, not a silent pilot change.
- Unknown non-baseline method names fall into the DR branch; manifests must be validated against the intended supported set. The reviewed development manifest uses intended names.
- The finite queue's child bootstrap is intentionally Windows-specific (`NORMAL_PRIORITY_CLASS`), so it is not a portable local macOS launcher. It was not run locally or modified.
- Whole-process CPU includes import/setup time. Individual LU, SVD and inherited material inner solves are non-preemptible; late candidates are excluded, but overrun compute remains charged. No early-termination or budget-equality claim should discard that overhead.
- Paid physical selection in the new adapter uses all eligible periodic audits, not just stage ends. This is a declared audit policy change requiring matched controls. No PR restart is introduced.

## Equal-area classic source selection

`full_current.py` now accepts explicit `config['physics_family']='equal_area_disk'`, coordinated with the physics worker. It does not shrink equal-area radius `h/sqrt(pi)` to satisfy the original non-overlap guard.

- Cold classic CIE/Dual/A+ selects the actual `paper_cie_fft` loop, its EqualAreaGreen/PaperOneHop actions and historical stage coefficient thresholds.
- G19 requires its existing physical warm start, which the historical paper loop does not implement. It therefore keeps the actual R8 warm-start loop and replaces only its Green/OneHop factories with the compatible equal-area classes. It retains its R8 global stopping convention; supplied per-stage threshold fields are inactive and the receipt explicitly says so.
- All three G19 current-space controls share that same loop and stopping policy. This is an explicitly adapted G19 control, not a claim of a published equal-area G19 baseline.

Receipts record physics family, actual inherited loop entrypoint, effective coefficient thresholds and stopping policy. Legacy sources remain untouched.

The new equal-area selection was validated using the actual existing Austria `PUBLIC_DATA.npz` (Y and acquisition metadata only), unchanged radius, and at most three deliberately shortened updates. A deliberately high stage threshold stopped cold CIE after two updates in the paper loop; all three G19 spaces ran three updates with global threshold zero in the original R8 warm loop. G19 warm projection was observed for each space and full current retained it to roundoff. This directly tests the new selection/stopping plumbing and is not a reconstruction experiment or repeated operator identity check. Receipts: `equal_area_wiring_v2/VALIDATION.json` and per-family `EXECUTION.json`.

## Fixed-current profile diagnostic

New `src/profile_diagnostic.py` evaluates one fixed scalar-CIE current objective at beta=1 with the actual public G/E/S/Y and inherited per-source normalization. It computes stable exact weighted least-squares profiles using a pivoted thin QR of `C(m)U`, then a small receiver SVD. It does not form normal matrices, modify the material prior, or assert a common objective for dual-beta G19.

The three nested current families are original mf3 projected Fourier, full APCS with fixed DPCS, and full current. Their penalty split is

`Fourier = psi_restricted - psi_fullAPCS`,

`fixed_DPCS = psi_fullAPCS - psi_fullCurrent`.

At three already exposed public initializers:

| Case | Public DPCS rank | Fourier penalty | Fixed-DPCS penalty | Total |
|---|---:|---:|---:|---:|
| A r10_twod_0004 | 13 | 0.854071 | 0.775631 | 1.629703 |
| B r12_confirm24_0010 | 12 | 4.987261 | 0.199347 | 5.186608 |
| C r8_test_0001 | 13 | 0.292078 | 0.095878 | 0.387955 |

Five deterministic public material directions are negative original physical Phi gradient, negative full-profile gradient, negative Fourier-penalty gradient, negative DPCS-penalty gradient, and an unpropagated linear-port backprojection. Each is normalized in material infinity norm and probed at ±0.02 with the same legal box retraction. No direction or step is selected by truth or material recovery.

Total penalty varies over those probes: A `[1.594256,1.629703]`, B `[5.093511,5.187126]`, C `[0.382878,0.392973]`. This demonstrates material dependence in this fixed-block diagnostic; it does not establish a favorable or harmful correlation with unknown truth. Every base-profile tangent stationarity norm is approximately `7e-16`–`5e-15`. A new one-sided envelope-derivative wiring check has absolute first-order remainders below `1.3e-13` at material step `1e-6`.

Each case performs twelve material-profile evaluations (base, ten signed probes, one derivative-wiring point), each across three families; 36 thin QR/small-SVD profile solves. There was no rank-deficient fallback. Physical cost is only one LU and one adjoint per case (24/24/12 total complex RHS for A/B/C). Local measured total is approximately 2.81 CPU seconds, BLAS1, setup included and final output serialization excluded.

Artifacts: `profile_diagnostic_v2/{case}/PROFILE.json` and `PROFILE_MATERIALS_AND_BASE_CURRENTS.npz`. They contain source-weight/rank/DOF details, actual material directions, base currents, profile and penalty gradients, signed material probes, linear predictions, observed changes and complete cost counts.

## Early AP-pair forecast diagnostic

New `src/event_forecast_diagnostic.py` continues AP only to checkpoints 50, 128 and 512, then independently follows actual DR for 30 steps from each saved state. Forecast trajectories never feed back into AP continuation. It retains the actual checkpoint `c=P_C(z)`, `d=z-c`, material/ray traces, Gram coefficients and observed local deviation radii.

The ideal normal-ray first-posttie prediction, first angular-half branch flip, and first large material jump are separately reported. A large jump is descriptively defined as one quarter of the legal material interval (4.995), not a success gate. An angular-half crossing can be a smaller continuous movement and is not identified with the finite-threshold eigen-branch event.

| Case/AP | Predicted future step | Actual first angular-half flip | Actual first large jump | Initial relative pair defect |
|---|---:|---:|---:|---:|
| A/50 | 5 | 5 | 5 | 0.01442 |
| A/128 | 5 | 5 | 5 | 0.00468 |
| A/512 | 5 | 5 | 5 | 0.00274 |
| B/50 | 5 | 6 | 6 | 0.03900 |
| B/128 | 7 | 7 | 8 | 0.01853 |
| B/512 | 10 | 7 | 10 | 0.00861 |
| C/50 | 3 | 1 | 3 | 0.09761 |
| C/128 | 7 | 2 | 7 | 0.06170 |
| C/512 | 10 | 5 | 10 | 0.01549 |

The finite threshold matches first large jump on seven of these nine exposed traces, with the other two one step later; this is a descriptive local result, not independent statistical validation or a future-trajectory theorem. Angular-half crossings differ substantially on C and B/512. By DR30 actual trajectory deviation is 22.7–31.5 times the initial normal gap, so the exact stagnant-ray interpretation is not a trajectory approximation over this whole horizon.

`angular_branch_margin` uses the **observed actual per-cell deviation from the ideal ray** as its state-ball radius. Across all saved points it produces zero disagreements on certified cells, but it is an a posterior certificate: that radius is known only after running DR. It becomes nearly uninformative at step 30 (0–4 certified cells out of 576); this cannot be advertised as a future prediction guarantee. At step 5 the corresponding certified-cell counts range from 125 to 570.

Cost: 512 AP plus three independent 30-step DR traces per case, 2 physical LUs per case (original initial material and final longest-checkpoint material), 24/24/12 physical complex RHS. Total approximately 4.11 CPU seconds locally, BLAS1, including setup and trace serialization. A stdout-formatting error occurred after A's completed JSON/NPZ were written; it was fixed and A was reused without rerunning its science, then B/C completed. No result was lost or suppressed. No GT or recovery scoring was used.

Artifacts: `forecast_diagnostic_v2/{case}/FORECAST.json` and `AP_{50,128,512}_STATE_AND_TRACE.npz`. These are early-pair sensitivity/forecast-limit diagnostics, not repeated frozen operator identities and not an extension of the projection pilot's scientific budget.

## NFFT qualification

`NFFT_PRIMARY_AUDIT.md` records the actual 2018 original-source comparison. The state normalization, initialization, block update specification, material adaptation and stopping do not justify relabeling inherited CSI or full-current G19 as a mature original NFFT-SOM reproduction. That benchmark remains an explicit gap.
