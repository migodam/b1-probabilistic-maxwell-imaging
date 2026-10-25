# E5 sensing runner: code ready, execution NOT_RUN

`tools/e5_sensing.py --bank results/e2e3/finite_bank --reference-gate results/e2e3/inference/gate.json --output results/e5/finite_sensing [--ledger PATH]` must be executed by the parent's serialized E5 wrapper. Existing output is refused, preserving partial paid evidence. Syntax parsed only; worker did not import or run this tool.

The gate must state `decision=PASS`, `prior_id=FINITE_PRIOR_ENUM4_V1`, `reliable_actual_reference=true`. This is conditional finite-prior reference qualification only. The frozen supplement is recorded. No continuous-prior or scientific performance claim is emitted.

For each of 12 parents × 2 independent initial-noise repeats, t1 is recomputed exactly from one TX0/RX8 likelihood with its existing unique ID. All seven policies choose a new independent full16RX source action from indices0–3. Initial TX0/RX8 remains in the posterior even if source0 is chosen again: the newly drawn observation has a distinct ID and independent noise. Uniform uses `(2*parent+repeat) mod4`; Random uses seed2026102506 with case-specific substreams.

EIG and TaskEVSI use 32 current-posterior predictive outcomes with identical latent draws and standardized noise across candidate actions; no actual future observation or truth is passed to selection. EIG uses exact finite-grid posterior KL against current weights; TaskEVSI uses reduction of target-pose Euclidean squared-loss Bayes risk. These are 32-draw MC utilities without global certificates. Shared outcome computation is fully attributed to each policy as a standalone alternative.

At each current early posterior MAP, a fresh E5 Evaluator pays full4TX shared J, with no prior LU presumed available. Planned total: 24 FE + 480 RHS. Physical columns are epsinf, ell, px, py; singular-value comparisons use explicitly declared scales `[4.5,.25,2,1]` and noise whitening .01. Angle/retention use physical A (two wall columns) and B (two pose columns), SVD response-range projections, no pseudo-inverse and no jitter. Absolute wall information remains in physical coordinates. Numerical ranks and thresholds, singular values, principal angles, retention eigenvalues, and zero-response exclusions are saved. All maximizers resolve exact ties by lowest source index.

Durable action selections are saved **before** offline truth indices are opened or actual future noise generated. Actual full4TX future measurements/noise use seed2026102507 with separate case substreams and cached bank means. They are offline evaluation artifacts, never chooser inputs. No extra truth solves occur.

Outputs under the chosen directory:

- `MANIFEST.json`: configuration, reference gate, costs, parent-unit warning, partial failure or completion.
- Per `parent_XX_repeat_R/`: `PAID_LOCAL_J.npz`, `ACTION_SELECTION.json`, `POSTERIORS.npz` (early plus each policy), `ONLINE_RESULTS.json` (pose/occupancy/CDF, risks, independent factor IDs, CPU/wall timing, standalone and shared costs).
- `offline/FUTURE_parent_XX_repeat_R.npz`: actual independent future measurements/noise, factor IDs and truth index.
- `offline/GT_parent_XX_repeat_R.json`: ground-truth pose and squared errors, separate from online inference.
- `offline/PARENT_GROUPED_GT_SCORES.json`: each parent's two noise errors averaged per policy; 24 noise cases are not treated as 24 independent parents.
- `ONLINE_RESULTS.json`: aggregated online diagnostics.
- `READY.json`: created **only after actual completion**, status RAW_RESULTS_READY, scientific_claim NOT_RUN.

The bank's 11800 RHS / 2950 FE preparation is recorded separately as E2E3 shared setup, including its raw receipt when present. Every policy's standalone common finite-posterior attribution retains that setup; Random/Uniform selection itself has no Jacobian charge. Local controls have their real standalone 20 RHS / 1 FE attribution, and the same single paid J is shared between three controls per case in actual accounting. EIG/TaskEVSI never claim cheap online inference after hiding bank cost. Standalone policy alternatives must not be summed as actual paid charges.
