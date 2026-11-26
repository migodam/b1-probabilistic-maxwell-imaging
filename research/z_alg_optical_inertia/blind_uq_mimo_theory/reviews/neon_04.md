# Independent audit of the final numerical results

This is the reviewer’s reading of the supplied numerical snapshot. The follow-up on strong coupling and the two different residual upper bounds are detailed in [the response](RESPONSE_ZH.md). Local workspace addresses have been removed.

**1. Audit: Fixed Pair vs. Pool vs. Class (Verified)**
BENCHMARK_REPORT_ZH.md mostly separates these scopes across the 28 ROIs (24/28 exact `[0,4.5]`, 7 continuous; region_intervals.csv) and 21 `aperture_joint` endpoints (`1/7` zero, `3/7` patch, `3/7` unrestricted; 3 certified in verified.json). Two verified flaws remain:
- **Stale residuals:** Line 49 reports `< 0.813r`, `< 0.668r`, whereas n96_c01 and n96_c02 record `base_residual_ratio_upper = 0.8229`, `0.6770` (endpoints `0.9113/0.9579`, `0.7543/0.8817`).
- **Pool bias:** In roi_coupling.json, all 35 witnesses come from 48 Ω (`16` at 6×6, `17` at 8×8, `2` at 12×12, **0 at 16×16**), and every pass reuses only its source `q`. Listing `16×16` (`0/8`) beside `6×6`/`8×8` conflates empty pool with physical separation.

**2. Nonnegative Norm-Bound Fix (Confirmed)**
In attainable.py, `q0.nonnegative_part().sqrt()` and `q1.nonnegative_part().sqrt()` are valid: exact squared-modulus sums are $\ge 0$, so intersecting outward Arb enclosures with $[0,\infty)$ preserves certified upper bounds and removes round-3 `NaN`s (`n96_c01` ROI 2 becomes `0.1537`). (`quadratic_ratio_upper` takes midpoint `float(worst/r2)`, though `full_segment_strict` strictly checks `worst < 0`.)

**3. Physical Correction for the Report**
**Verified:** In roi_coupling_cells.csv, `noise_radius` is **identical across 16/48/144 Ω** per `(case, q)` (scaling strictly as `q/8`). At fixed 1 W power and 50 Ω termination, raising $\Omega_{\text{rad}}$ increases coupling $g^2=\Omega_{\text{rad}}$ and scattered amplitude (`q=12,16` witness distances scale `0.40–0.66×` at 16 Ω, `1.72–2.66×` at 144 Ω vs 48 Ω). Thus 144 Ω's `0/8` drop confounds mutual coupling with an SNR jump under a fixed noise floor. **Correction:** state this fixed noise floor explicitly (revising line 80); **conjecture:** normalizing relative SNR or refitting at 144 Ω would largely restore ambiguity.
