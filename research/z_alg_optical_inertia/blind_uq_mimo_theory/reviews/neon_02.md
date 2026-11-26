# Independent review 02

The following text is the reviewer’s dated snapshot. Numerical updates and the authors’ responses appear in [the response](RESPONSE_ZH.md). Local workspace addresses have been removed.

### Independent Technical Review of `blind-uq-mimo-27/round2`

**Personally inspected data files:** mechanism_class_separation.json, frequency_coherence_rank.json (plus material_scores32.csv, one_sided_witness.json, mechanism_witness.npz, and mechanism_witness_arb.json).

---

#### 1. 8-Channel "Whole-Class" Proof: Quantifiers, Selection Bias, and Assumptions (Verified)
class_separation.py and mechanism_class_separation.json do **not** prove whole-class separation ($\forall m_a \in C_+,\ \forall m_b \in C_0$ for presence vs. absence of a 15-cell ROI). They certify support-function separation only between two tiny $\ell_\infty$ boxes $B_\infty(m,w)\cap[0,4.5]^{1024}$ centered on a **single pre-selected witness pair** from mechanism_witness.npz (`with_region`: 283 nonzero cells, 15 ROI cells $=0.75$; `without_region`: 268 nonzero cells, ROI $=0$, outside-ROI max difference $1.963$). Because that pair was optimized to collide specifically at $6\times 6$, 400 MHz (`separation_over_2r = 0.0570`), its separation at $8\times 8$ (`2.4053`) reflects **inner-quantifier selection bias** ($\exists$ pair, locally) rather than outer-quantifier whole-class separation. Even locally, prior-model calibration fails by $w=0.01$ (`-0.000113`, just $0.22\%$ of $[0,4.5]$), nominal and directional-calibrated separation fail by $w=0.02$ (`-0.001676` and `-0.003488`), and masked absolute error fails at every tested width down to $w=0.0001$ (`-0.03575`).

Three physical/matrix assumptions further narrow scope:
1. **Same-operator noise & power scaling:** class_separation.py tests deterministic ball disjointness $\|F_a-F_b\|_F > 2r$ under identical operators with $r=\frac{3}{2}q\sigma_{\text{ref}}=0.001517$ ($q=8$), while model.py fixes total power at $1\,\text{W}$ ($\propto 1/\sqrt{q}$ per shot). Because $r\propto q$ while per-shot signal dilutes as $1/\sqrt{q}$, mechanism_witness_arb.json shows center separation/$2r$ *falling* from $2.405$ ($q=8$) to $1.662$ ($q=16$).
2. **Terminal-only calibration model:** `calibration_bounds` and `projected_calibration` model uncertainty solely as a terminal perturbation $(I+KD)^{-1}K$ ($\|D\|_2\le\eta=0.005$) and rely on synthetic readout noise in `finite_full_probe`, omitting element gain/phase/position errors that alter `E0_bare` and `S_bare` directly (unlike `variation["gain"]` in model.py).
3. **2D lossless, real, nondispersive TM disk-pulse VIE** on transmission-only arcs (`[95°,265°]` to `[-85°,85°]`).

---

#### 2. Two-Band Lift & 80-Cell Obstruction (Verified)
frequency_lift.py and frequency_coherence_rank.json rigorously verify that $D=G_{400}-G_{460}$ on 80 scattered cells has full rank 80 (inverse defect $\le 3.4\times 10^{-10}$; $\sigma_{\min}\ge 0.00255$ on $32^2$, $\ge 0.000549$ on $96^2$), forcing `negative_inertia_lower` $=80-2q \ge 2q$ (`tail_forced_zero = true`). Crucially, this is an obstruction to a **specific uniform-weight bilinear Gram lift** (`C = mu*gamma*(I - cap*(G2+G1^*)/2)` with `"full_field_independent_4p5": true` and `"uniform_common_sector_weight": true`), assuming strictly nondispersive real $\chi$. It does **not** prove two-band physics fails: mechanism_witness_arb.json proves the exact $6\times 6$ colliding pair separates cleanly at 460 MHz (`5.710`) and equal-power 400+460 MHz (`2.855`).

---

#### 3. Measured Image Metrics & Valid Practical Inference (Verified)
material_scores32.csv and baseline.json expose two decisive operational facts:
- All CIE baselines (`baseline`, `tighter_dual`, `b2_*`) sit far outside the physical noise ball (`residual_over_radius` $= 7.2\text{–}17.1$; e.g., `7.328` on `n32_c02` stationary because `state_term = 4.96e-4` maps to `22.0%` physical mismatch). Every dark-current variant (`b2_*`) changes `material_relative_L2` and IoU by $<0.001$—no imaging improvement.
- Forcing physical ball feasibility (`b1_wall_control`) **degrades `material_relative_L2` on all 8 cases** (`n32_c01` `0.8379→0.8652`, `n32_c08` `0.7410→0.8078`), while `held_relative_residual` remains `0.47–0.86` (`17–23×` held-out radius).

**Valid practical inference:** Without spatial priors beyond pixelwise $[0,4.5]$, this framework validly certifies *global FoV non-emptiness* (90–97% union in RESULTS_ZH.md) and *data non-identifiability*—**not** safe blind artifact deletion or tight local ROI bounds.

---

#### 4. Concrete Counterexample to Overstated Claims (Verified)
On `n32_c02` (stationary, $8\times 8$, 400 MHz), one_sided_witness.json certifies two full-wave materials inside the **same** $8\times 8$ noise ball: one zeroing `roi_02` (`residual_upper_over_radius = 0.9180 < 1`) and an unconstrained solution (`0.00128 < 1`) with `roi_02` mean `2.526`. Yet material_scores32.json shows `roi_02` has **83.3% ground-truth overlap** (`zero_truth_overlap_components = 0`), so `b1` pruning deletes a true scatterer and worsens `material_relative_L2` (`0.8346→0.8621`). Together with n32_c01.json (43-cell patch $\ge 2.0$ inside 172-cell `roi_01` at `0.5568r` vs. a zero-ROI witness), this directly refutes both whole-class $8\times 8$ ROI separation and treating "zero-feasible" as "artifact."

---

#### 5. One Creative Physical Test to Change the Conclusion (Conjecture)
**Proposal:** Convert the existing fixed two-arc array in model.py from transmission-only $q\times q$ readout to **full-aperture reciprocal $2q\times 2q$ reflection + transmission acquisition with in-situ multi-load termination switching ($Z_L\in\{50\,\Omega,\text{open},\text{short}\}$ on idle/hold ports) at fixed per-shot SNR**.

**Why this can alter the conclusion (conjecture, requiring verification):** First, adding monostatic/bistatic backscatter fills the missing $2k$ Ewald-sphere wedge responsible for `0.52–1.04` held-out residuals and doubles independent column rank from $q$ to $2q$. Second, empty-chamber multi-load switching physically identifies true terminal/coupling errors $\Delta M$ in situ, replacing the unverified terminal-only $\eta=0.005$ scalar and synthetic `finite_full_probe`. Caveat: even full-angle data cannot eliminate subwavelength dark currents under unconstrained pixelwise $[0,4.5]$, so testing this hardware protocol jointly with a physically grounded spatial-smoothness/TV prior is necessary to achieve nontrivial local ROI bounds.
