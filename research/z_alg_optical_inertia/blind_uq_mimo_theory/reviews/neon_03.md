# Independent review 03

The following text is the reviewer’s dated snapshot. Numerical updates and the authors’ responses appear in [the response](RESPONSE_ZH.md). Local workspace addresses have been removed.

### Targeted Critique of `round3`

**1. Closed interval vs. continuous path.**
*Verified:* In region_intervals.csv and region_summary.json, 22/28 ROIs attain both endpoints $0$ and $4.5$ (exact closed hull $[0,4.5]$), whereas only 5/28 certify a continuous feasible path across $[0,4.5]$—1/22 on $32^2$ versus 4/6 on $96^2$. The theoretical distinction in THEORY_REFINEMENT.md is sound.
*Critique:* The 4/6 continuous-path rate on $96^2$ reflects mesh degeneracy of the peak functional $p_R=\max_{i\in R}\chi_i$: spiking one sub-wavelength pixel ($1/9216$ of the domain) to $4.5$ barely perturbs 400 MHz scattering, whereas $32^2$ single-pixel paths stall at $0.75\text{–}3.75$ in 17/18 exact-hull ROIs. Also, roi_02_inner.json logs `zero_or_base_ratio_upper: NaN` because attainable.py applies Frobenius error `eF` entrywise so `q0` straddles zero before `sqrt()` (the Bernstein check and base ratio $0.107$ still validate the path—a reporting bug, not a proof flaw).

**2. Shared-impedance $w=0.05$ collision vs. all-pixel $w=0.02$ disjointness.**
*Verified:* On the supplied $6\times 6$ witness pair, mechanism_class_separation.json certifies all-pixel $8\times 8$ disjointness at $w=0.02$ under nominal and directional $0.25\,\Omega$ calibration error. Conversely, mechanism_class_collision.json and mechanism_class_collision.npz verify a $w=0.05$ pair colliding ($\le 0.9903r$) under one shared, passive, purely imaginary symmetric perturbation ($D_a=D_b$, $\|D\|_2<0.005$); nominal collision stays unresolved ($1.0033r$ at the 100-iteration cap).
*Critique:* First, $w\in[0.02,0.05]$ is $\le 1.1\%$ of $[0,4.5]$ around one pre-supplied pair, not blind UQ. Second, class_separation.py bounds calibration independently per class ($\|H_+\|_*+\|H_-\|_*$), whereas class_collision.py shows shared-error shrinkage ($0.013r$) is $\sim 9.5\times$ smaller than independent-error shrinkage ($0.124r$). *Conjecture:* Bounding shared-error difference sensitivity $\|H_+-H_-\|_*$ will materially widen the certifiable calibrated separation radius beyond $0.02$.

**3. Alternative-array refits under grid mismatch and held-out errors.**
*Verified:* In summary.csv and cross_audit.csv, fixed A-compatible pairs fail on swapped-arc layouts B and C (0/7), refits satisfy at most one branch on B (3/7, patch-only) and none on C (0/7), and 28/28 refits fail both A ($4.1\text{–}14.9r$) and held-out ports ($11.5\text{–}26.5r$). Meanwhile, grid_transfer.csv shows $32^2\leftrightarrow 96^2$ grid mismatch is $9.1\text{–}15.7r$ on truth ($2.0\text{–}8.1r$ on witnesses).
*Critique:* Refit failures ($1.1\text{–}6.6r$) are local L-BFGS-B outcomes without an unconstrained-refit control or global infeasibility certificate—even `patch` fails 4/7 on B and 7/7 on C (*conjecture:* blind-ROI mismatch or optimizer stagnation). Because grid mismatch exceeds refit residuals and B-feasible refits still miss held-out ports by $19\text{–}22r$, arc-switching breaks the specific A-trained pair inside a fixed discrete operator rather than proving continuum identifiability.

### Practical Implication
Stop certifying $1r$ discrete single-pixel peaks on $32^2$; pivot to **resolution-stable macroscopic functionals** (e.g., coarse-block averages) under a **combined noise-plus-discretization radius**, pairing **shared-calibration difference bounds ($\|H_+-H_-\|_*$)** with **unconstrained swapped-arc refit controls** to separate physical identifiability from mesh artifacts and local-optimizer failure.
