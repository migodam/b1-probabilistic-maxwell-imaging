### Independent Second Research Review (`fixed-operator-26`)

> **Scope & epistemic separation:** All workspace files were left untouched. Below, **[Verified]** marks statements directly checked from source code, stored JSON/CSV/NPZ data, or independent recalculation; **[Conjecture / Proposal]** marks untested hypotheses and proposed algorithmic or experimental direction changes.

---

### 1. Independent Numerical Verification [Verified]

I independently loaded raw binary64 `G, S, B, Y, radii, boxes` from [acquisition.npz] and recomputed the Hermitian quadratic assembly $\mathcal{A}$, negative inertia $\nu = n_-(K)$, Schur complement $H$, and $jt$-tail lower bound from scratch (pure-Python complex Jacobi eigensolver, no code reused from the repo):

* **Pair 15, condition `[0, 2, 6, 0]`** (`bg = 0.0625`, `width = 5.0`, `q = 4`, radius $r = 0.0604312790$, parameters $(\mu, \text{local}, \text{optical}) = (0.0078125, 0.25, 0.0)$):
  * **21 active cells** (dropping exact vacuum rows `[7, 8, 16, 17]`, where $\ell_i = u_i = 0$): recomputed $\nu = 1$, $\sqrt{jt}/r = \mathbf{2.422004693899}$ (full tail $2.499812$). This matches stored `E_strict_added` (`2.422004693897`) in [A/results/summary.json] and stored `known_zero_control` (`2.422004693890`) in [known_zero_control/pair_15.json] to **11–12 digits**.
  * **Full 25 cells** (Methods C/D at the exact same weights): recomputed $\nu = 1$, $\sqrt{jt}/r = \mathbf{0.518147368152} < 1$ (full tail $0.669324 < 1$), explaining why C/D left it `UNRESOLVED`.
* **Pair 8, condition `[1, 0, 5, 1]`** (`bg = 0.0`, `width = 3.0`, `q = 2`): across all 72 `(mu, local, optical)` grid points, full 25-cell C/D peaks at **$0.574609 < 1$**; dropping only E's 2 vacuum cells `[6, 11]` gives **$1.705689$** at C's best weights (stored E: `1.667010`); dropping all 19 known vacuum cells yields **$16.342277$** (stored `known_zero_control`: `16.336275`).
* **Native source-index consistency:** Across all **36/36 rows** of [native/source_index.json], `32 - deficit_source_only + nu - old_inertia[0] == direct_receiver_inertia[0]` holds with **0 mismatches**.

---

### 2. Unfair Comparisons & Experimental Bottlenecks [Verified]

1. **Major unfair comparison — Experiment A (`E` vs `C`/`D`):**
   * In [A/results/summary.json], `E` adds **+8** strict exclusions over `C`/`D` (`1124` vs `1116`, out of `1224` opposite non-cover boxes).
   * However, in [A/run.py], Method `E` sorts shared rows by `(u - l) / max(1, 0.5*(l + u))`, placing **zero-width known vacuum cells (`l = u = 0`) first** and eliminating $j_i = 0$ identically, whereas `C` and `D` keep those known vacuum cells in the state vector with only a finite quadratic penalty $w_i |j_i|^2$.
   * **Audit of all 8/8 `E_strict_added` boxes** (Pairs 8, 9, 15): **every single selected row** (`[6, 11]`, `[13]`, `[7, 8, 16, 17]`) had exact bounds `(0.0, 0.0)`. **0 / 2688 boxes** were newly excluded by genuinely nonzero-width interval subdivision (while `3469` interval leaves failed `uniform_positive_gap_nonpositive`, [A/results/summary.json]).
   * [A/results/known_zero_control.json] confirms the diagnosis: applying ordinary `C` after eliminating `upp == 0` rows yields **1144** exclusions (**+28** over `C`/`D`, **+20** over `E`), reproduces **8/8** of `E`'s additions with `0` false exclusions, and runs in `27.8 s` total. The reported `E > C/D` gain is therefore an artifact of asymmetric vacuum elimination, not evidence for nonzero-width shared-row charts.

2. **Secondary comparator asymmetry — Experiment A (`B` vs `C`):**
   * In [A/baseline_psd.py], `B` seeds its coordinate search *only* from `C`'s 72-point grid entries where `nu == 0`. On `structure` width ratios `3` and `5`, `B` scores **0/16** exclusions while `C` scores **14/16** and **12/16** ([A/results/summary.json]). While PSD duals genuinely fail when no $w \ge 0$ makes $K \succ 0$, gating `B`'s initialization strictly on `C`'s discrete grid conflates true PSD infeasibility with grid-seeding starvation.

3. **Ceiling effect — Experiment C (`diagnostic`):**
   * Even after widening the box 4× and raising noise to 2% (reducing strictly separating menu actions from `3180/3180` in [C/results.json] to `2016/3180` in [C/diagnostic/results.json]), **post-acquisition wrong-class exclusion saturates at the ceiling**: across all 36 policy-contract-batch targets per policy, `random` achieves **35/36** vs **36/36** for `uniform_gap`, `fisher_nuisance`, and `conditional_ellipsoid`, and `minimum_successful_tested_batch = 1` for **all 4 policies on all 6 targets** under both energy contracts ([C/diagnostic/results.json]). Only the *pre-acquisition certified edge count* separates policies (`conditional_ellipsoid` 34/36, `uniform_gap` 33/36, `fisher_nuisance` 31/36, `random` 25/36).

4. **Scope caveat — Experiment B amortization:**
   * The `27.6×` (100 boxes) and `38.0×` (1000 boxes) float speedups and `3.55 ms/box` vs `324.7 ms/box` sampled direct-Arb timing in [B/results/summary.json] ([B/benchmark.py]) are measured at a **single fixed scalar weight `w` and `optical = 0`**, so multi-weight screening scales setup cost linearly with the number of distinct background weights. Also note that on the three 12×12 objects, `arbitrary_high` gives **0** exclusions and `global_low_peak_excluded = false`, and the 96×96 native run hits `21` `node_budget_deferrals`.

---

### 3. Blind New Paired Study vs. Targeted Mechanism Replay [Verified]

These two evidentiary tiers must stay explicitly separated in reporting:

| Study / Artifact | Category | Valid Inference | Not Valid For |
| :--- | :--- | :--- | :--- |
| Exp A primary arms `A–E` (24 sealed clusters), Exp B primary maps, Exp C main (`results.json`) | **Blind new paired study** (preregistered, sealed truth before decision) | Unbiased prospective performance comparison under declared contracts | Generalizing beyond discrete box/noise contracts |
| Exp A `known_zero_control.json`, Exp A `historical_controls`, Exp C `diagnostic/results.json` (same 6 sealed media with post-hoc 4× width / 2% noise), `native/source_index.json` (archival 96×96 data, conservative $r_{32}^2 I_s$) | **Targeted mechanism replay** (post-hoc / archival diagnostics) | Algebraic identity verification, confounding isolation (e.g., proving `E`'s +8 came from zero rows), failure-mode anatomy | Claiming prospective policy superiority or primary-endpoint rankings |

---

### 4. Theoretical Audit of the Two Mathematical Questions

#### (a) Does $h = n - n_+(\mathcal{P}_{\text{material,old}})$ really give the output index at saturation? [Verified]

* **Algebraic truth (exact under invertible $K, P_o$):** In [response.tex] (`prop:sourceindex`), comparing the two congruences of the $(n+s+p)$ form $q(x,v,y) = [x;v]^*\mathcal{P}_o[x;v] + \|Sx - Yv - y\|^2$ yields, via Haynsworth inertia additivity (requiring $K, P_o$ nonsingular, **not** $K_{\text{mat}}$ or $R_y$):
  $$n_+(R_y) = p - h_{\text{src}} + \nu - n_+(P_o), \qquad n_0(R_y) = n_0(\mathcal{P}_o), \qquad n_-(R_y) = n_-(\mathcal{P}_o) - \nu - n_-(P_o),$$
  where $h_{\text{src}} = n - n_+(\mathcal{P}_o)$. I verified this both across all 36 rows of [native/source_index.json] and on a synthetic test with **singular $K_{\text{mat}}$** (where the Schur formula \eqref{eq:indexsmall} is undefined, yet the full block $\mathcal{P}_o$ in \eqref{eq:sourceform} predicts $n_+(R_y) = 2$ exactly).
* **Exact meaning at saturation ($n_+(P_o) = \nu$):** At saturation, $n_+(R_y) = p - h_{\text{src}}$, so **$h_{\text{src}}$ equals the exact non-positive output index $p - n_+(R_y) = n_-(R_y) + n_0(R_y)$** (strictly $n_-(R_y)$ only when $\mathcal{P}_o$ is nonsingular; if $\mathcal{P}_o$ has a kernel, $n_0(R_y) = n_0(\mathcal{P}_o) > 0$ also contributes to $h_{\text{src}}$).
* **Three essential caveats verified in [native/source_index.json]:**
  1. **Over-saturation ($n_+(P_o) > \nu$, 5 rows):** e.g., `solid_half_1.5`, `disk`, `s = 16, 32` have $h_{\text{src}} = 0$ yet $R_y$ has inertia `[30, 2]` (2 negative eigenvalues) because the class is already excluded by old data.
  2. **Under-saturation ($n_+(P_o) < \nu$, 7 rows):** e.g., `cap_4`, `disk`, `s = 8, 16, 32` have $R_y \succ 0$ (`[32, 0]`) *despite* $h_{\text{src}} \in \{5, 8, 11\} > 0$, because $n_+(P_o) < \nu$ cancels the deficit—here $R_y \succ 0$ is spurious because Theorem 1 does not hold without saturation.
  3. **Singular $P_o$:** If $P_o$ has a nullspace $U_o$, Theorem 1 restricts to $U_o^* Z^* Q Z_\star = 0$ with $P_o^\dagger$, and \eqref{eq:index} does not apply unmodified.

#### (b) Does $Q_s$ with $q$ positive generalized eigenvalues $\beta_i$ and $\Omega = \sigma^2 I$ make $\sigma^2 \sum(1/\beta_i)$ the minimum incident energy infimum to remove $h$? [Verified + Qualification]

Let $Q_s = C_{\mathcal{B}} - N^* K_{\text{mat}}^{-1} N$ ($K_{\text{mat}}$ invertible), $W_{\text{in}} \succ 0$ the incident-energy Gram matrix, $\beta_1 \ge \dots \ge \beta_q > 0$ the $q$ positive generalized eigenvalues of $(Q_s, W_{\text{in}})$, and $k = n_-(K_{\text{mat}})$ the number of negative material modes that must be cancelled to achieve $h_{\text{src}} = 0$. Transform to whitened variables $\widetilde{V} = W_{\text{in}}^{1/2} V$, $\widetilde{Q}_s = W_{\text{in}}^{-1/2} Q_s W_{\text{in}}^{-1/2}$, so total incident energy is $\mathcal{E}(V) = \operatorname{tr}(\widetilde{V}^* \widetilde{V})$.

* **Exact theorem (for fixed $\Omega = \sigma^2 I_s$ independent of $s$):**
  1. **When $k \le q$:** Achieving $n_+(\widetilde{V}^* \widetilde{Q}_s \widetilde{V} - \sigma^2 I_s) \ge k$ requires $s \ge k$ and $\lambda_k^\downarrow(B^{1/2} X B^{1/2}) > \sigma^2$, where $B = \operatorname{diag}(\beta_1,\dots,\beta_q)$ and $X = U_+^* \widetilde{V}\widetilde{V}^* U_+ \succeq 0$. Hence $B^{1/2} X B^{1/2} \succ \sigma^2 P_k$ for a rank-$k$ orthogonal projector $P_k$, so by conjugation and Ky Fan's minimum principle:
     $$\mathcal{E}(V) \ge \operatorname{tr}(X) > \sigma^2 \operatorname{tr}(P_k B^{-1}) \ge \sigma^2 \sum_{i=1}^{k} \frac{1}{\beta_i}.$$
     Conversely, $k$ $W_{\text{in}}$-orthogonal eigen-beams $\tilde{v}_i = \sqrt{(\sigma^2+\varepsilon)/\beta_i}\,\tilde{u}_i$ ($i=1,\dots,k$) achieve $\mathcal{E} \to \sigma^2 \sum_{i=1}^k 1/\beta_i$ as $\varepsilon \downarrow 0$. (Verified numerically on synthetic tests: 4,000 random $V$ stayed $\ge 2.20\times$ above $\sigma^2 \sum_{i=1}^k 1/\beta_i$ for $k=2$ and $\ge 5.53\times$ for $k=3$.)
  2. **Two exact precision catches:**
     * The sum runs over the **top $k = n_-(K_{\text{mat}})$** positive generalized eigenvalues $\beta_1 \ge \dots \ge \beta_k$, **not all $q$** when $q > k$ (matching [native/energy.py], which sums `range(raw)`—note that `energy.json` has not yet been generated).
     * When **$q < k$ (dictionary blindness)**, removing $h_{\text{src}}$ to $0$ is **impossible at any finite energy** ($\inf = +\infty$); $\sigma^2 \sum_{i=1}^q 1/\beta_i$ is only the infimum to reach the irreducible floor $h_{\text{src}} = k - q$.
* **Critical physical caveat on $\Omega = \sigma^2 I_s$:**
  * The formula holds *only* when $\sigma^2$ is **fixed independent of $s$** (a fixed total joint Frobenius noise budget across all $s$ beams).
  * Under **constant per-beam coherent noise** $\|e_j\|_2 \le r_1$ ([response.tex]), the valid deterministic majorant scales with the number of beams: $\Omega = s r_1^2 I_s$. Minimizing over $s \ge k$ then sets $s = k$ and multiplies the infimum by **$k$**: $\mathcal{E}_{\inf} = k r_1^2 \sum_{i=1}^{k} (1/\beta_i)$. Furthermore, [native/energy.py] sets `noise = disk.radius**2` (the archival **32-beam** joint radius squared $r_{32}^2$), which over-penalizes an optimal $k$-beam acquisition by $32/k$ relative to a fresh $k$-beam joint bound $k r_1^2 I_k$.

---

### 5. Proposed Physical Interpretation & Direction-Changing Correction [Conjecture / Proposal]

1. **Physical correction (make exact vacuum elimination mandatory, then pivot shared-row charts to low-rank physical ties):**
   * *Interpretation:* A box cell with $\ell_i = u_i = 0$ is a **known vacuum Dirichlet constraint** ($j_i \equiv 0$, zero scatterer support), not a constitutive interval to be relaxed with finite sector weight $w_i$. Leaving it in $K$ allows unphysical virtual currents in known vacuum and produced the entire `+8` apparent win of `E`.
   * *Correction:* Eliminate all `upp == 0` rows from $(G, S, B)$ **before** building `TraceBank` or running any of `B, C, D, E`. Then, if `E` is pursued further, replace independent per-pixel bisection (which fails exponentially at `3469` uniform-gap leaves) with **physically tied low-rank shared-row constraints** (e.g., $\chi_i = \chi_0$ across a homogeneous inclusion or piecewise-constant/radial shell), where $q$ cells share **1** real scalar parameter instead of $q$ independent intervals.

2. **Direction-changing experiment — Source-Eigenbeam Active Wave-Shaping at the $\mathcal{E}_{\inf}$ Phase Transition:**
   * *Diagnosis:* Experiment `C` currently operates far above the saturation/energy threshold (`random` wins 35/36), whereas `native/source_index.json` (`solid_half_2.5`, `cap_4`) operates *below* it ($h_{\text{src}} \in \{5, 7\}$ even with all 32 plane waves).
   * *Proposal:* Run `native/energy.py` to separate dictionary blindness ($q < k$) from energy-floor deficit ($q \ge k$), and construct a **new blind paired study** right across the exact source-index energy cliff $\mathcal{E} \in [0.5, 4] \times \sigma^2 \sum_{i=1}^k 1/\beta_i$ (using the physically matched $\Omega = s r_1^2 I_s$ scaling). Compare single-angle plane-wave selection (`random`, `uniform_gap`, `fisher_nuisance`) against **coherent eigen-beam synthesis** $v_i \propto \beta_i^{-1/2} u_i$ (the top $k$ generalized eigenvectors of $(Q_s, W_{\text{in}})$).
   * *Expected shift:* Physically, coherent multi-angle wave-shaping injects incident fluence directly into the $k$ unstable modes of $K_{\text{mat}}$, unlocking $h_{\text{src}} = 0$ (and hence bounded conditional ellipsoids $R_y \succ 0$) at substantially lower total incident energy than any single-angle plane-wave menu—creating a non-trivial regime where active physical illumination design genuinely beats random and uniform-angle baselines on post-acquisition survivors, not just pre-acquisition edges.
