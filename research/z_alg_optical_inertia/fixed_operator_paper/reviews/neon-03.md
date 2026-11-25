# Referee Report: *Material-set inversion with fixed constitutive traces and conditional scattering operators*

**Editorial Decision: Major Revision**

The manuscript addresses a genuine, well-formulated question in nonlinear finite-aperture inverse scattering: certifying exclusion of entire dielectric material classes and predicting prospective scattering bounds while retaining multiple scattering. The core algebraic identities are mathematically sound, the interval-arithmetic verification discipline is unusually thorough, and the authors deserve credit for honestly reporting negative outcomes (shared unknown-cell rows, isolated single-cell peaks, and laboratory loaded-system margins). However, direct audit of the compiled manuscript ([paper-text.txt]), its LaTeX sources, and the accompanying code and JSON/CSV outputs reveals **four concrete scientific and methodological flaws** that must be corrected before publication in *Inverse Problems*.

---

## 1. Independently Verified Scientific Claims

1. **Exact fixed-trace inertia invariance (Theorem 2.3) and known-vacuum mechanism (Section 2.3, Section 5.1).**
   - **Algebraic verification** ([material.tex]): Because $\operatorname{ran} W \perp \operatorname{ran} T_R^* + \mathcal{E}_{\le 0}(K_b)$, we have $T_R W = 0$ and $C_+ = W^* K_b W \succ 0$; eliminating $C_+$ commutes with all updates $D_\theta$, port substitutions $Z$, and thresholds $t$, shifting inertia by exactly $(n-d,0,0)$ without invertibility assumptions on $K_b$, $K_\theta$, or $D_\theta$.
   - **Data verification**: In Experiment A ([summary.json], [compressed.py]), full inertia (C) and fixed-trace reduction (D) agree on all 2688/2688 hypotheses (`EXCLUDED: 1116`, `WITNESS_FOUND: 1488`, `UNRESOLVED: 84`, zero false exclusions). The mechanism check ([known_zero_control.json]) confirms that exact global vacuum elimination reproduces all 8 added exclusions from Method E and yields 20 more (`1144` excluded, `56` unresolved).

2. **Source-side receiver index identity (Proposition 3.2) and generalized-eigenbeam energy infimum (Corollary 3.3).**
   - **Algebraic verification** ([response.tex], [energy.tex]): Evaluating the $(n+s+p)$-dimensional form on $(x,v,y)$ by translating $y$ gives $p + n_+(\mathcal{P}_o)$, while eliminating invertible $x$ then $v$ gives $(n-\nu) + n_+(P_o) + n_+(R_y)$, proving (3.6) and (3.7) exactly. Assuming a material in the class admits a physical solution across the full dictionary $\mathcal{B}$ (giving $n_+(Q_s) \le q$), the trace/interlacing proof of $E_{\rm inf} = \sigma^2 \sum_{i=1}^q \beta_i^{-1}$ is valid for fixed $\Omega = \sigma^2 I_s$. *(Minor precision fix: state "physical solution on the full dictionary" explicitly in Corollary 3.3 rather than only "physically nonempty".)*
   - **Data verification**: All 36 prefix rows in [source_index.json] match $n_+(R_y) = p - h_{\rm src} + \nu - n_+(P_o)$ against direct interval inertia of $R_y$. In [energy.json], `solid_half_1.5` certifies $q = 11$ positive dictionary modes and brackets $E_{\rm inf}$ in `[0.590827932, 0.590827945]` (disk) and `[3.744299996, 3.744300072]` (ring), while `solid_half_2.5` certifies $q = 13$ with at most $11$ positive dictionary modes (`[6, 11]`), ruling out $h_{\rm src} = 0$ at any energy. *(Reproducibility note: [energy.py] on disk does not match the schema/contents of [energy.json]—e.g., it outputs a bare list without `count_probes`, `physical_source_graph`, or the ring radius; the archived script must be synchronized with the reported JSON.)*

---

## 2. Specific Corrections Required (Five Focal Areas)

### A. Rank and Fixed-Trace Reduction Mismatch in Section 5.2 (Factual & Algorithmic Error)
- **What the paper claims** ([experiments.tex], [paper-text.txt]):
  > *"For a separate timing comparison, we used 100 and 1000 different material boxes with the same set of changed cells. The full and reduced threshold matrices had orders 151 and 25."*
- **What the data and code actually do**:
  - Direct inspection of [benchmark_100.json], [benchmark_1000.json], and [benchmark.py] shows the boxes do **not** share the same set of changed cells: they span **53 distinct 9-cell supports** (out of 100) and **64 distinct 9-cell supports** (out of 1000).
  - More importantly, `FixedKernel` in [certify.py] does **not** implement Theorem 2.3 (which fixes $R$, permits indefinite/singular $K_b$ and singular $D_\theta$, and shifts inertia by $(n-d,0,0)$). Instead, `FixedKernel` implements a **positive-definite dual-border Woodbury reduction**: it requires $K_b \succ 0$ (via Frobenius floor $w(1-b\|G\|_F) > 0$), precomputes the full $2n\times 2n$ Schur Gram $T K_b^{-1} T^*$, slices $2|R|\times 2|R|$ blocks for *moving* supports $R$, and adds $-D_\theta^{-1}$, which **requires invertible $D_\theta$** ($a = -w(\ell+u-b)/2 \neq 0$, raising `Singular local D difference` whenever $\ell+u=b$) and shifts inertia by $(n-|R|, -|R|, 0)$.
- **Required correction**: Fix the misstatement in Section 5.2, explicitly distinguish this moving-support PD dual-border algorithm (with its $K_b \succ 0$ and invertible-$D_\theta$ restrictions) from Theorem 2.3, and reconcile Section 4.1 ([space.tex]), which asserts that fixed-trace cost depends on the *union* of changed cells across tests (true for Theorem 2.3, false for `FixedKernel`).

### B. Clean-vs-Noisy Ellipsoids and Absence of Indefinite Prospective Demonstrations (Section 3.1, 3.4, 5.3)
- **Verified sound theory**: Theorem 3.1 ([response.tex]) rightly bounds the **clean** operator $\mathcal{F}_X$ without double-counting old noise, and Section 3.4 rightly treats the noisy future set as the Minkowski sum $\mathcal{E}_i^{\rm clean}(a) + B_2(0, r_a)$ (adding $r_a \|w\|_2$ to the support function rather than falsely asserting $\diag(\Omega_{\rm old}, \Omega_{\rm new})$ majorizes cross-terms $E_{\rm old}^* E_{\rm new}$).
- **Verified experimental gaps**:
  1. **Zero indefinite ($\nu > 0$) full-wave prospective demonstrations**: In both Experiment C runs ([results.json], [diagnostic/results.json]), **all 192/192** candidate weight tests used $K \succ 0$ ($\nu = 0$, $P_o \prec 0$) on perturbative boxes of half-width $0.0008$ or $0.0032$ ([diagnostic/registration.json]), while all 96/96 indefinite trials on $[0,3]^{16}$ suffered an index deficit. The *only* $\nu > 0$ saturated conditional ellipsoid instantiated in the manuscript is the 3-cell algebraic control ([virtual/result.json], $\nu=1$, $G = 4G_{\rm orig}$), which a uniform $\nu=0$ weight already excludes at $1.96r$. Consequently, Experiment C only exercises classical positive-definite S-procedure conditioning, not the indefinite Pick-extension mechanism of Theorem 3.1.
  2. **Suboptimal trial direction and conservative batch noise norm**: Section 3.4 ([response.tex]) emphasizes that maximizing (3.9) over $\|w\|\le 1$ is a receiver-space SOCP attaining the exact Euclidean distance between ellipsoids. However, `split_query` in [design.py] **never solves the SOCP**—it only tests the single center-to-center vector $W \propto C_j A_\star - C_i A_\star$. Furthermore, for multi-beam batches ($m_\star = 2, 4$) under column-wise noise $r_{\rm col}$, `report_support` ([design.py]) penalizes with $\sqrt{m_\star}\,r_{\rm col}\|W\|_F$ instead of the exact column-wise dual norm $r_{\rm col}\sum_{k=1}^{m_\star}\|w_k\|_2 = r_{\rm col}\|W\|_{2,1}$. Thus the reported 4-beam failure ("Only 345 of 1260 four-beam fixed-total actions", 0 edges on 2 of 6 targets) reflects both energy splitting and an unoptimized $W$ / loose Frobenius majorant; Section 5.3 must state this explicitly (or solve the SOCP with $\|W\|_{2,1}$).

### C. The Occupied Native Receiver Counterexample Is Algebraically Trivial (Abstract, Section 5.3)
- **Verified discrepancy between abstract and occupied example**:
  - The abstract ([abstract.tex]) highlights that the conditional bound *"can remain unbounded in seven receiving directions even after saturation on a native full-wave array."* That 7-direction case comes from the `solid_half_2.5` disk prefix test ([source_index.json], $s=16, 32$, $\nu = k_o = 13$, $h_{\rm src} = 7$), where the measured low-peak disk is **incompatible** with the solid class ([native_replay.json]: excluded at $>1.7284r$).
  - To remove the non-occupancy ambiguity, [native.tex] ([compatible_receiver.json]) switches to **`solid_half_1.5`** with $\mu = 2^{-13}$, $\nu = n_-(K) = 0$, and **only $s = 2$ beams** (`source_indices = [0, 16]`), obtaining 23 positive and **9** negative receiver directions.
- **Why this does not resolve the nontrivial obstruction**:
  - Because $\nu = 0$, saturation ($n_+(P_o) = 0 = \nu$) holds **automatically** for any compatible data ($P_o \prec 0$).
  - Because $q = n_-(K_{\rm mat}) = 11$ and $s = 2$, the deficit $h_{\rm src} \ge q - s = 9$ is **trivially forced by underdetermined beam count $s < q$ alone** (even at zero noise; indeed at $s = 0$, $R_y = I_p - S K^{-1} S^*$ already has $q - \nu = 11$ negative modes simply because subtracting $S^* S$ removes the sole stabilization of $K_{\rm mat}$).
  - Meanwhile, [energy.json] and [source_index.json] prove that `solid_half_1.5` has **zero dictionary blindness** (`positive_dictionary_modes = 11 = q`, $h_{\rm src} = 0$ at $s = 16, 32$).
- **Required correction**: Tone down the claim in [native.tex] and the abstract, or—far better—replace it with an **occupied** test inside `solid_half_2.5` at $s = 16, 32$ ($s \ge q = 13$, genuine dictionary blindness $h_{\rm src} \ge 2$) to establish whether an actually compatible interior material can saturate $n_+(P_o) = \nu = 13$ while $R_y$ remains indefinite.

### D. Scope of Finite-Multiplier Comparisons in Experiment A (Table 1, Section 5.1)
- **Verified flaw in Method B ("Positive state form") vs Method C ("Indefinite state form")**:
  - Cross-tabulating [all_nodes.csv] with all 24 verified inertia logs in `experiments/A/results/replay_inertia/pair_01.json`–`pair_24.json` shows that Method C beats Method B on **191 boxes** (`1116` vs `925`, `B&!C = 0`).
  - **28 of those 191 boxes (14.7% overall, and 16 of 32 = 50% in the amplitude pairs!) were excluded by Method C using a strictly positive-definite state matrix $\nu = 0$ ($K \succ 0$)!** (Breakdown of `C&!B`: amplitude `16` at $\nu=0$, `16` at $\nu>0$; structure `12` at $\nu=0$, `147` at $\nu>0$.)
  - Inspecting [pair_01.json] (all 4 `C&!B` rows have $\nu = 0$, weights `[0.0078125, 0.25, 0.0]`, C ratio $> 2.02r$) and [baseline_psd.py] reveals why: for $\nu = 0$, Method B only tests $\operatorname{tr}(H)$ (including negative eigenvalues!) and diagonal sums $\sum \max(0, (U^* H U)_{kk})$ in two fixed bases ($U = I_s$ and Hadamard), whereas Method C at $\nu = 0$ computes $\mathcal{T}_0(H) = \sum_{i=1}^s \max(\lambda_i(H), 0) = \operatorname{tr}(H_+)$—which is **still a strictly positive-definite ($K \succ 0$) dual bound** ($D_{\rm res} \succeq H, D_{\rm res} \succeq 0 \implies \operatorname{tr} D_{\rm res} \ge \operatorname{tr}(P_+ H)$)!
  - Even on the remaining **163 boxes** where Method C used $\nu > 0$, Method B only does a 1-step coordinate search (`factor in (0.2, 5.0)`) rather than solving a convex SDP over $(w_i, \lambda_i)$ subject to $K(w,\lambda) \succeq 0$.
- **Required correction**: Upgrade Method B to include $\operatorname{tr}(H_+)$ for every $\nu = 0$ trial so that the B-vs-C gap cleanly isolates indefinite state blocks ($\nu > 0$) rather than confounding them with an artificial basis restriction (which currently doubles the apparent amplitude-pair gain of C over B). Also benchmark at least a representative subset of the 163 $\nu>0$ boxes against a real convex SDP over $K(w,\lambda) \succeq 0$ to distinguish genuine indefinite-state power from finite-bank suboptimality. *(Minor counting fix in [experiments.tex]: "In 119 distinct tests" should read "In 119 distinct material boxes (192 data–box tests)" to match `certified_neumann_divergence_examples: 192` in [summary.json].)*

### E. Original Contribution vs Known Results
- **Genuine contributions**:
  1. Integrating indefinite-state local conservation/sector forms (extending $K \succ 0$ electromagnetic bounds [5–8]) with exact singular-safe fixed-trace condensation (Theorem 2.3) and template-based core-cover set inversion on full-wave 2D TM models (including the 96×96 finite-aperture map).
  2. The **source-side receiver index formula** (Proposition 3.2) and **generalized-eigenbeam energy infimum** (Corollary 3.3), cleanly separating data-dependent saturation from a priori source/dictionary blindness.
- **Known/incremental elements requiring clearer delimitation**:
  - Propositions 2.1–2.2 rest directly on classical Schur inertia and Weyl/Fan interlacing [10–12]; Theorem 2.3 generalizes static condensation [25] specifically via retaining $\mathcal{E}_{\le 0}(K_b)$.
  - Theorem 3.1 is a finite-dimensional Pontryagin/Pick negative-square Schur extension [26], and Section 3.4 mirrors guaranteed active diagnosis via convex set separation [27]. Because Experiment C never leaves $\nu = 0$, the paper has **not yet demonstrated** that the indefinite ($\nu > 0$) extension provides an active-illumination advantage over standard positive-definite S-procedure conditioning.

---

## 3. Concrete Algorithmic & Experimental Proposals (Conjectures / Research Directions)

Clearly distinguished from the verified findings above, two focused additions could substantially strengthen both the scientific claim and algorithmic utility of the revised paper:

1. **Test genuine occupied $\nu > 0$ conditional ellipsoids (decisive experiment, conjecture to test):**
   - On the native 96×96 grid, test an interior witness of `solid_half_2.5` at $s = 16, 32$ to determine whether occupied saturation ($n_+(P_o) = \nu = 13$) and dictionary blindness ($h_{\rm src} \ge 2$) coexist.
   - On a compact 16–25 cell model with wider contrast intervals (where $K$ is indefinite on the tightest multipliers), test whether a $\nu > 0$ saturated conditional ellipsoid ($n_+(P_o) = \nu > 0$, $h_{\rm src} = 0$) separates competing classes when every $\nu = 0$ weight fails. If yes, Theorem 3.1's indefinite machinery becomes empirically load-bearing; if not, the paper should honestly conclude that $\nu > 0$ is primarily an exclusion tool while $\nu = 0$ dominates prospective bounding.
2. **Algorithmic upgrades grounded in the audit (proposed algorithms):**
   - **SOCP + exact $\|W\|_{2,1}$ + multi-multiplier ellipsoid intersection for active design**: Replace the single center-to-center vector and Frobenius noise bound in [design.py] with a small receiver-space SOCP maximizing $\operatorname{Re}\operatorname{tr}(W^*(C_j-C_i)A_\star) - \sum_{c\in\{i,j\}}\|R_c^{-1/2}W(A_\star^*\Theta_c A_\star)^{1/2}\|_* - 2r_{\rm col}\sum_{k=1}^{m_\star}\|w_k\|_2$, and optionally intersect clean ellipsoids across multiple weights $(w^{(r)},\lambda^{(r)})$ simultaneously (turning Section 3.4's virtual-beam intersection into a practical multi-multiplier separation cone).
   - **Indefinite-safe moving-support reduction**: Combine `FixedKernel`'s moving-support efficiency with Theorem 2.3's generality by eliminating only the positive subspace orthogonal to $\mathcal{E}_{\le 0}(K_b)$ once globally, retaining the small global $\mathcal{E}_{\le 0}(K_b)$ block, and fetching moving $2|R|\times 2|R|$ trace subblocks via a singular-safe $2\times 2$ factored border (avoiding $-D_\theta^{-1}$ when $\ell+u=b$).

---

## Summary of Decision Conditions (**Major Revision**)
To warrant acceptance in *Inverse Problems*, the revision should:
1. Correct Section 5.2 and Section 4.1 to accurately report `FixedKernel`'s moving-support PD dual-border algorithm (53/64 distinct supports) separately from Theorem 2.3.
2. Fix Method B in Experiment A to include $\operatorname{tr}(H_+)$ for $\nu = 0$ trials (removing the 28-box positive-definite confounding) and add an SDP spot-check on the remaining $\nu > 0$ wins.
3. Replace or properly qualify the trivial $s=2, \nu=0$ occupied native counterexample (`solid_half_1.5`) and align the abstract's 7-direction claim with occupancy status.
4. Explicitly disclose in Section 5.3 that Experiment C uses $\nu = 0$ exclusively, uses a naive center-to-center trial direction rather than the SOCP, and uses a conservative Frobenius batch noise bound—or rerun Experiment C with the SOCP and exact $\|W\|_{2,1}$ support, plus sync [energy.py] with [energy.json].
