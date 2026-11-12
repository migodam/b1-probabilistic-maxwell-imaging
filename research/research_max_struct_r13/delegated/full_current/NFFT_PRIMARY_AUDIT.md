# Original NFFT-SOM qualification audit

Verdict: **NOT_ORIGINAL_NFFT_SOM**. The implemented full-current G19 control is a matched current-space transplant. No independent mature original NFFT-SOM benchmark has been implemented or claimed.

Primary source inspected directly: [Wei and Chen, 2018, two-page APS/URSI paper](https://usncursi.org/archive/aps-ursi/2018/papers/1927.pdf). Its equations (1)–(6) establish ordinary LSIE, a determined-current offset plus complete Fourier correction, and the normalized objective. The text specifies conjugate-gradient optimization and displays results at iteration 300, but does not supply a full algorithm listing. It also distinguishes complete one-dimensional Fourier corrections from earlier two-dimensional TSOM columns.

| Item | Primary-paper content | Actual inherited implementation | Qualification consequence |
|---|---|---|---|
| Current family | Eq(3): `Jp=Jp_s+F alpha_p`, complete Fourier correction | full-current G19 retains JD and complete Fourier correction | Same admissible current family; this alone does not replicate optimization |
| Current residual | Eq(5): ordinary LSIE `J-m(E+GJ)` | G19 current block uses scalar CIE contraction at its declared beta | Different objective; setting a label to NFFT cannot remove the difference |
| Source normalization | Eq(6) divides the state residual by determined-current energy `||Jp_s||²`; the data denominator is printed `|V_p|²` | inherited CIE/CSI default state weights are `1/||Einc_p||²`, data weights `1/||Y_p||²` | State weights definitively differ; the conference text does not explicitly define the printed V symbol in the inspected formulation |
| Optimization details | Text names conjugate gradient | inherited code uses specific alternating one-step sourcewise/global PR+, exact current line search, shared-real material readout | The two-page paper does not establish that this is its PR variant, block scheduling or material step |
| Initialization | JD offset specified; material/coefficient initialization is not fully described | inherited CSI sets JD=0, coefficient=current and starts current at zero; G19 physically warms through LS | Full-space equivalence does not make these initial trajectories or coefficient stopping equivalent |
| Stopping | reconstruction figure shows iteration 300 | inherited G19 frozen budgets/stages; historical CIE coefficient stopping | A displayed iteration number is insufficient to reconstruct original stopping |
| Material adaptation | Contrast linked to permittivity; no common TV prior or this project's box is stated | project uses declared shared-real material, prior and box | These must be named as adaptations, not attributed silently to the paper |

The existing dense `research/research_fbg_r3_experiments/src/cie_tsom_baseline.py::run_cie(method='csi')` provides an ordinary LSIE unrestricted-current implementation. It would be a useful explicitly named project CSI control. Its JD=0/direct-current coordinates, incident-field state normalization and material/stopping choices prevent calling it the original NFFT-SOM implementation. Further, `research/research_max_struct_r8/src/cie_fft.py::_run_cie` explicitly rejects `method!='cie'`, so passing `method='csi'` into the matched G19 adapter is not a valid implementation shortcut.

An Eq(3)–(6) adaptation could reuse existing G/E/S, Fourier actions, retained adjoints and sourcewise PR code without a new physical backend. It would still need explicitly selected normalization, offset initialization, material rule, CG variant, stopping and real/box/prior adaptations. Since those choices are not fully determined by this primary source, adding such a row would not satisfy a request for a mature original-method benchmark. Parent decision: leave the gap explicit rather than weaken that qualification.

For a future independently qualified implementation, obtain original author code or the detailed referenced SOM update source and verify the missing algorithmic choices. Do not treat that follow-up as completed by the present source audit.
