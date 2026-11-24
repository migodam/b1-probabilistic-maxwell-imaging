# Frozen weak-scattering controls

The two support functions below are fixed before computing their fields or running the exclusion test. Use the same stored \(96^2\) Green operator \(G\), full-ring receiver map \(S\), and incident matrix \(B\) as the primary study. The cells are centered on the 2 m square. No materials or observations from the primary six-object sample enter this test.

1. **Ultraweak ellipse.** Let \(q=((x+0.21)/0.37)^2+((y-0.16)/0.24)^2\). Set \(\chi=10^{-6}[1-q]_+\).
2. **Modest paired disks.** Let \(q_1=((x+0.29)^2+(y+0.12)^2)/0.32^2\) and \(q_2=((x-0.31)^2+(y-0.19)^2)/0.22^2\). Set \(\chi=0.04\max(\mathbf1_{q_1\le1},0.6\mathbf1_{q_2\le1})\).
3. **Exact first-Born arm.** Reuse the modest paired disks and generate the clean data with \(S D_\chi B\). This control is paired with the modest full-wave case and does not add an independent shape.

For each shape generate the full-VIE current with the native FFT. The exact Born arm has its own deliberately linear generator. Use complex iid Gaussian draws from NumPy PCG64 with \(\texttt{SeedSequence([20261124183, arm])}\), where arm is 1, 2, or 3, respectively. Normalize once per arm to the complex Frobenius clean-data norm. Use the same draw for 1% and 2% within each arm. The radius convention is the realized norm rounded upward plus \(10^{-6}\|Y\|_F\), as in the primary study.

Run the existing data-only Born cone test on each observation at its fixed iteration limit of 64, for the nonnegative unbounded cone and the \([0,4]\) box. Once those runs are finished, use the known generator shape only offline to compute an explicit Born model along a single nonnegative scalar ray. Compute the real scalar projection analytically and clip its coefficient to the legal box for the bounded test. An interval-verified residual upper bound within the noise ball establishes a concrete Born fit. A lower bound beyond the ball excludes that model class. When neither inequality is shown, the class stays undecided. No generic weak-scattering decision follows from three acquisition arms.
