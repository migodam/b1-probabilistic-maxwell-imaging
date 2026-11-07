# Saved material-column energy audit (zero FW, no GT)

The same ten exposed saved G19 states were compared with ordinary scalar denominators at identical J/Q/m/beta and the exact original saved source metric. No GT was read, no inverse/phase scan or physical solve was run. This tests the residual material-Jacobian column energy only; it does not include the prior Hessian or infer causality for complete trajectories. All ratios are Green / ordinary. The common rho*beta² multiplies both GN diagonals equally and is omitted from the stored diagonal energies.

|Saved state|Qeff energy median ratio|Inverse-column norm² median ratio|Total column energy q10 / median / q90|Lower / upper bound cells|Spearman(distance to box,total ratio)|Weighted excitation distance|
|---|---:|---:|---:|---:|---:|---:|
|r10_twod_0001_s0|0.998055|1.011136|0.995828 / 1.009196 / 1.035163|34 / 0|-0.134454|3.038312|
|r10_twod_0001_s3|1.000726|1.044205|1.021144 / 1.047472 / 1.071386|355 / 0|0.487136|1.988079|
|r10_twod_0004_s0|1.055115|1.003525|0.965381 / 1.060206 / 1.334806|210 / 9|-0.262174|1.294628|
|r10_twod_0004_s3|1.032367|1.020118|1.005050 / 1.056038 / 1.175565|494 / 43|0.010738|1.512621|
|r10_twod_0022_s0|1.000346|1.009783|0.984949 / 1.011292 / 1.051966|10 / 7|-0.086301|2.070733|
|r10_twod_0022_s3|1.006320|1.024863|1.019625 / 1.035477 / 1.099745|467 / 30|-0.123271|1.418640|
|r10_conditional_0000_s0|0.997614|1.011339|0.982087 / 1.009033 / 1.040895|13 / 52|0.122050|1.959314|
|r10_conditional_0000_s3|1.000476|1.058288|1.031792 / 1.056179 / 1.084204|336 / 0|0.539071|1.129961|
|r10_conditional_0003_s0|1.031335|1.002260|0.930968 / 1.037722 / 1.214317|326 / 0|-0.182134|1.388039|
|r10_conditional_0003_s3|1.061577|1.005278|1.010246 / 1.067730 / 1.181127|530 / 8|-0.240076|1.210081|

The shell does not show generally weaker columns relative to ordinary at either saved endpoint: its Qeff energy median ratios are1.0551/1.0324, inverse-column norm² ratios1.0035/1.0201, and total ratios1.0602/1.0560. Ring total ratios are1.0092/1.0475. Across all ten states, total median ratios range1.0090–1.0677. Individual reductions remain in the full distribution: the lowest total ratio is0.83749 (shell stage0), while its q10 is0.96538. These measurements do not support a broad relative material-excitation attenuation specifically at the failed shell states. They do not establish an alternative causal mechanism or rule out directional/correlated losses that a diagonal statistic cannot capture.

At shell stage3,494 cells are exactly at the lower bound and43 at the upper bound; total-energy-ratio medians there are1.05775 and1.03122, versus1.06240 on39 interior cells. Ring stage3 has355 lower-bound cells and221 interior cells, with medians1.03842 and1.06220. The table reports continuous distance-to-box rank correlations as well; no new success threshold or near-box success rule was introduced. Bound membership uses only1e-12 numerical tolerance around the original0.02/20 limits.

Weighted excitation distance is ||E+beta J||_Omega / ||E||_Omega. Its actual saved values range1.12996–3.03831. These are distances from the fullK=G auxiliary excitation-degeneracy condition, not a test of that identity. The actual saved operator here is projected K=PGP, so the fullK exact degeneracy statement is not being numerically asserted for PGP. The condition concerns an auxiliary state, may violate data/APCS/DPCS, and is not physical Maxwell nonuniqueness.

Neither ordinary nor Green has an exactly zero Qeff-energy column in these states (total zero count 0). Full min/q10/median/q90/max for each of the three ratios are in [QUANTILES.csv](QUANTILES.csv); complete bound-group and raw energy summaries are in [ROWS.json](ROWS.json). Per-state NPZ files store complex Qeff, energies, ratios, source weights, material and box distances. [Frozen specification/source paths](FROZEN_DIAGNOSTIC.json) and [reproducible script](../../../experiments/diagnose_material_column_energy.py) identify every input.

Measured complete-child CPU: 0.793096s. Paid auxiliary work:10 rank19 small LU (dimension sum190),5880 small triangular complex RHS (10×[12 saved defect+576 inverse columns]),10 PGP actions/120 action RHS,40 source-metric products. Physical LU/FW/adjoint:0. The576 identity columns are auxiliary denominator solves, never576 physical LU. A preceding header-only schema inspection was not externally timed (unknown CPU, not imputed0); no numerical algorithm check or identity experiment was added.
