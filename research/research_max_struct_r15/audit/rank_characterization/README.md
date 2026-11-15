# R15 conditional spectral rank characterization

Implemented `../../src/joint_spectral.py` only. Entry points `conditional_four_spaces(G,S,real_material=False)` and `joint_spectral(...)` return four bases, full pre-conditional G/S right frames and singular values, sector eigenframes/eigenvalues, and metadata/counters/timing. CLI accepts `--public FILE [FILE...] --output NEWDIR [--real-material]`; existing outputs are rejected. Uses fixed tau=.1, numerical rank rtol1e-12; no threshold/rank tuning.

Bits name S sector first and conditional G-energy second. Full S right frame is split at sigma>=.1 sigmaSmax. Each sector compresses H=G*G and is split at eigenvalue>=.01 sigmaGmax². Real material stacks real/imaginary operator rows before **all** steps; returned original frames are full pre-conditional frames in that chosen scalar domain. Raw eigenvalues determine membership and are saved; only a separate diagnostic copy clips tiny negative roundoff. Significant negative eigenvalues are retained and flagged. Empty spaces are not padded. Full numerical S-null/G-null dimensions use relative singular-value cutoff1e-12 and are separate from strong/weak thresholds. A zero operator follows the literal >=0 selection rule, explicitly documented in metadata.

This is an S-conditioned G-energy split, **not** simultaneous SVD or exact intersections of independent global strong spaces. Global projector noncommutation is reported as the Frobenius norm of Ps Pg−Pg Ps. Cross-sector coupling is ||Vs* H Vw||F. These diagnostics carry no novelty, contraction or recovery claim.

## Actual authorized public check

Used the existing Gaussian `.venv_nn` Python; all BLAS environment limits1. Only r14_new_0000,0009,0021 were characterized, once each scalar domain. Exact complete G/S equality via `array_equal` grouped all three into one operator group in each invocation: later inputs reuse the computed decomposition. No hashes, model construction, forward solve, optimizer or GT read.

|Domain|U11|U10|U01|U00|Numerical S-null|Numerical G-null|Orthonormality Frobenius|Completeness Frobenius|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|complex current|11|0|9|556|528|0|3.68844e-13|3.68853e-13|
|real material|11|0|11|554|528|0|4.32688e-13|4.32789e-13|

All four dimensions sum576 in both domains. Global strong-G rank19 is distinct from the conditional strong total20/22. Cross-sector G-Gram Frobenius is1.1096211130 complex / .4366000945 real; global projector commutator Frobenius is1.1031586677 /1.0195268211. No significant-negative-eigenvalue issues were reported. Empty U10 was preserved.

Complex batch CPU=.763005s/wall=.804633s, including NPZ reads, exact-equality comparisons, decomposition diagnostics and compressed exports; decomposition CPU=.171697s. Real batch CPU=.410880s/wall=.434637s; decomposition CPU=.119326s; whole-process CPU captured at receipt=.540250s including imports. **First invocation did not capture import-inclusive whole-process CPU**; its original receipt is retained rather than replaced or inferred. Recorded batch CPU total1.173885s; this is not a complete all-process cost. Initial public header inspection/editing CPU was not instrumented. All execution receipts, including per-input cache reuse, remain in each domain directory.

Artifacts: `current_complex_v1/` and `real_material_v1/` each contain RECEIPT, EXECUTION_00/01/02, and operator_group_00/{METADATA.json,SPACES.npz}. SUMMARY.json gives compact exact values. No online solver or prior study was modified. Final mathematical/scientific review remains with parent.
