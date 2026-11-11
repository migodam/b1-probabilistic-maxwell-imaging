# External measured Fresnel and Austria anchors

Actual inherited R8 measured source files and32² Austria observations were reused. All12 predictions were sealed before heldout/Austria truth scoring. No original measured raw files were copied into this deliverable.

Fresnel uses original2001 single/twin dielectric measurements at2GHz,24²/0.16m square,27 training sources with39 actual receiver channels each,711 heldout channels, source multipole illumination fitted only to empty-target measurements. Exact source-specific receive/adjoint masks are retained. Empty-field calibration validation medians are about.02515/.02408; unknown absolute units, receiver-horn response and noise covariance remain unresolved. The fixed .15||Ytrain|| ball is explicitly MODEL_DISCREPANCY_DIAGNOSTIC, not a known-noise feasibility certificate. No pixelGT or material superiority claim is available.

Austria reuses actual R8 austria_eps4_400MHz PUBLIC_DATA.npz:150² forward generation,32² inversion,20 plane-wave sources,40 receiver channels,10% fixed-norm complex noise. Material is real chi=epsilon_r−1=3 on the target, not epsilon4 used as chi. The predetermined .12||Y|| ball is a discrepancy diagnostic, not iid CN99. The unchanged R8 equal-area disk-pulse acquisition/Maxwell backend is used; no new backend or full-author-reproduction claim.

All arms use the same real[0,20] box, inherited smoothed-TV/squared-difference prior, and common scaled L-BFGS-B corrector max1000/60CPU. Fresnel starts.1, Austria starts0. Direct uses up to2000/180CPU physical auxiliary iterations. Dual reuses actual R8 mature sourcewise CIE-Dual with bounded300 total PR iterations; Fresnel beta1/.2/.2/.2,mf3/6/10/12, Austria beta8/2/.5,mf5/8/10,currentbeta8. AP/DR use groupedv2,scale6,max5000/audit50/180CPU. These are declared bounded-work comparisons, not equal CPU.

| Case | Arm | Train residual | Heldout residual / Austria NRMSE | Phi | Paid CPU s |
|---|---|---:|---:|---:|---:|
|dielTM_dec8f|direct|0.08764|0.09150|0.00384105|20.77|
|dielTM_dec8f|Dual|0.08764|0.09150|0.00384105|19.93|
|dielTM_dec8f|AP|0.08764|0.09150|0.00384108|42.65|
|dielTM_dec8f|DR|0.08767|0.09152|0.00384328|28.32|
|twodielTM_8f|direct|0.09137|0.09857|0.00417549|22.06|
|twodielTM_8f|Dual|0.09138|0.09856|0.00417583|14.66|
|twodielTM_8f|AP|0.09145|0.09868|0.00418255|28.18|
|twodielTM_8f|DR|0.09138|0.09855|0.00417578|32.74|
|austria_eps4_400MHz|direct|0.54568|0.98656|0.148885|14.83|
|austria_eps4_400MHz|Dual|0.10726|1.41964|0.00575969|31.44|
|austria_eps4_400MHz|AP|0.54567|0.98941|0.148882|60.91|
|austria_eps4_400MHz|DR|0.12267|2.13198|0.00753476|82.85|

All bad auxiliary and corrected images are retained in EXTERNAL_MATERIALS.png and per-arm derived arrays. FRESNEL_HELDOUT_PREDICTIONS.png shows prediction against heldout measurements without fitted gain. Receipt corrector stops/costs distinguish budget caps from numerical convergence. No ball radius, calibration, prior or material bounds were selected by heldout/GT. Physical affine reachability of a data ball alone is not Maxwell-material feasibility.

Source files/API: sibling r8_worktree/research/research_max_struct_r8/data_local/fresnel/{dielTM_dec8f,twodielTM_8f}.exp; sibling RESULTS/paper_benchmarks/austria_eps4_400MHz/PUBLIC_DATA.npz and DATA_RECEIPT.json; R8 src/fresnel_adapter.py, fresnel_solver.py, paper_physics.py, paper_cie_fft.py; current R12 src/field_current.py groupedv2. Context diagnostics read from R8 FRESNEL_REPORT.md and INCIDENT_FIT_AUDIT/REFERENCE_POSE records; no claim to have read all previous reports.

Metadata correction: actual inherited Austria forward grid is150² (inverse32²). See METADATA_AMENDMENT.json; original frozen metadata remains preserved as history.
