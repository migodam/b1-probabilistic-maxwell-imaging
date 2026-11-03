# R4 independent development scores

**POST_SEAL_PUBLIC_DEVELOPMENT**: all4 outputs were closed before this scoring pass. Parent authorized already-published R2 truth for scoring only. No new physics, GT-ranked selection, or algorithm change.

Unchanged R2 gates: RMSE/12≤.15; prediction-only256-bin Otsu IoU≥.7; original noisy whitened residual≤25.165422197; public material box.

| Endpoint | NRMSE/12 | Otsu IoU | Relative residual | White residual | Joint |
|---|---:|---:|---:|---:|---|
|R2_native (published reference)|0.510695259|0.245901639|0.100734862|in inherited report|FAIL|
|R3_scaled_original_restart (published reference)|0.510717418|0.254098361|0.100733687|in inherited report|FAIL|
|cooperative_clipped|0.510717418|0.254098361|0.100733687|2417.615|FAIL|
|weak_gn_clipped|0.510717418|0.254098361|0.100733687|2417.615|FAIL|
|box_cooperative|0.510717418|0.254098361|0.100733687|2417.615|FAIL|
|tmatrix_completion|0.524685470|0.250000000|0.095863919|2300.740|FAIL|

T-matrix endpoint versus R3 scaled-original: residual 0.100733687 → 0.095863919; NRMSE 0.510717418 → 0.524685470; IoU 0.254098361 → 0.250000000.
Residual decreased: True; material NRMSE decreased: False; IoU increased: False. Original joint gate: FAIL.

These metrics describe one disclosed development endpoint. They do not establish a new-blind recovery, superiority over converged CIE, or a converged optimum. Original R2 blind results remain unchanged. Scoring adds0 LU and0 RHS.
