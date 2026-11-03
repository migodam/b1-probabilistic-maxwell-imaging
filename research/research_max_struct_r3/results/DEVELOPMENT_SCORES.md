# R3 published-case development scores

All 11 predictions were closed before this scoring pass. The material truth is already public R2 development data. No new physical solve. N=RMSE/12; IoU is reconstruction-only256-bin Otsu. Clean-target discrepancy has no Gaussian gate.

|Control|N|IoU|Target relative residual|Original noisy white residual|Original-noisy joint|
|---|---:|---:|---:|---:|---|
|noisy_strong_release|0.490477|0.294118|0.10066194|2415.893|FAIL|
|clean_native|0.389673|0.422680|0.10782185|2587.686|FAIL|
|passive_path|0.817495|0.204545|0.05101706|1224.412|FAIL|
|extra_original|0.510695|0.245902|0.10073486|2417.643|FAIL|
|globaldual_restart|0.510695|0.245902|0.10073486|2417.643|FAIL|
|scaled_original_restart|0.510717|0.254098|0.10073369|2417.615|FAIL|
|onecell_roots|0.510717|0.254098|0.10073369|2417.615|FAIL|
|neighbor_swaps|0.510717|0.254098|0.10073369|2417.615|FAIL|
|all_pair_swaps|0.510717|0.254098|0.10073369|2417.615|FAIL|
|paired_neighbor_rays|0.510717|0.254098|0.10073369|2417.615|FAIL|
|dressed_neighbor_rays|0.510717|0.254098|0.10073369|2417.615|FAIL|

Published native reference: N=0.510695, IoU=0.245902, original-noisy residual=0.10073486.

Lower data residual alone is not a material improvement. These development controls cannot modify the original R2 blind table or establish a new R3 blind rate.
