# E0 cached likelihood audit

Read-only replay: 2835 candidate-array rows, 30068 distinct-material pairs under their identical observation and sigma. New states/RHS/LU: 0/0/0.

A29 pools store `L`, a background-linear response, not full Maxwell. Full likelihood remains null for 2100 rows; `POOL_LINEAR_DIAGNOSTICS.csv` is explicitly a different approximation. N3 exact and T endpoint predictions are read from paid full caches. No forward solve was used to fill a missing prediction.

Point likelihood is `ell=||prediction-y||²/sigma²`, equivalently half the 512-channel whitened squared residual. `230` same-y pairs have |Delta ell|≤log(10); only `0` pass necessary phase/range at both endpoints. This does not establish P0_EXACT membership or posterior basin mass.

The original TRAIN cache supplies 4032 distinct-parent pairs at amplitude 1.2. A single frozen median development sigma per background is applied to both endpoints; nearest d_pred=27.1054307. Structural qualification must use separately frozen TRAIN calibration. No pair is declared a competing posterior mode.

Historical truth noise outside tau: 4/32. These are offline noise realizations, not physically illegal materials.
