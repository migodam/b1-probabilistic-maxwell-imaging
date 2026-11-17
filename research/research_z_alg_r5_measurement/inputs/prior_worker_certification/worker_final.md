Implemented and verified the directed-interval local Krawczyk certifier.

- [certification.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_certification/certification.py): JSON API, 64/128/256-bit intervals, fixed radius schedule, strict inclusion and row-sum contraction checks, complete attempt and cost reporting.
- [test_certification.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_certification/test_certification.py): historical R2 4-cell and 6-cell fixtures, including archived root and Jacobian checks.
- [tests.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_certification/tests.json): all four tested roots certified at 64 bits; shifted-target control returned `UNKNOWN` with all nine failed attempts preserved.
- [SUMMARY.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/z_alg_r5_certification/SUMMARY.md): mathematical scope, costs, results, and exclusions.

Final test: `PASS`; 0.327445 CPU seconds and 0.331998 wall seconds. Historical R2 files were not modified.