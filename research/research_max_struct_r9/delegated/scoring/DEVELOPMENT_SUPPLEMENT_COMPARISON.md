# Three exposed development follow-ups

All three were publicly sealed before offline grading. Original43 outputs are untouched. Each proposal row includes the complete pre-existing diagnostic cost once; upstream original43 endpoint costs are separate.

| Task | ampNRMSE | OtsuIoU | Original Phi | Inclusive CPU(s) | LU | complex RHS | Stop |
|---|---:|---:|---:|---:|---:|---:|---|
| obs_40ac5b79e8a74ac3__profile_restarted | 0.457854 | 0.226277 | 0.0212010655413 | 15.735366 | 1090 | 26064 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
| obs_b20d5a55677d4eca__curvature_anchor | 0.491521 | 0.294118 | 0.00508886909947 | 11.607191 | 905 | 21720 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |
| obs_b20d5a55677d4eca__curvature_proposal | 0.491533 | 0.296610 | 0.0050888690831 | 11.784921 | 919 | 23148 | CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL |

All three fail the unchanged material recovery rule.40ac restarting a fixed penalty does not reproduce the original profile-continuation recovery (.05807/.96183). B20d negative-curvature proposal and ordinary continuation reach nearly identical erroneous materials; a small Phi reduction is not recovery. No global landscape conclusion is inferred.
