# 3D fixed-threshold support audit

Existing sealed/scored primary predictions only. No optimizer, forward, new scoring threshold, or recovery label was introduced. Pair distances use the original inherited amplitude and fixed chi>1 masks. Phi values are the already-paid common original vector objective, not new evaluations.

Bracket 00 versus11: amplitude-normalized pair RMSE=0.0114701112205, symmetric relativeL2=0.0912385580699, Linf=0.655276348311; Phi difference=1.84533317443e-07.
Its fixed-threshold IoU changes 0.675→0.769230769231 through 6 cell flips. Maximum distance of either flipping value from chi=1 is 0.3979426727729938; exact cell values and truth labels are in EXACT_SUPPORT_FLIPS.

For this bracket comparison, zero flipping cells have both values within.01 ofchi1, and two have both within.1. The exact margins exceed tiny floating tolerances; a pure precision explanation is unsupported here. Similar scalar NRMSE values do not imply small pairwise material distance.

The fixed support threshold discretizes the differing fields and may amplify their apparent support-score difference. These tables do not establish a different recovery class. No threshold/optimizer was retuned. All60 method pairs across the four conditions and every flipping voxel remain visible.
