# R4 blind confirmation and exposed-history development

All R3-DIGITS targets were already exposed. V1, factor separation, V2, GT paths, spectrum analysis, and historical Hessian analysis are DEVELOPMENT or OFFLINE_DIAGNOSTIC. None is a fresh blind result.

The final V3 inverse source and settings were frozen at the UTC recorded in `audit/V3_SOURCE_SNAPSHOT/FREEZE.json` and `configs/FRESH_CONTRACT.json`, before any new target was saved. The first generator process failed on an unavailable font before any forward or target file; that failed attempt is retained in `audit/GENERATOR_FONT_FAILURE.json`. Only the generator font was changed to an installed font before actual generation. The inverse source and settings were unchanged.

A separate generator created three new anonymous targets, using the same fixed 576-cell nonlinear VIE operators. The digit labels, fonts, masks, exact materials, transformations and generation entropy were absent from the inverse inputs. No extra frequencies, known support, or templates were supplied. Entropy was not saved. The inverse runner receives anonymous Y, public acquisition/operators, bounds, generic common regularizer, and public initialization permissions.

Each of the six final methods has 60 seconds, 2500 full LU and 60000 exposed complex RHS, stopping at the first resource cap; a final block can overshoot an RHS cap by its recorded block size. Two processes request single-thread BLAS. The same full-physics model is used at all 8/12/24 material parameterization stages. All L-BFGS arms use the same material scaling; TRF uses the analytic Jacobian and its native Jacobian scaling.

Uniform-branch and fixed-uniform-multistart arms first pay for a common public-m1 anchor, then share remaining resources among declared seeds. The uniform ray uses no hidden material shape. Four spread-out detected valleys are retained before correction, with all truncations recorded. The ablation corrects only the detected raw-score minimum. Both the ordinary multistart and branch method are finite search procedures, not probability models.

Online top1 minimizes the exact original-data objective plus identical smooth/TV terms. All predictions must be closed before the independent scorer is allowed to read the new materials. The subsequent release of seedless material arrays allows exact score replay. Best@K is explicitly posthoc GT ranking.

The predeclared joint success gate is relative full-data residual <=0.001, material RMSE <=1.5 and fixed-threshold IoU >=0.7. Data fit, material quality, finite candidate diversity and costs are also reported separately. A high-residual endpoint is never a distinct solution of the measured data. A small-residual but poor-material endpoint is not unique GT recovery.

Isolation is by process inputs and access discipline, not an operating-system security sandbox. No source hashes, repeated source comparisons, or formula precision experiments are part of this protocol.
