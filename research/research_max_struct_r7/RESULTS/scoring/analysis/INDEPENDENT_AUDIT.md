# Independent post-seal audit

|Arm|Same-grid material /24|Same-grid raw CN99 /24|Same-grid material+CN99 /24|Cross-grid material /6|Cross-grid raw CN99 /6|
|---|---:|---:|---:|---:|---:|
|cie_single|18|16|15|3|3|
|cie_dual|16|16|14|4|4|
|global_dual|16|16|14|4|4|
|a_material|17|17|14|4|4|
|a_plus|17|17|14|4|4|
|physical_lbfgs|10|11|9|0|0|

Cross-grid noise-only passes are not complete model-error acceptance. Six cross-grid GT feasible witnesses remain NOT_AVAILABLE_DIFFERENT_OPERATOR.

## Exact case-ID success contrasts

a_material vs global_dual: wins r7_obj_12, r7_obj_20; losses r7_obj_09.
a_plus vs a_material: wins none; losses none.
a_plus vs cie_dual: wins r7_obj_12, r7_obj_20; losses r7_obj_09.
a_material vs cie_single: wins r7_obj_12, r7_obj_24; losses r7_obj_09, r7_obj_17.

The globalDual comparison loses obj09 because IoU .696629 is below .7 while globalDual IoU .712644 passes. Obj17 is a different, single-CIE comparison and is not substituted into the global comparison. Numeric scores were correct throughout.

## Full576-variable physical endpoints: two large positive cases

|Case|Arm|Sources|NRMSE|IoU|Data residual|Phi|Linear relative residual|LU|Complex RHS|CPU s|
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
|r7_obj_12|cie_single|6|0.721134|0.208589|0.0247355|0.000336799|2.3e-15|3048|36552|79.232|
|r7_obj_12|cie_dual|6|0.735001|0.159763|0.0216655|0.000265597|2.55e-15|2988|35832|75.807|
|r7_obj_12|global_dual|6|0.735114|0.161677|0.0216722|0.000265649|2.58e-15|3077|36900|73.626|
|r7_obj_12|a_material|6|0.0396456|0.923077|0.000942055|4.54303e-06|1.67e-15|1710|20496|76.915|
|r7_obj_12|a_plus|6|0.0391449|0.953125|0.00094195|4.5407e-06|1.51e-15|1745|20916|51.529|
|r7_obj_12|physical_lbfgs|6|0.337123|0|0.832173|0.346258|1.86e-15|332|3984|7.422|
|r7_obj_20|cie_single|12|0.0702233|0.956835|0.00106648|1.31479e-05|2.53e-15|2097|50280|48.999|
|r7_obj_20|cie_dual|12|0.296361|0.622093|0.0130284|0.00010584|2.58e-15|2730|65472|58.206|
|r7_obj_20|global_dual|12|0.296552|0.611429|0.0130272|0.000105841|2.59e-15|3035|72792|62.085|
|r7_obj_20|a_material|12|0.0698092|0.971223|0.00107145|1.31478e-05|2.5e-15|2123|50892|95.793|
|r7_obj_20|a_plus|12|0.068721|0.964286|0.00105732|1.31518e-05|2.52e-15|2005|48060|63.258|
|r7_obj_20|physical_lbfgs|12|0.467981|0|1.01596|0.516092|1.52e-15|2|48|0.656|

All final maps contain576 unknowns. Physical matrices are576×576, source and adjoint costs counted; auxiliary FFT actions are not substituted for the full physical verification. These are selected positive examples from the disclosed complete30-object ledger, not the whole benchmark.

## Initialization and verification

cie_dual: host initialization successes 4/6; m=1 successes 4/6; same existing objects, no extra independent n.
a_material: host initialization successes 5/6; m=1 successes 5/6; same existing objects, no extra independent n.
a_plus: host initialization successes 5/6; m=1 successes 5/6; same existing objects, no extra independent n.

All180 primary and47 extra endpoints are sealed, legal and scored; no failed case was dropped. Residual recomputation agrees exactly; Phi max discrepancy1.11e−16. Float32 budget rounding and final serialization-tail limits are quantified in JSON. The initial six-group panels and all extra/rational maps were visually inspected.

Rational-action endpoints:0/12 material successes on3 disclosed old failed anchors. Finite action rank2/rank8 and continued-anchor controls did not escape to acceptable material reconstructions.

The Aplus-vs-mature-sourcewiseDual endpoint and equalCPU/LU/RHS comparisons are SECONDARY_USER_REQUESTED and retained separately. No post-hoc noninferiority margin or GT-selected candidate was added.
