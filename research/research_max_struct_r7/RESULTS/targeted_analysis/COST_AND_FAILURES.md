# Targeted20 cost-matching and retained failures

All80 tasks produced sealed public predictions with zero execution exceptions. All20 materials and all failures are retained. No optimizer rerun, new physical call, GT selection or sample replacement occurred.

|Comparator|Axis|Aplus /20 successes|Comparator /20 successes|Mean NRMSE delta|95% paired object bootstrap CI|
|---|---|---:|---:|---:|---|
|cie_dual|CPU|12|9|-0.0690573|[-0.147958, -0.0135993]|
|cie_dual|LU|12|9|-0.0689606|[-0.150144, -0.0132592]|
|cie_dual|complex_RHS|12|9|-0.0689606|[-0.150407, -0.013336]|
|global_dual|CPU|12|9|-0.07142|[-0.15084, -0.0154717]|
|global_dual|LU|12|9|-0.0716868|[-0.153567, -0.0151446]|
|global_dual|complex_RHS|12|9|-0.0716882|[-0.154361, -0.0145651]|
|cie_single|CPU|10|8|-0.0721711|[-0.165012, 0.00548403]|
|cie_single|LU|12|8|-0.0873688|[-0.173746, -0.0196061]|
|cie_single|complex_RHS|12|8|-0.0873688|[-0.172994, -0.0188811]|

The endpoint means cover all20 pairs. The comparison with cie_single loses endpoint mean-NRMSE CI exclusion of zero under equal CPU (CI includes zero), despite the unadjusted full-allowance endpoint difference. It is therefore necessary to distinguish the full allowance from equal CPU.

|Retained Aplus failure|Family|NRMSE|IoU|Residual|Corrector stop|Known better GT-Phi witness|
|---|---|---:|---:|---:|---|---|
|r7_target_02|curved_bar|0.72721|0.19298|0.019286|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|True|
|r7_target_10|slit_shell|0.64903|0.13158|0.041059|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|True|
|r7_target_12|slit_shell|0.57454|0.13901|0.052519|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|True|
|r7_target_13|slit_shell|0.69084|0.10204|0.025148|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|True|
|r7_target_14|slit_shell|0.67103|0.12411|0.021107|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|True|
|r7_target_15|slit_shell|0.54192|0.18009|0.061145|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|True|
|r7_target_16|slit_shell|0.63102|0.10117|0.040668|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|True|
|r7_target_19|slit_shell|0.63321|0.14176|0.041469|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|True|

Aplus fails one curved_bar and seven slit_shell objects. All eight fail the separate CN99 check. Two shell failures report optimizer numerical convergence; the remaining six stop at the iteration allowance. These are convergence/stopping diagnostics rather than global-optimality, uniqueness or basin certificates. The saved-GT feasible objective is a strictly better offline witness for every listed failure and is never used to select predictions.

The shell success set is identical for Aplus and cie_single: targets11,17,18. The two shell gains over cie_dual/global_dual (11,18) are therefore also achieved by the simple cie_single baseline. Curved_bar target03 succeeds only with Aplus among these four arms; targets00,01,08 succeed with both dual baselines as well as Aplus but fail with cie_single. This describes fixed observed outcomes, not a posthoc selection rule.

Figures all20_group_01.png through all20_group_04.png show all20 GT maps and all80 material predictions; pairwise_cost_frontiers.png shows all three cost axes and all three comparator pairs. Visual inspection of group04 confirms readable labels and retained shell failures.
