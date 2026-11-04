# R7 targeted selected-family secondary confirmation

Twenty fresh independent materials, all retained. Ten curved_bar at contrast12/six full sources; ten slit_shell at contrast18/12 half sources; relative complex noise .001. This conditional study follows selection of two successful family/acquisition strata in primary30; results are never pooled and cannot establish universal recovery. Primary material success requires amplitude-normalized NRMSE<=.15, Otsu IoU>=.7, and legal material box. CN99 is separate.

|Arm|Material successes /20|Material+CN99 /20|Mean NRMSE|Median NRMSE|CPU total s|LU total|Complex RHS total|
|---|---:|---:|---:|---:|---:|---:|---:|
|a_plus|12|9|0.28525|0.072796|1089.69|43594|849516|
|cie_dual|9|8|0.35356|0.30891|1061.57|46876|909312|
|global_dual|9|8|0.35604|0.30831|1027.61|45126|870528|
|cie_single|8|6|0.37375|0.38023|855.19|43013|836184|

|Family|Arm|Material successes /10|Mean NRMSE|Median NRMSE|
|---|---|---:|---:|---:|
|curved_bar|a_plus|9|0.10931|0.040379|
|curved_bar|cie_dual|8|0.17844|0.040822|
|curved_bar|global_dual|8|0.17945|0.040424|
|slit_shell|a_plus|3|0.46118|0.60278|
|slit_shell|cie_dual|1|0.52868|0.6269|
|slit_shell|global_dual|1|0.53263|0.62683|
|curved_bar|cie_single|5|0.27665|0.19198|
|slit_shell|cie_single|3|0.47084|0.59187|

|Comparison|Metric|Mean left-minus-right|95% paired-object bootstrap CI|
|---|---|---:|---|
|a_plus vs cie_dual|material_nrmse|-0.068311|[-0.142132, -0.0124211]|
|a_plus vs cie_dual|material_image_success|0.15|[0, 0.3]|
|a_plus vs cie_dual|cpu_seconds|1.40574|[-1.64537, 3.86227]|
|a_plus vs cie_dual|LU|-164.1|[-382.551, 12.255]|
|a_plus vs cie_dual|actual_complex_RHS|-2989.8|[-6655.85, 100.26]|
|a_plus vs global_dual|material_nrmse|-0.0707908|[-0.149351, -0.0153613]|
|a_plus vs global_dual|material_image_success|0.15|[0, 0.3]|
|a_plus vs global_dual|cpu_seconds|3.10385|[-0.755664, 7.24814]|
|a_plus vs global_dual|LU|-76.6|[-336.309, 196.406]|
|a_plus vs global_dual|actual_complex_RHS|-1050.6|[-6107.19, 4990.98]|
|a_plus vs cie_single|material_nrmse|-0.0885004|[-0.177438, -0.0198268]|
|a_plus vs cie_single|material_image_success|0.2|[0.05, 0.4]|
|a_plus vs cie_single|cpu_seconds|11.7251|[3.64437, 19.5203]|
|a_plus vs cie_single|LU|29.05|[-228.45, 259.002]|
|a_plus vs cie_single|actual_complex_RHS|666.6|[-3378.09, 4426.83]|

Ten-thousand paired20-object bootstrap resamples, seed2026110429, no multiplicity adjustment. Family10 tables are conditional descriptive evidence. All20 material realizations are independent; arm/acquisition rows are not new objects. Failures remain in success denominators; continuous metrics disclose finite pairs.

Actual CPU/LU/RHS comparisons use each object/pair observed common window [max(first saved cost),min(last saved cost)] at11 fixed fractions, taking the last paid original-Phi incumbent at or below budget. No interpolation, extrapolation or GT iterate choice. Sparse auxiliary and every-ten-observation checkpoint resolution and float32 storage can leave unused cost; exact budgets and tails are recorded. CPU includes setup and auxiliary work; equal LU/RHS is not equal CPU. CN99 reflects same-grid circular complex Gaussian noise. No new forward, LU or RHS calls occur in scoring.
