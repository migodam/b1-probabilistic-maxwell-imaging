# Main48:24 mixed-object phase summary

Only FINAL_CAP_ELIGIBLE rows from the already-authorized scored CSV are used. CPU_480 rows are not pooled or duplicated. There are24 independent mixed materials,144 method/object rows, plus24high materials. No GT array, prediction array, private file, forward, scoring or bootstrap is used.

|Method|Mixed materialPASS/24|Median globalNRMSE|Low-regionRMSE mean / median|High-regionRMSE mean / median|BackgroundRMSE mean / median|Estimated low mean|
|---|---:|---:|---:|---:|---:|---|
|cie_single|12|0.097790|1.554253 / 1.337725|3.346260 / 3.060160|0.983327 / 0.714148|MISSING|
|cie_dual|12|0.093614|1.370716 / 1.321013|3.081854 / 3.032359|0.854048 / 0.715189|MISSING|
|global_dual|13|0.092210|1.372524 / 1.336200|3.095793 / 3.010508|0.851979 / 0.713483|MISSING|
|fullG19|13|0.092404|1.368456 / 1.347356|3.057980 / 2.999296|0.849843 / 0.701503|MISSING|
|freezeA|10|0.105573|1.918503 / 1.416000|5.602082 / 3.282503|0.985828 / 0.735217|MISSING|
|contractedT|13|0.091786|1.380157 / 1.333061|3.068087 / 3.045891|0.849461 / 0.699786|MISSING|

RegionRMSE is in absolute chi units. The low truth phase ischi2.5; each per-object regional RMSE uses only that object’s truth-region cell count from the source table. Summary means/medians give each material equal weight. The CSV also reports a separately labelled pooled-cell RMSE=sqrt(sum(n_cells*RMSE_region²)/sum(n_cells)); those cells are not independent scientific samples.

GlobalNRMSE is dimensionless: globalRMSE over576cells divided byhigh_chi−.1, i.e.7.9/11.9/17.9 for the three declared contrasts. Thus global and regional numbers have different denominators/units. Regional errors do not change the inherited material gate.

Estimated low-phase mean is absent from the original scored CSV. It remainsnull/MISSING for all144mixed rows; it cannot be reconstructed from RMSE alone. No GT/prediction arrays were opened to supply it, and the earlier16-object low means must not be substituted into this48-object cohort.

Each high/mixed×high_chi8/12/18 group has8 independent objects per method. Continuous NRMSE and material-success counts for all36method/group rows are in REGIME_CONTRAST_GROUP_SUMMARY36.csv; all288per-object rows and the144mixed rows remain available. This is mechanical descriptive reporting only, with no new threshold or inferential claim.

## High/mixed × contrast groups

Each cell shows mean dimensionless NRMSE (material successes/8). Each row contains8independent materials, reused across methods; contrast amplitude denominators are7.9/11.9/17.9.

|Regime / highchi|Objects per method|CIE|Dual|Global|G19|freezeA|contractedT|
|---|---:|---:|---:|---:|---:|---:|---:|
|high / 8|8|0.218152 (3/8)|0.153040 (4/8)|0.152744 (4/8)|0.152122 (4/8)|0.265756 (3/8)|0.151956 (4/8)|
|high / 12|8|0.127486 (6/8)|0.128089 (6/8)|0.127479 (6/8)|0.127888 (6/8)|0.265192 (3/8)|0.128268 (6/8)|
|high / 18|8|0.272756 (0/8)|0.236388 (0/8)|0.236704 (0/8)|0.290218 (0/8)|0.367534 (0/8)|0.219263 (0/8)|
|mixed / 8|8|0.116786 (6/8)|0.116603 (5/8)|0.116649 (6/8)|0.116031 (6/8)|0.115981 (6/8)|0.116143 (6/8)|
|mixed / 12|8|0.099443 (3/8)|0.099715 (3/8)|0.099818 (3/8)|0.098750 (3/8)|0.155134 (2/8)|0.099097 (3/8)|
|mixed / 18|8|0.119867 (3/8)|0.095769 (4/8)|0.095922 (4/8)|0.095503 (4/8)|0.186001 (2/8)|0.095589 (4/8)|
