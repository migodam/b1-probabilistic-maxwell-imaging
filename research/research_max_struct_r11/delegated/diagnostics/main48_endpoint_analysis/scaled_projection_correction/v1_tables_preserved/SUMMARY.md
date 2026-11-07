# Main48 endpoint table integration

All288 authorized final cap rows strictly joined1:1 by(case_id,method) with the corrected exact-material-match endpoint table. No scoring/bootstrap/solver rerun; no GT arrays read.

|Method|MaterialPASS/48|CN99/48|PG<=1e-8/48|PGmedian|Scientific childCPU total/mean/median s|LU/RHS|
|---|---:|---:|---:|---:|---:|---:|
|cie_single|21|42|27|9.89e-09|7664.688 / 159.681 / 153.898|83125 / 1992672|
|cie_dual|22|44|27|9.87e-09|7475.969 / 155.749 / 153.125|79622 / 1908588|
|global_dual|23|44|29|9.88e-09|7602.297 / 158.381 / 151.953|82122 / 1968588|
|fullG19|23|43|25|9.94e-09|9973.719 / 207.786 / 197.305|85426 / 2047824|
|freezeA|16|31|32|9.71e-09|7825.266 / 163.026 / 168.688|62544 / 1498536|
|contractedT|23|45|27|9.88e-09|9426.984 / 196.396 / 192.172|83885 / 2010684|

PGcriterion is the descriptive originalgtol=.01/scale1e6=1e-8 applied to the saved unscaled originalPhi gradient including prior. It is not a new material success gate, a second-order local-minimum certificate, or information nonidentifiability evidence. Strict bound-sign counts and near1e-8 counts remain in the joined table.

case0022 remains in the main48 table. Its four interrupted/retry arms have unknown original full end-to-endCPU; valid successful-childOSCPU is known, but cannot be relabeled complete task expenditure. Exact end-to-end sums are null when affected, and source lower/observed-upper bounds remain separate. The47 sensitivity excludes the entire case0022 for all methods, preserving matched object sets. No retry is a new material.

See METHOD_SUMMARY/STOP_COUNTS, FAIL_CN99_PG_CROSSCOUNTS, full288 object differences and their47-case sensitivity, SUCCESS_SETS_BY_FAMILY and SUCCESS_SET_EQUALITY. All successes, failures and missing values retained.

Global, fullG19 and contractedT have exactly the same23 successful case IDs, including within every family; no reversal is hidden by the equal totals. All-method case sets are retained in SUCCESS_SETS_BY_FAMILY and their15 pairwise equality checks. On the matched47-object sensitivity, success counts stay21/22/23/23/16/23.

MaterialFAIL cross-counts for all48 (PGcriterion is descriptive only):

|Method|CN99in,PG<=tol|CN99in,PG>tol|CN99out,PG<=tol|CN99out,PG>tol|
|---|---:|---:|---:|---:|
|cie_single|15|7|4|1|
|cie_dual|13|11|1|1|
|global_dual|13|10|1|1|
|fullG19|13|8|3|1|
|freezeA|9|6|15|2|
|contractedT|12|11|2|0|

All288 have an exactly matched saved gradient, so noPG is missing. Corrector stop strings and this separately computed final-endpoint PGcriterion are both retained; the latter does not rewrite the original stop reason. The source helper uses unscaled originalPhi gradient including the same prior. No label here identifies a local minimum or information nonidentifiability.

Engineering processing costs: successful external child0.140931CPU seconds; two preserved startup failures0.035169+0.065941 seconds; recorded child total0.242041s. No physical work or GT arrays in any attempt. No additional study was run.
