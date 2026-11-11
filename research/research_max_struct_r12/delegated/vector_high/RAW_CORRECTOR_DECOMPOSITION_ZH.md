# 高对比度 true3D：已保存状态的 corrector 分解

本补充不改变算法、不跑新反演、不做新 Maxwell 正演/伴随，也不重新读取材料truth。只使用已保存材料、Phi、source-normalized residual、原quadraticR和已有离线材料指标。AP/DR原始物理current保存在audit里，因此可用原接收算子S读取raw残差；CIE/Dual的raw物理current未保存，精确raw全局CN99状态如实标缺。

原目标为 Phi=D+R，D=(1/P)Σ ||r_p||²/||Y_p||²（无1/2），R=2e-4 mean_edge(chi_i-chi_j)²+1e-6 mean_cell chi²。

|对象|方法|raw D|corrected D|raw R|corrected R|raw→corrected relL2|分解|
|---|---|---:|---:|---:|---:|---|---|
|random_high8|CIE_sourcewise|0.025810576|0.00013576631|0.00016134558|0.00059158292|0.457843→0.129785|data_improves_prior_worsens|
|random_high8|Dual_sourcewise|0.027210833|0.00013576434|0.00016278517|0.00059158489|0.454820→0.129783|data_improves_prior_worsens|
|random_high8|DR|0.00010210215|0.00013576744|0.00071731734|0.0005915818|0.061631→0.129783|data_worsens_prior_improves|
|random_high8|AP|0.0001252536|0.00013576599|0.0006516546|0.00059158325|0.103885→0.129786|data_worsens_prior_improves|
|random_mixed2_8|CIE_sourcewise|0.017657269|0.00013032712|4.1715501e-05|0.00024504379|0.618636→0.227237|data_improves_prior_worsens|
|random_mixed2_8|Dual_sourcewise|0.01786649|0.00013032539|4.3466567e-05|0.00024504552|0.608627→0.227232|data_improves_prior_worsens|
|random_mixed2_8|DR|0.00010034028|0.00013032556|0.00035855561|0.00024504535|0.078315→0.227239|data_worsens_prior_improves|
|random_mixed2_8|AP|0.0001941365|0.00013032811|0.00026823625|0.0002450428|0.216416→0.227260|both_data_and_prior_improve|
|random_high8_crossgrid8|CIE_sourcewise|0.055152674|0.002160323|0.00023697772|0.0012555988|0.397404→0.432218|data_improves_prior_worsens|
|random_high8_crossgrid8|Dual_sourcewise|0.058700806|0.0021603263|0.0002378937|0.0012555956|0.395368→0.432220|data_improves_prior_worsens|
|random_high8_crossgrid8|DR|0.0035434609|0.0021541298|0.0015323913|0.0012609406|0.491405→0.434481|both_data_and_prior_improve|
|random_high8_crossgrid8|AP|0.0026371795|0.0021603336|0.0013266345|0.0012555883|0.424022→0.432217|both_data_and_prior_improve|

|对象|方法|raw Phi|corrected Phi|raw CN99|corrected CN99|
|---|---|---:|---:|---|---|
|random_high8|CIE_sourcewise|0.0259719219|0.000727349236|False (norm界认证；exact残差MISSING)|False, ratio=1.14389|
|random_high8|Dual_sourcewise|0.027373618|0.000727349236|False (norm界认证；exact残差MISSING)|False, ratio=1.14388|
|random_high8|DR|0.000819419496|0.000727349236|True, ratio=0.99017|False, ratio=1.14389|
|random_high8|AP|0.000776908203|0.000727349235|False, ratio=1.09776|False, ratio=1.14389|
|random_mixed2_8|CIE_sourcewise|0.017698985|0.000375370909|False (norm界认证；exact残差MISSING)|False, ratio=1.12571|
|random_mixed2_8|Dual_sourcewise|0.0179099569|0.00037537091|False (norm界认证；exact残差MISSING)|False, ratio=1.12570|
|random_mixed2_8|DR|0.000458895883|0.00037537091|True, ratio=0.98395|False, ratio=1.12570|
|random_mixed2_8|AP|0.000462372751|0.000375370911|False, ratio=1.37192|False, ratio=1.12571|
|random_high8_crossgrid8|CIE_sourcewise|0.055389652|0.00341592185|False (norm界认证；exact残差MISSING)|False, ratio=4.57623|
|random_high8_crossgrid8|Dual_sourcewise|0.0589386993|0.00341592185|False (norm界认证；exact残差MISSING)|False, ratio=4.57623|
|random_high8_crossgrid8|DR|0.00507585225|0.00341507045|False, ratio=5.86882|False, ratio=4.56969|
|random_high8_crossgrid8|AP|0.00396381393|0.00341592185|False, ratio=5.04883|False, ratio=4.57624|

|对象|方法|raw选择迭代/周期|raw候选CPU|完整proposal实际CPU|corrected选择eval / 总acceptednit|完整CPU|
|---|---|---|---:|---:|---|---:|
|random_high8|CIE_sourcewise|auxcycle120|2.603|2.603|184 / 182|4.268|
|random_high8|Dual_sourcewise|auxcycle120|2.708|2.708|187 / 186|4.398|
|random_high8|DR|projection2000|14.272|34.734|179 / 176|38.064|
|random_high8|AP|projection1900|20.007|41.671|185 / 182|44.011|
|random_mixed2_8|CIE_sourcewise|auxcycle120|3.440|3.441|317 / 315|7.606|
|random_mixed2_8|Dual_sourcewise|auxcycle120|3.333|3.334|310 / 309|7.544|
|random_mixed2_8|DR|projection2400|24.908|40.745|270 / 269|45.797|
|random_mixed2_8|AP|projection3000|30.071|39.646|284 / 283|42.864|
|random_high8_crossgrid8|CIE_sourcewise|auxcycle120|4.186|4.187|237 / 236|7.780|
|random_high8_crossgrid8|Dual_sourcewise|auxcycle96|3.451|4.491|254 / 253|8.255|
|random_high8_crossgrid8|DR|projection50|0.574|35.972|202 / 201|37.880|
|random_high8_crossgrid8|AP|projection1450|8.847|26.477|212 / 211|28.394|

CIE/Dual的exact raw global residual缺失，但其CN99不可行可由已知source-normalized D严格推出：P D min_p||Y_p||² ≤ ||r||² ≤ P D max_p||Y_p||²。六个baseline raw下界均已大于delta，故无需新forward即可认证球外；精确残差仍不补造。

selected candidate出现时CPU与跑完整个固定proposal实际CPU是不同数字，均保留。corrector所选paid evaluation序号准确；该evaluation对应accepted iteration未保存，不能臆造。每项已付费物理/aux计数及AP/DR fixed-H计数见JSON。

同CPU前缀的raw-only比较（只考虑proposal，按原Phi选，不按材料误差选）：

|对象|CPU ceiling|方法|已付费候选CPU|Phi|prior|材料relL2|
|---|---:|---|---:|---:|---:|---:|
|random_high8|10|CIE_sourcewise|2.603|0.025971922|0.00016134558|0.457843|
|random_high8|30|CIE_sourcewise|2.603|0.025971922|0.00016134558|0.457843|
|random_high8|60|CIE_sourcewise|2.603|0.025971922|0.00016134558|0.457843|
|random_high8|120|CIE_sourcewise|2.603|0.025971922|0.00016134558|0.457843|
|random_high8|10|Dual_sourcewise|2.708|0.027373618|0.00016278517|0.454820|
|random_high8|30|Dual_sourcewise|2.708|0.027373618|0.00016278517|0.454820|
|random_high8|60|Dual_sourcewise|2.708|0.027373618|0.00016278517|0.454820|
|random_high8|120|Dual_sourcewise|2.708|0.027373618|0.00016278517|0.454820|
|random_high8|10|DR|9.719|0.00081970723|0.00071760035|0.062917|
|random_high8|30|DR|14.272|0.0008194195|0.00071731734|0.061631|
|random_high8|60|DR|14.272|0.0008194195|0.00071731734|0.061631|
|random_high8|120|DR|14.272|0.0008194195|0.00071731734|0.061631|
|random_high8|10|AP|9.581|0.00078932837|0.00062447567|0.127390|
|random_high8|30|AP|20.007|0.0007769082|0.0006516546|0.103885|
|random_high8|60|AP|20.007|0.0007769082|0.0006516546|0.103885|
|random_high8|120|AP|20.007|0.0007769082|0.0006516546|0.103885|
|random_mixed2_8|10|CIE_sourcewise|3.440|0.017698985|4.1715501e-05|0.618636|
|random_mixed2_8|30|CIE_sourcewise|3.440|0.017698985|4.1715501e-05|0.618636|
|random_mixed2_8|60|CIE_sourcewise|3.440|0.017698985|4.1715501e-05|0.618636|
|random_mixed2_8|120|CIE_sourcewise|3.440|0.017698985|4.1715501e-05|0.618636|
|random_mixed2_8|10|Dual_sourcewise|3.333|0.017909957|4.3466567e-05|0.608627|
|random_mixed2_8|30|Dual_sourcewise|3.333|0.017909957|4.3466567e-05|0.608627|
|random_mixed2_8|60|Dual_sourcewise|3.333|0.017909957|4.3466567e-05|0.608627|
|random_mixed2_8|120|Dual_sourcewise|3.333|0.017909957|4.3466567e-05|0.608627|
|random_mixed2_8|10|DR|9.617|0.00046535749|0.00036491149|0.094606|
|random_mixed2_8|30|DR|24.908|0.00045889588|0.00035855561|0.078315|
|random_mixed2_8|60|DR|24.908|0.00045889588|0.00035855561|0.078315|
|random_mixed2_8|120|DR|24.908|0.00045889588|0.00035855561|0.078315|
|random_mixed2_8|10|AP|9.986|0.00076592935|0.00019157112|0.317466|
|random_mixed2_8|30|AP|29.677|0.00046333895|0.00026752001|0.217625|
|random_mixed2_8|60|AP|30.071|0.00046237275|0.00026823625|0.216416|
|random_mixed2_8|120|AP|30.071|0.00046237275|0.00026823625|0.216416|
|random_high8_crossgrid8|10|CIE_sourcewise|4.186|0.055389652|0.00023697772|0.397404|
|random_high8_crossgrid8|30|CIE_sourcewise|4.186|0.055389652|0.00023697772|0.397404|
|random_high8_crossgrid8|60|CIE_sourcewise|4.186|0.055389652|0.00023697772|0.397404|
|random_high8_crossgrid8|120|CIE_sourcewise|4.186|0.055389652|0.00023697772|0.397404|
|random_high8_crossgrid8|10|Dual_sourcewise|3.451|0.058938699|0.0002378937|0.395368|
|random_high8_crossgrid8|30|Dual_sourcewise|3.451|0.058938699|0.0002378937|0.395368|
|random_high8_crossgrid8|60|Dual_sourcewise|3.451|0.058938699|0.0002378937|0.395368|
|random_high8_crossgrid8|120|Dual_sourcewise|3.451|0.058938699|0.0002378937|0.395368|
|random_high8_crossgrid8|10|DR|0.574|0.0050758522|0.0015323913|0.491405|
|random_high8_crossgrid8|30|DR|0.574|0.0050758522|0.0015323913|0.491405|
|random_high8_crossgrid8|60|DR|0.574|0.0050758522|0.0015323913|0.491405|
|random_high8_crossgrid8|120|DR|0.574|0.0050758522|0.0015323913|0.491405|
|random_high8_crossgrid8|10|AP|8.847|0.0039638139|0.0013266345|0.424022|
|random_high8_crossgrid8|30|AP|8.847|0.0039638139|0.0013266345|0.424022|
|random_high8_crossgrid8|60|AP|8.847|0.0039638139|0.0013266345|0.424022|
|random_high8_crossgrid8|120|AP|8.847|0.0039638139|0.0013266345|0.424022|

上述是相同CPU上限下已可获得的固定proposal状态，未花完的时间不虚计为开销；已完成CIE/Dual后其raw状态仍可取得。完整proposal费用与此前包含corrector的预封存前缀另列。该后续诊断不产生新的独立确认或新GT选择。

新增成本：Maxwell LU/sourceRHS/adjointRHS全为0；仅原外部接收算子装配1次和对已保存current做360个接收列动作。原始数据、预测与seal未改写。
