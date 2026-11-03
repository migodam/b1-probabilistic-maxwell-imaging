# R5 material-only 扩展：离线评分

新增3个已公开开发对象各1次material-only；历史R4引用为同3例的global/onehop/native三个原端点，共9条。新盲测对象0。历史参考来自另一平台，不称同运行或等算力对照。

独立封存时间：2026-11-03T12:08:17.935316+00:00；本次首次GT读取：2026-11-03T12:08:17.935928+00:00。原STAGE_I_SEAL未修改；评分新增物理LU/RHS为0。

历史合同：RMSE/12≤.15，预测自身Otsu256 IoU≥.7，公共99%噪声球与材料盒；未增加门或改变参数/选择。

| 对象 | 来源 | 臂 | RMSE/12 | IoU | 相对残差 | Phi | 历史联合门 |
|---|---|---|---:|---:|---:|---:|---|
|obs_42cb233e4f224bb3|new_R5_extension|material_only|0.766322569|0.106280193|0.066004220|0.00220906641|FAIL|
|obs_42cb233e4f224bb3|historical_R4_reference|global|0.741539073|0.121827411|0.090322616|0.00410523163|FAIL|
|obs_42cb233e4f224bb3|historical_R4_reference|onehop|0.761367699|0.110047847|0.065757057|0.00219301446|FAIL|
|obs_42cb233e4f224bb3|historical_R4_reference|native_extended|0.743118035|0.130653266|0.090835767|0.00415147687|FAIL|
|obs_63abafa750fc483c|new_R5_extension|material_only|0.404823109|0.403846154|0.111383239|0.0062238167|FAIL|
|obs_63abafa750fc483c|historical_R4_reference|global|0.575754689|0.225352113|0.097670798|0.00479410222|FAIL|
|obs_63abafa750fc483c|historical_R4_reference|onehop|0.576876052|0.239436620|0.095880909|0.00462148321|FAIL|
|obs_63abafa750fc483c|historical_R4_reference|native_extended|0.577783509|0.251798561|0.095780989|0.00461172046|FAIL|
|obs_100bd65498644322|new_R5_extension|material_only|0.072525170|0.917355372|0.001009799|6.25954263e-06|PASS|
|obs_100bd65498644322|historical_R4_reference|global|0.073369280|0.924369748|0.000999205|6.33450284e-06|PASS|
|obs_100bd65498644322|historical_R4_reference|onehop|0.071620592|0.917355372|0.001010864|6.24289412e-06|PASS|
|obs_100bd65498644322|historical_R4_reference|native_extended|0.085988378|0.875000000|0.001037607|6.64932294e-06|PASS|

GT只用于封存后的离线评分，不能用于在线择根或后续参数调整。三个开发对象的描述性结果不能解释为总体成功率或新盲测优势。
