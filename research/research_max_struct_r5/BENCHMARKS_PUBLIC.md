# R5 公开 benchmark 资格调查（不是最终运行裁决）

核对日期：2026-11-03。仅浏览官方入口、原作者文章/仓库与小型metadata；没有下载数据集、没有运行反演或读取GT数组。`NOT_RUN`，无新physics/LU/RHS。最终是否纳入本轮由主线程决定。

## Fresnel：来源确认，实际下载资格尚未闭环

[Institut Fresnel官方入口](https://www.fresnel.fr/3Ddatabase/)同时列出2001二维均匀、2005二维非均匀（TE/TM）及2009三维数据，声明科研免费使用并链接IOP。三个落点本轮均受限： [2001](https://iopscience.iop.org/article/10.1088/0266-5611/17/6/301/)、[2005](https://iopscience.iop.org/article/10.1088/0266-5611/21/6/S09/)、[2009](https://iopscience.iop.org/article/10.1088/0266-5611/25/2/024001/)。`https://www.fresnel.fr/2Ddatabase/`亦未能访问。官方入口存在不代表本机已取得data；当前为 **ACCESS_PARTIAL / DOWNLOAD_NOT_VERIFIED**。[官方说明](https://www.fresnel.fr/3Ddatabase/)

|数据|已核实合同（原作者全文）|格式、标定、GT及适配缺口|
|---|---|---|
|2005二维非均匀|2–10 GHz、1 GHz步长；FoamMetExt至18 GHz；TE/TM。源/收半径1.67m；接收1°步长，源附近±60°不可测；源8或18个，按target而定。|ASCII，10行header；7列：source index、receiver index、GHz、total实/虚、incident实/虚；散射由两者相减。尺寸/材料由文章给出，非现成逐cell sealed GT。论文比较只使用一个复标定系数，由反向接收incident参考给出；不能用未知物材料反拟合增益而不记为额外信息。|
|2009三维|3–8 GHz、0.25步长；source半径1.796m，9个azimuth×9个polar；receiver在赤道，10°步长并有遮挡。两种source极化，receiver为φ分量。|每target两份`PP.exp`/`TP.exp`，10列：θs、φs、θr、freq、incident实/虚、total实/虚、scattered实/虚。scattered已漂移纠正/标定；饱和点以0替代，必须mask，不能当真零。几何及材料在文中，Myster后来已披露，不能叫本轮fresh blind。|

二维来源为Geffrin–Sabouroux–Eyraud，*Free space experimental scattering database continuation: experimental set-up and measurement precision*，2005，DOI [10.1088/0266-5611/21/6/S09](https://doi.org/10.1088/0266-5611/21/6/S09)，[作者上传全文](https://www.researchgate.net/publication/231069726_Free_space_experimental_scattering_database_continuation_Experimental_set-up_and_measurement_precision)，§4及§5，ppS124–S125；几何与材料见§2/§3。三维来源为Geffrin–Sabouroux，*Continuing with the Fresnel database: experimental setup and improvements in 3D scattering measurements*，2009，DOI [10.1088/0266-5611/25/2/024001](https://doi.org/10.1088/0266-5611/25/2/024001)，[作者上传全文](https://www.researchgate.net/publication/230990377_Continuing_with_the_Fresnel_database_Experimental_setup_and_improvements_in_3D_scattering_measurements)，§2/§3/§5。官方站上的[2008原作者短文](https://www.fresnel.fr/perso/belkebir/Articles/Geffrin_APL_08.pdf)独立解释金属球Mie标定、平面波照明与漂移处理；不可替代2009文件的实际header核对。

## 与当前disk-pulse VIE的最短适配路径（工程推断）

已读本地 `audit/BACKEND_REUSE.md`：当前FBG九数字是exp(−iωt)、k=2π、二维scalar-TM disk-pulse VIE，公开G/E/S固定。Fresnel **TM纯介质target** 最接近此模型；TE、金属混合与三维target都不能直接塞入这套scalar数据合同。

最短路径是先取得一份2005纯介质TM原始`.exp`及原header；解析source/receiver角度、频率和有效mask，固定一个频点；依据测得incident参考制定一次性的相位/幅度和长度无量纲映射；用既有disk-pulse核为真实几何重建G/E/S（不可复用九数字原G/E/S）；用独立标定或已知reference检验传播、符号、单位与模型差异，再冻结预算/初始化/正则化与残差评分。几何标注可供最终离线图形评分，但材料公称值/边界存在制造及测量误差，不能假装精确pixel GT。

实际下载与header尚未验证、标定协议尚未冻结，当前不能报告“Fresnel已可直接运行”。若只画标称target生成同算子synthetic观测，那是新的合成实验，不是Fresnel实测benchmark。2009需要已有vector DDA的真实单位、入射方向/极化、有限距离接收分量及mask适配，成本与模型变化更大。

## eispy2d：公开工具框架，不等于实测数据库

已浏览[官方仓库](https://github.com/andre-batista/eispy2d)与[原论文v2](https://arxiv.org/html/2111.02185)，§2、§3.1/§3.4。身份由[arXiv页面](https://arxiv.org/abs/2111.02185v2)及Crossref核对：Costa Batista, Adriano & Batista；正式IEEE Access 2025 DOI [10.1109/ACCESS.2025.3573679](https://doi.org/10.1109/ACCESS.2025.3573679)。仓库README仍有“under review”旧语句，应以DOI身份区分。仓库为GPL-3.0。

论文模型为二维TMz、均匀背景、线性各向同性非磁非色散材料，单频；采用exp(+iωt)/Hankel第二类。`InputData`散射场`NM×NS` complex V/m，材料`Ny×Nx`，总场`Ny*Nx×NS`按C顺序；Configuration定义frequency/wavelength、domain与观测圆。TestSet可合成随机targets、Benchmark可作算法统计比较。已读论文全文未找到“Fresnel”字符串；**本轮没有证据证明它附带可直接导入的官方Fresnel实测集**。不能将合成工具资格写成实测资格。[论文](https://arxiv.org/html/2111.02185)

最短适配是做数组/元数据导出包装及明确时谐符号、contrast虚部、单元积分与归一化映射，继续使用既有physics；若采用其forward solver生成data，应披露跨离散或同离散关系。没有安装库、克隆/执行仓库或替换backend；现成importer/API未逐个核实，不能承诺开箱即用。该候选是否优于整理既有九数字、是否值得本轮接入，留主线程裁决。

## 状态

- 官方来源/论文身份：CONFIRMED。
- 实测数据原文件下载/解析：NOT_RUN / ACCESS_PARTIAL。
- GT下载/读取：NOT_RUN。
- 新反演、benchmark数值结果：NOT_RUN。
- 结论权限：仅资格与缺口提取；没有发布、推送或最终是否运行判断。
