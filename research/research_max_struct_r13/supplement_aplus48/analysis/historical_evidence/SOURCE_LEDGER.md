# R4–R13历史证据账本（提取材料，非新科学结论）

本账本逐cohort保留原报告口径。所有来源是本工作树现存文件；未读取任何GT NPZ、未运行反演、未改旧源码/结果。`E###`与CSV完全对应。数字是已报告值；`independent_object_count`若为历史组或多组，必须结合exposure与cohort读取。R13补充32个运行发生在对象已公开后；不能称新盲48。

## 必须区分的轴

- 材料归一化：二维amplitude NRMSE、三维relative chi L2、boxrange NRMSE不同。
- 成功标签：R4/R6部分联合数据资格；R7材料标签与CN99分列；R10后固定chi>1 support取代Otsu主口径；R13不新增合取成功门。
- 支出：固定3000步、body CPU、whole-child CPU、队列wall、共同付费prefix不能互换。
- 初值：uniform材料、NN材料+冷current、截断physical current、固定affine anchor、G19合法热重启、DR完整J/Q不同。
- A-material、R7最终A+ arc/guard、R13冻结A+ quintic/no-guard和G19不是同一版本。

## 分组来源与数字

### E001 · R4 · old b20d difficult bar
来源：`research/research_max_struct_r4/RESEARCH_REPORT_ZH.md:48-72`；## 2. 从公式得到的 A：传播本构缺陷，而非近似正演
对象/暴露：1；public historical development。
版本/表示：sourcewise CIE8000 / global / self / exact4cell / onehop / anti-phase。初值：same original init。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;12sources48receivers。预算/时钟：global/onehop mf3/6/10 x500; native8000; common corrector1000; no equal CPU claim。
估计器：corrected。指标：RMSE/12, prediction Otsu IoU, physical residual。
正面及报告数字：onehop .041672/.939759/.0009679 vs global .389183/.422680/.107821; raw onehop .08327/.77451。
负面及边界：all correctors hit1000;32grid onehop/global .053068/.053102 both recover; advantage disappears。

### E002 · R4 · fresh6 ring/bar/island
来源：`research/research_max_struct_r4/RESEARCH_REPORT_ZH.md:74-95`；## 4. 六个新鲜对象：没有超越成熟 CIE
对象/暴露：6；independent-at-original-scoring。
版本/表示：onehop / global / mature sourcewise native8000。初值：same init。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;CN independent complex noise; matched24。预算/时钟：same full frozen schedules; corrector1000; no wall truncation。
估计器：corrected primary minPhi。指标：joint descriptive .15/.7 + data/legality。
正面及报告数字：all three4/6 on same four objects。
负面及边界：onehop NRMSE1lower5higher vs global; minPhi portfolio still4/6; two failed bars pick worst material arm。

### E003 · R4 · fresh6 compute
来源：`research/research_max_struct_r4/RESEARCH_REPORT_ZH.md:133-145`；## 7. 计算成本与可复现性
对象/暴露：6；independent-at-original-scoring。
版本/表示：onehop/global/native8000。初值：same。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：CPU3804.34s; foreground elapsed1179.03s; summed workerwall3900.01s; LU18653/RHS445932。
估计器：cost。指标：not material success。
正面及报告数字：complete18 runs; no dropped failures。
负面及边界：onehop wall1147.09 vs global577.95 vs native2174.97; clocks cannot pool。

### E004 · R5 · four public mechanisms including old b20d
来源：`research/research_max_struct_r5/RESEARCH_REPORT_ZH.md:13-48`；### 1. 为什么R4在历史弯条上有效？
对象/暴露：4；public exposed;no independent confirmation。
版本/表示：global/global;onehop/onehop;onehop/ordinary;ordinary/onehop。初值：same original init; stage handoff m,J,Fourier coeff。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：mf3/6/10x500+common1000corrector;7 finite interventions。
估计器：corrected。指标：RMSE/12, IoU, residual。
正面及报告数字：ordinary current+onehop material .041480/.939759; localizes b20d benefit。
负面及边界：current-only .558514; stage switches .388613/.558738;42cb .766323 and63ab .404823 fail。

### E005 · R6 · historical9digits
来源：`research/research_max_struct_r6/RESEARCH_REPORT_ZH.md:7-29`；## 最终判断
对象/暴露：9；public historical regression。
版本/表示：DPCS-locked Fourier/APCS mf3/6/10 then12; CIE/Dual sourcewise PR; Global global PR; A/A+ global PR unless sourcewise named; A-material nonlinear onehop。初值：uniform host .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;12x48; historical clean contract。预算/时钟：3000aux +3000corrector;direct not included。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint;historical clean gate。
正面及报告数字：Dual9/9。
负面及边界：CIE8/9;Global8/9;A8/9;A digit5 .1142 but IoU .6667。

### E006 · R6 · historical difficult6
来源：`research/research_max_struct_r6/RESEARCH_REPORT_ZH.md:33-48`；## 2. 历史困难对象：正面结果确实扩大，但范围有限
对象/暴露：6；public historical development。
版本/表示：DPCS-locked Fourier/APCS mf3/6/10 then12; CIE/Dual sourcewise PR; Global global PR; A/A+ global PR unless sourcewise named; A-material。初值：host .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：3000aux+3000corrector;actual costs differ。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：A3/6 vs CIE/Dual/Global1/6;7b8d .04686 andb20d .04128。
负面及边界：42cb .70840;63ab .40475;e499 .42484 fail。

### E007 · R6 · fresh6 independent
来源：`research/research_max_struct_r6/RESEARCH_REPORT_ZH.md:52-75`；## 3. 六个新对象：没有新增确认优势
对象/暴露：6；independent-at-original-scoring。
版本/表示：DPCS-locked Fourier/APCS mf3/6/10 then12; CIE/Dual sourcewise PR; Global global PR; A/A+ global PR unless sourcewise named;A-material。初值：host .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;matched24;.1%complex RMS;12x48。预算/时钟：3000aux+3000corrector;thinring all hit cap;open shell reachesPGtol。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint;CN99 separate。
正面及报告数字：all four4/6。
负面及边界：A vs Dual3lower3higher;vs Global1lower5higher;thin ring andopen shell fail。

### E008 · R6 · 21 main task compute andFFT2 engineering
来源：`research/research_max_struct_r6/RESEARCH_REPORT_ZH.md:93-114`；## 5. 真正FFT接入：已运行，成本收益有界
对象/暴露：21；9historical+6historical+6fresh kept separate。
版本/表示：DPCS-locked Fourier/APCS mf3/6/10 then12; CIE/Dual sourcewise PR; Global global PR; A/A+ global PR unless sourcewise named;dense Fourier vs native FFT。初值：same host。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：main total CPU CIE897.22/Dual932.81/Global956.57/A1504.32; main elapsed1140.37;FFT55.36s。
估计器：cost+corrected。指标：CPU not wall。
正面及报告数字：A FFT CPU59.02→49.76 and70.58→54.44。
负面及边界：A still46%/43% dearer than FFTDual;A main1.61xDual;two FFT objects not new cohort。

### E009 · R7 · primary30
来源：`research/research_max_struct_r7/RESEARCH_REPORT_ZH.md:15-49`；## 1. 实际完成与分母
对象/暴露：30；independent-at-original-scoring。
版本/表示：DPCS-locked Fourier/APCS mf3/6/10 then12; CIE/Dual sourcewise PR; Global global PR; A/A+ global PR unless sourcewise named;final A+ arc quadratic proposals+priorMM+gtolguard。初值：uniform .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;24matched+6cross36→24;6/12sources;full/half;noise.1/1/.1/.1/.5%;bundled conditions。预算/时钟：3000aux+3000corrector;direct6000; no wall cap;CPU CIE1543.60/Dual1539.15/Global1512.97/A2495.14/A+1711.34/direct838.18。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint;support area>=.5;prediction Otsu。
正面及报告数字：CIE21/Dual20/Global20/A21/A+21/direct10 outof30; A+−Dual mean−.03778。
负面及边界：A+ CPU11.2%aboveDual;mean gains dominated few objects;A+IoU .69663 boundary loss;no broad successrate proof。

### E010 · R7 · primary30 paid prefixes
来源：`research/research_max_struct_r7/RESEARCH_REPORT_ZH.md:53-61`；## 3. 相同实际成本改变了什么判断
对象/暴露：30；independent-at-original-scoring。
版本/表示：same six versions。初值：uniform .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：per-pair common actualCPU/LU/RHS; float32 snapshots every10 physicaleval; sparse stage points。
估计器：paid minPhi prefix。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：A+vsDual common maxCPU/LU/RHS22/30vs20/30。
负面及边界：A vs Global commonCPU17/30vs20/30;A+ final21 differs prefix22;cannot promote prefix intoendpoint。

### E011 · R7 · postselected family20
来源：`research/research_max_struct_r7/CHAT_HANDOFF.md:18-20`；## 最后的家族确认
对象/暴露：20；new objects but family selection post30。
版本/表示：A+ final arcguard /Dual/Global/CIE。初值：uniform .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;new objects after selecting bar/shell families from30。预算/时钟：fullCPU1089.69/1061.57/1027.61/855.19;commonCPU separate。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：A+/Dual/Global/CIE12/9/9/8;bar9/8/8/5;shell3/1/1/3。
负面及边界：postselected families cannot merge into prospective50;8 A+ failures retained;single matches all A+ shell successes。

### E012 · R8 · 42 NN confirmation:ID24/OOD12/cross6
来源：`research/research_max_struct_r8/RESEARCH_REPORT_ZH.md:57-80`；## 3. 神经初始化与独立确认
对象/暴露：42；independent-at-original-scoring。
版本/表示：CIE/Dual/Global/A/A+;sourcewise/global A+ final arcguard;direct。初值：uniform vs frozen NN; optional projected physical current;fixed affine anchor separately。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;6/12sources;full/half;.1/1/3%noise;cross36→24。预算/时钟：aux3000+corrector3000;direct6000; hardware/cohort CPU stratified。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint;prediction Otsu;mask>.1+.5amplitude。
正面及报告数字：NNdirect means .0965/.1624/.0843 vs uniformdirect .2497/.4047/.1959。
负面及边界：NNglobalA+ .2961/.5939/.3293;NNsourcewiseA+ .3218/.5916/.3281;uniformDual .3684/.6383/.3282;no universal NN→aux benefit。

### E013 · R8 · 42 NN initializer anddirect
来源：`research/research_max_struct_r8/NEURAL_WARM_STARTS.md:5-57`；## 模型、数据和公开选择
对象/暴露：42；independent original test; now public。
版本/表示：24737param two-scale CNN; physics10features; fixed checkpoint。初值：1536train192val;120epochs;no testGT selection。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;trainingadditionalprior。预算/时钟：training+generation402.558889 Mac cumulativeCPU;not remeasuredR13;deploymentfeatures/inference charged。
估计器：NNraw→directfinal。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：NNraw .14690/.21898/.12052→direct .09647/.16245/.08431。
负面及边界：same NN→A+ often loses benefit;NNdirect finalCPU1.06–2.74xuniformdirect。

### E014 · R8 · fixed OOD r8_test_0027
来源：`research/research_max_struct_r8/RESEARCH_REPORT_ZH.md:82-86`；### 3.1 一个具体的路径反例
对象/暴露：1；predeclared illustration within42。
版本/表示：direct vs global/sourcewise A+ arcguard。初值：same NN; material init vsaux interface。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：direct accepted1814 PGTOL;aux commoncorrector3000 cap。
估计器：precorrector vs corrected。指标：NRMSE andPhi。
正面及报告数字：NNdirect Phi.00042866 NRMSE.12082。
负面及边界：A+ precursor Phi.01565<NNinitial.30965 butfinal.00046350 NRMSE1.01370;not explained simply bydirect6000cap。

### E015 · R8 · 42 fixed affine anchor andportfolio
来源：`research/research_max_struct_r8/CHAT_HANDOFF.md:23-25`；## 这轮真正发现了什么
对象/暴露：42；development mechanism+frozen confirmation extension。
版本/表示：NNfixedcurrentanchor;two complete branches minfinalPhi。初值：same frozenNN;Jd+(I-Pd)jphys anchor。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：both branchescharged;portfolioCPU2.38–3.17xNNdirect。
估计器：corrected/portfolio primary。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：anchor partly repairs OOD/cross;portfolioID mean.07754。
负面及边界：NNdirect+A+ portfolio OOD.34498/cross.22882 worse;8/8 initialdefect improvement notpredictive。

### E016 · R9 · R8 highcontrast audit andR7 reference
来源：`research/research_max_struct_r9/RESEARCH_REPORT_ZH.md:17-19`；## 1. 继承与对照资格
对象/暴露：26；reanalysis alreadypublic;26and18 neverpool。
版本/表示：R8 uniform CIE/Dual/Global/A/globalA+/sourcewiseA+ vsR7 highcontrast18。初值：uniform。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;different acquisition/noise/object contracts。预算/时钟：audit630R8 receipts;no rerun。
估计器：inherited final。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：R7 highcontrastA/A+11/18 vsDual/Global9/18。
负面及边界：R8 high26 counts6/3/5/5/4/5;A+−Global mean−.01505 median−.0000713;no A+exclusive success。

### E017 · R9 · new6 noise adaptation
来源：`research/research_max_struct_r9/RESEARCH_REPORT_ZH.md:23-43`；## 2. 物理、数据和评分合同
对象/暴露：6；independent-at-original-scoring。
版本/表示：old fixed rtol.001 DPCS rank17 vsCN99 adaptive rank;all four CIEvariants。初值：uniform。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;contrast2.5/12/18;noise.1/1/3%;12x48 full/half;matched24。预算/时钟：aux3000+corrector3000;direct6000。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：public noise rank releases shared currentdirections。
负面及边界：changescurrent feasible space;not A+exclusive;new6 small。

### E018 · R9 · new6 quantitative
来源：`research/research_max_struct_r9/RESEARCH_REPORT_ZH.md:76-98`；## 5. 新六材料：修复有效，A+独有优势没有出现
对象/暴露：6；independent-at-original-scoring。
版本/表示：CIE/Dual/Global/A+ old andnoiseadapted;relaxation;direct。初值：uniform。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;seeE016。预算/时钟：mean fullCPU old48.43/46.06/47.96/52.67;adapted38.84/40.81/42.98/46.37。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：allfour old2/6→adapted4/6;adaptedmeans .18666/.18882/.18766/.18552。
负面及边界：A+ no extrasuccess;19%moreCPUthansingle;relaxation2/6;portfolio4/6。

### E019 · R10 · main24
来源：`research/research_max_struct_r10/RESEARCH_REPORT_ZH.md:38-65`；## 2. 新24材料主确认：全部结果
对象/暴露：24；independent-at-original-scoring。
版本/表示：CIE/Dual/Global;A+global/sourcewise;G19;twofold30。初值：uniform;inheritedcurrentinterface。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;12x48;1%noise;8familiesxlow/high/mixed;same24。预算/时钟：3000aux+3000correctoraccepted/12000eval;fullchildCPU157.71/153.98/150.13/169.49/174.98/187.66/192.53。
估计器：primary corrected。指标：NRMSE<=.15 ANDfixedchi>1 IoU>=.7 ANDlegal;Otsusecondary。
正面及报告数字：CIE11/.16287;Dual10/.16346;Global10/.16369。
负面及边界：A+global8/.18044;sourcewise10/.18045;G198/.17669;Twofold8/.18099;not universaldenominatorwin。

### E020 · R10 · three main24 decisive cases
来源：`research/research_max_struct_r10/RESEARCH_REPORT_ZH.md:67-73`；### 三个决定性的新正反例
对象/暴露：3；subsetof24;notnew3。
版本/表示：CIE/Dual/A+global/G19/twofold。初值：uniform。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;sameasmain24。预算/时钟：samefullschedules;actualcostdifferent。
估计器：primary corrected。指标：amplitudeNRMSE。
正面及报告数字：ringCIE.42288/A+.42275→G19.09889。
负面及边界：shellCIE.11109→A+/G19.51514;bridgeCIE.09741/A+.09608→G19.33700。

### E021 · R10 · new4 conditional deployment
来源：`research/research_max_struct_r10/RESEARCH_REPORT_ZH.md:119-137`；### 4.3 补偿CIE，而不是强行替换所有材料读出
对象/暴露：4；independent-at-original-scoring。
版本/表示：CIE→G19 vsCIE→Dual; triggeroutsideCN99;retainminPhi。初值：CIEuniform run thenfallback fromprotocol;do not infer sameendpoint restart。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;fournewhighobjects。预算/时钟：deploymentwholechildCPU692.84/1015.55/1013.44;sharedCIEprefixnotadditive。
估计器：primary policy。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：CIE2/4;CIE→Dual2/4;CIE→G193/4;means.43222/.43222/.19711。
负面及边界：only2triggers;coatedbodyfails;noGlobalfallback;cannotuniquelyattributetoGreen。

### E022 · R11 · independent16
来源：`research/research_max_struct_r11/RESEARCH_REPORT_ZH.md:13-15`；## 最终判断
对象/暴露：16；independent-at-original-scoring。
版本/表示：CIE/Dual/Global/G19;samecurrentcapacity;G19rank19notfullG。初值：sameinit;inheritedinterface。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;high/lowmixed。预算/时钟：fullCPU G19/CIE1.269;commonLUdifferentfromCPU。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：CIE4/.23134;Dual5/.21068;Global6;G197/.14485;newringexclusivevsGlobal。
负面及边界：3newsuccessvsCIE notallGreen-specific;only1vsGlobal。

### E023 · R11 · independent48
来源：`research/research_max_struct_r11/RESEARCH_REPORT_ZH.md:7-11`；## 最终判断
对象/暴露：48；independent-at-original-scoring。
版本/表示：CIE/Dual/Global/G19/freezeA/contractedT。初值：sameinit;inheritedinterface。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;8families;high/mixedchi8/12/18。预算/时钟：effectivechildCPU49968.922;G19/Dual1.334;interruptionCPUunknowncase22。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：counts21/22/23/23/16/23;contractedTfixesbridge/combdamage。
负面及边界：Global/G19/contractedTidentical23successset;noG19exclusive;cannotpoolwith16。

### E024 · R11 · independent48 actualcostprefix
来源：`research/research_max_struct_r11/RESEARCH_REPORT_ZH.md:53-59`；## 3. 同成本、停止与材料误差的分开记账
对象/暴露：48；independent-at-original-scoring。
版本/表示：DualvsG19/contractedT;conditionalCIEportfolios。初值：sameinit。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：B=mincompletepaidCPU per pair;case22endtoendunknown;47sensitivity。
估计器：paid minPhi prefix。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：conditionalCIE→G19/Dual/Globalall23/48。
负面及边界：Dual/G19commonCPU22/21;Dual/contractedT22/20;endpointbestnotcostwin。

### E025 · R11 · paired8 materials two grids
来源：`research/research_max_struct_r11/RESEARCH_REPORT_ZH.md:11-11`；## 最终判断
对象/暴露：8；independent-at-original-scoring。
版本/表示：CIE/Dual/Global/G19。初值：sameinit。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;matched24andcross48→24 paired;targetfine2x2mean。预算/时钟：G19CPU20.8%/34.7%moreDual;total10921.40625s。
估计器：full primary +commonCPU。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint。
正面及报告数字：bothgridsCIE/Dual/Global4/8;G196/8;ringandshellrestoredfullendpoint。
负面及边界：commonCPUcrossincrementchanges: ring+mixedisland;not shell;not16independentmaterials。

### E026 · R11 · true3D8matched+2crossconditions
来源：`research/research_max_struct_r11/RESEARCH_REPORT_ZH.md:73-79`；## 4. 扩大确认与真正三维
对象/暴露：8；independent-at-original-scoring。
版本/表示：CIE/Dual/Global/fullGreenrank120/contracted;notG19。初值：same perobjecthost。
合同：12cube1728realchi;5184complexcurrent/source;60dipoles;CM+RR;3%noise;cross18→12。预算/时钟：all50reach320acceptedcorrectorcap;Mac/GPUclocknotpool。
估计器：primary。指标：relativechiL2 not2DamplitudeNRMSE。
正面及报告数字：matchedmeans .22844/.22608/.23127/.22780/.22774;fullGreen5/8lowerDual。
负面及边界：fullGreenoverallmeanhigherDual;cross .46483/.45439/.43356/.46714/.46319;allcrossoutsideball。

### E027 · R12 · independent24
来源：`research/research_max_struct_r12/RESEARCH_REPORT_ZH.md:66-101`；## 4. 独立24确认：收益没有成为通用优势
对象/暴露：24；independent-at-original-scoring。
版本/表示：CIE/Dual/Global/A+/G19/DR/physicalsmoothmultistart/NNdirect。初值：uniform exceptfrozenNNdirect;G19warmphysicalcurrent。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;1%CN;12x48;levelset/branch/island/continuous;chi8/18/mixed。预算/时钟：wholechildcap240;144proposal96corrector;actualcostunequal;prefix30/60/120/240。
估计器：primary corrected。指标：amplitude NRMSE=RMSE/declared contrast; material label <=.15 AND IoU>=.7 AND box legal; CN99 separate unless explicitly joint;supportfixedchi>1 whereinherited; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：means .202932/.185456/.184007/.171692/.143891/.184631/.619361/.149514;G19bestmean。
负面及边界：oldlabels9/9/9/9/9/7/4/9;DR11ballinsidefailuresIoUonly;notallnoimage。

### E028 · R12 · 24 CPU-prefix
来源：`research/research_max_struct_r12/RESEARCH_REPORT_ZH.md:87-101`；### 同成本前缀与配对不确定性
对象/暴露：24；independent-at-original-scoring。
版本/表示：CIE/Dual/A+/G19/DR/NNdirect。初值：asmain24。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：wholechild30/60/120/240 prefixes;notbodyCPU。
估计器：paid minPhi prefix。指标：all24 meanamplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：G19 .21718/.14988/.14411/.14389;NNdirect .15044/.14881/.14951/.14951。
负面及边界：DR .18551/.18396/.18408/.18463;DR−G19+.040740CI[+.001934,+.095095];NNadditionalprior。

### E029 · R12 · oldR11bad3 originalbodybudget
来源：`research/research_max_struct_r12/RESEARCH_REPORT_ZH.md:48-56`；## 3. 最清楚的正例：从已陷入的坏材料重新出发
对象/暴露：3；public exposed。
版本/表示：physicalcontinuation/smoothmultistart/AP/DR。初值：same savedbadmaterial fullphysicalJ/Q。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：paid240bodyCPU;notwholechild。
估计器：primary prefix。指标：amplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：DRring.16585/shell.11377/bridge.10617。
负面及边界：original lacks matchedG19restart;usepairedsupplementforG19comparison。

### E030 · R12 · NN4 originalinterface
来源：`research/research_max_struct_r12/RESEARCH_REPORT_ZH.md:120-122`；## 6. NN、实测及经典对象没有被跳过
对象/暴露：4；publicR8objects。
版本/表示：uniformdirect/NNdirect/NNDR/NNCIE。初值：frozenR8NN fullstate forDR;coldCIEinterface asrecorded。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：originallocalCPU;notpairedsupplementcost。
估计器：primary corrected。指标：amplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：means .4470/.1419/.0862/.3089。
负面及边界：noNNAPsinglefactororiginalarm;cannotattributeallgainstoreflection。

### E031 · R12 · 3Dcorrectorcounterexample
来源：`research/research_max_struct_r12/RESEARCH_REPORT_ZH.md:108-118`；## 5. 材料选择的问题：共同corrector为何会毁掉更好的结果
对象/暴露：2；development。
版本/表示：DRrawvscommoncorrector;CIE/Dual/APfinalnear。初值：sameDRcandidate。
合同：highchi8/mixed2and8;6cubevectorDDA;cross8→6cleanmismatch26.16%。预算/时钟：finiteoriginalcorrector;notgeneralconvergencetheorem。
估计器：raw→corrected。指标：relativechiL2;CNballqualification。
正面及报告数字：raw .06163/.07831 inball。
负面及边界：corrected.12978/.22724 outball;D increaseswhileRdecreases;moderategroupdifferencesvanishwithlongercorrector。

### E032 · R12 · pairedA oldbad3
来源：`research/research_max_struct_r12/supplement_paired/RESEARCH_REPORT_ZH.md:23-39`；## 2. A：同一个坏终点，谁能改善材料
对象/暴露：3；public exposed。
版本/表示：G19legalrestart/DR/AP/direct。初值：sameoldbadmaterial;G19physicalcurrentLSprojectedthenPRreset;APDRsamefullJ/Q。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：240wholechildcap;144/96;CPUs518.234/426.859/234.109/237.203。
估计器：primary originalPhi。指标：amplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：G19/DRring .114232/.140651;DRshell.112699vsG19.515137;DRbridge.097269vsG19.371690。
负面及边界：notallbadendpointsDRbeatsG19;G19restartnotunbrokencontinuation。

### E033 · R12 · pairedB full24policy fourtriggers
来源：`research/research_max_struct_r12/supplement_paired/RESEARCH_REPORT_ZH.md:39-58`；## 3. B：把救援连接到同一批24对象
对象/暴露：24；public exposed。
版本/表示：G19→G19/DR/AP/direct。初值：sameG19primary;20bypass;4outsideCN99restartfullphysicalstate。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：120wholechildcap72/48;originalG19CPU4290.1875;DR+336.875=4627.0625。
估计器：policy primary。指标：amplitudeNRMSE;oldjointcountsseparate; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：G19mean.143891→DR.106261;3/4triggeredimproveandenterball;oldjoint9→11。
负面及边界：G19/AP/direct .143883;0010stillfails;20bypassesnotDRrunnodamage;7.85%extraCPU。

### E034 · R12 · pairedC NN4
来源：`research/research_max_struct_r12/supplement_paired/RESEARCH_REPORT_ZH.md:58-78`；## 4. C：NN、完整场接口与反射的贡献
对象/暴露：4；public exposed。
版本/表示：NNG19/NNDR/NNAP/NNdirect。初值：samefrozenNN;G19warmLSphysicalcurrent;PRreset;APDRidenticalfullstate。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：240wholechildcap;CPUs790.359/544.328/550.375/374.297;sharedNNprep.763496separate。
估计器：primary;secondaryminRseparate。指标：amplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：means.086768/.086094/.187543/.142623;DR31%lessCPUthanG19。
负面及边界：AP0001.487912;AP0022.062564matchesDRgain;secondaryG19.080104vsDR.082901reversesmeanordering。

### E035 · R12 · pairedcorrector
来源：`research/research_max_struct_r12/supplement_paired/RESEARCH_REPORT_ZH.md:78-86`；## 5. Corrector：作用取决于到达的材料区域
对象/暴露：3；public exposed;A3/B24/C4notpooled。
版本/表示：DRraw/correctedA3/B24/C4。初值：samecandidate。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：paidcorrectorwindow;no newcorrector。
估计器：raw→corrected。指标：amplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：A.200448→.116873;Bfull24.110176→.106261;C.100916→.086094。
负面及边界：C0004.074037→.078824;C0022.055806→.067208despiteDandRdown。

### E036 · R13 · fresh48main
来源：`research/research_max_struct_r13/RESEARCH_REPORT_ZH.md:46-80`；## 4. 新 48 对象：同实际 CPU 的主比较
对象/暴露：48；independent-at-original-scoring。
版本/表示：CIE/Dual/G19/fullAPCS/fullcurrent/DR/event/random/direct;G19warmcurrent;CIE/Dualcold。初值：uniform .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;36matched12cross;6/12sources24/48receivers;1/2%noise;bundledconditions。预算/时钟：wholeprocessCPUcap120;proposal72corrector48;direct120;meanCPU106.15/106.26/119.08/117.24/115.86/84.50/86.33/84.50/59.46。
估计器：primary minimumeligiblepaidPhi。指标：continuousamplitudeNRMSE;fixedchi>1IoU;no newjointsuccessgate;ballnotrecovery; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：meansCIE.19454/Dual.17985/G19.16758/DR.17672;G19NRMSEbestofcore。
负面及边界：G19meanIoU.68577vsDual.71592;fullAPCS.29235/fullcurrent.30985;G19nativeonly1/48completes;crosscoarsetruthalloutsideball。

### E037 · R13 · samefresh48 stratum
来源：`research/research_max_struct_r13/RESEARCH_REPORT_ZH.md:64-80`；## 4. 新 48 对象：同实际 CPU 的主比较
对象/暴露：48；independent-at-original-scoring。
版本/表示：CIE/Dual/G19/DR。初值：uniform .1。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：same120cap。
估计器：primary。指标：amplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：matched36 .21391/.19437/.17334/.19979;cross12 .13643/.13628/.15029/.10751。
负面及边界：differentgrid/noise/acquisitionbundles;notpuregridcausality;cannotturnmeansintooneranking。

### E038 · R13 · fixed16mechanism NN/A+
来源：`research/research_max_struct_r13/RESEARCH_REPORT_ZH.md:82-100`；## 5. 事件、NN、A+：固定 16 对象机制审查
对象/暴露：16；independent-at-original-scoring。
版本/表示：CIE/Dual/G19/A+quinticnoguard/DR/event/NNDR/NNG19/NNdirect。初值：uniformvsfrozenR8NN;16everythirdID。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：same120wholeprocesscap;A+meanactual120.39overrunpaid。
估计器：primary。指标：amplitudeNRMSE;nonewjointgate; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：CIE.14920/Dual.12786/G19.15008/A+.10966/DR.12313/event.12189;NNDR.11936/NNG19.12656/NNdirect.27223。
负面及边界：16not48;NNDR−uniform−.00377CIincludes0;NNG19−.02352CIincludes0;NNtrainingextraadditionalprior。

### E039 · R13 · M4localreference16
来源：`research/research_max_struct_r13/RESEARCH_REPORT_ZH.md:76-80`；### 基线充分性与完整 current
对象/暴露：16；independent-at-original-scoring。
版本/表示：CIE/Dual/G19/DR/direct/physicalGlobalDual/FullJTRF。初值：uniformsame16IDs;notNN。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：360CPUonM4;nativeCIE/Dual/G19complete;hardwarelibrariesdifferfromi9。
估计器：primary。指标：amplitudeNRMSE; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：means .14521/.12464/.13737/.11509/.41081/.36171/.37670。
负面及边界：M4vs i9 notpurebudgetablation;physicalGlobalDualnotglobalPRDual。

### E040 · R13 · G19badendpoint rescue48policy
来源：`research/research_max_struct_r13/RESEARCH_REPORT_ZH.md:102-123`；## 6. G19 坏终点救援：同起点，完整 48 分母
对象/暴露：48；independent-at-original-scoring。
版本/表示：G19restart/fullcurrentrestart/directcontinuation/DR/event;posthocuniformDRportfoliosseparate。初值：9badG19primaries+39bypass;fullphysicalrestart;PRreset。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：extra120CPUpertrigger;fullCPU G195716/DR6385.66/event6349.20。
估计器：policy primary。指标：amplitudeNRMSE;ballnotrecovery; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：G19.16758→DR.15135/event.14947;twoobjectlargerecoveries。
负面及边界：G19restart.16762;fullcurrent.16763;direct.16769;posthocCIE/Dual→uniformDR .13377/.13317notnewconditionalqueue。

### E041 · R13 · core48corrector
来源：`research/research_max_struct_r13/RESEARCH_REPORT_ZH.md:125-133`；## 7. Corrector、驻点与候选选择
对象/暴露：48；independent-at-original-scoring。
版本/表示：CIE/Dual/G19/DR/event。初值：samepaidrawcandidates。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4。预算/时钟：commoncorrector48CPUwithin120。
估计器：raw→corrected。指标：amplitudeNRMSE;minimumPhi primary; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：raw→corrected .19148→.19454/.18295→.17985/.20908→.16758/.20067→.17672/.21822→.18149。
负面及边界：CIEmeanworsens;errorimprove/worsen33/15,35/13,39/9,40/8,41/7;nearestrawballrulesworse。

### E042 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv:4-4`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：cie;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：primary。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.1945359164780224;IoU_mean=0.6945573018135299;ball=35/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E043 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv:9-9`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：dual_cie;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：primary。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.17984533048443338;IoU_mean=0.7159174283411804;ball=37/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E044 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv:14-14`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：g19;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：primary。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.16758203911504158;IoU_mean=0.6857688472283993;ball=39/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E045 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv:19-19`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：dr;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：primary。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.17671933834812928;IoU_mean=0.687260322602408;ball=39/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E046 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv:24-24`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/METHOD_VIEW_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：a_plus;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：primary。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.15210017110259402;IoU_mean=0.7506265910285571;ball=40/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E047 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv:2-2`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：cie;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：common_actual_CPU。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.19468445986195548;IoU_mean=0.6947807585920515;ball=34/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E048 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv:4-4`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：dual_cie;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：common_actual_CPU。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.17978262252064484;IoU_mean=0.7165681857022953;ball=36/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E049 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv:6-6`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：g19;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：common_actual_CPU。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.18221555528147837;IoU_mean=0.6551056339386284;ball=34/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E050 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv:8-8`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：dr;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：common_actual_CPU。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.17671944219211944;IoU_mean=0.687260322602408;ball=39/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。

### E051 · R13 · supplement samefresh48
来源：`research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv:10-10`；research/research_max_struct_r13/supplement_aplus48/analysis/comparison48/COMMON_PUBLIC_BUDGET_SUMMARY.csv
对象/暴露：48；extension32 newruns onalreadypublicobjects+16existing;notnewblind48。
版本/表示：a_plus;a_plus=R13quinticnoguardfrozenvariant。初值：uniform .1;coldA+/CIE/Dual;warmG19。
合同：scalar-TM disk-pulse VIE; inverse24x24=576 real chi; box[.02,20]; original smooth/TV1e-4;samecore36matched12cross。预算/时钟：120wholeprocesscap72/48;commonactualCPUlimitperobjectminmethods whereviewnamed。
估计器：common_actual_CPU。指标：amplitudeNRMSE=RMSE/(max(coarse truth)-.1);ballcountisdataqualification; for fresh scalar confirmation RMSE/(max reference chi-.1), R12 confirmation report line7 and R13 score_results.py:49,67-82; do not substitute a nominal contrast。
正面及报告数字：NRMSE_mean=0.1574596474373644;IoU_mean=0.7326757284750928;ball=37/48。
负面及边界：no R7arcguard48experiment;no universalA+versionclaim。
